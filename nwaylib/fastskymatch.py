"""
Functions for finding pairs within a certain search distance "err".

Very fast method based on hashing.
"""
from __future__ import print_function, division
import numpy
import itertools
from collections import defaultdict
import os
import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import astropy.units as u
from numpy import sin, cos, arctan2, hypot, arccos, arcsin, pi, exp, log
import healpy
import tqdm
import joblib
cachedir = 'cache'
if not os.path.isdir(cachedir): os.mkdir(cachedir)
mem = joblib.Memory(cachedir=cachedir, verbose=False)

def dist(apos, bpos):
	"""
	Angular separation (in degrees) between two points on a sphere.
	http://en.wikipedia.org/wiki/Great-circle_distance
	"""
	(a_ra, a_dec), (b_ra, b_dec) = apos, bpos
	lon1 = a_ra / 180 * pi
	lat1 = a_dec / 180 * pi
	lon2 = b_ra / 180 * pi
	lat2 = b_dec / 180 * pi
	sdlon = sin(lon2 - lon1)
	cdlon = cos(lon2 - lon1)
	slat1 = sin(lat1)
	slat2 = sin(lat2)
	clat1 = cos(lat1)
	clat2 = cos(lat2)

	num1 = clat2 * sdlon
	num2 = clat1 * slat2 - slat1 * clat2 * cdlon
	denominator = slat1 * slat2 + clat1 * clat2 * cdlon

	return arctan2(hypot(num1, num2), denominator) * 180 / pi

def dist3d(apos, bpos):
	"""
	Angular separation in ra & dec between two points on a sphere.
	"""
	(a_ra, a_dec), (b_ra, b_dec) = apos, bpos
	a_ra  = numpy.where(a_ra  == -99, numpy.nan, a_ra)
	a_dec = numpy.where(a_dec == -99, numpy.nan, a_dec)
	b_ra  = numpy.where(b_ra  == -99, numpy.nan, b_ra)
	b_dec = numpy.where(b_dec == -99, numpy.nan, b_dec)
	with numpy.errstate(divide='ignore', invalid='ignore'):
		a = SkyCoord(a_ra, a_dec, frame="icrs", unit="deg")
		b = SkyCoord(b_ra, b_dec, frame="icrs", unit="deg")
	
		localframe = SkyOffsetFrame(origin=a)
		na = a.transform_to(localframe)
		nb = b.transform_to(localframe)
		dra  = na.lon - nb.lon
		ddec = na.lat - nb.lat

		separation = a.separation(b).to(u.degree).value
	
		dra  = dra.to(u.degree).value
		ddec = ddec.to(u.degree).value
	
		return separation, dra, ddec

def get_tablekeys(table, name, tablename=''):
	keys = sorted(table.dtype.names, key=lambda k: 0 if k.upper() == name else 1 if k.upper().startswith(name) else 2)
	assert len(keys) > 0 and name in keys[0].upper(), 'ERROR: No "%s" column found in input catalogue "%s". Only have: %s' % (name, tablename, ', '.join(table.dtype.names))
	return keys[0]


def get_healpix_resolution_degrees(nside):
	resol = healpy.pixelfunc.nside2resol(nside) / pi * 180
	# according to monte carlo simulations, distances up to this factor
	# are completely contained within the pixel and its neighbors
	resfactor = 0.7
	return resfactor * resol

@mem.cache(ignore=['logger'])
def crossproduct(radectables, err, logger, pairwise_errs=[]):
	# check if away from the poles and RA=0
	use_flat_bins = True
	for ra, dec in radectables:
		if not(err < 1 and (ra > 10*err).all() and (ra < 360-10*err).all() and (numpy.abs(dec) < 45).all()):
			use_flat_bins = False
			break
	
	if use_flat_bins:
		logger.log('matching: using fast flat-sky approximation for this match')
	else:
		# choose appropriate nside for err (in deg)
		nside = 1
		for nside_next in range(30): 
			# largest distance still contained within pixels
			dist_neighbors_complete = get_healpix_resolution_degrees(2**nside_next)
			# we are looking for a pixel size which ensures bigger distances than the error radius
			# but we want the smallest pixels possible, to reduce the cartesian product
			if dist_neighbors_complete < err:
				# too small, do not accept
				# sources within err will be outside the neighbor pixels
				break
			nside = 2**nside_next
		resol = get_healpix_resolution_degrees(nside) * 60 * 60
		logger.log('matching: healpix hashing on pixel resolution ~ %f arcsec (nside=%d)' % (resol, nside))
	
	buckets = defaultdict(lambda : [[] for _ in range(len(radectables))])
	primary_cat_keys = None
	
	pbar = tqdm.tqdm(total=sum([len(t[0]) for t in radectables]))
	for ti, (ra_table, dec_table) in enumerate(radectables):
		if use_flat_bins:
			for ei, (ra, dec) in enumerate(zip(ra_table, dec_table)):
				i, j = int(ra / err), int(dec / err)
				
				# put in bucket, and neighbors
				for jj, ii in (j,i), (j,i+1), (j+1,i), (j+1, i+1):
					k = (ii, jj)
					# only primary catalogue is allowed to define new buckets
					if ti == 0 or k in buckets:
						buckets[k][ti].append(ei)
				pbar.update()
		else:
			# get healpixels
			ra, dec = ra_table, dec_table
			phi = ra / 180 * pi
			theta = dec / 180 * pi + pi/2.
			i = healpy.pixelfunc.ang2pix(nside, phi=phi, theta=theta, nest=True)
			j = healpy.pixelfunc.get_all_neighbours(nside, phi=phi, theta=theta, nest=True)
			# only consider four neighbours in one direction (N)
			# does not work, sometimes A is south of B, but B is east of A
			# so need to consider all neighbors, and deduplicate later
			neighbors = numpy.hstack((i.reshape((-1,1)), j.transpose()))
			
			# put in bucket, and neighbors
			if ti == 0:
				# only primary catalogue is allowed to define new buckets
				for ei, keys in enumerate(neighbors):
					for k in keys:
						buckets[k][ti].append(ei)
					pbar.update()
			else:
				for ei, keys in enumerate(neighbors):
					for k in keys:
						if k in primary_cat_keys:
							buckets[k][ti].append(ei)
					pbar.update()
			if ti == 0:
				primary_cat_keys = set(buckets.keys())
				
	pbar.close()
	
	# add no-counterpart options
	results = set()
	# now combine within buckets
	logger.log('matching: collecting from %d buckets, creating cartesian products ...' % len(buckets))
	#print('matching: %6d matches expected after hashing' % numpy.sum([
	#	len(lists[0]) * numpy.product([len(li) + 1 for li in lists[1:]]) 
	#		for lists in buckets.values()]))
	
	#pbar = logger.progress(ndigits=5, maxval=len(buckets)).start()
	pbar = tqdm.tqdm(total=len(buckets))
	while buckets:
		k, lists = buckets.popitem()
		pbar.update()
		# add for secondary catalogues the option of missing source
		for l in lists[1:]:
			l.insert(0, -1)
		# create the cartesian product
		local_results = itertools.product(*[sorted(l) for l in lists])
		
		# if pairwise filtering is requested, use it to trim down solutions
		if pairwise_errs:
			local_results = numpy.array(list(local_results))
			#nstart = len(local_results)
			for tablei, tablej, errij in pairwise_errs:
				indicesi = local_results[:,tablei]
				indicesj = local_results[:,tablej]
				# first find entries that actually have both entries
				mask_both = numpy.logical_and(indicesi >= 0, indicesj >= 0)
				#if not mask_both.any():
				#	continue
				# get the RA/Dec
				rai, deci = radectables[tablei]
				raj, decj = radectables[tablej]
				rai, deci = rai[indicesi[mask_both]], deci[indicesi[mask_both]]
				raj, decj = raj[indicesj[mask_both]], decj[indicesj[mask_both]]
				# compute distances
				mask_good = dist((rai, deci), (raj, decj)) < errij * 60 * 60
				# select the ones where one is missing, or those within errij
				mask_good2 = ~mask_both
				mask_good2[mask_both][mask_good] = True
				#print(mask_good2.sum(), mask_both.shape, mask_both.sum(), mask_good.sum())
				local_results = local_results[mask_good2,:]
			
			#print("compression:%d/%d" % (len(local_results), nstart))
			results.update([tuple(l) for l in local_results])
		else:
			results.update(local_results)
		del local_results
	pbar.close()

	n = len(results)
	logger.log('matching: %6d unique matches from cartesian product. sorting ...' % n)
	# now make results unique by sorting
	results = numpy.array(sorted(results))
	return results

# use preferred newer astropy command if available
if hasattr(pyfits.BinTableHDU, 'from_columns'):
	fits_from_columns = pyfits.BinTableHDU.from_columns
else:
	fits_from_columns = pyfits.new_table

def match_multiple(tables, table_names, err, fits_formats, logger, circular=True, pairwise_errs=[]):
	"""
	computes the cartesian product of all possible matches,
	limited to a maximum distance of err (in degrees).
	
	tables: input FITS table
	table_names: names of the tables
	fits_formats: FITS data type of the columns of each table
	
	returns 
	results: cartesian product of all possible matches (smaller than err)
	cat_columns: table with separation distances in arcsec
	header: which columns were used in each table for RA/DEC
	
	"""
	logger.log('')
	logger.log('matching with %f arcsec radius' % (err * 60 * 60))
	logger.log('matching: %6d naive possibilities' % numpy.product([len(t) for t in tables]))

	logger.log('matching: hashing')

	ra_keys = [get_tablekeys(table, 'RA', tablename=tablename) for table, tablename in zip(tables, table_names)]
	logger.log('    using RA  columns: %s' % ', '.join(ra_keys))
	dec_keys = [get_tablekeys(table, 'DEC', tablename=tablename) for table, tablename in zip(tables, table_names)]
	logger.log('    using DEC columns: %s' % ', '.join(dec_keys))

	ratables = [(t[ra_key], t[dec_key]) for t, ra_key, dec_key in zip(tables, ra_keys, dec_keys)]
	resultstable = crossproduct(ratables, err, logger=logger, pairwise_errs=pairwise_errs)
	results = resultstable.view(dtype=[(table_name, resultstable.dtype) for table_name in table_names]).reshape((-1,))

	keys = []
	for table_name, table in zip(table_names, tables):
		keys += ["%s_%s" % (table_name, n) for n in table.dtype.names]
	
	logger.log('merging in %d columns from input catalogues ...' % sum([1 + len(table.dtype.names) for table in tables]))
	cat_columns = []
	pbar = tqdm.tqdm(total=sum([1 + len(table.dtype.names) for table in tables]))
	for table, table_name, fits_format in zip(tables, table_names, fits_formats):
		tbl = table[results[table_name]]
		# set missing to nan
		mask_missing = results[table_name] == -1
		
		pbar.update()
		for n, format in zip(table.dtype.names, fits_format):
			k = "%s_%s" % (table_name, n)
			keys.append(k)
			col = tbl[n]
			#print('   setting "%s" to -99 (%d affected; column format "%s")' % (k, mask_missing.sum(), format))
			try:
				col[mask_missing] = -99
			except Exception as e:
				logger.log('   setting "%s" to -99 failed (%d affected; column format "%s"): %s' % (k, mask_missing.sum(), format, e))
			
			fitscol = pyfits.Column(name=k, format=format, array=col)
			cat_columns.append(fitscol)
			pbar.update()
	pbar.close()
	
	tbhdu = fits_from_columns(pyfits.ColDefs(cat_columns))
	header = dict(
		COLS_RA = ' '.join(["%s_%s" % (ti, ra_key) for ti, ra_key in zip(table_names, ra_keys)]),
		COLS_DEC = ' '.join(["%s_%s" % (ti, dec_key) for ti, dec_key in zip(table_names, dec_keys)])
	)
	
	logger.log('    adding angular separation columns')
	max_separation = numpy.zeros(len(results))
	for i in range(len(tables)):
		a_ra  = tbhdu.data["%s_%s" % (table_names[i], ra_keys[i])]
		a_dec = tbhdu.data["%s_%s" % (table_names[i], dec_keys[i])]
		for j in range(i):
			k = "Separation_%s_%s" % (table_names[i], table_names[j])
			k1 = k + "_ra"
			k2 = k + "_dec"
			if circular:
				keys += [k, k1, k2]
			else:
				keys += [k]
			
			b_ra  = tbhdu.data["%s_%s" % (table_names[j], ra_keys[j])]
			b_dec = tbhdu.data["%s_%s" % (table_names[j], dec_keys[j])]
			
			if circular:
				col = dist((a_ra, a_dec), (b_ra, b_dec))
			else:
				col, col_ra, col_dec = dist3d((a_ra, a_dec), (b_ra, b_dec))
			
			valid_input = numpy.logical_and(a_ra != -99, b_ra != -99)
			assert not numpy.isnan(col[valid_input]).any(), ['%d distances are nan' % numpy.isnan(col[valid_input]).sum(), 
				a_ra[numpy.isnan(col)], a_dec[numpy.isnan(col)], 
				b_ra[numpy.isnan(col)], b_dec[numpy.isnan(col)]]
			
			col[a_ra == -99] = numpy.nan
			col[b_ra == -99] = numpy.nan
			if not circular:
				col_ra[a_ra == -99] = numpy.nan
				col_ra[b_ra == -99] = numpy.nan
				col_dec[a_ra == -99] = numpy.nan
				col_dec[b_ra == -99] = numpy.nan
			max_separation = numpy.nanmax([col * 60 * 60, max_separation], axis=0)
			# store distance in arcsec 
			cat_columns.append(pyfits.Column(name=k, format='E', array=col * 60 * 60))
			if not circular:
				cat_columns.append(pyfits.Column(name=k1, format='E', array=col_ra * 60 * 60))
				cat_columns.append(pyfits.Column(name=k2, format='E', array=col_dec * 60 * 60))
	
	cat_columns.append(pyfits.Column(name="Separation_max", format='E', array=max_separation))
	cat_columns.append(pyfits.Column(name="ncat", format='I', array=(resultstable > -1).sum(axis=1)))
	keys.append("Separation_max")
	mask = max_separation < err * 60 * 60
	for c in cat_columns:
		c.array = c.array[mask]
	
	logger.log('matching: %6d matches after filtering by search radius' % mask.sum())
	logger.log('')
	return results[mask], cat_columns, header

def wraptable2fits(cat_columns, extname):
	tbhdu = fits_from_columns(pyfits.ColDefs(cat_columns))
	hdu = pyfits.PrimaryHDU()
	import datetime, time
	now = datetime.datetime.fromtimestamp(time.time())
	nowstr = now.isoformat()
	nowstr = nowstr[:nowstr.rfind('.')]
	hdu.header['DATE'] = nowstr
	hdu.header['ANALYSIS'] = 'NWAY matching'
	tbhdu.header['EXTNAME'] = extname
	hdulist = pyfits.HDUList([hdu, tbhdu])
	return hdulist

def array2fits(table, extname):
	cat_columns = pyfits.ColDefs([pyfits.Column(name=n, format='E',array=table[n]) 
		for n in table.dtype.names])
	return wraptable2fits(cat_columns, extname)


