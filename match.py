import numpy
import itertools
import os, sys
import pyfits
import progressbar
from numpy import sin, cos, arccos, pi, exp, log

# angular dist = arccos(  cos(90 - d1) + sin(90 - d1)*sin(90 - d2)*cos(a1 - a2) )
def dist((a_ra, a_dec), (b_ra, b_dec)):
	a1 = a_ra / 180 * pi
	d1 = a_dec / 180 * pi
	a2 = b_ra / 180 * pi
	d2 = b_dec / 180 * pi
	return arccos(sin(d1) * sin(d2) + cos(d1)*cos(d2)*cos(a1 - a2)) * 180 / pi

def get_tablekeys(table, name):
	keys = sorted(table.dtype.names, key=lambda k: 0 if k.upper() == name else 1 if k.upper().startswith(name) else 2)
	assert len(keys) > 0 and name in keys[0].upper(), 'ERROR: No "%s"  column found in input catalogue. Only have: %s' % (name, ', '.join(table.dtype.names))
	return keys[0]

def match_multiple(tables, table_names, err):
	buckets = {}

	print 'matching with %f arcsec radius' % (err * 60 * 60)
	print 'matching: %6d naive possibilities' % numpy.product([len(t) for t in tables])

	print 'matching: hashing'

	ra_keys = [get_tablekeys(table, 'RA') for table in tables]
	print 'using RA  columns: %s' % ', '.join(ra_keys)
	dec_keys = [get_tablekeys(table, 'DEC') for table in tables]
	print 'using DEC columns: %s' % ', '.join(dec_keys)
	#ra_min = min([table[k].min() for k, table in zip(ra_keys, tables)])
	#ra_max = max([table[k].max() for k, table in zip(ra_keys, tables)])
	#dec_min = min([table[k].min() for k, table in zip(dec_keys, tables)])
	#dec_max = max([table[k].max() for k, table in zip(dec_keys, tables)])
	#print ra_min, ra_max, dec_min, dec_max

	pbar = progressbar.ProgressBar(widgets=[
		progressbar.Percentage(), progressbar.Counter('%3d'),
		progressbar.Bar(), progressbar.ETA()], maxval=sum([len(t) for t in tables])).start()
	for ti, table in enumerate(tables):
		ra_key = ra_keys[ti]
		dec_key = dec_keys[ti]
		for ei, e in enumerate(table):
			ra, dec = e[ra_key], e[dec_key]
			i, j = int(ra / err), int(dec / err)
			
			#for jj, ii in (j-1,i-1), (j-1,i), (j-1,i+1), (j,i-1), (j,i), (j,i+1), (j+1,i-1), (j+1,i), (j+1, i+1):
			for jj, ii in (j,i), (j,i+1), (j+1,i), (j+1, i+1):
				#k = int(ii + 10000*jj)
				k = (ii, jj)
				if k not in buckets:
					buckets[k] = [[] for _ in range(len(tables))]
				#print ' bucket', k, 'table', ti, e
				buckets[k][ti].append(ei)
			pbar.update(pbar.currval + 1)
	pbar.finish()

	results = []
	def are_close(*l):
		# all pairs
		for i, a in enumerate(l):
			for j, b in enumerate(l[:i]):
				if numpy.abs(a[ra_keys[i]] - b[ra_keys[j]]) > err or \
					numpy.abs(a[dec_keys[i]] - b[dec_keys[j]]) > err:
					return False
				#if dist((a[ra_keys[i]], a[dec_keys[i]]),
				#	(b[ra_keys[j]], b[dec_keys[j]])) > err:
				#	return False
				
		return True
	
	# now combine in buckets
	print 'matching: %6d matches after hashing' % numpy.sum([numpy.product([len(li) for li in lists]) for lists in buckets.values()])

	print 'matching: collecting from buckets'
	pbar = progressbar.ProgressBar(widgets=[
		progressbar.Percentage(), progressbar.Counter('%3d'), 
		progressbar.Bar(), progressbar.ETA()], maxval=len(buckets)).start()
	for lists in buckets.values():
		if any([len(li) == 0 for li in lists]):
			continue
		#results += [tuple(pair) for pair in itertools.product(*lists) if are_close(*pair)]
		results += itertools.product(*lists)
		pbar.update(pbar.currval + 1)
	pbar.finish()

	# now make results unique by sorting
	print 'matching: %6d matches from hashing' % len(results)
	results = list(set(results))
	print 'matching: %6d unique matches' % len(results)

	keys = []
	for table_name, table in zip(table_names, tables):
		keys += ["%s_%s" % (table_name, n) for n in table.dtype.names]
	
	results = numpy.array(results, dtype=[(table_name, 'i') for table_name in table_names])
	
	print 'merging...'
	cat_columns = []
	pbar = progressbar.ProgressBar(widgets=[
		progressbar.Percentage(), progressbar.Counter('%3d'), 
		progressbar.Bar(), progressbar.ETA()], 
		maxval=sum([1 + len(table.dtype.names) for table in tables])).start()
	for table, table_name in zip(tables, table_names):
		tbl = table[results[table_name]]
		pbar.update(pbar.currval + 1)
		for n in table.dtype.names:
			k = "%s_%s" % (table_name, n)
			keys.append(k)
			cat_columns.append(pyfits.Column(name=k, format='E', array=tbl[n]))
			pbar.update(pbar.currval + 1)
	pbar.finish()
	
	tbhdu = pyfits.new_table(pyfits.ColDefs(cat_columns))
	
	print 'merging: adding angular separation columns'
	cols = []
	mask = numpy.zeros(len(results)) == 0
	for i in range(len(tables)):
		a_ra  = tbhdu.data["%s_%s" % (table_names[i], ra_keys[i])]
		a_dec = tbhdu.data["%s_%s" % (table_names[i], dec_keys[i])]
		for j in range(i):
			k = "Separation_%s_%s" % (table_names[i], table_names[j])
			keys.append(k)
			
			b_ra  = tbhdu.data["%s_%s" % (table_names[j], ra_keys[j])]
			b_dec = tbhdu.data["%s_%s" % (table_names[j], dec_keys[j])]
			
			col = dist((a_ra, a_dec), (b_ra, b_dec))
			
			mask = numpy.logical_and(mask, col < err)
			cat_columns.append(pyfits.Column(name=k, format='E', array=col))
	
	for c in cat_columns:
		c.array = c.array[mask]
	
	print 'matching: %6d matches after filtering' % mask.sum()
	
	return results, cat_columns

def wraptable2fits(cat_columns, extname):
	tbhdu = pyfits.new_table(pyfits.ColDefs(cat_columns))
	hdu = pyfits.PrimaryHDU()
	import datetime, time
	now = datetime.datetime.fromtimestamp(time.time())
	nowstr = now.isoformat()
	nowstr = nowstr[:nowstr.rfind('.')]
	hdu.header.update('CREATOR', """Johannes Buchner <jbuchner@mpe.mpg.de>""")
	hdu.header.update('DATE', nowstr)
	hdu.header.update('ANALYSIS', 'match test tables')
	tbhdu.header.update('EXTNAME', extname)
	hdulist = pyfits.HDUList([hdu, tbhdu])
	return hdulist

def array2fits(table, extname):
	cat_columns = pyfits.ColDefs([pyfits.Column(name=n, format='E',array=table[n]) 
		for n in table.dtype.names])
	return wraptable2fits(cat_columns, extname)

if __name__ == '__main__':
	filenames = sys.argv[1:]
	
	numpy.random.seed(0)

	def gen():
		ra  = numpy.random.uniform(size=1000)
		dec = numpy.random.uniform(size=1000)
		return numpy.array(zip(ra, dec), dtype=[('ra', 'f'), ('dec', 'f')])

	for fitsname in filenames:
		if os.path.exists(fitsname):
			continue
		print 'generating toy data for %s!' % fitsname
		hdulist = array2fits(gen(), fitsname.replace('.fits', ''))
		hdulist[0].header.update('GENERAT', 'match test table, random')
		hdulist.writeto(fitsname, clobber=False)
	
	tables = [pyfits.open(fitsname)[1] for fitsname in filenames]
	table_names = [t.name for t in tables]
	tables = [t.data for t in tables]
	
	err = 0.03
	
	results, columns = match_multiple(tables, table_names, err)
	tbhdu = pyfits.new_table(pyfits.ColDefs(columns))
	
	hdulist = wraptable2fits(tbhdu, 'MATCH')
	hdulist[0].header.update('ANALYSIS', 'match table from' + ', '.join(table_names))
	hdulist[0].header.update('INPUT', ', '.join(filenames))
	hdulist.writeto('match.fits', clobber=True)

	print 'plotting'
	import matplotlib.pyplot as plt
	for table_name in table_names:
		plt.plot(tbhdu.data['%s_RA'  % table_name], 
			 tbhdu.data['%s_DEC' % table_name], '.')
	
	for r in tbhdu.data:
		plt.plot(r['%s_RA'  % table_name],
			 r['%s_DEC' % table_name],
			'-', color='b', alpha=0.3)
	
	plt.savefig('matchtest.pdf')


