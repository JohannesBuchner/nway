"""
Functions for finding pairs within a certain search distance "err".

Very fast method based on hashing.
"""
from __future__ import print_function, division
from nwaylib.fastskymatch import dist, array2fits, match_multiple, wraptable2fits, fits_from_columns
import astropy.io.fits as pyfits
import numpy
import nwaylib.progress as progress
import nwaylib.logger as logger
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import numpy.testing as test

def test_dist():
	
	d = dist((53.15964508, -27.92927742), (53.15953445, -27.9313736))
	print('distance', d)
	assert not numpy.isnan(d)
	
	ra = numpy.array([ 53.14784241,  53.14784241,  53.14749908, 53.16559982,     53.19423676,  53.1336441 ])
	dec = numpy.array([-27.79363823, -27.79363823, -27.81790352, -27.79622459,      -27.70860672, -27.76327515])
	ra2 = numpy.array([ 53.14907837,  53.14907837,  53.1498642 , 53.16150284,     53.19681549,  53.13626862])
	dec2 = numpy.array([-27.79297447, -27.79297447, -27.81404877, -27.79223251,  -27.71365929, -27.76314735])

	d = dist((ra, dec), (ra2, dec2))
	print('distance', d)
	assert not numpy.isnan(d).any()

def run_match(nfiles, ngen=40):
	numpy.random.seed(0)

	def gen():
		ra  = numpy.random.uniform(size=ngen)
		dec = numpy.random.uniform(size=ngen)
		return numpy.array(list(zip(ra, dec)), dtype=[('ra', 'f'), ('dec', 'f')])
	
	filenames = ['test_input_%d.fits' % i for i in range(nfiles)]

	for fitsname in filenames:
		print('generating toy data for %s!' % fitsname)
		hdulist = array2fits(gen(), fitsname.replace('.fits', ''))
		hdulist[0].header['GENERAT'] = 'match test table, random'
		hdulist.writeto(fitsname, **progress.kwargs_overwrite_true)

	tables = [pyfits.open(fitsname)[1] for fitsname in filenames]
	table_names = [t.name for t in tables]
	tables = [t.data for t in tables]

	err = 0.03
	fits_formats = []
	for table in tables:
		fits_formats.append([c.format for c in table.columns])

	results, columns, header = match_multiple(tables, table_names, err, fits_formats, logger=logger.NormalLogger())
	tbhdu = fits_from_columns(pyfits.ColDefs(columns))

	hdulist = wraptable2fits(tbhdu, 'MATCH')
	hdulist[0].header['ANALYSIS'] = 'match table from' + ', '.join(table_names)
	hdulist[0].header['INPUT'] = ', '.join(filenames)
	hdulist.writeto('test_match%d.fits' % nfiles, **progress.kwargs_overwrite_true)

	print('plotting')
	tbhdu = pyfits.open('test_match%d.fits' % nfiles)[1]
	for table_name in table_names:
		ra = tbhdu.data['%s_RA'  % table_name]
		dec = tbhdu.data['%s_DEC'  % table_name]
		
		print(ra.shape, dec.shape)
		assert len(ra)>20
		assert len(dec)==len(ra)

def test_dist_astropyconsistent():
	# positions:
	dec = numpy.linspace(0, 90, 4000)[1:-1]
	ra = 0
	shift = 1e-5
	offset_ra  = shift
	offset_dec = 0 * shift

	# precision:
	sigma_ra  = 1 * u.arcsec
	sigma_dec = 1 * u.arcsec
	sigma_ra2  = sigma_ra / 10
	sigma_dec2 = sigma_dec / 10
	sigma_ra_tot = (sigma_ra**-2 + sigma_ra2**-2)**-0.5
	sigma_dec_tot = (sigma_dec**-2 + sigma_dec2**-2)**-0.5

	a = SkyCoord(ra, dec, frame="icrs", unit="deg")
	b = SkyCoord(ra+offset_ra, dec+offset_dec, frame="icrs", unit="deg")
	localframe = SkyOffsetFrame(origin=a)
	na = a.transform_to(localframe)
	nb = b.transform_to(localframe)
	
	# compute projected distance on the sky using our function
	d = dist((ra, dec), (ra+offset_ra, dec+offset_dec))
	# compute comparison with astropy
	diff = b.separation(a).to(u.deg)

	# compute simple flat version
	dra  = na.lon - nb.lon
	ddec = na.lat - nb.lat
	simplelocaldiff = (dra**2 + ddec**2)**0.5
	test.assert_almost_equal(diff / u.deg, simplelocaldiff / u.deg)

	

def test_match_2():
	run_match(2, ngen=1000)

def test_match_3():
	run_match(3, ngen=400)

def test_match_4():
	run_match(4, ngen=100)

def test_match_5():
	run_match(5, ngen=40)
