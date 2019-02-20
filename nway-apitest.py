#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

__doc__ = """Multiway association between astrometric catalogue. Use --help for usage.

Example: nway.py --radius 10 --prior-completeness 0.95 --mag GOODS:mag_H auto --mag IRAC:mag_irac1 auto cdfs4Ms_srclist_v3.fits :Pos_error CANDELS_irac1.fits 0.5 gs_short.fits 0.1 --out=out.fits
"""

import sys
import numpy
from numpy import log10, pi, exp, logical_and
import astropy.io.fits as pyfits
import argparse
import nwaylib.logger as logger
import nwaylib.fastskymatch as match
import nwaylib.bayesdistance as bayesdist
import nwaylib.magnitudeweights as magnitudeweights
import nwaylib

# set up program arguments

def table_from_fits(fitsname, poserr_value=None, area=None, magnitude_columns=[]):
	fits_table = pyfits.open(fitsname)[1]
	table_name = fits_table.name
	ra = fits_table.data['RA']
	dec = fits_table.data['DEC']
	if 'pos_err' in fits_table.data.columns.names:
		poserr = fits_table.data['pos_err']
	else:
		assert poserr_value is not None, ('"pos_err" column not found in file "%s", and no poserr_value passed' % fitsname)
		poserr = poserr_value * numpy.ones(len(ra))
	if area is None:
		area = fits_table.header['SKYAREA'] * 1.0 
	
	# magnitude columns
	mags = []
	maghists = []
	magnames = []
	#for mag in magnitude_columns:
	for col_name, magfile in magnitude_columns:
		assert col_name in fits_table.data.dtype.names
		
		mag_all = fits_table.data[col_name]
		# mark -99 as undefined
		mag_all[mag_all == -99] = numpy.nan
		
		mags.append(mag_all)
		magnames.append(col_name)
		if magfile == 'auto':
			maghists.append(None)
		else:
			bins_lo, bins_hi, hist_sel, hist_all = numpy.loadtxt(magfile).transpose()
			maghists.append((bins_lo, bins_hi, hist_sel, hist_all))
	
	return dict(name=table_name, ra=ra, dec=dec, error=poserr, area=area, mags=mags, maghists=maghists, magnames=magnames)
	# area in square degrees
	# error in arcsec
	# ra/dec in degrees
	# mag: column of something
	# maghists: either (bin, sel, all) tuple or None (for auto)

result = nwaylib.nway_match(
	[
	table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
	table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0),
	],
	match_radius = 20, # in arcsec
	prior_completeness = 0.9,
)
assert len(result) == 37836, (len(result), result)
min_output_columns = ['Separation_max', 'ncat', 'dist_bayesfactor', 'dist_post', 'p_single', 'prob_has_match', 'prob_this_match']
extra_columns = ['Separation_XMM_OPT']
for col in min_output_columns + extra_columns:
	assert col in result.columns, ('looking for', col, 'in', result.columns)

result = nwaylib.nway_match(
	[
	table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
	table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0, magnitude_columns=[('MAG', 'auto')]),
	],
	match_radius = 20, # in arcsec
	prior_completeness = 0.9,
	store_mag_hists = False,
	mag_include_radius = 4.0, # in arcsec
)
assert len(result) == 37836, (len(result), result)
extra_columns = ['Separation_XMM_OPT', 'bias_OPT_MAG']
for col in min_output_columns + extra_columns:
	assert col in result.columns, ('looking for', col, 'in', result.columns)

result = nwaylib.nway_match(
	[
	table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
	table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0, magnitude_columns=[('MAG', 'auto')]),
	],
	match_radius = 20, # in arcsec
	prior_completeness = 0.9,
	mag_include_radius = 4.0, # in arcsec
)
assert len(result) == 37836, (len(result), result)
for col in min_output_columns + extra_columns:
	assert col in result.columns, ('looking for', col, 'in', result.columns)

result = nwaylib.nway_match(
	[
	table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
	table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0, magnitude_columns=[('MAG', 'auto')]),
	table_from_fits('doc/COSMOS_IRAC.fits', poserr_value=0.5, area=2.0, magnitude_columns=[('mag_ch1', 'auto')]),
	],
	match_radius = 20,
	prior_completeness = 0.9,
)
assert len(result) == 387601, (len(result), result)
extra_columns = ['Separation_XMM_OPT', 'Separation_OPT_IRAC', 'Separation_XMM_IRAC', 'bias_OPT_MAG', 'bias_IRAC_mag_ch1']
for col in min_output_columns + extra_columns:
	assert col in result.columns, ('looking for', col, 'in', result.columns)

result = nwaylib.nway_match(
	[
	table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
	table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0, magnitude_columns=[('MAG', 'OPT_MAG_fit.txt')]),
	table_from_fits('doc/COSMOS_IRAC.fits', poserr_value=0.5, area=2.0, magnitude_columns=[('mag_ch1', 'IRAC_mag_ch1_fit.txt')]),
	],
	match_radius = 20,
	prior_completeness = 0.9,
	logger=logger.NullOutputLogger(),
)
assert len(result) == 387601, (len(result), result)
for col in min_output_columns + extra_columns:
	assert col in result.columns, ('looking for', col, 'in', result.columns)


"""
if not filenames[0].endswith('shifted.fits'):
	print()
	print()
	print('  You can calibrate a p_any cut-off with the following steps:')
	print('   1) Create a offset catalogue to simulate random sky positions:')
	shiftfile = filenames[0].replace('.fits', '').replace('.FITS', '') + '-fake.fits'
	shiftoutfile = outfile + '-fake.fits'
	print('      nway-create-fake-catalogue.py --radius %d %s %s' % (args.radius*2, filenames[0], shiftfile))
	print('   2) Match the offset catalogue in the same way as this run:')
	newargv = []
	i = 0
	while i < len(sys.argv):
		v = sys.argv[i]
		if v == filenames[0]:
			newargv.append(shiftfile)
		elif v == '--mag':
			newargv.append(v)
			v = sys.argv[i+1]
			newargv.append(v)
			if sys.argv[i+2] == 'auto':
				newargv.append(v.replace(':', '_') + '_fit.txt')
			else:
				newargv.append(sys.argv[i+2])
			i = i + 2
		elif v == '--out':
			newargv.append(v)
			i = i + 1
			newargv.append(shiftoutfile)
		elif v.startswith('--out='):
			newargv.append('--out=' + shiftoutfile)
		else:
			newargv.append(v)
		i = i + 1
	print('      ' + ' '.join(newargv))
	print('   3) determining the p_any cutoff that corresponds to a false-detection rate')
	print('      nway-calibrate-cutoff.py %s %s' % (outfile, shiftoutfile))
	print()
	
# write out fits file
print()
print('creating output FITS file ...')
tbhdu = match.fits_from_columns(pyfits.ColDefs(columns))

hdulist = match.wraptable2fits(tbhdu, 'MULTIMATCH')
hdulist[0].header['METHOD'] = 'multi-way matching'
hdulist[0].header['INPUT'] = ', '.join(filenames)
hdulist[0].header['TABLES'] = ', '.join(table_names)
hdulist[0].header['BIASING'] =  ', '.join(biases.keys())
hdulist[0].header['NWAYCMD'] = ' '.join(sys.argv)
for k, v in args.__dict__.items():
	hdulist[0].header.add_comment("argument %s: %s" % (k, v))
hdulist[0].header.update(match_header)
print('    writing "%s" (%d rows, %d columns) ...' % (outfile, len(tbhdu.data), len(columns)))
hdulist.writeto(outfile, **progress.kwargs_overwrite_true)

import nwaylib.checkupdates
nwaylib.checkupdates.checkupdates()


"""
