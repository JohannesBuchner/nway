#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Demonstration of how to use the Python API of NWAY."""

from __future__ import print_function, division

import numpy
import astropy.io.fits as pyfits
import nwaylib.logger as logger
import nwaylib


# set up program arguments

def table_from_fits(fitsname, poserr_value=None, area=None, magnitude_columns=[]):
	"""Convert FITS input to the way the API can understand."""
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
	
	# area in square degrees
	# error in arcsec
	# ra/dec in degrees
	# mag: column of something
	# maghists: either (bin, sel, all) tuple or None (for auto)
	return dict(name=table_name, ra=ra, dec=dec, error=poserr, area=area, mags=mags, maghists=maghists, magnames=magnames)


# run example 1 from the manual

result = nwaylib.nway_match(
	[
		table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
		table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0),
	],
	match_radius=20,  # in arcsec
	prior_completeness=0.9,
)
assert len(result) == 37836, (len(result), result)
min_output_columns = ['Separation_max', 'ncat', 'dist_bayesfactor', 'dist_post', 'p_single', 'prob_has_match', 'prob_this_match']
extra_columns = ['Separation_XMM_OPT']
for col in min_output_columns + extra_columns:
	assert col in result.columns, ('looking for', col, 'in', result.columns)

# run example 2 from the manual

result = nwaylib.nway_match(
	[
		table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
		table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0, magnitude_columns=[('MAG', 'auto')]),
	],
	match_radius=20,  # in arcsec
	prior_completeness=0.9,
	store_mag_hists=False,
	mag_include_radius=4.0,  # in arcsec
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
	match_radius=20,  # in arcsec
	prior_completeness=0.9,
	mag_include_radius=4.0,  # in arcsec
)
assert len(result) == 37836, (len(result), result)
for col in min_output_columns + extra_columns:
	assert col in result.columns, ('looking for', col, 'in', result.columns)

# run example 3 from the manual

result = nwaylib.nway_match(
	[
		table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
		table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0, magnitude_columns=[('MAG', 'auto')]),
		table_from_fits('doc/COSMOS_IRAC.fits', poserr_value=0.5, area=2.0, magnitude_columns=[('mag_ch1', 'auto')]),
	],
	match_radius=20,  # in arcsec
	prior_completeness=0.9,
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
	match_radius=20,  # in arcsec
	prior_completeness=0.9,
	logger=logger.NullOutputLogger(),
)
assert len(result) == 387601, (len(result), result)
for col in min_output_columns + extra_columns:
	assert col in result.columns, ('looking for', col, 'in', result.columns)
