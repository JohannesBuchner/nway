#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

__doc__ = """Multiway association between astrometric catalogue. Use --help for usage.

Example: nway.py --radius 10 --prior-completeness 0.95 --mag GOODS:mag_H auto --mag IRAC:mag_irac1 auto cdfs4Ms_srclist_v3.fits :Pos_error CANDELS_irac1.fits 0.5 gs_short.fits 0.1 --out=out.fits
"""

import numpy
import numpy as np
#from numpy import log10, pi, exp, logical_and
import astropy.io.fits as pyfits
#import argparse
#import nwaylib.logger as logger
#import nwaylib.fastskymatch as match
#import nwaylib.bayesdistance as bayesdist
#import nwaylib.magnitudeweights as magnitudeweights
import nwaylib
import corner
import matplotlib.pyplot as plt

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

c2 = pyfits.getdata('doc/COSMOS_OPTICAL.fits', 1)
c3 = pyfits.getdata('doc/COSMOS_IRAC.fits', 1)

# generate many random combinations
N = 100000
c2indices = np.random.randint(len(c2), size=N)
c3indices = np.random.randint(len(c3), size=N)

X2 = np.vstack((c2['MAG'][c2indices], c3['mag_ch1'][c3indices]))
W2 = np.ones(N)
Y2 = np.zeros(N, dtype=bool)


def compute_ml_biases(match_tables, table, mag_include_radius, mag_exclude_radius, magauto_post_single_minvalue, store_mag_hists, logger):
	if mag_include_radius is not None:
		selection = table['Separation_max'].values < mag_include_radius
		selection_weights = numpy.ones(len(selection))
		selection_possible = table['Separation_max'].values < mag_exclude_radius
	else:
		selection = (table['dist_post'] > magauto_post_single_minvalue).values
		selection_weights = table['dist_post'].values[selection]
		selection_possible = (table['dist_post'] > 0.01).values

	secure_result = table[selection]

	# the candidates for which we will make predictions 
	# this is just the complement -- those we did not use for training
	insecure_result = table[~selection]
	Xp = np.vstack((c2['MAG'][insecure_result.OPT], c3['mag_ch1'][insecure_result.IRAC])).transpose()

	# good target training sample (class 1)
	X1 = np.vstack((c2['MAG'][secure_result.OPT], c3['mag_ch1'][secure_result.IRAC]))
	# add field training sample (class 0)
	X = np.hstack((X1, X2)).transpose()
	# define classes
	Y1 = np.ones(len(secure_result), dtype=bool)
	Y = np.hstack((Y1, Y2))
	assert len(X) == len(Y), (X.shape, Y.shape, X1.shape, X2.shape, Y1.shape, Y2.shape)
	# weights
	W = np.hstack((selection_weights, W2))
	logger.log("visualising training inputs ...")
	corner.corner(X1.transpose())
	plt.savefig("nwayml-train-target.pdf", bbox_inches='tight')
	plt.close()
	corner.corner(X2.transpose())
	plt.savefig("nwayml-train-field.pdf", bbox_inches='tight')
	plt.close()

	logger.log("training random forest using %d samples ..." % len(Y))
	# apply random forest classifier to distinguish target from field
	import sklearn.ensemble
	clf = sklearn.ensemble.RandomForestClassifier(n_estimators=400)
	clf.fit(X, Y, sample_weight=W)
	logger.log("predicting for %d samples..." % len(Xp))
	proba = clf.predict_proba(Xp)
	newbf = np.log(proba[:,1] + 1e-10) - np.log(proba[:,0] + 1e-10)

	insecure_interesting = selection_possible[~selection]
	newbf_p = newbf[insecure_interesting]
	Xp_p = Xp[insecure_interesting,:]

	# show the most extreme cases
	logger.log("visualising predictions ...")
	thresh_lo, thresh_hi = np.percentile(newbf_p, [10, 90])
	print(thresh_lo, thresh_hi, (newbf_p > thresh_hi).sum(), (newbf_p < thresh_lo).sum())
	#
	corner.corner(Xp_p[newbf_p >= thresh_hi,:])
	plt.suptitle("P>%.3f%%" % (np.exp(thresh_hi)*100))
	plt.savefig("nwayml-result-target.pdf", bbox_inches='tight')
	plt.close()
	corner.corner(Xp_p[newbf_p <= thresh_lo,:])
	plt.suptitle("P<%.3f%%" % (np.exp(thresh_lo)*100))
	plt.savefig("nwayml-result-field.pdf", bbox_inches='tight')
	plt.close()
	return newbf, ~selection


def apply_ml_biasing(match_tables, table, mag_include_radius, mag_exclude_radius, magauto_post_single_minvalue, store_mag_hists, logger):
	newbf, to_update = compute_ml_biases(match_tables, table, mag_include_radius, mag_exclude_radius, magauto_post_single_minvalue, store_mag_hists, logger)

	log_bf = table['dist_bayesfactor'].values
	total = log_bf.copy()
	total[to_update] = newbf
	
	return table, total

result = nwaylib.nway_match(
	[
	table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
	table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0, magnitude_columns=[('MAG', 'auto')]),
	table_from_fits('doc/COSMOS_IRAC.fits', poserr_value=0.5, area=2.0, magnitude_columns=[('mag_ch1', 'auto')]),
	],
	match_radius = 20,
	prior_completeness = 0.9,
	biasing_function=apply_ml_biasing,
)

"""

percentile = 50
p_any_threshold = np.percentile(result.prob_has_match[result.match_flag == 1], percentile)
mask_unrelated = result.prob_this_match < 0.01
mask = np.logical_and(result.prob_has_match > p_any_threshold, result.prob_this_match > 0.5)
mask_insecure = ~mask
secure_result = result[mask]
insecure_result = result[mask_insecure]

# generate many random combinations
N = 100000
c2indices = np.random.randint(len(c2), size=N)
c3indices = np.random.randint(len(c3), size=N)

X1 = np.vstack((c2['MAG'][secure_result.OPT], c3['mag_ch1'][secure_result.IRAC]))
W1 = secure_result.prob_has_match * secure_result.prob_this_match
X2 = np.vstack((c2['MAG'][c2indices], c3['mag_ch1'][c3indices]))
W2 = np.ones(N, dtype=bool)
X = np.hstack((X1, X2)).transpose()
W = np.hstack((W1, W2))
Y = np.hstack((np.ones(len(W1), dtype=bool),np.zeros(len(W2), dtype=bool)))
Xp = np.vstack((c2['MAG'][insecure_result.OPT], c3['mag_ch1'][insecure_result.IRAC])).transpose()

# apply random forest classifier to distinguish target from field
import corner
import matplotlib.pyplot as plt
corner.corner(X1.transpose())
plt.savefig("nwayml-train-target.pdf", bbox_inches='tight')
plt.close()
corner.corner(X2.transpose())
plt.savefig("nwayml-train-field.pdf", bbox_inches='tight')
plt.close()

print("training random forest...")
import sklearn.ensemble
rfc = sklearn.ensemble.RandomForestClassifier(n_estimators=400)
rfc.fit(X, Y, sample_weight=W)
Yp = rfc.predict_log_proba(Xp)[:,1]

corner.corner(Xp[Yp > np.log(0.9),:])
plt.savefig("nwayml-result-target.pdf", bbox_inches='tight')
plt.close()
corner.corner(Xp[Yp < np.log(0.1),:])
plt.savefig("nwayml-result-field.pdf", bbox_inches='tight')
plt.close()

# add Yp as a bayes factor weight
result.bf[mask_insecure] += Yp

nwaylib.nway_match(
	[
	table_from_fits('doc/COSMOS_XMM.fits', area=2.0),
	table_from_fits('doc/COSMOS_OPTICAL.fits', poserr_value=0.1, area=2.0, magnitude_columns=[('MAG', 'auto')]),
	table_from_fits('doc/COSMOS_IRAC.fits', poserr_value=0.5, area=2.0, magnitude_columns=[('mag_ch1', 'auto')]),
	],
	match_radius = 20,
	prior_completeness = 0.9)

ncats = len(match_tables)

table, resultstable, separations, errors = _create_match_table(match_tables, match_radius, logger=logger)

if not len(table) > 0:
	raise EmptyResultException('No matches.')

source_densities, source_densities_plus = _compute_source_densities(match_tables, logger=logger)

# first pass: find secure matches and secure non-matches

prior, log_bf = _compute_single_log_bf(match_tables, source_densities, source_densities_plus, table, separations, errors, prior_completeness, logger=logger)
table = table.assign(
	dist_bayesfactor_uncorrected=log_bf,
	dist_bayesfactor=log_bf
)

if consider_unrelated_associations:
	_correct_unrelated_associations(table, separations, errors, ncats, source_densities, source_densities_plus, logger=logger)
	log_bf = table['dist_bayesfactor']

# add the additional columns
post = bayesdist.posterior(prior, log_bf)
table = table.assign(dist_post=post)

# find magnitude biasing functions
table, total = _apply_magnitude_biasing(match_tables, table, mag_include_radius, mag_exclude_radius, magauto_post_single_minvalue, store_mag_hists, logger=logger)

table = _compute_final_probabilities(match_tables, table, prob_ratio_secondary, prior, total, logger=logger)

table = _truncate_table(table, min_prob, logger=logger)
"""
