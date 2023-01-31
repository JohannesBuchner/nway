"""Multiway association between astrometric catalogues"""
from __future__ import print_function, division

import numpy
from numpy import log10, pi
import pandas
from collections import OrderedDict

from .logger import NullOutputLogger, NormalLogger
from . import fastskymatch as match
from . import bayesdistance as bayesdist
from . import magnitudeweights as magnitudeweights

__author__ = """Johannes Buchner"""
__email__ = 'johannes.buchner.acad@gmx.com'
__version__ = '4.5.4'


class UndersampledException(Exception):
	pass
class EmptyResultException(Exception):
	pass

def nway_match(match_tables, match_radius, prior_completeness,
	mag_include_radius=None, mag_exclude_radius=None, magauto_post_single_minvalue=0.9,
	prob_ratio_secondary = 0.5,
	min_prob=0., consider_unrelated_associations=True, 
	store_mag_hists=True,
	logger=NormalLogger()):
	"""
	match_tables: list of catalogues, each a dict with entries:
		- name (short catalog name, no spaces, used in output columns)
		- ra (RA in degrees)
		- dec (dec in degrees)
		- error (positional error in arcsec)
		- area (sky area covered by the catalogue in square degrees)
		- mags (list of additional columns to consider as priors)
		- magnames (short name for each entry in mags, used in output columns)
		- maghists (list of prior information for each entry in mags)
		  use None to automatically build a histogram of target/non-target sources.
		  Otherwise, supply the histogram manually: bins_lo, bins_hi, hist_sel, hist_all.
	
	match_radius: maximum radius in arcsec to consider.
		More distant counterparts are cut off. 
		Set to a very large value, e.g. 5 times the largest positional error.
		Setting to larger values does not change the result, but
		setting to smaller values improves performance.
	
	prior_completeness: expected fraction of sources in the primary catalogue 
		that are expected to have a counterpart (e.g., 90%).
		If an array is passed, completeness for each catalog. First
		entry has to be 1.
	
	mag_include_radius:
		search radius for building the magnitude histogram of target sources. If None (default), the Bayesian posterior is used
	
	mag_exclude_radius:
		exclusion radius for building the magnitude histogram of field sources. If None (default), mag_include_radius is used
	
	magauto_post_single_minvalue:
		minimum posterior probability (default: 0.9) for the magnitude histogram of secure target sources. Used in the Bayesian procedure.
	
	prob_ratio_secondary: for each primary catalogue entry, the most probable associations is marked with match_flag=1. Comparably good solutions receive match_flag=2 (all others are match_flag=0). The prob_ratio_secondary parameter controls the lowest probability ratio p(secondary)/p(primary) still considered comparable.
	
	min_prob: Minimum probability for truncating the output table.
		default: 0 (no truncation)
	
	consider_unrelated_associations: Correction for missing sources in 
		partial assocations. Keep to True.
	
	store_mag_hists: Write constructed mag hists to file (filename based on table name and mags).
	
	logger: NormalLogger for stderr output and progress bars, NullOutputLogger if silent
	"""
	if mag_exclude_radius is None:
		mag_exclude_radius = mag_include_radius
	if mag_include_radius is not None:
		if mag_include_radius >= match_radius:
			logger.warn('WARNING: magnitude radius is very large (>= matching radius). Consider using a smaller value.')

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
	
	return table


def _create_match_table(match_tables, match_radius, logger):
	# first match input catalogues, compute possible combinations in match_radius
	ratables = [(t['ra'], t['dec']) for t in match_tables]
	table_names = [t['name'] for t in match_tables]

	resultstable = match.crossproduct(ratables, match_radius / 60. / 60, logger=logger)
	#results = resultstable.view(dtype=[(t['name'], resultstable.dtype) for t in match_tables]).reshape((-1,))
	nresults = len(resultstable)
	keys = []
	columns = []
	for i, t in enumerate(match_tables):
		keys.append(t['name'])
		columns.append(resultstable[:,i])

	# matrix of separations
	separations = []
	errors = []
	invalid_separations = numpy.ones(nresults) * numpy.nan
	logger.log('    adding angular separation columns')
	max_separation = numpy.zeros(nresults)
	for i in range(len(match_tables)):
		a_ra  = ratables[i][0][resultstable[:,i]]
		a_dec = ratables[i][1][resultstable[:,i]]
		errors.append(match_tables[i]['error'][resultstable[:,i]])
		row = []
		for j in range(len(match_tables)):
			if i < j:
				b_ra  = ratables[j][0][resultstable[:,j]]
				b_dec = ratables[j][1][resultstable[:,j]]
				col = match.dist((a_ra, a_dec), (b_ra, b_dec))
				assert not numpy.isnan(col).any(), ['%d distances are nan' % numpy.isnan(col).sum(), 
					a_ra[numpy.isnan(col)], a_dec[numpy.isnan(col)], 
					b_ra[numpy.isnan(col)], b_dec[numpy.isnan(col)]]
				
				assert not (a_ra == b_ra).all()
				# store distance in arcsec
				col[resultstable[:,i] == -1] = numpy.nan
				col[resultstable[:,j] == -1] = numpy.nan
				col[a_ra == -99] = numpy.nan
				col[b_ra == -99] = numpy.nan
				col_arcsec = col * 60 * 60
				keys.append("Separation_%s_%s" % (table_names[i], table_names[j]))
				columns.append(col_arcsec)
				max_separation = numpy.nanmax([col_arcsec, max_separation], axis=0)
				row.append(col_arcsec)
			else:
				# mark the upper right triangle of the matrix as invalid data
				row.append(invalid_separations)
		separations.append(row)
		del row

	keys.append('Separation_max')
	columns.append(max_separation)
	keys.append('ncat')
	columns.append((resultstable > -1).sum(axis=1))

	# now truncate all columns:
	mask = max_separation < match_radius
	columns = [c[mask] for c in columns]
	errors = [e[mask] for e in errors]
	separations = [[cell[mask] for cell in row] for row in separations]
	#results = results[mask]
	resultstable = resultstable[mask,:]
	nresults = len(resultstable)

	logger.log('matching: %6d matches after filtering by search radius' % mask.sum())

	# now we have columns, which contains the distance information.
	assert len(columns) == len(keys), (len(columns), len(keys),  keys)

	table = pandas.DataFrame(OrderedDict(zip(keys, columns)))
	assert len(table) == nresults, (len(table), nresults,  mask.sum())

	return table, resultstable, separations, errors

def _compute_source_densities(match_tables, logger):
	source_densities = []
	source_densities_plus = []
	for i, match_table in enumerate(match_tables):
		n = len(match_table['ra'])
		area = match_table['area'] * 1.0 # in square degrees
		area_total = (4 * pi * (180 / pi)**2)
		density = n / area * area_total
		logger.log('%s "%s" (%d), density gives %.2e objects on entire sky' % ('Primary catalogue' if i == 0 else 'Catalogue', match_table['name'], n, density))
		# this takes into account that the source may be absent
		density_plus = (n + 1) / area * area_total
		source_densities.append(density)
		source_densities_plus.append(density_plus)

	# source can not be absent in primary catalogue
	source_densities_plus[0] = source_densities[0]
	source_densities_plus = numpy.array(source_densities_plus)
	source_densities = numpy.array(source_densities)
	return source_densities, source_densities_plus

def _compute_single_log_bf(match_tables, source_densities, source_densities_plus, table, separations, errors, prior_completeness, logger):
	logger.log('Computing distance-based probabilities ...')
	ncats = len(match_tables)

	if numpy.shape(prior_completeness) == ():
		prior_completeness = numpy.array([1.0] + [float(prior_completeness)**(1./(ncats-1)) for i in range(1, ncats)])
	
	if len(prior_completeness) != ncats:
		raise Exception('Prior completeness needs one value per catalog. Received "%s".' % prior_completeness)
	assert prior_completeness[0] == 1.0

	log_bf = numpy.zeros(len(table)) * numpy.nan
	prior = numpy.zeros(len(table)) * numpy.nan
	# handle all cases (also those with missing counterparts in some catalogues)
	for case in range(2**(ncats-1)):
		table_mask = numpy.array([True] + [(case // 2**(ti)) % 2 == 0 for ti in range(len(match_tables)-1)])
		# select those cases
		mask = True
		for i in range(1, len(match_tables)):
			if table_mask[i]: # require not nan
				mask = numpy.logical_and(mask, ~numpy.isnan(separations[0][i]))
			else:
				mask = numpy.logical_and(mask, numpy.isnan(separations[0][i]))
		# select errors
		errors_selected = [e[mask] for e, m in zip(errors, table_mask) if m]
		separations_selected = [[cell[mask] for cell, m in zip(row, table_mask) if m] 
			for row, m2 in zip(separations, table_mask) if m2]
		# here we should call the elliptical error variant if errors is a 2d array
		r = bayesdist.log_bf(separations_selected, errors_selected)
		assert r.shape == separations_selected[0][0].shape, (r.shape, separations_selected[0][0].shape)
		assert r.shape == errors_selected[0].shape, (r.shape, errors_selected[0].shape)
		assert r.shape == mask.sum(), (r.shape, mask.sum(), mask.shape)
		log_bf[mask] = r

		prior[mask] = source_densities[0] * numpy.product(prior_completeness[table_mask]) / numpy.product(source_densities_plus[table_mask])
		assert numpy.isfinite(prior[mask]).all(), (source_densities, prior_completeness[table_mask], numpy.product(source_densities_plus[table_mask]))


	assert numpy.isfinite(prior).all(), (prior, log_bf)
	assert numpy.isfinite(log_bf).all(), (prior, log_bf)
	return prior, log_bf

def _correct_unrelated_associations(table, separations, errors, ncats, source_densities, source_densities_plus, logger):
	logger.log('    correcting for unrelated associations ...')
	# correct for unrelated associations
	# identify those in need of correction
	# two unconsidered catalogues are needed for an unrelated association
	primary_id = table[table.columns[0]]
	#candidates = numpy.unique(primary_id[(ncat <= ncats - 2).values])
	pbar = logger.progress()
	for primary_id_key, group in pbar(table.groupby(primary_id, sort=False)):
		ncat_here = group['ncat']
		# list which ones we are missing
		best_logpost = 0
		# go through more complex associations
		for j in group[ncat_here > 2].index.values:
			missing_cats = [k for k, sep in enumerate(separations[0]) if numpy.isnan(sep[j])]
			# check if this association has sufficient overlap with the one we are looking for
			# it must contain at least two of the catalogues we are missing
			augmented_cats = []
			for k in missing_cats:
				if not numpy.isnan(separations[0][k][j]):
					augmented_cats.append(k)
			n_augmented_cats = len(augmented_cats)
			if n_augmented_cats >= 2:
				# ok, this is helpful.
				# identify the separations and errors
				separations_selected = [[[separations[k][k2][j]] for k2 in augmented_cats] for k in augmented_cats]
				errors_selected = [[errors[k][j]] for k in augmented_cats]
				# identify the prior
				prior_j = source_densities[augmented_cats[0]] / numpy.product(source_densities_plus[augmented_cats])
				# compute a log_bf
				log_bf_j = bayesdist.log_bf(numpy.array(separations_selected), numpy.array(errors_selected))
				logpost_j = bayesdist.unnormalised_log_posterior(prior_j, log_bf_j, n_augmented_cats)
				if logpost_j > best_logpost:
					#print('post:', logpost_j, log_bf_j, prior_j)
					best_logpost = logpost_j

		# ok, we have our correction factor, best_logpost
		# lets multiply it onto log_bf
		if best_logpost > 0:
			group['dist_bayesfactor'] += best_logpost

def _apply_magnitude_biasing(match_tables, table, mag_include_radius, mag_exclude_radius, magauto_post_single_minvalue, store_mag_hists, logger):
	biases = {}
	for i, t in enumerate(match_tables):
		table_name = t['name']
		for magvals, maghist, magname in zip(t['mags'], t['maghists'], t['magnames']):
			col_name = magname
			col = "%s_%s" % (table_name, col_name)
			mag = "%s:%s" % (table_name, col_name)
			logger.log('Incorporating bias "%s" ...' % mag)
			
			res = table[table.columns[i]].values
			res_defined = res != -1
			# get magnitudes of all
			# mark -99 as undefined
			mag_all = magvals
			mag_all[mag_all == -99] = numpy.nan
		
			# get magnitudes of selected
			mask_all = numpy.isfinite(mag_all)
			
		
			if maghist is None:
				if mag_include_radius is not None:
					selection = table['Separation_max'].values < mag_include_radius
					selection_possible = table['Separation_max'].values < mag_exclude_radius
					selection_weights = numpy.ones(len(selection))
				else:
					selection = (table['dist_post'] > magauto_post_single_minvalue).values
					selection_weights = table['dist_post'].values
					selection_possible = (table['dist_post'] > 0.01).values
				
				# ignore cases where counterpart is missing
				assert res_defined.shape == selection.shape, (res_defined.shape, selection.shape)
				selection = numpy.logical_and(selection, res_defined)
				selection_weights = selection_weights[res_defined]
				selection_possible = numpy.logical_and(selection_possible, res_defined)
				
				#print '   selection', selection.sum(), selection_possible.sum(), (-selection_possible).sum()
				
				#rows = results[table_name][selection].tolist()
				assert selection.shape == res.shape, (selection.shape, res.shape)
				rows, unique_indices = numpy.unique(res[selection], return_index=True)
				rows_weights = selection_weights[unique_indices]
				
				assert len(rows) > 0, 'No magnitude values within radius for "%s".' % mag
				mag_sel = magvals[rows]
				mag_sel_weights = rows_weights
				
				# remove vaguely possible options from alternative histogram
				assert selection_possible.shape == res.shape, (selection_possible.shape, res.shape)
				#print(res, selection, selection_possible)
				rows_possible = numpy.unique(res[selection_possible])
				mask_others = mask_all.copy()
				mask_others[rows_possible] = False
				
				# all options in the total (field+target sources) histogram
				mask_sel = ~numpy.logical_or(numpy.isnan(mag_sel), numpy.isinf(mag_sel))

				#print '      non-nans: ', mask_sel.sum(), mask_others.sum()

				logger.log('magnitude histogram of column "%s": %d secure matches, %d insecure matches and %d secure non-matches of %d total entries (%d valid)' % (col, mask_sel.sum(), len(rows_possible), mask_others.sum(), len(mag_all), mask_all.sum()))
				
				# make function fitting to ratio shape
				bins, hist_sel, hist_all = magnitudeweights.adaptive_histograms(mag_all[mask_others], mag_sel[mask_sel], weights=mag_sel_weights[mask_sel])
				if store_mag_hists:
					logger.log('magnitude histogram stored to "%s".' % (mag.replace(':', '_') + '_fit.txt'))
					with open(mag.replace(':', '_') + '_fit.txt', 'wb') as f:
						f.write(b'# lo hi selected others\n')
						numpy.savetxt(f,
							numpy.transpose([bins[:-1], bins[1:], hist_sel, hist_all]), 
							fmt = ["%10.5f"]*4)
				if mask_sel.sum() < 100:
					raise UndersampledException('ERROR: too few secure matches (%d) to make a good histogram. If you are sure you want to use this poorly sampled histogram, replace "auto" with the filename. You can also decrease the mag-auto-minprob parameter.' % mask_sel.sum())
			else:
				logger.log('magnitude histogramming: using user-supplied histogram for "%s"' % (col))
				bins_lo, bins_hi, hist_sel, hist_all = maghist #numpy.loadtxt(magfile).transpose()
				bins = numpy.array(list(bins_lo) + [bins_hi[-1]])
			func = magnitudeweights.fitfunc_histogram(bins, hist_sel, hist_all)
			if store_mag_hists:
				magnitudeweights.plot_fit(bins, hist_sel, hist_all, func, mag)
			magcol = magvals[res]
			magcol[~numpy.logical_and(res_defined, numpy.isfinite(magcol))] = -99
			#magcol[~res_defined] = -99
			weights = log10(func(magcol))
			# undefined magnitudes do not contribute
			weights[numpy.isnan(weights)] = 0
			biases[col] = weights

	# add the bias columns
	table = table.assign(**{'bias_%s' % col:10**weights for col, weights in biases.items()})
	log_bf = table['dist_bayesfactor'].values
	total = log_bf + sum(biases.values())
	
	return table, total

def _compute_final_probabilities(match_tables, table, prob_ratio_secondary, prior, total, logger):
	logger.log('')
	logger.log('Computing final probabilities ...')

	# add the posterior column
	post = bayesdist.posterior(prior, total)
	table = table.assign(p_single=post)

	# compute weights for group posteriors
	# 4pi comes from Eq. 
	ncat = table['ncat'].values
	log_post_weight = bayesdist.unnormalised_log_posterior(prior, total, ncat)
	table = table.assign(log_post_weight=log_post_weight)

	# flagging of solutions. Go through groups by primary id (IDs in first catalogue)

	table = table.assign(
		match_flag = numpy.zeros(len(table), dtype=int),
		prob_has_match = numpy.zeros(len(table)),
		prob_this_match = numpy.zeros(len(table)),
	)

	logger.log('    grouping by primary catalogue ID and flagging ...')

	def compute_group_statistics(group):
		# group
		# compute no-match probability
		values = group['log_post_weight'].values
		#values = grpvalues.values
		offset = values.max()
		bfsum = log10((10**(values - offset)).sum()) + offset
		if len(values) > 1:
			offset = values[1:].max()
			bfsum1 = log10((10**(values[1:] - offset)).sum()) + offset
		else:
			bfsum1 = 0
		assert group['ncat'].values[0] == 1, group['ncat']
		#print(group['ncat'].values[0])
		# for p_any, find the one without counterparts
		p_none = values[0]
		p_any = 1 - 10**(p_none - bfsum)
		# this avoids overflows in the no-counterpart solution, 
		# which we want to set to 0
		values[0] = bfsum1
		p_i = 10**(values - bfsum1)
		p_i[0] = 0
		
		best_val = p_i.max()
		
		# flag best & second best
		# ignore very poor solutions
		match_flag = numpy.where(best_val == p_i, 1, 
			numpy.where(p_i > prob_ratio_secondary * best_val, 2, 0))
		
		group['prob_has_match'] = p_any
		group['prob_this_match'] = p_i
		group['match_flag'] = match_flag
		return group

	table = table.groupby(table.columns[0], sort=False).apply(compute_group_statistics)
	del table['log_post_weight']
	return table

def _truncate_table(table, min_prob, logger):
	# cut away poor posteriors if requested
	if min_prob > 0:
		mask = ~(table['prob_this_match'] < min_prob)
		logger.log('    cutting away %d (below p_i minimum)' % (len(mask) - mask.sum()))
		table = table[mask]

	return table
