#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

__doc__ = """Multiway association between astrometric catalogue. Use --help for usage.

Example: nway.py --radius 10 --prior-completeness 0.95 --mag GOODS:mag_H auto --mag IRAC:mag_irac1 auto cdfs4Ms_srclist_v3.fits :Pos_error CANDELS_irac1.fits 0.5 gs_short.fits 0.1 --out=out.fits
"""

import sys
import numpy
from numpy import log10, pi, exp, logical_and
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import argparse
from nwaylib import progress
import nwaylib.fastskymatch as match
import nwaylib.bayesdistance as bayesdist
import nwaylib.magnitudeweights as magnitudeweights

# set up program arguments

class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Johannes Buchner (C) 2013-2017 <johannes.buchner.acad@gmx.com>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--radius', type=float, required=True,
	help='exclusive search radius in arcsec for initial matching')

parser.add_argument('--mag-radius', default=None, type=float,
	help='search radius for building the magnitude histogram of target sources. If not set, the Bayesian posterior is used.')

parser.add_argument('--mag-exclude-radius', default=None, type=float,
	help='exclusion radius for building the magnitude histogram of field sources. If not set, --mag-radius is used.')

parser.add_argument('--prior-completeness', metavar='COMPLETENESS', default=1, type=float,
	help='expected matching completeness of sources (prior)')

parser.add_argument('--ignore-unrelated-associations', dest='consider_unrelated_associations', action='store_false',
	help='Ignore in the calculation source pairings unrelated to the primary source (not recommended)')

parser.set_defaults(consider_unrelated_associations=True)

parser.add_argument('--mag', metavar='MAGCOLUMN+MAGFILE', type=str, nargs=2, action='append', default=[],
	help="""name of <table>:<column> for magnitude biasing, and filename for magnitude histogram 
	(use auto for auto-computation within mag-radius).
	Example: --mag GOODS:mag_H auto --mag IRAC:mag_irac1 irac_histogram.txt""")

parser.add_argument('--acceptable-prob', metavar='PROB', type=float, default=0.5,
	help='ratio limit up to which secondary solutions are flagged')

parser.add_argument('--min-prob', type=float, default=0,
	help='lowest probability allowed in final catalogue. If 0, no trimming is performed.')

parser.add_argument('--out', metavar='OUTFILE', help='output file name', required=True)

parser.add_argument('catalogues', type=str, nargs='+',
	help="""input catalogue fits files and position errors.
	Example: cdfs4Ms_srclist_v3.fits :Pos_error CANDELS_irac1.fits 0.5 gs_short.fits 0.1
	""")


# parsing arguments
args = parser.parse_args()

print('NWAY arguments:')

diff_secondary = args.acceptable_prob
outfile = args.out

filenames = args.catalogues[::2]
print('    catalogues: ', ', '.join(filenames))
pos_errors = args.catalogues[1::2]
print('    position errors/columns: ', ', '.join(pos_errors))

fits_tables = []
table_names = []
tables = []
source_densities = []
source_densities_plus = []
fits_formats = []
for fitsname in filenames:
	fits_table = pyfits.open(fitsname)[1]
	fits_tables.append(fits_table)
	table_name = fits_table.name
	table_names.append(table_name)
	table = fits_table.data
	fits_formats.append([c.format for c in fits_table.columns])
	tables.append(table)

	n = len(table)
	assert 'SKYAREA' in fits_table.header, 'file "%s", table "%s" does not have a field "SKYAREA", which should contain the area of the catalogue in square degrees' % (fitsname, table_name)
	area = fits_table.header['SKYAREA'] * 1.0 # in square degrees
	area_total = (4 * pi * (180 / pi)**2)
	density = n / area * area_total
	print('      from catalogue "%s" (%d), density gives %.2e on entire sky' % (table_name, n, density))
	# this takes into account that the source may be absent
	density_plus = (n + 1) / area * area_total
	source_densities.append(density)
	source_densities_plus.append(density_plus)

# source can not be absent in primary catalogue
source_densities_plus[0] = source_densities[0]
source_densities_plus = numpy.array(source_densities_plus)
min_prob = args.min_prob

match_radius = args.radius / 60. / 60 # in degrees
#if args.mag_radius is not None:
#	mag_radius = match_radius # in arc sec
#else:
mag_include_radius = args.mag_radius # in arc sec
mag_exclude_radius = args.mag_exclude_radius # in arc sec
if mag_exclude_radius is None:
	mag_exclude_radius = mag_include_radius

magnitude_columns = args.mag
print('    magnitude columns: ', ', '.join([c for c, _ in magnitude_columns]))

for mag, magfile in magnitude_columns:
	table_name, col_name = mag.split(':', 1)
	assert table_name in table_names, 'table name specified for magnitude ("%s") unknown. Known tables: %s' % (table_name, ', '.join(table_names))
	ti = table_names.index(table_name)
	col_names = tables[ti].dtype.names
	assert col_name in col_names, 'column name specified for magnitude ("%s") unknown. Known columns in table "%s": %s' % (mag, table_name, ', '.join(col_names))


# first match input catalogues, compute possible combinations in match_radius
results, columns, match_header = match.match_multiple(tables, table_names, match_radius, fits_formats)
table = match.fits_from_columns(pyfits.ColDefs(columns)).data

assert len(table) > 0, 'No matches.'

# first pass: find secure matches and secure non-matches
print('Computing distance-based probabilities ...')

print('  finding position error columns ...')
# get the separation and error columns for the bayesian weighting
errors    = []
for ti, (table_name, pos_error) in enumerate(zip(table_names, pos_errors)):
	if pos_error[0] == ':':
		# get column
		k = "%s_%s" % (table_name, pos_error[1:])
		assert k in table.dtype.names, 'ERROR: Position error column for "%s" not in table "%s". Have these columns: %s' % (k, table_name, ', '.join(table.dtype.names))
		table_errors = tables[ti][pos_error[1:]]
		print('    Position error for "%s": found column %s: Values are [%f..%f]' % (table_name, k, table_errors.min(), table_errors.max()))
		if table_errors.min() <= 0:
			print('WARNING: Some separation errors in "%s" are 0! This will give invalid results (%d rows).' % (k, (table_errors <= 0).sum()))
		if table_errors.max() > match_radius * 60 * 60:
			print('WARNING: Some separation errors in "%s" are larger than the match radius! Increase --radius to >> %s' % (k, table_errors.max()))
		errors.append(table[k])
	else:
		print('    Position error for "%s": using fixed value %f' % (table_name, float(pos_error)))
		if float(pos_error) > match_radius * 60 * 60:
			print('WARNING: Given separation error for "%s" is larger than the match radius! Increase --radius to >> %s' % (k, float(pos_error)))
		errors.append(float(pos_error) * numpy.ones(len(table)))

print('  finding position columns ...')
# table is in arcsec, and therefore separations is in arcsec
separations = []
for ti, a in enumerate(table_names):
	row = []
	for tj, b in enumerate(table_names):
		if ti < tj:
			k = 'Separation_%s_%s' % (b, a)
			assert k in table.dtype.names, 'ERROR: Separation column for "%s" not in merged table. Have columns: %s' % (k, ', '.join(table.dtype.names))
			row.append(table[k])
		else:
			row.append(numpy.ones(len(table)) * numpy.nan)
	separations.append(row)

print('  building primary_id index ...')
primary_id_key = match.get_tablekeys(tables[0], 'ID')
primary_id_key = '%s_%s' % (table_names[0], primary_id_key)

primary_ids = []
primary_id_start = []
last_primary_id = None
primary_id_column = table[primary_id_key]
for i, pid in enumerate(primary_id_column):
	if pid != last_primary_id:
		last_primary_id = pid
		primary_ids.append(pid)
		primary_id_start.append(i)

primary_id_end = primary_id_start[1:] + [len(primary_id_column)]

# compute n-way position evidence
print('  computing probabilities ...')

log_bf = numpy.zeros(len(table)) * numpy.nan
prior = numpy.zeros(len(table)) * numpy.nan
# handle all cases (also those with missing counterparts in some catalogues)
for case in range(2**(len(table_names)-1)):
	table_mask = numpy.array([True] + [(case // 2**(ti)) % 2 == 0 for ti in range(len(tables)-1)])
	ncat = table_mask.sum()
	# select those cases
	mask = True
	for i in range(1, len(tables)):
		if table_mask[i]: # require not nan
			mask = numpy.logical_and(mask, ~numpy.isnan(separations[0][i]))
		else:
			mask = numpy.logical_and(mask, numpy.isnan(separations[0][i]))
	# select errors
	errors_selected = [e[mask] for e, m in zip(errors, table_mask) if m]
	separations_selected = [[cell[mask] for cell, m in zip(row, table_mask) if m] 
		for row, m in zip(separations, table_mask) if m]
	log_bf[mask] = bayesdist.log_bf(separations_selected, errors_selected)

	prior[mask] = source_densities[0] * args.prior_completeness / numpy.product(source_densities_plus[table_mask])
	assert numpy.isfinite(prior[mask]).all(), (source_densities, args.prior_completeness, numpy.product(source_densities_plus[table_mask]))

assert numpy.isfinite(prior).all(), (prior, log_bf)
assert numpy.isfinite(log_bf).all(), (prior, log_bf)
columns.append(pyfits.Column(name='dist_bayesfactor', format='E', array=log_bf))

ncat = table['ncat']
ncats = len(tables)

if args.consider_unrelated_associations:
	candidates = numpy.where(ncat <= ncats - 2)[0]
	if len(candidates) > 0:
		print('    correcting for unrelated associations ...')
		# correct for unrelated associations
		# identify those in need of correction
		# two unconsidered catalogues are needed for an unrelated association
		pbar = progress.bar(ndigits=6)
		for i in pbar(candidates):
			# list which ones we are missing
			missing_cats = [k for k, sep in enumerate(separations[0]) if numpy.isnan(sep[i])]
			pid = table[primary_id_key][i]
			pid_index = primary_ids.index(pid)
			best_logpost = 0
			# go through more complex associations
			for j in range(primary_id_start[pid_index], primary_id_end[pid_index]):
				if not (ncat[j] > 2): continue
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
				log_bf[i] += best_logpost
		columns.append(pyfits.Column(name='dist_bayesfactor_corrected', format='E', array=log_bf))
	else:
		print('      correcting for unrelated associations ... not necessary')

# add the additional columns
post = bayesdist.posterior(prior, log_bf)
columns.append(pyfits.Column(name='dist_post', format='E', array=post))

# find magnitude biasing functions
if magnitude_columns:
	print()
	print('Incorporating magnitude biases ...')
biases = {}
for mag, magfile in magnitude_columns:
	print('    magnitude bias "%s" ...' % mag)
	table_name, col_name = mag.split(':', 1)
	assert table_name in table_names, 'table name specified for magnitude ("%s") unknown. Known tables: %s' % (table_name, ', '.join(table_names))
	ti = table_names.index(table_name)
	col_names = tables[ti].dtype.names
	assert col_name in col_names, 'column name specified for magnitude ("%s") unknown. Known columns in table "%s": %s' % (mag, table_name, ', '.join(col_names))
	ci = col_names.index(col_name)
	
	res = results[table_name]
	res_defined = results[table_name] != -1
	
	# get magnitudes of all
	mag_all = tables[ti][col_name]
	# mark -99 as undefined
	mag_all[mag_all == -99] = numpy.nan
	
	# get magnitudes of selected
	mask_all = ~numpy.logical_or(numpy.isnan(mag_all), numpy.isinf(mag_all))

	col = "%s_%s" % (table_name, col_name)
	
	if magfile == 'auto':
		if mag_include_radius is not None:
			if mag_include_radius >= match_radius * 60 * 60:
				print('WARNING: magnitude radius is very large (>= matching radius). Consider using a smaller value.')
			selection = table['Separation_max'] < mag_include_radius
			selection_possible = table['Separation_max'] < mag_exclude_radius
		else:
			selection = post > 0.9
			selection_possible = post > 0.01
		
		# ignore cases where counterpart is missing
		assert res_defined.shape == selection.shape, (res_defined.shape, selection.shape)
		selection = numpy.logical_and(selection, res_defined)
		selection_possible = numpy.logical_and(selection_possible, res_defined)
		
		#print '   selection', selection.sum(), selection_possible.sum(), (-selection_possible).sum()
		
		#rows = results[table_name][selection].tolist()
		rows = list(set(results[table_name][selection]))
		
		assert len(rows) > 1, 'No magnitude values within radius for "%s".' % mag
		mag_sel = mag_all[rows]
		
		# remove vaguely possible options from alternative histogram
		rows_possible = list(set(results[table_name][selection_possible]))
		mask_others = mask_all.copy()
		mask_others[rows_possible] = False
		
		# all options in the total (field+target sources) histogram
		mask_sel = ~numpy.logical_or(numpy.isnan(mag_sel), numpy.isinf(mag_sel))

		#print '      non-nans: ', mask_sel.sum(), mask_others.sum()

		print('magnitude histogram of column "%s": %d secure matches, %d insecure matches and %d secure non-matches of %d total entries (%d valid)' % (col, mask_sel.sum(), len(rows_possible), mask_others.sum(), len(mag_all), mask_all.sum()))
		
		# make function fitting to ratio shape
		bins, hist_sel, hist_all = magnitudeweights.adaptive_histograms(mag_all[mask_others], mag_sel[mask_sel])
		print('magnitude histogram stored to "%s".' % (mag.replace(':', '_') + '_fit.txt'))
		with open(mag.replace(':', '_') + '_fit.txt', 'wb') as f:
			f.write(b'# lo hi selected others\n')
			numpy.savetxt(f,
				numpy.transpose([bins[:-1], bins[1:], hist_sel, hist_all]), 
				fmt = ["%10.5f"]*4)
		if mask_sel.sum() < 100:
			print('ERROR: too few secure matches to make a good histogram. If you are sure you want to use this poorly sampled histogram, replace "auto" with the filename.')
			sys.exit(1)
	else:
		print('magnitude histogramming: using histogram from "%s" for column "%s"' % (magfile, col))
		bins_lo, bins_hi, hist_sel, hist_all = numpy.loadtxt(magfile).transpose()
		bins = numpy.array(list(bins_lo) + [bins_hi[-1]])
	func = magnitudeweights.fitfunc_histogram(bins, hist_sel, hist_all)
	magnitudeweights.plot_fit(bins, hist_sel, hist_all, func, mag)
	weights = log10(func(table[col]))
	# undefined magnitudes do not contribute
	weights[numpy.isnan(weights)] = 0
	biases[col] = weights


# add the bias columns
for col, weights in biases.items():
	columns.append(pyfits.Column(name='bias_%s' % col, format='E', array=10**weights))

print()
print('Computing final probabilities ...')

# add the posterior column
total = log_bf + sum(biases.values())
post = bayesdist.posterior(prior, total)
columns.append(pyfits.Column(name='p_single', format='E', array=post))

# compute weights for group posteriors
# 4pi comes from Eq. 
log_post_weight = bayesdist.unnormalised_log_posterior(prior, total, ncat)

# flagging of solutions. Go through groups by primary id (IDs in first catalogue)
index = numpy.zeros_like(post)
prob_has_match = numpy.zeros_like(post)
prob_this_match = numpy.zeros_like(post)

match_header['COL_PRIM'] = primary_id_key
match_header['COLS_ERR'] = ' '.join(['%s_%s' % (ti, poscol) for ti, poscol in zip(table_names, pos_errors)])
print('    grouping by column "%s" and flagging ...' % (primary_id_key))

pbar = progress.bar(ndigits=6)
pid_index = primary_ids.index(pid)
best_log_bf = 0
# go through more complex associations
for primary_id, ilo, ihi in pbar(list(zip(primary_ids, primary_id_start, primary_id_end))):
	# group
	mask = slice(ilo, ihi)
	
	# compute no-match probability
	values = log_post_weight[mask]
	offset = values.max()
	bfsum = log10((10**(values - offset)).sum()) + offset
	if len(values) > 1:
		offset = values[1:].max()
		bfsum1 = log10((10**(values[1:] - offset)).sum()) + offset
	else:
		bfsum1 = 0
	
	# for p_any, find the one without counterparts
	p_none = float(values[0])
	p_any = 1 - 10**(p_none - bfsum)
	# this avoids overflows in the no-counterpart solution, 
	# which we want to set to 0
	values[0] = bfsum1
	p_i = 10**(values - bfsum1)
	p_i[0] = 0
	
	prob_has_match[mask] = p_any
	prob_this_match[mask] = p_i
	
	best_val = p_i.max()
	
	# flag best & second best
	# ignore very poor solutions
	index[mask] = numpy.where(best_val == p_i, 1, 
		numpy.where(p_i > diff_secondary * best_val, 2, 0))
	
	
columns.append(pyfits.Column(name='p_any', format='E', array=prob_has_match))
columns.append(pyfits.Column(name='p_i', format='E', array=prob_this_match))

#index[ncat == 1] == -1
# add the flagging column
columns.append(pyfits.Column(name='match_flag', format='I', array=index))

# cut away poor posteriors if requested
if min_prob > 0:
	mask = -(prob_this_match < min_prob)
	print('    cutting away %d (below p_i minimum)' % (len(mask) - mask.sum()))

	for c in columns:
		c.array = c.array[mask]

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



