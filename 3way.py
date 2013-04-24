#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

__doc__ = """Probabilistic Cross-Identification of Astronomical Sources

The program performs multiway matching between catalogue. Use --help for usage."""


import sys
import scipy, scipy.optimize, scipy.interpolate, scipy.signal, scipy.integrate
import numpy
from numpy import log, pi, exp, logical_and
import matplotlib.pyplot as plt
import pyfits
import argparse
import match
import bayesdist

parser = argparse.ArgumentParser(description=__doc__,
	epilog="Johannes Buchner (C) 2013 <jbuchner@mpe.mpg.de>",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--radius', default=10, type=float,
	help='exclusive search radius in arcsec for initial matching')

parser.add_argument('--mag-radius', default=5, type=float,
	help='search radius for magnitude histograms')

#nx/(1887*15e+18)
parser.add_argument('--prior', default=0.95/6.4e+20, type=float,
	help='prior: density of sources expected')

parser.add_argument('--mag', type=str, action='append',
	help='name of <table>_<column> for magnitude biasing', default=[])

parser.add_argument('--acceptable-prob', type=float, default=0.005,
	help='limit up to which secondary solutions are flagged')

parser.add_argument('--min-prob', type=float, default=0,
	help='lowest probability allowed')

parser.add_argument('--out', help='output file name', required=True)

parser.add_argument('catalogues', type=str, nargs='+',
	help='input catalogue fits files')

args = parser.parse_args()

diff_secondary = args.acceptable_prob
outfile = args.out

filenames = args.catalogues[::2] # args.catalogues # sys.argv[1:]
pos_errors = args.catalogues[1::2]

tables = [pyfits.open(fitsname)[1] for fitsname in filenames]
table_names = [t.name for t in tables]
tables = [t.data for t in tables]
min_prob = args.min_prob

match_radius = args.radius / 60. / 60 # in degrees
mag_radius = args.mag_radius / 60. / 60 # in degrees

# match input catalogues

results, columns = match.match_multiple(tables, table_names, match_radius)
table = pyfits.new_table(pyfits.ColDefs(columns)).data


# compute magnitude distributions
# use adaptive binning for that


def plot_fit(bin_mag, hist_sel, hist_all, func, name):
	plt.figure()
	hist_n = (hist_sel + 1e-3) / (hist_all + 1e-3)
	plt.plot(bin_mag[:-1], hist_all / hist_all.sum(), '--', 
		drawstyle='steps-pre', label='all')
	plt.plot(bin_mag[:-1], hist_sel / hist_sel.sum(), '--', 
		drawstyle='steps-pre', label='selected')
	plt.plot(bin_mag[:-1], hist_n / hist_n.sum(), '--', 
		drawstyle='steps-pre', label='ratio histogram')
	mags = numpy.linspace(bin_mag.min(), bin_mag.max(), 400)
	plt.plot(mags, exp(func(mags)), '-', drawstyle='steps-pre', label='fit')
	plt.legend(loc='best')
	plt.ylabel('normalized weight')
	plt.xlabel('magnitude')
	plt.savefig(name.replace(':', '_') + '_fit.pdf', bbox_inches='tight')
	plt.close()

def fitfunc_histogram(bin_mag, hist_sel, hist_all):
	bin_n = (hist_sel + 1e-3) / (hist_all + 1e-3)
	w = scipy.signal.flattop(4) #, 1)
	w /= w.sum()
	bin_n_smooth = scipy.signal.convolve(bin_n, w, mode='same')
	interpfunc = scipy.interpolate.interp1d(bin_mag[:-1], 
		bin_n_smooth, bounds_error=False, fill_value=bin_n.min(), kind='quadratic')
	norm, err = scipy.integrate.quadrature(interpfunc, bin_mag.min(), bin_mag.max())
	return lambda mag: log(interpfunc(mag) / norm)

def fitfunc(mag_all, mag_sel):
	# make histogram
	func_sel = scipy.interpolate.interp1d(
		numpy.linspace(0, 1, len(mag_sel)), 
		sorted(mag_sel))
	#x = numpy.linspace(numpy.nanmin(mag_all), numpy.nanmax(mag_all), 10)
	x = func_sel(numpy.linspace(0, 1, 20))
	hist_sel, bins = numpy.histogram(mag_sel, bins=x,    normed=True)
	hist_all, bins = numpy.histogram(mag_all, bins=bins, normed=True)
	return bins, hist_sel, hist_all


biases = {}

for mag in args.mag:
	table_name, col_name = mag.split(':', 1)
	assert table_name in table_names, 'table name specified for magnitude (%s) unknown. Known tables: %s' % (table_name, ', '.join(table_names))
	ti = table_names.index(table_name)
	col_names = tables[ti].dtype.names
	assert col_name in col_names, 'column name specified for magnitude (%s) unknown. Known columns in table %s: %s' % (mag, table_name, ', '.join(col_names))
	ci = col_names.index(col_name)
	
	# get magnitudes of all
	mag_all = tables[ti][col_name]
	
	# get magnitudes of selected
	mask_all = -numpy.logical_or(numpy.isnan(mag_all), numpy.isinf(mag_all))
	
	rows = list(set(results[table_name][table['Separation_max'] < mag_radius]))
	mag_sel = mag_all[rows]
	mask_radius = table['Separation_max'] < mag_radius
	mask_sel = -numpy.logical_or(numpy.isnan(mag_sel), numpy.isinf(mag_sel))
	col = "%s_%s" % (table_name, col_name)
	print 'selected %d matches for magnitude histogramming: %d magnitude values for %s' % (mask_radius.sum(), len(mag_sel), col)
	
	# make function fitting to ratio shape
	bins, hist_sel, hist_all = fitfunc(mag_all[mask_all], mag_sel[mask_sel])
	func = fitfunc_histogram(bins, hist_sel, hist_all)
	plot_fit(bins, hist_sel, hist_all, func, mag)
	weights = func(table[col])
	biases[col] = weights

# add the additional columns

# compute n-way position evidence of merged catalogue

errors    = []
for table_name, pos_error in zip(table_names, pos_errors):
	if pos_error[0] == ':':
		# get column
		k = "%s_%s" % (table_name, pos_error[1:])
		assert k in table.dtype.names, 'ERROR: Position error column for %s not in table %s. Have columns: %s' % (k, table_name, ', '.join(table.dtype.names))
		print 'using column', (table[k].min(), table[k].max())
		errors.append(table[k])
	else:
		errors.append(float(pos_error[1:]) * numpy.ones(len(table)))

separations = []
for ti, a in enumerate(table_names):
	row = []
	for tj, b in enumerate(table_names):
		if ti < tj:
			k = 'Separation_%s_%s' % (b, a)
			assert k in table.dtype.names, 'ERROR: Separation column for %s not in merged table. Have columns: %s' % (k, ', '.join(table.dtype.names))
			row.append(table[k])
		else:
			row.append(numpy.ones(len(table)) * numpy.nan)
	separations.append(row)

ln_bf = bayesdist.log_bf(separations, errors)

columns.append(pyfits.Column(name='bf', format='E', array=ln_bf/log(10)))
columns.append(pyfits.Column(name='bfpost', format='E', array=bayesdist.posterior(prior, ln_bf)))

for col, weights in biases.iteritems():
	columns.append(pyfits.Column(name='bias_%s' % col, format='E', array=weights))

total = ln_bf + sum(biases.values())
columns.append(pyfits.Column(name='post', format='E', array=bayesdist.posterior(prior, total)))

post = bayesdist.posterior(prior, total)
index = numpy.zeros_like(post)

primary_id_key = match.get_tablekeys(tables[0], 'ID')
print 'grouping by %s from %s' % (primary_id_key, table_names[0])
primary_id_key = '%s_%s' % (table_names[0], primary_id_key)

for primary_id in set(table[primary_id_key]):
	# group
	mask = table[primary_id_key] == primary_id
	group_posterior = post[mask]
	best_index = group_posterior.argmax()
	best_val = group_posterior[best_index]
	
	# flag second best ( max('post') - post < diff_post )
	mask1 = mask.copy()
	mask1[best_val != group_posterior] = False
	mask2 = mask.copy()
	mask2[best_val - group_posterior > diff_secondary] = False
	index[mask2] = 2
	# flag best
	index[mask1] = 1
	maskbad = mask.copy()
	maskbad[group_posterior < min_prob] = False
	index[maskbad] = -1

columns.append(pyfits.Column(name='match_flag', format='I', array=index))

tbhdu = pyfits.new_table(pyfits.ColDefs(columns))
hdulist = match.wraptable2fits(tbhdu, 'MULTIMATCH')
hdulist[0].header.update('METHOD', 'multi-way matching')
hdulist[0].header.update('INPUT', ', '.join(filenames))
hdulist[0].header.update('TABLES', ', '.join(table_names))
hdulist[0].header.update('BIASING', ', '.join(biases.keys()))
for k, v in args.__dict__.iteritems():
	hdulist[0].header.add_comment("argument %s: %s" % (k, v))
hdulist.writeto(outfile, clobber=True)

print 'wrote %s (%d rows, %d columns)' % (outfile, len(table), len(columns))




