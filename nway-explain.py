#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Explain a multiway association of astrometric catalogue. Use --help for usage.

Example: nway-explain.py out.fits 179
"""

from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy
import astropy.io.fits as pyfits
import sys
import argparse
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from matplotlib.backends.backend_pdf import PdfPages

class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Johannes Buchner (C) 2013-2016 <johannes.buchner.acad@gmx.com>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('matchcatalogue', type=str, 
	help="""nway output catalogue""")

parser.add_argument('id', type=str,
	help='ID to explain (from primary catalogue)')



# parsing arguments
args = parser.parse_args()

#print('loading catalogue %s' % args.matchcatalogue)
f = pyfits.open(args.matchcatalogue)
header = f[0].header
data = f[1].data
primary_id_col = header['COL_PRIM']
#print('    searching for %s == %s' % (primary_id_col, args.id))
if issubclass(data.dtype[primary_id_col].type, numpy.integer):
	mask = data[primary_id_col] == int(args.id)
elif issubclass(data.dtype[primary_id_col].type, float):
	mask = data[primary_id_col] == float(args.id)
else:
	mask = data[primary_id_col] == args.id
#print('    %d rows found' % (mask.sum()))
if mask.sum() == 0:
	print('ERROR: ID not found. Was searching for %s == %s' % (primary_id_col, args.id))
	sys.exit(1)
#print()
# make a plot of the positions

plt.figure(figsize=(12,12))
plt.axis('equal')
cols_ra = header['COLS_RA'].split(' ')
cols_dec = header['COLS_DEC'].split(' ')
cols_err = header['COLS_ERR'].split(' ')
tablenames = header['TABLES'].split(', ')
center_ra = data[cols_ra[0]][mask][0]
center_dec = data[cols_dec[0]][mask][0]
p_any = data['p_any'][mask][0]

print('NWAY results for Source %s:' % args.id)
print()
if p_any > 0.8:
	print('This source probably has a counterpart (p_any=%.2f)' % p_any)
elif p_any < 0.1:
	print('This source probably does not a counterpart (p_any=%.2f)' % p_any)
else:
	print('It is uncertain if this source has a counterpart (p_any=%.2f)' % p_any)
print()
print("Assuming it has a counterpart, we have the following possible associations:")
print()

j_option = 0
def print_option(name, i):
	global j_option
	j_option += 1
	matchflag = data['match_flag'][mask][i]
	matchflagstars = '**' if matchflag == 1 else ('*' if matchflag==2 else '')
	if p_any < 0.1:
		matchflagstars = ''
	if matchflag == 0:
		print('Association %d: probability p_i=%.2f ' % (j_option, data['p_i'][mask][i]))
	else:
		print('Association %d%s[match_flag==%d]: probability p_i=%.2f ' % (j_option, matchflagstars, matchflag, data['p_i'][mask][i]))
	print('     Involved catalogues:  %s ' % (name))
	for col in header['BIASING'].split(', '):
		if col.strip() == '': continue
		bias = data['bias_' + col][mask][i]
		if bias >= 2:
			print('     prior %-15s increased the probability (bias_%s=%.2f)' % (col, col, bias))
		elif bias <= 0.5:
			print('     prior %-15s decreased the probability (bias_%s=%.2f)' % (col, col, bias))
	print()

def convx(ra):
	return (ra - center_ra) * 3600
def convy(dec):
	return (dec - center_dec) * 3600
def converr(err):
	return err * 3600

markers = ['x', '+', '^', '<', '>', 'v', 'p'] * 10
colors = ['b', 'c', 'g', 'r', 'k', 'brown'] * 10
all_options = []
ii = numpy.argsort(data['p_i'][mask])[::-1]
for col_ra, col_dec, col_err, marker, color in zip(cols_ra, cols_dec, cols_err, markers, colors):
	tblname, err = col_err.split('_', 1)
	if err.startswith(':'):
		err_components = err[1:].split(':')
		if len(err_components) == 1:
			ra_err = data[tblname + '_' + err_components[0]][mask]
			dec_err = ra_err
			pa_err = numpy.zeros(ra_err.shape)
		elif len(err_components) == 2:
			ra_err = data[tblname + '_' + err_components[0]][mask]
			dec_err = data[tblname + '_' + err_components[1]][mask]
			pa_err = numpy.zeros(ra_err.shape)
		elif len(err_components) == 3:
			ra_err = data[tblname + '_' + err_components[0]][mask]
			dec_err = data[tblname + '_' + err_components[1]][mask]
			pa_err = data[tblname + '_' + err_components[2]][mask]
		else:
			assert len(err_components) == 3, err_components
	else:
		ra_err = numpy.ones(mask.sum()) * float(err)
		dec_err = numpy.ones(mask.sum()) * float(err)
		pa_err = numpy.zeros(ra_err.shape)
	pos = set(zip(data[col_ra][mask], data[col_dec][mask], ra_err, dec_err, pa_err))
	ras = numpy.array([ra for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99])
	decs = numpy.array([dec for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99])
	ra_errs = numpy.array([ra_err for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99]) / 60. / 60.
	dec_errs = numpy.array([dec_err for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99]) / 60. / 60.
	pa_errs = numpy.array([pa_err for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99])
	r,  = plt.plot(convx(ras), convy(decs), marker=marker, mec=color, mfc='None', ms=8, mew=2, ls=' ', label='%s %s' % (col_ra, col_dec))
	patches = [
		Ellipse((convx(ra), convy(dec)),
		2 * converr(ra_err),
		2 * converr(dec_err),
		angle=90 - pa_err)
			for ra, dec, ra_err, dec_err, pa_err in zip(ras, decs, ra_errs, dec_errs, pa_errs)]
	p = PatchCollection(patches)
	p.set_facecolor('None')
	p.set_edgecolor(color)
	plt.gca().add_collection(p)
	options = [(-99, -99)]
	for i in ii:
		ra, dec = data[col_ra][mask][i], data[col_dec][mask][i]
		if (ra, dec) not in options:
			options.append((ra, dec))
	all_options.append(options)

def graph_make(all_options, highlight=None):
	for j, options in enumerate(all_options):
		if j != 0:
			plt.plot(j, 0, marker='x', color='gray')
		for k in range(len(options)-1):
			plt.plot(j, -(k + 1), marker='o', color='k')

def graph_highlight(all_options, selected):
	x = numpy.arange(len(all_options))
	plt.plot(x, selected, '-', color='k')


out_options = []
outfilename = '%s_explain_%s_options.pdf' % (args.matchcatalogue, args.id)
pp = PdfPages(outfilename)
maxy = max([len(o)-1 for o in all_options])
maxx = len(all_options)-1
for i in ii:
	plt.figure(figsize=(3+maxx, 3))
	plt.axis('off')
	graph_make(all_options)
	j = []
	name = []
	for col_ra, col_dec, options, tablename in zip(cols_ra, cols_dec, all_options, tablenames):
		radec = data[col_ra][mask][i], data[col_dec][mask][i]
		plt.text(len(j), 0.1, col_ra + '\n' + col_dec, 
			rotation=90, size=6, ha='center', va='bottom')
		k = options.index(radec)
		j.append(-k)
		if k == 0:
			name.append("")
		elif k == 1:
			name.append("%s" % tablename)
	print_option('-'.join(name), i)
	graph_highlight(all_options, j)
	plt.text(-0.1, -1, 'p_i=%.2f' % data['p_i'][mask][i], ha='right', va='center')
	plt.text(maxx + 0.1, 0, '$\leftarrow$ absent', ha='left', va='center')
	plt.ylim(-maxy-0.5, 0.5)
	plt.xlim(-0.5, maxx+0.5)
	plt.savefig(pp, format='pdf', bbox_inches='tight')
	plt.close()
pp.close()
print('plotting to %s' % outfilename)

# go through each association and highlight
for j, i in enumerate(numpy.argsort(data['p_single'][mask])[::-1][:3]):
	ras = []
	decs = []
	for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
		if data[col_ra][mask][i] == -99:
			continue
		ra = data[col_ra][mask][i]
		dec = data[col_dec][mask][i]
		ras.append(ra)
		decs.append(dec)
	
	plt.plot(convx(ras), convy(decs), '-', lw=(3-j), label='top %s by distance (p_single=%.2f, %d cat.)' % (j+1, data['p_single'][mask][i], data['ncat'][mask][i]), color='y')

mask2 = numpy.logical_and(mask, data['match_flag'] == 1)
ras = []
decs = []
first = True
for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
	ra = float(data[col_ra][mask2])
	dec = float(data[col_dec][mask2])
	if ra == -99:
		continue
	if not first:
		plt.text(convx(ra), convy(dec), ' 1', va='top', ha='left', alpha=0.5, size=16, fontweight='bold')
	first = False
	if ra == -99:
		continue
	ras.append(ra)
	decs.append(dec)
plt.plot(convx(ras), convy(decs), '-', lw=1.7, label='p_i=%.2f (match_flag=1)' % (float(data['p_i'][mask2])), color='orange')

mask2 = numpy.logical_and(mask, data['match_flag'] == 2)
for i in numpy.where(mask2)[0]:
	ras = []
	decs = []
	first = True
	for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
		ra = data[col_ra][i]
		dec = data[col_dec][i]
		if ra == -99: 
			continue
		if not first:
			plt.text(convx(ra), convy(dec), 
				' 2', va='top', ha='left', alpha=0.5, 
				size=16, fontweight='bold')
		first = False
		ras.append(ra)
		decs.append(dec)
	plt.plot(convx(ras), convy(decs), '-', lw=0.5, label='p_i=%.2f (match_flag=2)' % (data['p_i'][i]), color='yellow')

plt.xlabel('$\Delta$RA [arcsec]')
plt.ylabel('$\Delta$DEC [arcsec]')
plt.title('Source %s, p_any=%.2f' % (args.id, p_any))
xlo, xhi = plt.xlim()
ylo, yhi = plt.ylim()
hi = max(-xlo, xhi, -ylo, yhi)
plt.ylim(-hi, hi) # DEC
plt.xlim(hi, -hi) # RA goes the other way
plt.legend(loc='best', numpoints=1, prop=dict(size=8))
outfilename = '%s_explain_%s.pdf' % (args.matchcatalogue, args.id)
print('plotting to %s' % outfilename)
print()
print("Disclaimer: These results assume that the input (sky densities, positional errors, and priors) are correct.")
print()
plt.savefig(outfilename, bbox_inches='tight')
plt.close()
