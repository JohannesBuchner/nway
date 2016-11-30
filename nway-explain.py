#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

__doc__ = """Explain a multiway association of astrometric catalogue. Use --help for usage.

Example: nway-explain.py out.fits 179
"""

import matplotlib.pyplot as plt
import numpy
import astropy.io.fits as pyfits
import sys
import argparse
from matplotlib.patches import Circle
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

print('loading catalogue %s' % args.matchcatalogue)
f = pyfits.open(args.matchcatalogue)
header = f[0].header
data = f[1].data
primary_id_col = header['COL_PRIM']
print('    searching for %s == %s' % (primary_id_col, args.id))
if issubclass(data.dtype[primary_id_col].type, numpy.integer):
	mask = data[primary_id_col] == int(args.id)
elif issubclass(data.dtype[primary_id_col].type, numpy.float):
	mask = data[primary_id_col] == float(args.id)
else:
	mask = data[primary_id_col] == args.id
print('    %d rows found' % (mask.sum()))
if mask.sum() == 0:
	print('ERROR: ID not found. Was searching for %s == %s' % (primary_id_col, args.id))
	sys.exit(1)

# make a plot of the positions

plt.figure(figsize=(12,12))
plt.axis('equal')
cols_ra = header['COLS_RA'].split(' ')
cols_dec = header['COLS_DEC'].split(' ')
cols_err = header['COLS_ERR'].split(' ')
center_ra = data[cols_ra[0]][mask][0]
center_dec = data[cols_dec[0]][mask][0]


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
		err = data[tblname + '_' + err[1:]][mask]
	else:
		err = numpy.ones(mask.sum()) * float(err)
	pos = set(zip(data[col_ra][mask], data[col_dec][mask], err))
	ras = numpy.array([ra for ra, dec, err in pos if ra != -99])
	decs = numpy.array([dec for ra, dec, err in pos if ra != -99])
	errs = numpy.array([err for ra, dec, err in pos if ra != -99]) / 60. / 60.
	r,  = plt.plot(convx(ras), convy(decs), marker=marker, mec=color, mfc='None', ms=8, mew=2, ls=' ', label='%s %s' % (col_ra, col_dec))
	#plt.circles(ras, decs, errs, facecolor='None', edgecolor=r.get_color())
	patches = [Circle((convx(ra), convy(dec)), converr(err)) for ra, dec, err in zip(ras, decs, errs)]
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


outfilename = '%s_explain_%s_options.pdf' % (args.matchcatalogue, args.id)
with PdfPages(outfilename) as pp:
	maxy = max([len(o)-1 for o in all_options])
	maxx = len(all_options)-1
	for i in ii:
		plt.figure(figsize=(3+maxx, 3))
		plt.axis('off')
		graph_make(all_options)
		j = []
		for col_ra, col_dec, options in zip(cols_ra, cols_dec, all_options):
			radec = data[col_ra][mask][i], data[col_dec][mask][i]
			plt.text(len(j), 0.1, col_ra + '\n' + col_dec, 
				rotation=90, size=6, ha='center', va='bottom')
			j.append(-options.index(radec))
		graph_highlight(all_options, j)
		plt.text(-0.1, -1, 'p_i=%.2f' % data['p_i'][mask][i], ha='right', va='center')
		plt.ylim(-maxy-0.5, 0.5)
		plt.xlim(-0.5, maxx+0.5)
		plt.savefig(pp, format='pdf', bbox_inches='tight')
		plt.close()

# go through each association and highlight
for j, i in enumerate(numpy.argsort(data['p_i'][mask])[::-1][:3]):
	ras = []
	decs = []
	for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
		if data[col_ra][mask][i] == -99:
			continue
		ra = data[col_ra][mask][i]
		dec = data[col_dec][mask][i]
		ras.append(ra)
		decs.append(dec)
	
	plt.plot(convx(ras), convy(decs), '-', lw=(3-j), label='top %s by distance (p_i=%.2f)' % (j+1, data['p_i'][mask][i]), color='y')

for col in header['BIASING'].split(', '):
	for col_ra, col_dec, color in zip(cols_ra, cols_dec, markers):
		if col.split('_', 2)[0] != col_ra.split('_', 2)[0]:
			continue
		bias = data['bias_' + col]
		mask1 = numpy.logical_and(mask, data[col_ra] != -99)
		mask2 = numpy.logical_and(mask1, bias > 1)
		if mask2.any():
			ra = data[col_ra][mask2]
			dec = data[col_dec][mask2]
			plt.plot(convx(ra), convy(dec), 
				's ', mew=4, ms=20, mec='g', mfc='None', label='%s strong boost > 1' % col, alpha=0.3)
		mask2 = numpy.logical_and(mask, bias < -0.5)
		if mask2.any():
			ra = data[col_ra][mask2]
			dec = data[col_dec][mask2]
			plt.plot(convx(ra), convy(dec), 'd ', mew=4, ms=20, mec='r', mfc='None', label='%s strong reject < -1' % col, alpha=0.3)

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
plt.plot(convx(ras), convy(decs), '-', lw=1.7, label='Most probable association (match_flag=1)', color='orange')

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
	plt.plot(convx(ras), convy(decs), '-', lw=0.5, label='Secondary, similarly probable association (match_flag=2)', color='yellow')

plt.xlabel('$\Delta$RA [arcsec]')
plt.ylabel('$\Delta$DEC [arcsec]')
xlo, xhi = plt.xlim()
ylo, yhi = plt.ylim()
hi = max(-xlo, xhi, -ylo, yhi)
plt.xlim(-hi, hi)
plt.ylim(-hi, hi)
plt.legend(loc='best', numpoints=1, prop=dict(size=8))
outfilename = '%s_explain_%s.pdf' % (args.matchcatalogue, args.id)
print('plotting to %s' % outfilename)
plt.savefig(outfilename, bbox_inches='tight')
plt.close()





