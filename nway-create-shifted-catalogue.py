#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

__doc__ = """Create a shifted catalogue for testing the false association rate.

Example: nway-create-shifted-catalogue.py --radius 20 --shift-ra 0 --shift-dec 60 COSMOS-XMM.fits shifted-COSMOS-XMM.fits
"""

import sys
import numpy
from numpy import log10, pi, exp, logical_and
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import argparse
import nwaylib.progress as progress

import nwaylib.fastskymatch as match

class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Johannes Buchner (C) 2013-2017 <johannes.buchner.acad@gmx.com>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--shift-dec', default=0, type=float,
	help='Shift to add in dec (arcsec)')

parser.add_argument('--shift-ra', default=0, type=float,
	help='Shift to add in ra (arcsec)')

parser.add_argument('--radius', type=float, required=True,
	help='Remove sources which are near original sources, within this radius (arcsec).')

parser.add_argument('inputfile', type=str, help="input catalogue fits file")
parser.add_argument('outputfile', help='output catalogue fits file')

# parsing arguments
args = parser.parse_args()

outfile = args.outputfile
filename = args.inputfile
radius = args.radius

print('opening', filename)
inputfitsfile = pyfits.open(filename)
header_hdu = inputfitsfile[0]
table = inputfitsfile[1].data

if args.shift_ra==0 and args.shift_dec==0:
	print('ERROR: You have to set either shift-ra or shift-dec to non-zero')
	sys.exit(1)

ra_key  = match.get_tablekeys(table, 'RA')
print('    using RA  column: %s' % ra_key)
dec_key = match.get_tablekeys(table, 'DEC')
print('    using DEC column: %s' % dec_key)

ra_orig = table[ra_key]
dec_orig = table[dec_key]

ra = ra_orig + args.shift_ra / 60. / 60
dec = dec_orig + args.shift_dec / 60. / 60

# for each of them, check that there is no collision
excluded = []

pbar = progress.bar(ndigits=6)
for i, (ra_i, dec_i) in pbar(list(enumerate(zip(ra, dec)))):
	d = match.dist((ra_i, dec_i), (ra_orig, dec_orig))
	excluded.append((d * 60 * 60 < radius).any())

excluded = numpy.array(excluded)
print('removed %d sources which collide with original positions' % (excluded.sum()))

table[ra_key] = ra
table[dec_key] = dec

newcolumns = []
for col in table.columns:
	newcolumns.append(pyfits.Column(name=col.name, format=col.format, array=col.array[~excluded]))

tbhdu = match.fits_from_columns(newcolumns)
print('writing "%s" (%d rows)' % (outfile, len(tbhdu.data)))
for k, v in inputfitsfile[1].header.items():
	if k not in tbhdu.header:
		tbhdu.header[k] = v

hdulist = pyfits.HDUList([header_hdu, tbhdu])
hdulist.writeto(outfile, **progress.kwargs_overwrite_true)



