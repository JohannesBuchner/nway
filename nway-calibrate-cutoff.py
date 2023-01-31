#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

__doc__ = """Characterise the false association rate and efficiency of a match
with a offset (fake) match and a real match.

Example: nway-calibrate-cutoff.py example2.fits example2-shifted-match.fits
"""

import sys
import numpy
from numpy import log10, pi, exp, logical_and
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import argparse

import nwaylib.fastskymatch as match

class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Johannes Buchner (C) 2013-2017 <johannes.buchner.acad@gmx.com>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('realfile', help="match output using real catalogue")
parser.add_argument('fakefile', help='match output using fake catalogue')

# parsing arguments
args = parser.parse_args()

f = pyfits.open(args.realfile)
data = f[1].data
mask = data['ncat'] == 1
p_any0 = data['p_any'][mask]
mask = data['match_flag'] == 1
p_any = data['p_any'][mask]
p_i = data['p_i'][mask]

f = pyfits.open(args.fakefile)
data = f[1].data
mask = data['ncat'] == 1
p_any0_offset = data['p_any'][mask]
mask = data['match_flag'] == 1
p_any_offset = data['p_any'][mask]
p_i_offset = data['p_i'][mask]

cutoffs = numpy.linspace(0, 1, 101)

efficiency = [(p_any0 > cutoff).mean() for cutoff in cutoffs]
error_rate = [(p_any0_offset > cutoff).mean() for cutoff in cutoffs]

plt.figure(figsize=(5,5))
plt.plot(cutoffs, efficiency, '-', color='k', label='selection efficiency', lw=2)
plt.plot(cutoffs, error_rate, '--', color='r', label='false selection rate')

plt.ylabel('fraction(>p_any)')
plt.xlabel('p_any')
plt.legend(loc='lower left', prop=dict(size=10))
plt.savefig(args.realfile + '_p_any_cutoffquality.pdf', bbox_inches='tight')
plt.savefig(args.realfile + '_p_any_cutoffquality.png', bbox_inches='tight')
plt.close()
print('created plot "%s_p_any_cutoffquality.pdf"' % args.realfile)

plt.figure(figsize=(5,5))
plt.plot(p_any, p_i, '. ', color='r', label='real')
plt.plot(p_any_offset, p_i_offset, '. ', color='gray', label='offset')

plt.ylabel('p_i')
plt.xlabel('p_any')
plt.legend(loc='lower left', prop=dict(size=10))
plt.savefig(args.realfile + '_p_any_p_i.pdf', bbox_inches='tight')
plt.savefig(args.realfile + '_p_any_p_i.png', bbox_inches='tight')
plt.close()
print('created plot "%s_p_any_p_i.pdf"' % args.realfile)

for rate in [0.01,0.03,0.05,0.1]:
	print()
	# find where error-rate is < rate
	mask = numpy.array(error_rate) < rate
	if not mask.any():
		print('A false detection rate of <%d%% is not possible.' % (rate*100))
	else:
		i = numpy.min(numpy.where(mask)[0])
		#print(i)
		print('For a false detection rate of <%d%%' % (rate*100))
		print('--> use only counterparts with p_any>%.2f (%.2f%% of matches)' % (cutoffs[i], efficiency[i] * 100))
