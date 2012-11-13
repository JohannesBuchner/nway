#!/usr/bin/env python
__doc__ = """
 Probabilistic Cross-Identification of Astronomical Sources
 [2012-11-08] Tamas Budavari <budavari@jhu.edu>

Reference: Budavari & Szalay (2008), ApJ, 679:301-309
Authors: Tamas Budavari (C) 2012
Authors: Johannes Buchner (C) 2012

The functions implement the static Bayes factor
formulas for high-accuracy detections, see Figure 1 and
Table 1 of Budavari & Szalay (2008)

The program performs 3-way matching
Usage: python 3way.py <association.fits> <hist_ratio_1> <hist_ratio_1>

association.fits is a fits file containing all possible cross matches (cross product)

The hist_ratio file is a (csv) table containing magnitude bins and the 
a-priori ratio in the last column.

The output is association.fits_out.csv, a csv file containing the cross matches.
"""

import numpy
from numpy import log, pi, exp

arcsec2rad = 3600 * 180 / pi
log_arcsec2rad = log(arcsec2rad)
log_arcsec2rad2 = 2 * log_arcsec2rad
log_arcsec2rad4 = 4 * log_arcsec2rad
log_2_arcsec2rad2 = log(2) + log_arcsec2rad2
log_4_arcsec2rad4 = log(4) + log_arcsec2rad4

nx = 1887
sum_prob = 0

def posterior(prior, bf):
	return 1 / ( 1 + (1-prior) / prior / bf)

# for numerical stability, the same in log:
def log_posterior(prior, log_bf):
	return -log( 1 + (1 - prior) * exp(-log_bf - log(prior)))

# Natural log of the 2-way Bayes factor, see eq.(16)
# psi=separation s1 and s2=accuracy of coordinates
def log_bf2(psi, s1, s2):
	s2 = s1*s1 + s2*s2;
	return log_2_arcsec2rad2 - log(s2) - psi*psi / 2 / s2

# Natural log of the 3-way Bayes factor, see eq.(17)
def log_bf3(p12,p23,p31, s1,s2,s3):
	ss1 = s1*s1
	ss2 = s2*s2
	ss3 = s3*s3
	s4 = ss1*ss2 + ss2*ss3 + ss3*ss1
	q = ss3 * p12*p12 + ss1 * p23*p23 + ss2 * p31*p31
	return log_4_arcsec2rad4 - log(s4) - q / 2 / s4


# read input files

import sys

import scipy, scipy.optimize
# fit functions to ratio function shape

def plot_fit(bin_mag, bin_n, func, name):
	try:
		import matplotlib.pyplot as plt
	except:
		print 'ERROR: plotting not possible, install matplotlib'
		return
	plt.figure()
	plt.plot(bin_mag, bin_n, '--', drawstyle='steps-pre', label='histogram')
	mags = numpy.linspace(bin_mag.min(), bin_mag.max(), 400)
	plt.plot(mags, exp(func(mags)), '-', drawstyle='steps-pre', label='fit')
	plt.legend(loc='best')
	plt.savefig(name + '_fit.pdf')
	plt.close()


def fitfunc(bin_mag, bin_n):
	a = 8 / 40.
	b = 17
	c = 2.2
	d = 0.8
	a,b,c,d = 0.157358, 17.000000, 1.668907, 0.767761
	
	def minfunc((a,b,c,d)): # function to minimize
		# simple chi^2
		predicted = a*(bin_mag-b)**c*exp(-d*(bin_mag-b)) 
		#plot_fit(bin_mag, bin_n, lambda bin_mag: a*(bin_mag-b)**c*exp(-d*(bin_mag-b)), 'opt')
		#print predicted, bin_n
		#print 'press return'
		#sys.stdin.readline()
		chi2 = ((predicted - bin_n)**2).sum()
		if numpy.isnan(chi2):
			return 1e300
		#print a,b,c,d,chi2
		
		return chi2
	
	a,b,c,d = scipy.optimize.fmin(minfunc, x0=[a,b,c,d], disp=0)
	
	print 'params: %f, %f, %f, %f' % (a,b,c,d)
	# return logarithm
	
	return lambda mag: log(a) + c * log(mag - b) + -d * (mag - b)

if len(sys.argv) != 3 + 1:
	print __doc__
	sys.exit(1)

print 'fitting magnitude function 1'
ratio1 = numpy.loadtxt(sys.argv[2])
func1 = fitfunc(ratio1[:,0], ratio1[:,-1]) # first and last column
plot_fit(ratio1[:,0], ratio1[:,-1], func1, sys.argv[2])

print 'fitting magnitude function 2'
ratio2 = numpy.loadtxt(sys.argv[3])
func2 = fitfunc(ratio2[:,0], ratio2[:,-1]) # first and last column
plot_fit(ratio2[:,0], ratio2[:,-1], func2, sys.argv[3])

# so far so good, now lets calculate probabilities

print 'loading catalogue'
import pyfits
xcat = pyfits.open(sys.argv[1])[1].data
outfilename = sys.argv[1] + '_prob.csv'

xpos = xcat['separation_X_irac']
pos1 = xcat['separation_irac_opt']
pos2 = xcat['separation_x_opt']
xacc = xcat['ERROR_RADIUS']

ln_bf = log_bf3(xpos, pos1, pos2, xacc, 0.5, 0.1)

assert not numpy.any(numpy.isnan(ln_bf))

mag1 = xcat['I'] # I-band
mag2 = xcat['I'] # TODO: which column 

ln_bf_i = ln_bf + numpy.where(mag1 > -99, func1(mag1), 0) + numpy.where(mag2 > -99, func2(mag2), 0)

assert not numpy.any(numpy.isnan(ln_bf_i))
assert not numpy.any(numpy.isinf(ln_bf_i))

# TODO: do the computation of 15e+18 explicitly here for documentation
log_prob = log_posterior(nx/(1887*15e+18), ln_bf)
log_post = log_posterior(1/15e+18, ln_bf_i)

print 'writing output'
fout = file(outfilename, 'w')
# TODO: give nicer names
fout.write('NAME,BF,prob,bf_i,posterior\n')

for row in zip(xcat['NAME'], ln_bf/log(10), log_prob/log(10), ln_bf_i/log(10), log_post/log(10)):
	fout.write('%s,%g,%g,%g,%g\n' % row)

print 'all done. see %s' % outfilename


