"""
Computes adaptively binned histograms of two magnitude distributions.

Provides smooth biasing functions and plots them.
"""

import scipy, scipy.interpolate, scipy.signal, scipy.integrate
import numpy
from numpy import log, pi, exp, logical_and
import matplotlib.pyplot as plt

# compute magnitude distributions
# use adaptive binning for that

"""
ratio between the two histograms. Minor offsets to avoid zero and inf.
"""
def ratio(hist_sel, hist_all):
	return (hist_sel + 1e-2) / (hist_all + 1e-2)

"""
Plotting
"""
def plot_fit(bin_mag, hist_sel, hist_all, func, name):
	plt.figure()
	hist_n = ratio(hist_sel, hist_all)
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

"""
creates the biasing functions
"""
def fitfunc_histogram(bin_mag, hist_sel, hist_all):
	bin_n = ratio(hist_sel, hist_all)
	w = scipy.signal.flattop(4) #, 1)
	w /= w.sum()
	bin_n_smooth = scipy.signal.convolve(bin_n, w, mode='same')
	interpfunc = scipy.interpolate.interp1d(bin_mag[:-1], 
		bin_n_smooth, bounds_error=False, fill_value=bin_n.min(), kind='quadratic')
	norm, err = scipy.integrate.quad(interpfunc, bin_mag.min(), bin_mag.max(),
		epsrel=1e-2)
	return lambda mag: log(interpfunc(mag) / norm)

"""
creates the histograms for the two columns in an adaptive way (based on mag_sel)
with the same binning.
"""
def adaptive_histograms(mag_all, mag_sel):
	# make histogram
	func_sel = scipy.interpolate.interp1d(
		numpy.linspace(0, 1, len(mag_sel)), 
		sorted(mag_sel))
	# choose bin borders based on cumulative distribution, using 20 points
	x = func_sel(numpy.linspace(0, 1, 20))
	hist_sel, bins = numpy.histogram(mag_sel, bins=x,    normed=True)
	hist_all, bins = numpy.histogram(mag_all, bins=bins, normed=True)
	return bins, hist_sel, hist_all


