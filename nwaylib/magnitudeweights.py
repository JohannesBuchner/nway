"""
Computes adaptively binned histograms of two magnitude distributions.

Provides smooth biasing functions and plots them.
"""

import scipy, scipy.interpolate, scipy.signal, scipy.integrate
import numpy
from numpy import log10, pi, exp, logical_and
import matplotlib.pyplot as plt

# compute magnitude distributions
# use adaptive binning for that

"""
ratio between the two histograms. Minor offsets to avoid zero and inf.
"""
def ratio(hist_sel, hist_all):
	with numpy.errstate(divide='ignore', invalid='ignore'):
		return numpy.where(hist_all == 0, 100, hist_sel / hist_all)

"""
fraction of selected
"""
def fraction(bin_mag, hist_sel, hist_all):
	# need a re-normalisation so that unknown values can be set to prob=1
	m = hist_all > 0
	hist_sel_m = hist_sel[m]
	hist_all_m = hist_all[m]
	ratio = hist_sel_m / hist_all_m
	delta = (bin_mag[1:] - bin_mag[:-1])[m]
	avg = (ratio * hist_all_m * delta).sum() / (delta * hist_all_m).sum()
	assert avg != 0
	
	r = numpy.ones(len(hist_all))
	r[m] = hist_sel_m / hist_all_m / avg
	assert numpy.isfinite(r).all(), r[~numpy.isfinite(r)]
	assert (r > 0).all(), r[~(r > 0)]
	return r

"""
Plotting
"""
def plot_fit(bin_mag, hist_sel, hist_all, func, name):
	mags = numpy.linspace(bin_mag.min(), bin_mag.max(), 400)
	plt.figure()
	plt.subplot(2, 1, 1)
	hist_n = ratio(hist_sel, hist_all)
	plt.plot(bin_mag[:-1], hist_all, '-', 
		drawstyle='steps-post', label='all')
	plt.plot(bin_mag[:-1], hist_sel, '-', 
		drawstyle='steps-post', label='selected')
	plt.legend(loc='best')
	plt.ylabel('normalized weight')
	plt.xlabel(name)
	plt.xlim(mags.min(), mags.max())
	plt.subplot(2, 1, 2)
	plt.plot(bin_mag[:-1], hist_n, '-',
		drawstyle='steps-post', label='ratio')
	plt.plot(mags, func(mags), '-', label='fit')
	plt.legend(loc='best')
	plt.ylabel('normalized weight')
	plt.xlabel(name)
	plt.xlim(mags.min(), mags.max())
	plt.savefig(name.replace(':', '_') + '_fit.pdf', bbox_inches='tight')
	plt.close()

"""
creates the biasing functions
"""
def fitfunc_histogram(bin_mag, hist_sel, hist_all):
	bin_n = ratio(hist_sel, hist_all)
	# w = scipy.signal.gaussian(5, 1)
	# w /= w.sum()
	# bin_n_smooth = scipy.signal.convolve(bin_n, w, mode='same')
	# no smoothing
	bin_n_smooth = bin_n
	interpfunc = scipy.interpolate.interp1d(bin_mag,
		list(bin_n_smooth) + [bin_n_smooth[-1]],
		bounds_error=False, kind='zero')
	return interpfunc

"""
creates the histograms for the two columns in an adaptive way (based on mag_sel)
with the same binning.
"""
def adaptive_histograms(mag_all, mag_sel):
	func_sel = scipy.interpolate.interp1d(
		numpy.linspace(0, 1, len(mag_sel)), 
		sorted(mag_sel))
	# choose bin borders based on cumulative distribution, using 20 points
	x = func_sel(numpy.linspace(0, 1, 15))
	x = numpy.asarray(list(x) + [mag_all.max() + 1])
	#print x
	# linear histogram (for no adaptiveness):
	##x = numpy.linspace(mag_all.min(), mag_all.max(), 20)
	hist_sel, bins = numpy.histogram(mag_sel, bins=x,    density=True)
	hist_all, bins = numpy.histogram(mag_all, bins=bins, density=True)
	return bins, hist_sel, hist_all


