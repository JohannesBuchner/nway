"""
Probabilistic Cross-Identification of Astronomical Sources

Reference: Budavari & Szalay (2008), ApJ, 679:301-309
Authors: Johannes Buchner (C) 2013-2016
Authors: Tamas Budavari (C) 2012
"""
from __future__ import print_function, division
import numpy
from numpy import log, log10, pi, exp, logical_and, where, e

# use base 10 everywhere

log_arcsec2rad = log(3600 * 180 / pi)

def log_posterior(prior, log_bf):
	"""
	Returns log10 posterior probability normalised against alternative hypothesis
	that the sources are unrelated (Budavari+08)
	"""
	return -log10( 1 + (1 - prior) * 10**(-log_bf - log10(prior)))

def posterior(prior, log_bf):
	"""
	Returns posterior probability normalised against alternative hypothesis
	that the sources are unrelated (Budavari+08)
	"""
	with numpy.errstate(over='ignore'):
		return 1. / (1 + (1 - prior) * 10**(-log_bf - log10(prior)))

def unnormalised_log_posterior(prior, log_bf, ncat):
	"""
	Returns posterior probability (without normalisation against alternatives)
	"""
	return log_bf + log10(prior)


def log_bf2(psi, s1, s2):
	"""
	log10 of the 2-way Bayes factor, see eq.(16)
	psi separation 
	s1 and s2=accuracy of coordinates
	"""
	s = s1*s1 + s2*s2;
	return (log(2) + 2 * log_arcsec2rad - log(s) - psi*psi / 2 / s) * log10(e)

def log_bf3(p12,p23,p31, s1,s2,s3):
	"""
	log10 of the 3-way Bayes factor, see eq.(17)
	"""
	ss1 = s1*s1
	ss2 = s2*s2
	ss3 = s3*s3
	s = ss1*ss2 + ss2*ss3 + ss3*ss1
	q = ss3 * p12**2 + ss1 * p23**2 + ss2 * p31**2
	return (log(4) + 4 * log_arcsec2rad - log(s) - q / 2 / s) * log10(e)

def log_bf(p, s):
	"""
	log10 of the multi-way Bayes factor, see eq.(18)

	p: separations matrix (NxN matrix of arrays)
	s: errors (list of N arrays)
	"""
	
	n = len(s)
	# precision parameter w = 1/sigma^2
	w = [numpy.asarray(si, dtype=numpy.float)**-2. for si in s]
	norm = (n - 1) * log(2) + 2 * (n - 1) * log_arcsec2rad
	
	wsum = numpy.sum(w, axis=0)
	s = numpy.sum(log(w), axis=0) - log(wsum)
	q = 0
	for i, wi in enumerate(w):
		for j, wj in enumerate(w):
			if i < j:
				q += wi * wj * p[i][j]**2
	exponent = - q / 2 / wsum
	return (norm + s + exponent) * log10(e)

# vectorized in the following means that many 2D matrices/vectors are going to be handled.
# i.e., each entry in the matrix or vector, is a vector of numbers.

def matrix_multiply(A, B):
	""" 2D matrix multiplication, vectorized """
	(a11, a12), (a21, a22) = A
	(b11, b12), (b21, b22) = B
	return (a11*b11+a12*b21, a11*b12+a12*b22), (a21*b11+a22*b21, a21*b12+a22*b22)

def matrix_det(A):
	""" 2D matrix determinant, vectorized """
	(a, b), (c, d) = A
	return a * d - b * c

def apply_vector_right(A, b):
	""" multiply 2D vector with 2D matrix, vectorized """
	(a11, a12), (a21, a22) = A
	(b1, b2) = b
	return a11*b1+a12*b2, a21*b1+a22*b2
	
def apply_vector_left(a, B):
	""" multiply 2D matrix with 2D vector, vectorized """
	a1, a2 = a
	(b11, b12), (b21, b22) = B
	return a1*b11+a2*b21, a1*b12+a2*b22

def vector_multiply(a, b):
	""" multiply two vectors, vectorized """
	a1, a2 = a
	b1, b2 = b
	return a1*b1+a2*b2

def apply_vABv(v, A, B):
	""" compute v^T A B v, vectorized """
	return vector_multiply(v, apply_vector_right(matrix_multiply(A, B), v))

def make_covmatrix(sigma_x, sigma_y, rho=0):
	""" create covariance matrix from given standard deviations and normalised correlation rho, vectorized """
	return (sigma_x**2, rho * sigma_x * sigma_y), (rho * sigma_x * sigma_y, sigma_y**2)

def make_invcovmatrix(sigma_x, sigma_y, rho=0):
	""" create inverse covariance matrix from given standard deviations and normalised correlation rho, vectorized """
	F = 1.0 / (sigma_x**2 * sigma_y**2 * (1 - rho**2))
	return (F * sigma_y**2, F * -rho * sigma_x * sigma_y), \
		(F * -rho * sigma_x * sigma_y, F * sigma_x**2)

def convert_from_ellipse(a, b, phi):
	""" create covariance parameters from ellipse major axis a, minor axis b and angle phi; vectorized """
	a2 = a**2
	b2 = b**2
	s = numpy.sin(phi)
	c = numpy.cos(phi)
	s2 = s**2
	c2 = c**2
	sigma_x = (a2 * s2 + b2 * c2)**0.5
	sigma_y = (a2 * c2 + b2 * s2)**0.5
	rho = c * s * (a2 - b2) / (sigma_x * sigma_y)
	return sigma_x, sigma_y, rho

def log_bf_elliptical(separations_ra, separations_dec, pos_errors):
	"""
	log10 of the multi-way Bayes factor, see eq.(18)

	separations_ra: RA separations matrix (NxN matrix of arrays)
	separations_dec: DEC separations matrix (NxN matrix of arrays)
	pos_errors: errors (list of (N,3) arrays) giving sigma_RA, sigma_DEC, rho
	"""
	
	# p ~ separations, s ~ pos_errors
	n = len(pos_errors)
	
	error_matrices = [make_invcovmatrix(si, sj, rho) 
		for si, sj, rho in pos_errors]
	
	# precision parameter w = 1/sigma^2
	w = [matrix_det(mi)**0.5 for mi in error_matrices]
	norm = (n - 1) * log(2) + 2 * (n - 1) * log_arcsec2rad
	
	wsum = numpy.sum(w, axis=0)
	s = numpy.sum(log(w), axis=0) - log(wsum)
	q = 0
	for i, Mi in enumerate(error_matrices):
		for j, Mj in enumerate(error_matrices):
			if i < j:
				v = (separations_ra[i][j], separations_dec[i][j])
				q += apply_vABv(v, Mi, Mj)
	exponent = - q / 2 / wsum
	return (norm + s + exponent) * log10(e)


def test_log_bf():
	import numpy.testing as test
	sep = numpy.array([0., 0.1, 0.2, 0.3, 0.4, 0.5])
	for psi in sep:
		print(psi)
		print('  ', log_bf2(psi, 0.1, 0.2), )
		poserrs = [0.1, 0.2]
		psi2 = [[None, psi]]
		print('  ', log_bf(psi2, poserrs))
		test.assert_almost_equal(log_bf2(psi, 0.1, 0.2), log_bf(psi2, poserrs))
		test.assert_almost_equal(log_bf(psi2, poserrs),
			log_bf_elliptical([[None, 0]], psi2, [convert_from_ellipse(0.1,0.1,0), convert_from_ellipse(0.2,0.2,0)]))
		test.assert_almost_equal(log_bf(psi2, poserrs),
			log_bf_elliptical(psi2, [[None, 0]], [convert_from_ellipse(0.1,0.1,0), convert_from_ellipse(0.2,0.2,0)]))
		
	for psi in sep:
		print(psi)
		bf3 = log_bf3(psi, psi, psi, 0.1, 0.2, 0.3)
		print('  ', bf3)
		g = log_bf([[None, psi, psi], [psi, None, psi], [psi, psi, None]], [0.1, 0.2, 0.3])
		print('  ', g)
		test.assert_almost_equal(bf3, g)
	q = numpy.zeros(len(sep))
	print(log_bf(numpy.array([[numpy.nan + sep, sep, sep], [sep, numpy.nan + sep, sep], [sep, sep, numpy.nan + sep]]), 
		[0.1 + q, 0.2 + q, 0.3 + q]))

	print("Elliptical:")
	print("-----------")
	print(log_bf_elliptical([[None, 0]], [[None, 0]], [convert_from_ellipse(0.1,0.1,0), convert_from_ellipse(0.2,0.2,0)]))
	print()
	print("Circular:")
	print("-----------")
	print(log_bf([[None, 0]], [0.1,0.2]))


