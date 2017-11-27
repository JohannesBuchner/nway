"""
Probabilistic Cross-Identification of Astronomical Sources

Reference: Budavari & Szalay (2008), ApJ, 679:301-309
Authors: Johannes Buchner (C) 2013-2016
Authors: Tamas Budavari (C) 2012
"""
from nwaylib.bayesdistance import *

def test_log_bf():
	import numpy.testing as test
	sep = numpy.array([0., 0.1, 0.2, 0.3, 0.4, 0.5])
	for psi in sep:
		print(psi)
		print('  ', log_bf2(psi, 0.1, 0.2), )
		print('  ', log_bf([[None, psi]], [0.1, 0.2]), )
		test.assert_almost_equal(log_bf2(psi, 0.1, 0.2), log_bf([[None, psi]], [0.1, 0.2]))
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



