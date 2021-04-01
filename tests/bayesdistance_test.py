"""
Probabilistic Cross-Identification of Astronomical Sources

Reference: Budavari & Szalay (2008), ApJ, 679:301-309
Authors: Johannes Buchner (C) 2013-2020
Authors: Tamas Budavari (C) 2012
"""
import numpy
import numpy.testing as test
from nwaylib.bayesdistance import log_bf, log_bf2, log_bf3, make_invcovmatrix, apply_vABv, convert_from_ellipse

def test_log_bf_consistent2():
	# test that log_bf function reduces to log_bf2 for 2 catalogs:
	sep = numpy.array([0., 0.1, 0.2, 0.3, 0.4, 0.5])
	for psi in sep:
		print(psi)
		print('  ', log_bf2(psi, 0.1, 0.2), )
		print('  ', log_bf([[None, psi]], [0.1, 0.2]), )
		test.assert_almost_equal(log_bf2(psi, 0.1, 0.2), log_bf([[None, psi]], [0.1, 0.2]))

def test_log_bf_consistent3():
	sep = numpy.array([0., 0.1, 0.2, 0.3, 0.4, 0.5])
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

def test_ellipse_conversion():
	rng = numpy.random.RandomState(1)
	sigma_ra = rng.uniform(1, 100, size=100)
	sigma_dec = rng.uniform(1, 100, size=100)
	angles = rng.uniform(0, 180, size=100)
	# circular case:
	sigma_x, sigma_y, rho = convert_from_ellipse(sigma_ra, sigma_ra, angles)
	test.assert_almost_equal(rho, 0.)
	test.assert_almost_equal(sigma_x, sigma_ra)
	test.assert_almost_equal(sigma_y, sigma_ra)
	del sigma_x, sigma_y, rho

	# aligned case:
	sigma_x, sigma_y, rho = convert_from_ellipse(sigma_ra, sigma_dec, 0)
	test.assert_almost_equal(rho, 0.)
	test.assert_almost_equal(sigma_y, sigma_ra)
	test.assert_almost_equal(sigma_x, sigma_dec)
	del sigma_x, sigma_y, rho

	# rotated by 90 degrees
	sigma_x, sigma_y, rho = convert_from_ellipse(sigma_ra, sigma_dec, numpy.pi / 2)
	test.assert_almost_equal(rho, 0.)
	test.assert_almost_equal(sigma_x, sigma_ra)
	test.assert_almost_equal(sigma_y, sigma_dec)
	del sigma_x, sigma_y, rho

	sigma_x, sigma_y, rho_corr = convert_from_ellipse(sigma_ra, sigma_dec, angles)
	assert numpy.all(numpy.abs(rho_corr) > 1e-6), rho_corr

def test_corrmatrix_distances_circular():
	# difference in positions
	N = 10
	shift = numpy.arange(N)
	dra  = shift
	ddec = 0 * shift
	dvec = (dra, ddec)

	# precision:
	sigma_ra  = 1 * numpy.ones(N)
	sigma_dec = 1 * numpy.ones(N)
	sigma_ra2  = sigma_ra / 10
	sigma_dec2 = sigma_dec / 10
	sigma_ra_tot = (sigma_ra**-2 + sigma_ra2**-2)**-0.5
	sigma_dec_tot = (sigma_dec**-2 + sigma_dec2**-2)**-0.5

	# compute distance in units of sigma^2 by hand:
	delta1 = dra**2 / sigma_ra**2 + ddec**2 / sigma_dec**2
	delta2 = dra**2 / sigma_ra2**2 + ddec**2 / sigma_dec2**2
	delta = delta1 + delta2
	assert delta[0] == 0, delta[0]
	assert numpy.isclose(delta[1], 101), delta[1]
	test.assert_almost_equal(delta, (dra / sigma_ra_tot)**2 + (ddec / sigma_dec_tot)**2)

	# compute distance in units of sigma using bayesdistance tools:
	symsigmaI1 = make_invcovmatrix(sigma_ra, sigma_dec)
	symsigmaI2 = make_invcovmatrix(sigma_ra2, sigma_dec2)
	delta_matrix = apply_vABv(dvec, symsigmaI1, symsigmaI2)
	# make sure it is the same as the computation by hand above
	test.assert_almost_equal(delta, delta_matrix)

	# make sure the order of catalogs does not matter:
	delta_matrix2 = apply_vABv(dvec, symsigmaI2, symsigmaI1)
	test.assert_almost_equal(delta_matrix2, delta_matrix)

def test_corrmatrix_distances_elliptical():
	# difference in positions
	N = 4
	dra  = numpy.array([0., 1., 1., 1.])
	ddec = numpy.array([0., 0., 1., -1.])
	dvec = (dra, ddec)

	sigma_ra  = 1 * numpy.ones(N)
	sigma_dec = 2 * numpy.ones(N)
	sigma_ra2  = sigma_ra * 10
	sigma_dec2 = sigma_dec * 10
	#sigma_ra_tot = (sigma_ra**-2 + sigma_ra2**-2)**-0.5
	#sigma_dec_tot = (sigma_dec**-2 + sigma_dec2**-2)**-0.5
	sigma_x_corr, sigma_y_corr, rho_corr = convert_from_ellipse(sigma_ra, sigma_dec, -45 * numpy.pi / 180)
	print('circular:', sigma_ra, sigma_dec)
	print('elliptical:', sigma_x_corr, sigma_y_corr, rho_corr)

	assert numpy.all(numpy.abs(rho_corr) > 0.001), rho_corr

	symsigmaI1 = make_invcovmatrix(sigma_ra, sigma_dec)
	symsigmaI2 = make_invcovmatrix(sigma_ra2, sigma_dec2)
	delta_symm = apply_vABv(dvec, symsigmaI1, symsigmaI2)

	# use correlated errors:
	asymsigmaI1 = make_invcovmatrix(sigma_x_corr, sigma_y_corr, rho_corr)
	asymsigmaI2 = make_invcovmatrix(sigma_ra2, sigma_dec2)
	delta_corr = apply_vABv(dvec, asymsigmaI1, asymsigmaI2)
	# this must be different
	print(delta_symm, delta_corr)
	test.assert_almost_equal(delta_corr[0], delta_symm[0])
	# step along covariance should be about sqrt(2) smaller!
	assert delta_corr[2] < delta_symm[2] / 1.4, (delta_corr[2], delta_symm[2])
	# step in the orthogonal should be about sqrt(2) larger
	assert delta_corr[3] > delta_symm[3] * 1.4, (delta_corr[3], delta_symm[3])
