"""
Tests for the prior_callback argument of nwaylib.nway_match.

The callback is given the match table and returns one probability per
candidate combination, which is ingested into p_i and p_any like the
per-catalogue magnitude priors.
"""
import numpy

import nwaylib
import nwaylib.logger as logger


def make_match_tables(seed=42, n_primary=30, n_secondary=100):
	rng = numpy.random.RandomState(seed)
	# small field around ra=250 deg, dec=10 deg (flat-sky approximation safe)
	ra0, dec0 = 250.0, 10.0
	half_side = 30.0 / 60 / 60 # 30 arcsec in degrees
	tables = [
		dict(name='primary',
			ra=ra0 + rng.uniform(-half_side, half_side, n_primary),
			dec=dec0 + rng.uniform(-half_side, half_side, n_primary),
			error=0.5 * numpy.ones(n_primary),
			area=1.0, mags=[], maghists=[], magnames=[]),
		dict(name='secondary',
			ra=ra0 + rng.uniform(-half_side, half_side, n_secondary),
			dec=dec0 + rng.uniform(-half_side, half_side, n_secondary),
			error=0.2 * numpy.ones(n_secondary),
			area=1.0, mags=[], maghists=[], magnames=[]),
	]
	return tables, rng


def make_type_prior_callback(types_primary, types_secondary):
	"""
	prior probability: 1 when the mock property is identical, 0 otherwise.
	Combinations without a secondary entry are neutral.
	"""
	def prior_callback(table):
		prob = numpy.ones(len(table))
		has_candidate = table['secondary'].values >= 0
		same = (types_primary[table['primary'].values[has_candidate]]
			== types_secondary[table['secondary'].values[has_candidate]])
		prob[has_candidate] = same.astype(float)
		return prob
	return prior_callback


def test_prior_callback_none():
	# without a callback, matching works as before
	tables, rng = make_match_tables()
	result = nwaylib.nway_match(tables, match_radius=5, prior_completeness=0.9,
		logger=logger.NullOutputLogger())
	assert len(result) > 0
	assert numpy.isfinite(result['prob_this_match']).all()
	assert numpy.isfinite(result['prob_has_match']).all()
	assert 'bias_prior_callback' not in result.columns


def test_prior_callback_selects_identical_property():
	tables, rng = make_match_tables()
	n1, n2 = len(tables[0]['ra']), len(tables[1]['ra'])
	types1 = rng.randint(0, 2, n1)
	types2 = rng.randint(0, 2, n2)
	result = nwaylib.nway_match(tables, match_radius=5, prior_completeness=0.9,
		prior_callback=make_type_prior_callback(types1, types2),
		logger=logger.NullOutputLogger())
	assert len(result) > 0
	assert 'bias_prior_callback' in result.columns
	# after the grouping step, the primary ids are stored in the table index
	prim = result.index.get_level_values('primary').values
	sec = result['secondary'].values
	has_candidate = sec >= 0
	same = types1[prim[has_candidate]] == types2[sec[has_candidate]]
	p_i = result['prob_this_match'].values
	# combinations with a different property are excluded (p_i == 0)
	assert (p_i[has_candidate][~same] == 0).all()
	# combinations with an identical property are selected
	assert p_i[has_candidate][same].max() > 0.5
	assert numpy.isfinite(p_i).all()
	# primaries without an identical-property candidate have no match
	for primary_id, group in result.groupby(level='primary', sort=False):
		sec_g = group['secondary'].values
		has = sec_g >= 0
		same_g = types1[group.index.get_level_values('primary').values[has]] == types2[sec_g[has]]
		if not same_g.any():
			assert (group['prob_has_match'] == 0).all()


def test_prior_callback_all_zero():
	# a callback that excludes every combination: no matches at all, no NaNs
	tables, rng = make_match_tables()

	def prior_callback(table):
		return numpy.zeros(len(table))

	result = nwaylib.nway_match(tables, match_radius=5, prior_completeness=0.9,
		prior_callback=prior_callback, logger=logger.NullOutputLogger())
	assert len(result) > 0
	assert numpy.isfinite(result['prob_this_match']).all()
	assert numpy.isfinite(result['prob_has_match']).all()
	assert (result['prob_this_match'].values == 0).all()
	assert (result['prob_has_match'].values == 0).all()
	assert (result['match_flag'].values == 0).all()
