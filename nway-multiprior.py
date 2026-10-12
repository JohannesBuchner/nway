#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Multiway association with a prior_callback function.

Demonstrates the prior_callback argument of nwaylib.nway_match: a
user-provided function is given the match table and computes a probability
for each candidate combination, which is ingested into p_i and p_any like
the per-catalogue magnitude priors.

This example creates two mock catalogues:
- a primary catalogue with a few sources with 10 arcsec positional errors,
  each with a mock "color" property (green, blue, orange, red, brown)
- a dense secondary catalogue (on average 2 sources within a 10 arcsec
  radius) with 1 arcsec positional errors and randomly assigned colors

The prior_callback gives a probability of 1 when the colors of the primary
and secondary entries are identical, and 0 otherwise.

Finally, a plot shows the primary (stars) and secondary (circles) entries
in their mock color, with lines between the high-probability matches
(p_i > 0.2).
"""

import astropy.io.fits as pyfits
import astropy.table
import matplotlib.pyplot as plt
import numpy
from matplotlib.lines import Line2D

import nwaylib
import nwaylib.logger as logger

# some helper functions ...


def write_mock_catalogue(filename, name, ra, dec, error, colors, area):
	"""
	Writes a mock catalogue as FITS file with the columns that nway
	expects (RA, DEC in degrees, pos_err in arcsec) plus a string-valued
	"color" column, and the SKYAREA (square degrees) header field.
	"""
	table = astropy.table.Table()
	table['RA'] = ra
	table['DEC'] = dec
	table['pos_err'] = error
	table['color'] = colors
	hdu = astropy.io.fits.table_to_hdu(table)
	hdu.name = name
	hdu.header['SKYAREA'] = area
	pyfits.HDUList([pyfits.PrimaryHDU(), hdu]).writeto(filename, overwrite=True)
	print('wrote "%s" (%d sources)' % (filename, len(table)))


def table_from_fits(fitsname):
	"""
	Builds a match table dict (see nwaylib.nway_match) from a mock
	catalogue FITS file.
	"""
	fits_table = pyfits.open(fitsname)[1]
	data = fits_table.data
	return dict(name=fits_table.name,
		ra=data['RA'], dec=data['DEC'], error=data['pos_err'],
		area=fits_table.header['SKYAREA'], # square degrees
		mags=[], maghists=[], magnames=[])


def read_colors(fitsname):
	"""
	Reads the color column of a mock catalogue as array of strings.
	"""
	colors = pyfits.open(fitsname)[1].data['color']
	return numpy.array([c.decode() if isinstance(c, bytes) else str(c) for c in colors])


# create the catalogues

# constants

# we assume each object has some property, here color
# this property has to match across catalogs
COLORS = ['green', 'blue', 'orange', 'red', 'brown']

RA0, DEC0 = 250.0, 10.0 # field centre in degrees
FIELD_SIDE = 120.0 # arcsec
# SKYAREA: area covered by the catalogues, in square degrees
field_area = FIELD_SIDE**2 / 3600.0**2
N_PRIMARY = 15 # a few primary sources
MATCH_RADIUS = 30.0 # arcsec, >> the largest positional error
PRIMARY_FILENAME = 'mock_primary_catalog.fits'
SECONDARY_FILENAME = 'mock_secondary_catalog.fits'

rng = numpy.random.RandomState(42)

# primary sources with 10 arcsec positional errors, placed away from the
# field edges (by more than the match radius), so that each has the full
# candidate coverage
half = FIELD_SIDE / 2.0
prim_ra_off = rng.uniform(-half + MATCH_RADIUS, half - MATCH_RADIUS, N_PRIMARY)
prim_dec_off = rng.uniform(-half + MATCH_RADIUS, half - MATCH_RADIUS, N_PRIMARY)
prim_colors = rng.choice(COLORS, N_PRIMARY)
write_mock_catalogue(PRIMARY_FILENAME, 'XRAY',
	RA0 + prim_ra_off / 3600.0, DEC0 + prim_dec_off / 3600.0,
	numpy.full(N_PRIMARY, 10.0), prim_colors, field_area)

# dense secondary catalogue (on average 2 sources within a 10 arcsec radius,
# i.e. a density of 2/(pi * 10^2) per square arcsec) with 1 arcsec
# positional errors and randomly assigned colors
n_secondary = rng.poisson(2.0 / (numpy.pi * 10.0**2) * FIELD_SIDE**2)
sec_ra_off = rng.uniform(-half, half, n_secondary)
sec_dec_off = rng.uniform(-half, half, n_secondary)
sec_colors = rng.choice(COLORS, n_secondary)
write_mock_catalogue(SECONDARY_FILENAME, 'OPT',
	RA0 + sec_ra_off / 3600.0, DEC0 + sec_dec_off / 3600.0,
	numpy.full(n_secondary, 1.0), sec_colors, field_area)

# load the mock catalogues and match them
prim_colors = read_colors(PRIMARY_FILENAME)
sec_colors = read_colors(SECONDARY_FILENAME)

# the catalogues are named after the FITS extensions (upper-cased
# per the FITS standard); the match table columns are named after them
primary_name = 'XRAY'
secondary_name = 'OPT'

match_tables = [table_from_fits(PRIMARY_FILENAME), table_from_fits(SECONDARY_FILENAME)]


# first: distance-only matching, without color information

result_nocallback = nwaylib.nway_match(match_tables,
	match_radius=MATCH_RADIUS, prior_completeness=0.9,
	logger=logger.NullOutputLogger())
sel_nocallback = result_nocallback['prob_this_match'] > 0.2
sec_nocb = result_nocallback[secondary_name].values[sel_nocallback]
has_nocb = sec_nocb >= 0
n_diff_nocallback = (prim_colors[result_nocallback.index.get_level_values(primary_name).values[sel_nocallback][has_nocb]]
	!= sec_colors[sec_nocb[has_nocb]]).sum()
print()
print('without prior_callback: %d matches with p_i > 0.2, %d of them connect different colors' % (
	sel_nocallback.sum(), n_diff_nocallback))


# now we enforce that the colors have to be consistent in the matches
# we define a function which looks at each match table row and returns a probability for the combination of colors:

def prior_callback(table):
	"""
	Prior probability for each candidate combination:
	1 when the colors of the primary and secondary entries are identical,
	0 otherwise. Combinations without a secondary entry (the
	no-counterpart solutions) are neutral.
	"""
	prob = numpy.ones(len(table))
	# only compute a probability for matched rows:
	has_candidate = table['OPT'].values >= 0
	# extract the colors here
	same = (prim_colors[table['XRAY'].values[has_candidate]]
		== sec_colors[table['OPT'].values[has_candidate]])
	prob[has_candidate] = same.astype(float)
	return prob


result = nwaylib.nway_match(match_tables,
	match_radius=MATCH_RADIUS, prior_completeness=0.9,
	prior_callback=prior_callback,
	logger=logger.NormalLogger())


# check that the color prior was ingested into p_i and p_any

sel = result['prob_this_match'] > 0.2
prim_sel = result.index.get_level_values(primary_name).values[sel]
sec_sel = result[secondary_name].values[sel]
has_sec = sec_sel >= 0
same = prim_colors[prim_sel[has_sec]] == sec_colors[sec_sel[has_sec]]
assert same.all(), 'high-probability matches with different colors found'

print()
print('with prior_callback: %d matches with p_i > 0.2, all of them connect identical colors' % sel.sum())
print('primaries with a highlighted match: %d of %d' % (len(numpy.unique(prim_sel)), N_PRIMARY))
print()
print(result.loc[sel, [secondary_name, 'Separation_%s_%s' % (primary_name, secondary_name),
	'p_single', 'prob_has_match', 'prob_this_match', 'match_flag']])


# plot the catalogues and the high-probability matches

plt.figure(figsize=(7, 7))
plt.subplot(1, 1, 1)
for j in range(len(prim_sel)):
	plt.plot([prim_ra_off[prim_sel[j]], sec_ra_off[sec_sel[j]]],
		[prim_dec_off[prim_sel[j]], sec_dec_off[sec_sel[j]]],
		'-', color=prim_colors[prim_sel[j]], lw=1.5, alpha=0.6, zorder=1)
plt.scatter(sec_ra_off, sec_dec_off, marker='o', c=sec_colors,
	s=45, alpha=0.75, edgecolors='k', linewidths=0.4, zorder=2)
plt.scatter(prim_ra_off, prim_dec_off, marker='*', c=prim_colors,
	s=400, edgecolors='k', linewidths=0.8, zorder=3)
plt.xlabel('RA offset [arcsec]')
plt.ylabel('DEC offset [arcsec]')
plt.gca().invert_xaxis()
plt.title('nway prior_callback example: lines connect matches with p_i > 0.2')
plt.legend(handles=[
		Line2D([], [], color='k', marker='*', ls='', markersize=14, label='primary (stars)'),
		Line2D([], [], color='k', marker='o', ls='', markersize=8,
			markerfacecolor='none', label='secondary (circles)'),
	] + [Line2D([], [], color=color, lw=2, label='color "%s"' % color) for color in COLORS],
	loc='upper left', fontsize=9)
plt.tight_layout()
plt.savefig('nway-multiprior.png', dpi=150)
print()
print('plot written to "nway-multiprior.png"')
