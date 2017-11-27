"""
Functions for finding pairs within a certain search distance "err".

Very fast method based on hashing.
"""
from __future__ import print_function, division
from nwaylib.fastskymatch import *

def test_dist():
	
	d = dist((53.15964508, -27.92927742), (53.15953445, -27.9313736))
	print('distance', d)
	assert not numpy.isnan(d)
	
	ra = numpy.array([ 53.14784241,  53.14784241,  53.14749908, 53.16559982,     53.19423676,  53.1336441 ])
	dec = numpy.array([-27.79363823, -27.79363823, -27.81790352, -27.79622459,      -27.70860672, -27.76327515])
	ra2 = numpy.array([ 53.14907837,  53.14907837,  53.1498642 , 53.16150284,     53.19681549,  53.13626862])
	dec2 = numpy.array([-27.79297447, -27.79297447, -27.81404877, -27.79223251,  -27.71365929, -27.76314735])

	d = dist((ra, dec), (ra2, dec2))
	print('distance', d)
	assert not numpy.isnan(d).any()
	
