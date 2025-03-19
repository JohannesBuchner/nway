#!/usr/bin/env python
from __future__ import division, print_function

import inspect
import sys

import astropy.io.fits as pyfits

if len(sys.argv) != 4:
	sys.stderr.write("""SYNOPSIS: %s <catalogue.fits> <tablename> <skyarea>

tablename: name of the catalogue
skyarea: catalogue area in square degrees

Author: Johannes Buchner (C) 2013-2025
""" % sys.argv[0])
	sys.exit(1)

f = pyfits.open(sys.argv[1])
print('current', f[1].name, 'SKYAREA:', f[1].header.get('SKYAREA', None))
assert '_' not in sys.argv[2], 'Table name must not contain underscore "_".'
f[1].name = sys.argv[2]
f[1].header['SKYAREA'] = float(sys.argv[3])
print('new    ', f[1].name, 'SKYAREA:', f[1].header.get('SKYAREA', None))

args = inspect.signature(pyfits.writeto).parameters if hasattr(inspect, 'signature') else inspect.getargspec(pyfits.writeto).args
if 'overwrite' in args:
	arg_overwrite = 'overwrite'
else:
	arg_overwrite = 'clobber'
f.writeto(sys.argv[1], **{arg_overwrite:True})
