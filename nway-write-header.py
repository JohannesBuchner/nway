import astropy.io.fits as pyfits
import sys
import inspect

if len(sys.argv) != 4:
	sys.stderr.write("""SYNOPSIS: %s <catalogue.fits> <tablename> <skyarea>

tablename: name of the catalogue
skyarea: catalogue area in square degrees

Author: Johannes Buchner (C) 2014-2017
""" % sys.argv[0])
	sys.exit(1)

f = pyfits.open(sys.argv[1])
print 'current', f[1].name, 'SKYAREA:', f[1].header.get('SKYAREA', None)
f[1].name = sys.argv[2]
f[1].header['SKYAREA'] = float(sys.argv[3])
print 'new    ', f[1].name, 'SKYAREA:', f[1].header.get('SKYAREA', None)

if 'overwrite' in inspect.getargspec(pyfits.writeto).args:
	arg_overwrite = 'overwrite'
else:
	arg_overwrite = 'clobber'
f.writeto(sys.argv[1], **{arg_overwrite:True})



