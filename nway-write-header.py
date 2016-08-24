import astropy.io.fits as pyfits
import sys

f = pyfits.open(sys.argv[1])
print 'current', f[1].name, 'SKYAREA:', f[1].header.get('SKYAREA', None)
f[1].name = sys.argv[2]
f[1].header.update('SKYAREA', float(sys.argv[3]))
print 'new    ', f[1].name, 'SKYAREA:', f[1].header.get('SKYAREA', None)
f.writeto(sys.argv[1], clobber=True)



