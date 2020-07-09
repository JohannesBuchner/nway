import numpy
from numpy import pi, sin, exp, cos, log, log10

fout = open('randomcatO.csv', 'wb')
fout.write(b'ID,RA,DEC\n')

rarange  = numpy.arange(150.0, 150.1, 1./60./60.)
decrange = numpy.arange(2.0, 2.1,     1./60./60.)

print('area: %.3f deg2' % ((rarange[-1] - rarange[0]) * (decrange[-1] - decrange[0])))

for i, ra in enumerate(rarange):
	ids = 1 + i*len(decrange) + numpy.arange(len(decrange))
	numpy.savetxt(fout, numpy.transpose([ids,numpy.random.normal(decrange, 0.01/60/60)]), fmt='%%d,%.6f,%%.6f' % ra)

fout.close()

fout = open('randomcatX.csv', 'w')
fout2 = open('randomcatR.csv', 'w')
fout.write('ID,RA,DEC,pos_err,a,b,phi\n')
fout2.write('ID,RA,DEC,pos_err,a,b,phi\n')

id = 1
ractr  = rarange.mean()
decctr = decrange.mean()

lastr = 0.
for i in range(10):
	r = (1. + i**2) / 60. / 60.
	
	if i % 2 == 0:
		a  = (r - lastr) / 10.
		b  = (r - lastr) / 10.
	else:
		a  = (r - lastr) / 10. * 2
		b  = (r - lastr) / 10. / 2
	symerror = (a**2 + b**2)**0.5
	
	for j, angle in enumerate(numpy.arange(0, 2*pi, pi/6)):
		
		ra  = ractr + r * cos(angle)
		dec = decctr + r * sin(angle)
		rotangle = angle if (i+1) % 4 == 0 else angle + pi/2
		
		fout.write("%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (id, ra, dec, symerror*60*60,  a*60*60, b*60*60, rotangle/pi*180))
		fout2.write("%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (id, numpy.random.normal(ra, symerror), numpy.random.normal(dec, symerror), symerror*60*60,  a*60*60, b*60*60, rotangle/pi*180))
		id += 1
	lastr = r
fout.close()
fout2.close()
