import matplotlib.pyplot as plt
import numpy
import astropy.io.fits as pyfits
import sys
import argparse
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Johannes Buchner (C) 2013-2016 <johannes.buchner.acad@gmx.com>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('matchcatalogue', type=str, 
	help="""3way output catalogue""")

parser.add_argument('id', type=str,
	help='ID to explain (from primary catalogue)')



# parsing arguments
args = parser.parse_args()

print 'loading catalogue %s' % args.matchcatalogue
f = pyfits.open(args.matchcatalogue)
header = f[0].header
data = f[1].data
primary_id_col = header['COL_PRIM']
print '    searching for %s == %s' % (primary_id_col, args.id)
if issubclass(data.dtype[primary_id_col].type, numpy.integer):
	print 'using int'
	mask = data[primary_id_col] == int(args.id)
elif issubclass(data.dtype[primary_id_col].type, numpy.float):
	mask = data[primary_id_col] == float(args.id)
else:
	mask = data[primary_id_col] == args.id
print '    %d rows found' % (mask.sum())
if mask.sum() == 0:
	print 'ERROR: ID not found. Was searching for %s == %s' % (primary_id_col, args.id)
	sys.exit(1)

# make a plot of the positions

plt.figure(figsize=(7,7))
plt.axis('equal')
cols_ra = header['COLS_RA'].split(' ')
cols_dec = header['COLS_DEC'].split(' ')
cols_err = header['COLS_ERR'].split(' ')
center_ra = data[cols_ra[0]][mask][0]
center_dec = data[cols_dec[0]][mask][0]
markers = ['x', '+', '^', '<', '>', 'v', 'p'] * 10
colors = ['b', 'c', 'g', 'r', 'k', 'brown'] * 10
for col_ra, col_dec, col_err, marker, color in zip(cols_ra, cols_dec, cols_err, markers, colors):
	tblname, err = col_err.split('_', 1)
	if err.startswith(':'):
		err = data[tblname + '_' + err[1:]][mask]
	else:
		err = numpy.ones(mask.sum()) * float(err)
	pos = set(zip(data[col_ra][mask], data[col_dec][mask], err))
	ras = numpy.array([ra for ra, dec, err in pos if ra != -99]) - center_ra
	decs = numpy.array([dec for ra, dec, err in pos if ra != -99]) - center_dec
	errs = numpy.array([err for ra, dec, err in pos if ra != -99]) / 60. / 60.
	r,  = plt.plot(ras, decs, marker=marker, mec=color, mfc='None', ms=8, mew=2, ls=' ', label='%s %s' % (col_ra, col_dec))
	#plt.circles(ras, decs, errs, facecolor='None', edgecolor=r.get_color())
	patches = [Circle((ra, dec), err) for ra, dec, err in zip(ras, decs, errs)]
	p = PatchCollection(patches)
	p.set_facecolor('None')
	p.set_edgecolor(color)
	plt.gca().add_collection(p)
	

# go through each association and highlight
for j, i in enumerate(numpy.argsort(data['bfpost'][mask])[::-1][:3]):
	ras = []
	decs = []
	for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
		if data[col_ra][mask][i] == -99:
			continue
		ra = data[col_ra][mask][i] - center_ra
		dec = data[col_dec][mask][i] - center_dec
		ras.append(ra)
		decs.append(dec)
	
	plt.plot(ras, decs, '-', lw=(3-j), label='top %s by distance (bfpost=%.2f)' % (j+1, data['bfpost'][mask][i]))

for col in header['BIASING'].split(', '):
	for col_ra, col_dec, color in zip(cols_ra, cols_dec, markers):
		if col.split('_', 2)[0] != col_ra.split('_', 2)[0]:
			continue
		bias = data['bias_' + col]
		print bias[mask]
		mask1 = numpy.logical_and(mask, data[col_ra] != -99)
		mask2 = numpy.logical_and(mask1, bias > 1)
		if mask2.any():
			ra = data[col_ra][mask2] - center_ra
			dec = data[col_dec][mask2] - center_dec
			plt.plot(ra, dec, 's ', mew=4, ms=20, mec='g', mfc='None', label='%s strong boost > 1' % col, alpha=0.3)
		mask2 = numpy.logical_and(mask, bias < -0.5)
		if mask2.any():
			ra = data[col_ra][mask2] - center_ra
			dec = data[col_dec][mask2] - center_dec
			plt.plot(ra, dec, 'd ', mew=4, ms=20, mec='r', mfc='None', label='%s strong reject < -1' % col, alpha=0.3)

#for col in data.columns:
#	if col.startswith('bias')
#for j, i in enumerate(numpy.argsort(data['bfpost'][mask])[::-1][:3]):

mask2 = numpy.logical_and(mask, data['match_flag'] == 1)
ras = []
decs = []
first = True
for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
	ra = data[col_ra][mask2] - center_ra
	dec = data[col_dec][mask2] - center_dec
	if not first and ra[0] != -99:
		plt.text(ra[0], dec[0], ' 1', va='top', ha='left', alpha=0.5, size=16, fontweight='bold')
	first = False
	ras.append(ra[ra !=-99])
	decs.append(dec[ra !=-99])
plt.plot(ras, decs, '--', lw=7, label='Most probable association (match_flag=1)', color='orange')

mask2 = numpy.logical_and(mask, data['match_flag'] == 2)
ras = []
decs = []
first = True
for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
	ra = data[col_ra][mask2] - center_ra
	dec = data[col_dec][mask2] - center_dec
	if not first:
		for rai, deci in zip(ra, dec):
			if rai != -99:
				plt.text(rai, deci, ' 2', va='top', ha='left', alpha=0.5, size=16, fontweight='bold')
	first = False
	ras.append(ra[ra != -99])
	decs.append(dec[ra != -99])
plt.plot(ras, decs, '--', lw=5, label='Socondary, similarly probable association (match_flag=2)', color='yellow')

plt.xlabel('$\Delta$RA')
plt.ylabel('$\Delta$DEC')
xlo, xhi = plt.xlim()
ylo, yhi = plt.ylim()
hi = max(-xlo, xhi, -ylo, yhi)
plt.xlim(-hi, hi)
plt.ylim(-hi, hi)
plt.legend(loc='best', numpoints=1, prop=dict(size=8))
outfilename = '%s_explain_%s.pdf' % (args.matchcatalogue, args.id)
print 'plotting to %s' % outfilename
plt.savefig(outfilename, bbox_inches='tight')
plt.close()





