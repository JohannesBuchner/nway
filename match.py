import numpy
import itertools
import os, sys
import pyfits

filenames = sys.argv[1:]

numpy.random.seed(0)
#tables = [numpy.loadtxt(t) for t in tables]

def gen():
	ra  = numpy.random.uniform(size=200)
	dec = numpy.random.uniform(size=200)
	return numpy.array(zip(ra, dec), dtype=[('ra', 'f'), ('dec', 'f')])

def array2fits(table, extname):
	cat_columns = pyfits.ColDefs([pyfits.Column(name=n, format='E',array=table[n]) 
		for n in table.dtype.names])
	tbhdu = pyfits.new_table(cat_columns)
	hdu = pyfits.PrimaryHDU()
	import datetime, time
	now = datetime.datetime.fromtimestamp(time.time())
	nowstr = now.isoformat()
	nowstr = nowstr[:nowstr.rfind('.')]
	hdu.header.update('CREATOR', """Johannes Buchner <jbuchner@mpe.mpg.de>""")
	hdu.header.update('DATE', nowstr)
	hdu.header.update('ANALYSIS', 'match test tables')
	tbhdu.header.update('EXTNAME', extname)
	hdulist = pyfits.HDUList([hdu, tbhdu])
	return hdulist

for fitsname in filenames:
	if os.path.exists(fitsname):
		continue
	print 'generating toy data for %s!' % fitsname
	hdulist = array2fits(gen(), fitsname.replace('.fits', ''))
	hdulist[0].header.update('GENERAT', 'match test table, random')
	hdulist.writeto(fitsname, clobber=False)

tables = [pyfits.open(fitsname)[1] for fitsname in filenames]
table_names = [t.name for t in tables]
tables = [t.data for t in tables]

err = 0.03

buckets = {}

print 'naive', numpy.product([len(t) for t in tables])

import progressbar
pbar = progressbar.ProgressBar(widgets=[
	progressbar.Percentage(), progressbar.Counter('%3d'),
	progressbar.Bar(), progressbar.ETA()], maxval=sum([len(t) for t in tables])).start()

for ti, table in enumerate(tables):
	for e in table:
		ra, dec = e['ra'], e['dec']
		i, j = int(ra / err), int(dec / err)
		
		#for jj, ii in (j-1,i-1), (j-1,i), (j-1,i+1), (j,i-1), (j,i), (j,i+1), (j+1,i-1), (j+1,i), (j+1, i+1):
		for jj, ii in (j,i), (j,i+1), (j+1,i), (j+1, i+1):
			#k = int(ii + 10000*jj)
			k = (ii, jj)
			if k not in buckets:
				buckets[k] = [[] for _ in range(len(tables))]
			#print ' bucket', k, 'table', ti, e
			buckets[k][ti].append(e)
		pbar.update(pbar.currval + 1)
pbar.finish()

"""
pbar = progressbar.ProgressBar(widgets=[
	progressbar.Percentage(), progressbar.Counter('%3d'),
	progressbar.Bar(), progressbar.ETA()], maxval=3).start()

hashes0 = []
hashes1 = []
hashes2 = []
hashes3 = []
print 'hashing ...'
for ti, table in enumerate(tables):
	ra, dec = table['ra'], table['dec']
	i = numpy.array(ra / err / 3,  dtype=numpy.int)
	j = numpy.array(dec / err / 3, dtype=numpy.int)
	hashes0.append(i + 10000 * j)
	hashes1.append(i + 1 + 10000 * j)
	hashes2.append(i + 10000 * (j + 1))
	hashes3.append(i + 1 + 10000 * (j + 1))

pbar.update(pbar.currval + 1)
print 'collecting hashes...'
all_hashes = set()
for hashes in [hashes0, hashes1, hashes2, hashes3]:
	for k in hashes:
		all_hashes.update(k)

pbar.update(pbar.currval + 1)
print 'bucketing...'
buckets = {}
for h in all_hashes:
	buckets[h] = []
	for ti, table in enumerate(tables):
		buckets[h].append(table[numpy.logical_or(
			numpy.logical_or(hashes0[ti] == h, hashes1[ti] == h),
			numpy.logical_or(hashes2[ti] == h, hashes3[ti] == h))])

pbar.update(pbar.currval + 1)
pbar.finish()
"""

results = []
def are_close(*l):
	# all pairs
	for i, a in enumerate(l):
		for b in l[:i]:
			if numpy.abs(a['ra'] - b['ra']) > err or numpy.abs(a['dec'] - b['dec']) > err:
				return False
	return True

print ' combining ...'
# now combine in buckets
print 'hashing', numpy.sum([numpy.product([len(li) for li in lists]) for lists in buckets.values()])

pbar = progressbar.ProgressBar(widgets=[
	progressbar.Percentage(), progressbar.Counter('%3d'), 
	progressbar.Bar(), progressbar.ETA()], maxval=len(buckets)).start()
for lists in buckets.values():
	if any([len(li) == 0 for li in lists]):
		continue
	#	results += [tuple(pair) for pair in itertools.product(*lists)]
	results += [tuple(pair) for pair in itertools.product(*lists) if are_close(*pair)]
	#print ' bucket', k, [len(li) for li in lists]
	#print '    ', lists
	#print '    ', product
	pbar.update(pbar.currval + 1)
pbar.finish()

# now make results unique by sorting
print 'filtered', len(results)
results = set(results)
print 'unique', len(results)

table = []

keys = []
for r in results:
	row = []
	for t in r:
		row += list(t)
	table.append(row)
for table_name, table in zip(table_names, tables):
	keys += ["%s_%s" % (table_name, n) for n in table.dtype.names]
print keys
table = numpy.array(table, dtype=[(k, 'f') for k in keys])
hdulist = array2fits(table, 'MATCH')
hdulist[0].header.update('ANALYSIS', 'match table from' + ', '.join(table_names))
hdulist[0].header.update('INPUT', ', '.join(filenames))
hdulist.writeto('match.fits', clobber=True)


print 'plotting'
import matplotlib.pyplot as plt

for ti in range(len(tables)):
	plt.plot([r[ti]['ra'] for r in results], [r[ti]['dec'] for r in results], '.')

for r in results:
	plt.plot([t['ra'] for t in r], [t['dec'] for t in r], '-', color='b', alpha=0.3)

plt.savefig('matchtest.pdf')

