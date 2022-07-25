print('This is NWAY Python API')


import sys
import numpy
from numpy import log10, pi, exp, logical_and
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import argparse
import nwaylib 
import nwaylib.progress as progress
print('nwaylib file', nwaylib.__file__)
import nwaylib.fastskymatch as match



def table_from_fits(fitsname, poserr_value=None, area=None, magnitude_columns=[]):
	fits_table = pyfits.open(fitsname)[1]
	table_name = fits_table.name
	ra = fits_table.data['RA']
	dec = fits_table.data['DEC']
	if 'pos_err' in fits_table.data.columns.names:
		poserr = fits_table.data['pos_err']
	else:
		assert poserr_value is not None, ('"pos_err" column not found in file "%s", and no poserr_value passed' % fitsname)
		poserr = poserr_value * numpy.ones(len(ra))
	if area is None:
		area = fits_table.header['SKYAREA'] * 1.0 
	
	# magnitude columns
	mags = []
	maghists = []
	magnames = []
	#for mag in magnitude_columns:
	for col_name, magfile in magnitude_columns:
		assert col_name in fits_table.data.dtype.names
		
		mag_all = fits_table.data[col_name]
		# mark -99 as undefined
		mag_all[mag_all == -99] = numpy.nan
		
		mags.append(mag_all)
		magnames.append(col_name)
		if magfile == 'auto':
			maghists.append(None)
		else:
			bins_lo, bins_hi, hist_sel, hist_all = numpy.loadtxt(magfile).transpose()
			maghists.append((bins_lo, bins_hi, hist_sel, hist_all))
	
	return dict(name=table_name, ra=ra, dec=dec, error=poserr, area=area, mags=mags, maghists=maghists, magnames=magnames)
	# area in square degrees
	# error in arcsec
	# ra/dec in degrees
	# mag: column of something
	# maghists: either (bin, sel, all) tuple or None (for auto)



def create_shifted_catalogue(inputfile: str, outputfile:str,
shift_dec: float = 0, shift_ra: float = 0,  radius: float = 0):


    filename = inputfile
    outfile = outputfile
    print('opening', filename)
    inputfitsfile = pyfits.open(filename)
    header_hdu = inputfitsfile[0]
    table = inputfitsfile[1].data

    if shift_ra==0 and shift_dec==0:
        raise ValueError('ERROR: You have to set either shift-ra or shift-dec to non-zero')

    ra_key  = match.get_tablekeys(table, 'RA')
    print('    using RA  column: %s' % ra_key)
    dec_key = match.get_tablekeys(table, 'DEC')
    print('    using DEC column: %s' % dec_key)

    ra_orig = table[ra_key]
    dec_orig = table[dec_key]

    ra = ra_orig + shift_ra / 60. / 60
    dec = dec_orig + shift_dec / 60. / 60

    # for each of them, check that there is no collision
    excluded = []

    pbar = progress.bar(ndigits=6)
    for i, (ra_i, dec_i) in pbar(list(enumerate(zip(ra, dec)))):
        d = match.dist((ra_i, dec_i), (ra_orig, dec_orig))
        excluded.append((d * 60 * 60 < radius).any())

    excluded = numpy.array(excluded)
    print('removed %d sources which collide with original positions' % (excluded.sum()))

    table[ra_key] = ra
    table[dec_key] = dec

    newcolumns = []
    for col in table.columns:
        newcolumns.append(pyfits.Column(name=col.name, format=col.format, array=col.array[~excluded]))

    tbhdu = match.fits_from_columns(newcolumns)
    print('writing "%s" (%d rows)' % (outfile, len(tbhdu.data)))
    for k, v in inputfitsfile[1].header.items():
        if k not in tbhdu.header:
            tbhdu.header[k] = v

    hdulist = pyfits.HDUList([header_hdu, tbhdu])
    hdulist.writeto(outfile, **progress.kwargs_overwrite_true)




def calibrate_cutoff(realfile_table, fakefile_table):


    data = realfile_table
    mask = data['ncat'] == 1
    p_any0 = data['prob_has_match'][mask]
    mask = data['match_flag'] == 1
    p_any = data['prob_has_match'][mask]
    p_i = data['prob_this_match'][mask]
    
    # mask = data['ncat'] == 1
    # p_any0 = data['p_any'][mask]
    # mask = data['match_flag'] == 1
    # p_any = data['p_any'][mask]
    # p_i = data['p_i'][mask]

    data = fakefile_table
    mask = data['ncat'] == 1
    p_any0_offset = data['prob_has_match'][mask]
    mask = data['match_flag'] == 1
    p_any_offset = data['prob_has_match'][mask]
    p_i_offset = data['prob_this_match'][mask]

    # mask = data['ncat'] == 1
    # p_any0_offset = data['p_any'][mask]
    # mask = data['match_flag'] == 1
    # p_any_offset = data['p_any'][mask]
    # p_i_offset = data['p_i'][mask]

    cutoffs = numpy.linspace(0, 1, 101)

    efficiency = [(p_any0 > cutoff).mean() for cutoff in cutoffs]
    error_rate = [(p_any0_offset > cutoff).mean() for cutoff in cutoffs]



    for rate in [0.01,0.03,0.05,0.1]:
        print()
        # find where error-rate is < rate
        mask = numpy.array(error_rate) < rate
        if not mask.any():
            print('A false detection rate of <%d%% is not possible.' % (rate*100))
        else:
            i = numpy.min(numpy.where(mask)[0])
            #print(i)
            print('For a false detection rate of <%d%%' % (rate*100))
            print('--> use only counterparts with p_any>%.2f (%.2f%% of matches)' % (cutoffs[i], efficiency[i]*100))



    plt.figure(figsize=(5,5))
    plt.plot(cutoffs, efficiency, '-', color='k', label='completeness (selection efficiency)', lw=2)
    plt.plot(cutoffs, purity:=1-numpy.asanyarray(error_rate), '--', color='r', label='purity (1-false selection rate)')
    #find the intersection of the two curves: efficiency and purity
    cutoff_intersection = numpy.argmin(numpy.abs(efficiency - purity))
    print('The efficiency is %.2f%%' % (tmp_val:=efficiency[cutoff_intersection]*100))
    print('The purity is  %.2f%%' % (purity[cutoff_intersection]*100))
    plt.axvline(cutoffs[cutoff_intersection], color='k', ls='--', label=f'purity=completeness={tmp_val:.2g}%')


    plt.ylabel('fraction(>p_any)')
    plt.xlabel('p_any')
    plt.legend(loc='lower left', prop=dict(size=10))
    #plt.savefig(realfile + '_p_any_cutoffquality.pdf', bbox_inches='tight')
    #plt.savefig(realfile + '_p_any_cutoffquality.png', bbox_inches='tight')
    plt.show()
    #print('created plot "%s_p_any_cutoffquality.pdf"' % realfile)



    plt.figure(figsize=(5,5))
    plt.plot(cutoffs, efficiency, '-', color='k', label='selection efficiency', lw=2)
    plt.plot(cutoffs, error_rate, '--', color='r', label='false selection rate')

    plt.ylabel('fraction(>p_any)')
    plt.xlabel('p_any')
    plt.legend(loc='lower left', prop=dict(size=10))
    #plt.savefig(realfile + '_p_any_cutoffquality.pdf', bbox_inches='tight')
    #plt.savefig(realfile + '_p_any_cutoffquality.png', bbox_inches='tight')
    plt.show()
    #print('created plot "%s_p_any_cutoffquality.pdf"' % realfile)

    plt.figure(figsize=(5,5))
    plt.plot(p_any, p_i, '. ', color='r', label='real')
    plt.plot(p_any_offset, p_i_offset, '. ', color='gray', label='offset')

    plt.ylabel('p_i')
    plt.xlabel('p_any')
    plt.legend(loc='lower left', prop=dict(size=10))
    #plt.savefig(realfile + '_p_any_p_i.pdf', bbox_inches='tight')
    #plt.savefig(realfile + '_p_any_p_i.png', bbox_inches='tight')
    plt.show()
    #print('created plot "%s_p_any_p_i.pdf"' % args.realfile)

    return cutoffs, efficiency, error_rate
