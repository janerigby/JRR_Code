''' Deredden the individual MagE spectra for Milky Way extinction, as per suggestion of 
referee of Stacked paper.  jrigby, July 2017
Want to walk through every comb1.txt or combwC1.txt file, run it through MW dereddening, and spit out the same file
but corrected for MW reddening.'''

import jrr
import pandas
import matplotlib.pyplot as plt
import numpy as np
import re
import to_precision  

mage_mode = 'reduction'
zchoice = 'stars'

(spec_path, line_path) = jrr.mage.getpath(mage_mode)
speclist = jrr.mage.wrap_getlist(mage_mode, which_list="all", zchoice=zchoice)
for label in speclist.index:  # For each file in spectra-filenames-redshifts.txt
    infile = speclist.filename[label]
    EBV = speclist['EBV_MW'][label]
    print "Looking at ", label, infile, EBV
    header = jrr.util.read_header_from_file(infile, comment="#")  # save the header
    sp = pandas.read_table(spec_path + infile, delim_whitespace=True, comment="#", header=0)
    jrr.mage.deredden_MW_extinction(sp, EBV, colwave='wave', colf='fnu', colfu='noise', colcont='cont_fnu', colcontu='cont_uncert')
    sp.replace([np.inf, -np.inf], np.nan, inplace=True)
    # Man. this to_precision() stuff is clunky, but I can't figure out how to do this in pandas.
    sp['wave']   = sp['wave'].map(  lambda x: to_precision.to_precision(x, 8, notation='standard'), na_action='ignore')
    sp['obswave'] = sp['obswave'].map(lambda x: to_precision.to_precision(x, 8, notation='standard'), na_action='ignore')
    sp['fnu'] = sp['fnu'].map(lambda x: to_precision.to_precision(x, 4, notation='sci'), na_action='ignore')
    sp['noise'] = sp['noise'].map(lambda x: to_precision.to_precision(x, 4, notation='sci'), na_action='ignore')
    if 'cont_fnu' in sp.keys() :
        sp['cont_fnu']    = sp['cont_fnu'].map(   lambda x: to_precision.to_precision(x, 4, notation='sci'), na_action='ignore')
        sp['cont_uncert'] = sp['cont_uncert'].map(lambda x: to_precision.to_precision(x, 4, notation='sci'), na_action='ignore')
    outfile = re.sub('.txt', '_MWdr.txt', infile)
    sp.to_csv('temp', sep='\t', index=False)
    jrr.util.put_header_on_file('temp', header, outfile)
     
