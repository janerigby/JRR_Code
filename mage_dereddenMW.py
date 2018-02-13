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

labels = ('planckarc_pos1', 'planckarc_slit4a', 'planckarc_slit4bc', 'planckarc', 'PSZ0441_slitA', 'PSZ0441_slitB', 'PSZ0441', 'SPT0310_slitA', 'SPT0310_slitB', 'SPT0310', 'SPT2325')
# Adding "Friends of Megasaura" sample

(spec_path, line_path) = jrr.mage.getpath(mage_mode)
#speclist = jrr.mage.wrap_getlist(mage_mode, which_list="all", zchoice=zchoice, MWdr=False)  # MWdr=False to read the files that are not corrected for MWreddening
speclist = jrr.mage.wrap_getlist(mage_mode, which_list="labels", labels=labels, zchoice=zchoice, MWdr=False)  # MWdr=False to read the files that are not corrected for MWreddening
for label in speclist.index:  # For each file in spectra-filenames-redshifts.txt
    infile = speclist.filename[label]
    EBV = speclist['EBV_MW'][label]
    print "Looking at ", label, infile, EBV, 
    header = jrr.util.read_header_from_file(infile, comment="#")  # save the header
    header += ("# MW reddening correction of E(B-V)="+str(EBV)+" , Rv=3.1, CCM has been applied\n")
    sp = pandas.read_table(spec_path + infile, delim_whitespace=True, comment="#", header=0, dtype=np.float64)
    jrr.util.check_df_for_Object(sp)   # Check if Pandas has erroneously read column as Object instead of float64
    jrr.mage.deredden_MW_extinction(sp, EBV, colwave='wave', colf='fnu', colfu='noise', colcont='cont_fnu', colcontu='cont_uncert')
    sp.replace([np.inf, -np.inf], np.nan, inplace=True)
    sp.replace(["Infinity", "-Infinity"], np.nan, inplace=True)
    outfile = re.sub('.txt', '_MWdr.txt', infile)
    jrr.util.check_df_for_Object(sp)   # Check if Pandas has erroneously read column as Object instead of float64
    sp.to_csv('temp', sep='\t', index=False, float_format="% 1.8E", na_rep=' NaN          ')
    jrr.util.put_header_on_file('temp', header, outfile)
     
