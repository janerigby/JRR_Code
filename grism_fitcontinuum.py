from __future__ import print_function
import jrr
import pandas
import numpy as np
from scipy.interpolate import CubicSpline
from astropy.io import fits
from astropy.table import Table
import glob
import re
import os
from matplotlib import pyplot as plt

''' This is a simple script to fit continuua to the WFC3 G102 and G41 grism spectra of S1723 and S2340.
The inputs are the spectra reduced using Grizli by Michael Florian.  The outputs are the spectra with 
an extra column, the continuum fit.  Continuum fitting is done by XIDL x_continuum b/c I couldn't find
an equivalent package within Python.  In principle, all 3 steps can be run at once, but in practice it
was easier to run them one at a time bc it took several tries to get the IDL continuum fits right.  
It shouldn't matter what dir this script gets run in.
jrigby, may 2018
'''

#### Setup  #########   By hand RUN THIS FOR EACH DIRECTORY, BOTH GRISMS
which_grism = "102" # "141" # 
which_dirnum = 0 # Step through the dirs by hand
Make_IDL_script = False
Run_IDL_script  = False
Paste_spectra_continuua = True
################

home = os.path.expanduser("~")
wdir1 = home + '/Dropbox/Grism_S1723/Grizli_redux/'
wdir2 = home + '/Dropbox/Grism_S2340/Grizli_redux/'
#              0X                      1                       2                                 3X                 4
alldirs = (wdir1 + '1Dsum/', wdir1 + '1Dbyclumps/', wdir1 + '1D_complete_images_A2_A3/', wdir2 + '1Dsum/', wdir2 + '1Dbyclumps/') # first run
#alldirs = (wdir2 + '1Dbyclumps_individuals/', wdir2 + '1Dbyclumps_stacks/') # catch-up for Michael's new S2340 extractions
grismdir = alldirs[which_dirnum]  # Step through
contdir  = grismdir + "Wcont/"

thefiles =  [ os.path.basename(x) for x in glob.glob(grismdir + "*"+which_grism+"*txt") ]   
idlscript = 'fit_continuum_script' + which_grism + '.idl'
grism_header = "# Wavelength (A)  Flux (ergs/sec/cm2/A)   Flux Uncertainty (1 sigma)  JRR fitted continuum\n" # From Florian's reductions
extra_header = "# Milky Way dereddening (MWdr) has been applied by grism_fitcontinuum.py\n"
EBV = jrr.grism.get_MWreddening_S1723()         # Get the Milky Way reddening value



## PART 1:  Write the IDL script to fit the continuum.
if Make_IDL_script :
    print("Keystroke S then DONE -- not q -- when finished, or IDL will go nuts.")
    f = open(idlscript, 'w')  # Over-ride the IDL script to fit the continuum
    if not os.path.exists(contdir):  os.makedirs(contdir)  # Make the continuum dir to hold output, if it does not already exist
    for specinfile in thefiles :
        root = re.sub('.txt', '', specinfile)
        f.write("contfile   = \"" + contdir + root + '_cont.fits\"\n')
        f.write("specinfile = \"" + grismdir + specinfile + "\"\n")
        f.write("readcol, specinfile,  f=\'d,d,d\', wave, fnu, dfnu\n")
        f.write("x_continuum, fnu, dfnu, wave=wave, lsig=2, inflg=4, xsize=1800, ysize=600, OUTFIL=contfile\n;\n")
    f.close()  # Close the cloudy_script

    
## PART 2:  RUN the idl script   idl> @fit_continuum_script.idl
if Run_IDL_script :
    os.system("/Applications/exelis/idl85/bin/idl " + idlscript)

## PART 3:  Combine the spectra w their continuum fits into one file.
if Paste_spectra_continuua :
    for specinfile in thefiles :
        root = re.sub('.txt', '', specinfile)
        colnames = ('wave', 'flam', 'flam_u')
        sp_grism = pandas.read_table(grismdir + specinfile, delim_whitespace=True, comment="#", skiprows=1, names=colnames)
        contfile =  re.sub('.txt', '_cont.fits', specinfile)
        wcontfile =  re.sub('.txt', '_wcontMWdr.txt', specinfile)
        data_in, header = fits.getdata(contdir + contfile, header=True)
        if data_in.shape[0] != sp_grism.shape[0] : raise Exception("ERROR: continuum file and parent spectrum file have different # of pixels.")
        #t = Table(np.vstack(data_in))  # Temp step, into Astropy tables to get rid of bigendian littlendian
        sp_cont = Table(np.vstack(data_in)).to_pandas()
        if sp_grism.shape[0] != sp_cont.shape[0] : raise Exception("ERROR: continuum file and parent spectrum data frames have different lengths.")
        sp_grism['cont'] = sp_cont['col0']   #  Insert continuum column into spectrum
        jrr.spec.deredden_MW_extinction(sp_grism, EBV, colwave='wave', colf='flam', colfu='flam_u', colcont='cont')  # Apply Milky Way dereddening (MWdr)
        sp_grism.to_csv('temp', index=False, na_rep='NaN')
        jrr.util.put_header_on_file('temp', grism_header+extra_header,  contdir + wcontfile)
