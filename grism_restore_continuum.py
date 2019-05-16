''' Inexplicably, some of Jane's continuum fits (*_cont.fits) have gone missing for S2340, G102.  
Rather than remake them from scratch, will attempt to retrieve them from the *wcontMWdr.txt  files.
Let's not overthink this.
jrigby, May 2019''' 

import glob
from os.path import exists, basename
from os import makedirs
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas
from matplotlib import pyplot as plt

wrongMWdir = 'WrongMWdr/'
newcontdir = 'Restorecont/'

wcontMWdr_files =  [ basename(x) for x in glob.glob(wrongMWdir + "*_wcontMWdr.txt") ]    # Find the files. 

for wcontMWdr in wcontMWdr_files :
    #wcontMWdr = 'sgas2340_image2_roll043_G141_wcontMWdr.txt'  # testing
    newcontfile  = wcontMWdr.replace('_wcontMWdr.txt', '_cont.fits')
    if not exists(newcontdir):   makedirs(newcontdir)
    df = pandas.read_csv(wrongMWdir + wcontMWdr, comment="#")
    df2 = df.copy(deep=True)
    for thiskey in ('wave', 'flam', 'flam_u', 'cont', 'MWredcor') :
        df2[thiskey] /= df2['MWredcor']  # undo the MW dereddening correction
    hdu = fits.PrimaryHDU(df2['cont'].values)
    hdu.writeto(newcontdir + newcontfile, overwrite=True)
    print("Wrote new file", newcontdir + newcontfile)
    
#
##  Let's figure out how a continuum file is supposed to look.
#data_in, header = fits.getdata(contfile, header=True)
#sp_cont = Table(np.vstack(data_in)).to_pandas()
## Oh, it's just 1D array of numbers.
