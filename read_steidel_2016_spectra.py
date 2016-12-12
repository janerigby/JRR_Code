from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import pandas
from re import sub
import jrr

#infile = "KBSS-LM1.uv.fnu.fits"
#uncert_file = "KBSS-LM1.uv.fnu.sig.fits"
df = jrr.mage.convert_chuck_UVspec()

mosfire = ("KBSS-LM1.H.flam.fits", "KBSS-LM1.J.flam.fits", "KBSS-LM1.K.flam.fits")
for infile in mosfire :
    df2 = jrr.mage.convert_chuck_mosfire(infile)
    txtfile = sub(".fits", ".csv", infile)
    df2.to_csv(txtfile, sep='\t')    
