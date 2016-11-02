from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import pandas
from re import sub

def read_chuck_UV_spec(infile, uncert_file=None, outfile=None) :
    ''' Takes a 1D spectrum with the wavelength stuck in the WCS header,          
    gets the wavelength array out, and packages the spectrum into a nice          
    pandas data frame.  Returns DF, and dumps a pickle file and csv file too.                  
    This was written to reach Chuck Steidel's stacked spectrum, and is not
    general enough to make into a function.  Example of munging: the wcs.dropaxis business.'''
    sp = fits.open(infile)
    header = sp[0].header
    wcs = WCS(header)
    wcs2 = wcs.dropaxis(1)  # Kill the WCS's dummy 2nd dimension                  
    index = np.arange(header['NAXIS1'])
    temp =  (np.array(wcs2.wcs_pix2world(index, 0))).T
    wavelength = (10**temp)[:,0]
    fnu = sp[0].data
    if uncert_file :
        sp2 = fits.open(uncert_file)      # Get the uncertainty                   
        fnu_u = sp2[0].data
    else : fnu_u = np.zeros_like(fnu)

    # Make a pandas data frame                                                    
    foo = np.array((wavelength, fnu, fnu_u))
    print "DEBUG", foo.shape, wavelength.shape, fnu.shape, fnu_u.shape
    df = pandas.DataFrame(foo.T, columns=("wave", "fnu", "fnu_u"))
    if not outfile :
        outfile = sub(".fits", ".p", infile)
    df.to_pickle(outfile)
    txtfile = sub(".fits", ".csv", infile)
    df.to_csv(txtfile)    
    return(df)


def read_mosfire_temp(infile, outfile=None) :
    sp = fits.open(infile)
    header = sp[0].header
    wcs = WCS(header)
    wcs2 = wcs.dropaxis(1)  # Kill the WCS's dummy 2nd dimension
    index = np.arange(header['NAXIS1'])
    temp =  (np.array(wcs2.wcs_pix2world(index, 0))).T
    wavelength = (temp)[:,0]
    (flam, flam_u) = sp[0].data

    # Make a pandas data frame
    foo = np.array((wavelength, flam, flam_u))
    print "DEBUG", foo.shape, wavelength.shape, flam.shape, flam_u.shape
    df = pandas.DataFrame(foo.T, columns=("wave", "flam", "flam_u"))
    if not outfile :
        outfile = sub(".fits", ".p", infile)
    df.to_pickle(outfile)
    return(df)

infile = "KBSS-LM1.uv.fnu.fits"
uncert_file = "KBSS-LM1.uv.fnu.sig.fits"
df = read_chuck_UV_spec(infile, uncert_file)
txtfile = sub(".fits", ".csv", infile)
df.to_csv(txtfile, sep='\t')    

mosfire = ("KBSS-LM1.H.flam.fits", "KBSS-LM1.J.flam.fits", "KBSS-LM1.K.flam.fits")
for infile in mosfire :
    df2 = read_mosfire_temp(infile)
    txtfile = sub(".fits", ".csv", infile)
    df2.to_csv(txtfile, sep='\t')    



