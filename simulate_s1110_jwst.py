from astropy.io import fits
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_fwhm_to_sigma   # a useful constant
from astropy.stats import mad_std
from astropy.io import fits
import numpy as np
import subprocess #  from comment by Adam on Python me dammit.
import os
import glob
from os.path import expanduser, basename
import pyds9
np.set_printoptions(precision=3)

def get_filenames(psf_dir) :
    myfiles = [ basename(x) for x in glob.glob(psf_dir + "*.fits") ]
    return(myfiles)


# Before starting ipython, need to type this in bash term:
#>>> export XPA_METHOD="local"

d = pyds9.DS9('foo1')  # start ds9.  'd' is the way to call ds9

outdir = 'Output/'

# Here is Mike's fitting to image plane image A2, subsampled by 4x from 0.03"/pix to 0.0075"/pixel.
# This is very close to to Marshall's JWST subsampled SW PSFs that I have, of 0.0078"/pixel.
image, header_in = fits.getdata('606out4x4_sub2.fits', 2, header=True) # Grab model extension

psf_dir = '/Users/jrrigby1/Texts/Proposals/JWST/C1/S1110/Simulations/JWST_PSFs/NIRCam/SW/'

d.set_np2arr(image) # sending ndarray im directly to ds9
d.set("colorbar no")   # example of manipulating the ds9 window
d.set("scale limits 0 10")  # example of manipulating the ds9 window
d.set("zoom to 1 1")
d.set('frame 2')


#Loop this over all NIRCam filters.
# Careful, the pixelscale is different for the LW channels (0.00777) compared to SW (0.01575)
# For now, just use SW.
psf_files = get_filenames(psf_dir)
for psf_file in psf_files :
    psf, psf_header =  fits.getdata(psf_dir + psf_file, 0, header=True)
    print("Using PSF ", psf_file, psf.shape)
    # Making a small PSF cutout to speed convolution.  shape needs to be odd for convolve
    midpt = 644 # the middle pixel of the PSF
    psfwin = 50 
    small_psf = psf[midpt - psfwin : midpt + psfwin +1,    midpt - psfwin : midpt + psfwin +1]
    convolved_image = convolve(image, small_psf)
    outfile = "S1110_simulated_" + psf_header['FILTER'] + '.fits'
    fits.writeto(outdir + outfile, convolved_image, overwrite=True)
    d.set('frame 3')
    d.set_np2arr(convolved_image)
    d.set("scale limits 0 10")  # example of manipulating the ds9 window

# The imput model is the same for all 
# Should I resample with actual pixels?  Or dirzzling?

