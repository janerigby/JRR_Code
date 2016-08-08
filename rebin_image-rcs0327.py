'''This code takes the source-plane reconstruction of lensed galaxy RCS0327 and bins
it down to the native (unlensed) resolution of HST>  It's a quick way to fake
the spatial resolution (but not the noise level) of CANDELS.  Written to
answer the question, 'What would RCS0327 look like to HST were it not lensed?'
jrigby, late 2015
'''

from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import Box2DKernel

import numpy as np
from skimage.measure import block_reduce
import subprocess #  from comment by Adam on Python me dammit.

thisdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/RCS0327/HST/Reconstructed"

# Calculate binning factors for RCS0327.
pixscale = 0.2 * 0.04  #(arcsec/pixel)   The scale of the input image.

# Super kludgy!  hardcoding waves here, but getting filenames from command line.  Order may change.
waves = np.array((0.98, 1.25, 1.32, 1.60, 0.39, 0.606, 0.814))
diff_lim = waves * 1E-6 / 2.5 * 206265.
print "Waves (um)", waves
print "Diff lim (\")", diff_lim

print "Diffraction limit is this many pixels:", diff_lim / pixscale

ksize = np.sqrt( (diff_lim/pixscale)**2 - 2**2)  # kernel is quad diff of final * initial resolns
# assumes 2 pixels per resoln element in reconstucted image, to be vaguely nyquist
print "Gaussian kernel size: ", ksize

# Grab the files.
dothis = "ls " + thisdir + "/*crop*fits"
print "Want to ", dothis
files = subprocess.check_output(dothis, shell=True).split()

# set up final HST pixel scale
final_scale = 0.04 #arcsec/pixel
binby = int(final_scale / pixscale)
print "Will bin to HST pixels of", final_scale, "binning by", binby 

for ii, thisfile in enumerate(files) :
    data_in, header = fits.getdata(thisfile, header=True)
    print "image should have wave", thisfile, waves[ii]

    # convolve with a Gaussian
    kernel = Gaussian2DKernel(ksize[ii])
    data_out = convolve(data_in, kernel)
    newname = thisdir + "/FG" + str(waves[ii]) + ".fits"   #output are FG[wave].fits, blurred but not rebinned
    fits.writeto(newname, data_out, header, clobber=True)

    # rebin to HST pixel sizes
    rebinned = block_reduce(data_out, block_size=(binby,binby), func=np.sum)  # in skimage
    newname = thisdir + "/FGrebin" + str(waves[ii]) + ".fits"  # blurred and rebinned output images
    fits.writeto(newname, rebinned, header, clobber=True)
