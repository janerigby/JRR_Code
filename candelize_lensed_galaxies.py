'''candelize_lensed_galaxies.py
This code takes the source-plane reconstruction of a lensed galaxy, convolves it to 
the normal diffraction limit of a telescope without lensing, then bins it to a more
normal pixel scale.  Running with HST resoln and pixel scale to fake 
the spatial resolution (but not the noise level) of CANDELS.  Written to
answer the question, 'What would RCS0327 look like to HST were it not lensed?'
jrigby, late 2015.  Improved 2/2016 to run on s1110.  Now uses real HST PSFs.
'''
from astropy.io import fits
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import Box2DKernel
from astropy.stats import gaussian_fwhm_to_sigma   # a useful constant
from astropy.stats import mad_std
from scipy.misc import imresize
import numpy as np
from skimage.measure import block_reduce
import subprocess #  from comment by Adam on Python me dammit.
import os


def candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam) :
    '''Inputs:
    thisdir :     what directory to look for the input src plane reconstructed images.
    outdir  :     what directory to write the output
    waves:        np array containing wavelengths (in micron) of input images
    in_pixscale:  the input pixel scale of the reconstructed src plane images, in arcsec/pixel
    out_pixscale: the desired pixel scale of the output 'candelized images', in arcsec/pixel
    tel_diam:     the diameter of the telescope in meters.  Assumes diffraction-limited performance
    Outputs:      creates output images in outdir
    '''

    np.set_printoptions(precision=3)
    diff_lim = waves * 1E-6 / tel_diam * 206265.  # Diffraction limit FWHM (")
    print "Waves (um)\t\t", waves
    print "Diff lim (\")\t\t", diff_lim
    dif_in_pix = diff_lim / in_pixscale    # still a fwhm
    print "Diff limit (pix)\t", dif_in_pix 

    ksize  = np.sqrt( (dif_in_pix)**2 - 2**2)  # kernel is quad diff of final & initial resolns.
    print "Gaussian fwhm (pix):\t", ksize
    ksize *=  gaussian_fwhm_to_sigma           # convert from fwhm to stdev needed by Gaussian2DKernel
    # assumes 2 pixels per resoln element in reconstucted image, to be vaguely nyquist
    print "Gaussian sigma (pix):\t", ksize

    # set up final pixel scale
    print "Input scale (\"/pix)\t", in_pixscale
    print "Output scale (\"/pix)\t", out_pixscale
    binby = int(out_pixscale / in_pixscale)  # this is the bin factor needed by rebinned, for each wave
    print "Binning by\t\t", binby 

    # Grab the files.
    dothis = "ls " + thisdir + "*out*fits"
    print "DEBUG, Want to ", dothis
    files = subprocess.check_output(dothis, shell=True).split()

    if not os.path.exists(outdir):  # Make output directory if it doesn't already exist.
        os.makedirs(outdir)

    print "#image  wave (check that wave matches file, since it's kludged!)"
    for ii, thisfile in enumerate(files) :
        data_in, header = fits.getdata(thisfile, header=True)
        print thisfile, waves[ii]
                
        # convolve with a Gaussian
        kernel = Gaussian2DKernel(ksize[ii], factor=10)  #ksize is a stdev, not fwhm
        kernel.normalize()
        kernname = outdir + "kern" +  str(waves[ii]) + ".fits"
        fits.writeto(kernname, kernel.array, clobber=True)  # Write the kernel to a file, for galfit later
    
        data_out = convolve(data_in, kernel)
        newname = outdir + "/FGgaus" + str(waves[ii]) + ".fits"   #output are FG[wave].fits, blurred but not rebinned
        fits.writeto(newname, data_out, header, clobber=True)                
        # rebin
        rebinned = block_reduce(data_out, block_size=(binby,binby), func=np.sum)  # in skimage
        newname = outdir + "/FGgausrebin" + str(waves[ii]) + ".fits"  # blurred and rebinned output images
        fits.writeto(newname, rebinned, header, clobber=True)

        # Instead of a gaussian, use Kate's measured PSFs as kernel
        psf_dir = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/PSF/"
        psf_file = psf_dir + "F" + str(waves[ii]) + "_psf.fits"
        psf_in =  fits.getdata(psf_file)
        psf_subsampled = imresize(psf_in, 0.03/in_pixscale)  # Kate's PSFs are 0.03"/pix.  Input src plane is finer
        #print "DEBUG, psf was resampled to ", psf_subsampled.shape, " from ", psf_in.shape
        data_out = convolve_fft(data_in, psf_subsampled)
        newname = outdir + "/FGpsf" + str(waves[ii]) + ".fits"   #output are FGpsf[wave].fits, blurred by psf but not rebinned
        fits.writeto(newname, data_out, header, clobber=True)
        # rebin
        rebinned = block_reduce(data_out, block_size=(binby,binby), func=np.sum)  # in skimage
        newname = outdir + "/FGpsfrebin" + str(waves[ii]) + ".fits"  # blurred by psf and rebinned output images
        fits.writeto(newname, rebinned, header, clobber=True)
        
        # Write 1-sigma uncertainty image, for galfit
        #
        #print "DEBUGGING, mad_std was", mad_std(rebinned)
        sig_out      = np.ones(shape=rebinned.shape) * mad_std(rebinned) # kludge, testing, crap!
        sig_out_name = outdir + "FGrebin" + str(waves[ii]) + "_fakesig.fits"  # fake input sigma images for galfit.
        fits.writeto(sig_out_name, sig_out, clobber=True)
        #print "DEBUGGING, wrote sigma file ", sig_out_name 
        
    print "Wrote CANDELIZED images (convolved to normal diffraction limit of HST, binned to normal HST pixel scale)."
    print "They are the FGrebin*fits images in"
    print "\t", outdir
    return(0)  # success

do_s1110   = True
do_rcs0327 = False
do_s1438   = False

# Let's candelize s1110 
if do_s1110 :
    tel_diam = 2.4003  # telescope diameter (m).  Currently set for HST
    in_pixscale = 0.1 * 0.03  #(arcsec/pixel)   The pixel scale of the input source plane reconstructions
    # From Keren: pixel scale = 0.1 of image plane pixel scale (0.03") = 0.003"/pix
    out_pixscale = 0.03   #arcsec/pixel   A compromise for CANDELS, which used 0.03" UVIS, 0.06" IR.
#
# Super kludgy!  hardcoding waves here, but get filenames from command line.  Order may change.  Check output,
# that waves correspond to filename
    waves = np.array((1.05, 1.60, 0.39, 0.606))  #units are micron
    thisdir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/drizzled_source_images_23Feb2016/Aa/"
    outdir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/Candelize/Aa/"
    candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam)  # image Aa

    waves = np.array((1.05, 1.60))  #units are micron
    thisdir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/drizzled_source_images_23Feb2016/Aa_minus_cluster_gal/"
    outdir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/Candelize/Aa_minus_cluster_gal/"
    candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam)  # image Aa minus cluster galaxy

    waves = np.array((1.05, 1.60, 0.39, 0.606))  #units are micron
    thisdir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/drizzled_source_images_23Feb2016/Ab/"
    outdir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/Candelize/Ab/"
    candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam)  # image Ab
#
    thisdir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/drizzled_source_images_23Feb2016/Ac/"
    outdir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/Candelize/Ac/"
    candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam)  # image Ac

# Next, candelize RCS0327
if do_rcs0327 :
    thisdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/RCS0327/HST/Reconstructed/"
    outdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/RCS0327/HST/Reconstructed/Candelize/"
    tel_diam = 2.5  # telescope diameter (m).  Currently set for HST
    in_pixscale = 0.2 * 0.04  #(arcsec/pixel)
    out_pixscale =  0.04 #arcsec/pixel
    waves = np.array((1.25, 1.60, 0.39, 0.606, 0.814))
    candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam)

# Next, candelize SDSS1438
if do_s1438 :
    # Warning, b/c Kate Whitaker didn't create an F140W PSF, we don't have one.  I used the F125W PSF instead.  KLUDGE.
    thisdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/s1438p14/HST/sourceplane/Aa/"
    outdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/s1438p14/HST/Candelized/Aa/"
    tel_diam = 2.5  # telescope diameter (m).  Currently set for HST
    in_pixscale = 0.03 * 0.2  # From README.  "scale = 0.2 (of the original pixel scale, which was 0.03)"
    out_pixscale =  0.03 #arcsec/pixel
    waves = np.array((1.40, 0.606, 0.814))
    candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam)

    thisdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/s1438p14/HST/sourceplane/Ab/"
    outdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/s1438p14/HST/Candelized/Ab/"
    candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam)

    thisdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/s1438p14/HST/sourceplane/Ac/"
    outdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/s1438p14/HST/Candelized/Ac/"
    candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam)

    thisdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/s1438p14/HST/sourceplane/CI/"
    outdir  = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/s1438p14/HST/Candelized/CI/"
    candelize_srcplane(thisdir, outdir, waves, in_pixscale, out_pixscale, tel_diam)

print "ALL DONE."
