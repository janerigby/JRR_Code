'''jwstify_lensed_galaxies.py
This code takes the source-plane reconstruction of a lensed galaxy, convolves it to 
the normal diffraction limit of a telescope without lensing, then bins it to a more
normal pixel scale.  Now running with JWST PSFs from Webbpsf.
jrigby, 7/2016, based on candelized_lensed_galaxies.py
'''
from __future__ import print_function
from astropy.io import fits
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_fwhm_to_sigma   # a useful constant
from astropy.stats import mad_std
from scipy.misc import imresize
import numpy as np
from skimage.measure import block_reduce
import subprocess #  from comment by Adam on Python me dammit.
import os


def srcplane_to_jwst(indir, in_images, JWST_inst, JWST_filters, in_pixscale, outdir) :
    '''Take a source-plane reconstruction of a lensed galaxy, and predict what JWST would see.
    Inputs:
    indir :       directory to find the src plane reconstructed images.
    in_images:    list of input images, in indir
    JWST_inst:    which JWST instrument to get the PSFs for
    JWST_filters: list of JWST filters to get PSF for.  Must be same shape as in_images
    in_pixscale:  the input pixel scale of the reconstructed src plane images, in arcsec/pixel
    outdir  :     what directory to write the output
    Outputs:      creates output images in outdir
    '''

    # Here is where to put Marshall Perrin's WebbPSF PSFs
    psf_dir = "/Volumes/Apps_and_Docs/WORK/JWST/WebbPSF/"

    if len(in_images) != len(JWST_filters) :
        raise Exception("in_images and JWST_filters must be same shape.")
    
    np.set_printoptions(precision=3)

    if not os.path.exists(outdir):  # Make output directory if it doesn't already exist.
        os.makedirs(outdir)

    print("Pixel scales:  [PSF, subsampled] [src plane input]  [output_desired]   [output_got].  All arcsec\pix")
                
    for ii, thisfile in enumerate(input_images) :#
        print(thisfile, JWST_inst, JWST_filters[ii], end=' ')
        data_in, header_in = fits.getdata(indir + thisfile, header=True)
        psf_file = psf_dir + "PSF_" + JWST_inst + "_" + JWST_filters[ii] + "_revV-1.fits"
        #print "DEBUGGING, psf file is", psf_file

        # Load Marshall's subsampled PSF.  It was subsampled by 4x, according to the WebbPSF documentation.
        (psf_in, header) =  fits.getdata(psf_file, ext=0, header=True)
        psf_pixscale = header['pixelscl']

        #I don't think we want to use Marshall's un-subsampled PSF. But it gives the native pixel scale.
        (psf_not_subsamp, psf_not_subsamp_header) =  fits.getdata(psf_file, ext=1, header=True)  # not subsampled
        out_pixscale = psf_not_subsamp_header['pixelscl']  # the native size of the instrument pixels, from WebbPSF header.
                
        psf_resize_factor = psf_pixscale / in_pixscale 
        psf_subsampled = imresize(psf_in, psf_resize_factor)  # Kate's PSFs are 0.03"/pix.  Input src plane is finer
       #print "DEBUG, psf was resampled to ", psf_subsampled.shape, " from ", psf_in.shape
        data_out = convolve_fft(data_in, psf_subsampled)
        newname = outdir + JWST_inst + JWST_filters[ii]  + "_conv.fits"    # output is convolved by PSF but not rebinned
        header_in['pix_scale'] = in_pixscale
        fits.writeto(newname, data_out, header_in, clobber=True)
        # rebin
        binby = int(out_pixscale / in_pixscale)  # this is the bin factor needed by rebinned
        # THIS IS A KLUDGE, because block_reduce can only handle integer downsampling factors.  Should do something more sophisticated.
        #print "Binning by", binby, ", but should be by", out_pixscale/in_pixscale
        rebinned = block_reduce(data_out, block_size=(binby,binby), func=np.sum)  # in skimage
        newname = outdir + JWST_inst + JWST_filters[ii]  + "_conv_rebin.fits"   # output is convolved by PSF and rebinned
        header_in['pix_scale'] = in_pixscale * binby
        fits.writeto(newname, rebinned, header_in, clobber=True)
        print("Pixel scales:  ", psf_pixscale, in_pixscale, out_pixscale,  in_pixscale * binby)

        
    print("Wrote JWST-ified images (convolved to JWST PSF, binned to normal JWST pixel scale.)")
    print("They are the _conv_rebin.fits images in ")
    print("\t", outdir)
    return(0)  # success


instrument = ("NIRCam", "NIRISS", "NIRSpec", "MIRI")  # JWST instrument names        
do_s1110   = True

# Let's candelize s1110

if do_s1110 :
    indir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/drizzled_source_images_23Feb2016/Ab/"
    outdir  = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/s1110-paper2/JWSTify/Ab/"
    input_images = ("sdssj1110p6459_F606W_0.03g0.8_drc_sci_core_Ab_s01_out.fits", "sdssj1110p6459_F105W_0.03g0.5_drz_sci_core_Ab_s01_out.fits", "sdssj1110p6459_F160W_0.03g0.5_drz_sci_core_Ab_s01_out.fits", "sdssj1110p6459_F160W_0.03g0.5_drz_sci_core_Ab_s01_out.fits")
    HST_filters  = ('F606W', 'F105W', 'F160W') # HST filters of input data.
    JWST_filters = ('F070W', 'F115W', 'F150W', 'F460M') # Filters I want to map this to in JWST NIRCam. May want non-obvious mapping
    in_pixscale = 0.1 * 0.03  #(arcsec/pixel)   The pixel scale of the input source plane reconstructions
    # From Keren: pixel scale = 0.1 of image plane pixel scale (0.03") = 0.003"/pix
    srcplane_to_jwst(indir, input_images, "NIRCam", JWST_filters, in_pixscale, outdir)
#
print("ALL DONE.")
