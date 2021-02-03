#!/usr/bin/env /Users/jrrigby1/anaconda3/envs/astropy_stable3/bin/python -W ignore
from re import sub
from sys import argv
from os import remove
from ccdproc import ImageFileCollection
from ccdproc import Combiner, combine
from astropy.nddata import CCDData
from astropy.stats import sigma_clip

''' The purpose of this script is to create the residual bias ("delta dark") frame
    from IRAC Spitzer images, and then subtract if off each frame.  
    Creates bias-subtracted bcbcd images.  Was originally an iraf template script
    by jrigby 8/2014 (delta_bias.cl).  Converted to python delta_bias.py 1/2021.'''

if len(argv) != 3 :
    raise Exception("Required arguments for delta_bias.py are AOR (e.g. r15102976) and channel (e.g. ch1)")
aor = str(argv[1])
ch  = str(argv[2])

use_ds9=False

thisdir = '../../Downloaded' + '/' + aor + '/' + ch + '/bcd/'

ccd = {}  # dictionaries to hold the files
clipped = {}
clipped_ccd = {}
image_collection = ImageFileCollection(thisdir, glob_include='*cbcd.fits')  # All the files in this dir
for thisfile in image_collection.files:    
    # use astropy's sigma clip to mask with grow function.  output of sigma_clip is np masked array.  Convert to CCDData
    ccd[thisfile] = CCDData.read(thisdir + thisfile)
    clipped[thisfile] = sigma_clip(ccd[thisfile], sigma_lower=10, sigma_upper=3, maxiters=None, grow=2)
    clipped_ccd[thisfile] = CCDData(clipped[thisfile].data, mask=clipped[thisfile].mask, unit='MJy / sr')

# Median combine with extra sigma clipping, and write to a file
combiner = Combiner(list(clipped_ccd.values()))  
combiner.sigma_clipping(low_thresh=3, high_thresh=3)
combined_image = combiner.median_combine()
outdir = 'Dbias_' + ch + '_' + aor + '/'
outfile = ch + '_' + aor + '_dbias.fits'
combined_image.write(outdir + outfile, overwrite=True)

if use_ds9 :
    import pyds9
    d = pyds9.DS9('foo1')
    d.set("file " + thisdir + image_collection.files[0])
    d.set("frame 2")
    d.set_np2arr(clipped[image_collection.files[0]].mask.astype('float64'))   # This isn't a file, but a numpy array
    d.set("frame 3")
    d.set("file " + outfile)

# Now, subtract the dbias.fits file from each individual file, and write to Dbias_ch1_r15102976/ folder
# Had to fix this so it doesn't toss the header, which will break mopex downstream
for hdu, fname in image_collection.hdus(return_fname=True) :
    header = hdu.header
    data = hdu.data
    data_ccd = CCDData(hdu.data, unit='MJy/sr', header=header) # clunky.  Needed to iterate all of header, data, filename
    subtracted = data_ccd.subtract(combined_image, handle_meta='first_found') # This saves the header from data_ccd
    newname=sub('_cbcd.fits', '_bcbcd.fits', fname)
    subtracted.write(outdir + newname, overwrite=True)   
print("Subtracted delta bias from each cbcd frame, to make bcdbd frames.")

print("All done delta bias python script.")
