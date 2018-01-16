import jrr
import numpy as np
from astropy.io import fits
import pyds9
from matplotlib import pyplot as plt

def adjust_ds9(d):
    d.set("colorbar no")   # example of manipulating the ds9 window
    d.set("scale zscale")  # example of manipulating the ds9 window
    d.set("zoom to fit")
    d.set("zoom to 4 4")
    return(0)

imdir = "/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Mage/Redux/Mage_Nov2015/"
infiles = ("N2/s0033-2d-210/s0033-2d-210-015sum.fits", "N2/s0033-2d-211/s0033-2d-211-015sum.fits", "N3/s0033-2d-297/s0033-2d-297-015sum.fits", "N3/s0033-2d-298/s0033-2d-298-015sum.fits")

d = pyds9.DS9('foo1')  # start ds9.  'd' is the way to call ds9

# Mike wants to see profile in summed image
sum_im, sum_header = fits.getdata(imdir + "s0033_sum.fits", header=True)
sum_im = sum_im.astype(np.float64)
d.set_np2arr(sum_im) # sending ndarray im directly to ds9
adjust_ds9(d)
dslice = 4
for ii in range(13,25,dslice) :
    specslice  = np.nansum(sum_im[ii : ii+dslice, :], axis=0)
    sum_in_Lya = np.nansum(sum_im[ii : ii+dslice, 1066:1120])
    scaled_specslice = specslice  / sum_in_Lya
    plt.plot(scaled_specslice, label=ii)
plt.title("Summed image")
plt.xlim(1010,1130)
plt.ylim(-0.01,0.07)
plt.legend()
plt.show()

for thisfile in infiles:  # temp, just do first one
    print imdir + thisfile
    data_in, header = fits.getdata(imdir + thisfile, header=True)
    im = data_in.astype(np.float64)
    d.set_np2arr(im) # sending ndarray im directly to ds9
    adjust_ds9(d)
    dslice = 4
    for ii in range(13,25,dslice) :
        specslice = np.nansum(im[ii : ii+dslice, :], axis=0)
        sum_in_Lya = np.nansum(im[ii : ii+dslice, 1066:1120])
        scaled_specslice = specslice  / sum_in_Lya
        plt.plot(scaled_specslice, label=ii)
    plt.title(thisfile)
    plt.xlim(1010,1130)
    plt.ylim(-0.01,0.07)
    plt.legend()
    plt.show()



    
    
