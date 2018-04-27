import jrr
import numpy as np
import pandas
#import seaborn as sns
from astropy.io import fits
import pyds9
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from os.path import expanduser

def adjust_ds9(d):
    d.set("colorbar no")   # example of manipulating the ds9 window
    d.set("scale zscale")  # example of manipulating the ds9 window
    d.set("zoom to fit")
    d.set("zoom to 4 4")
    return(0)

# for running on satchmo
#imdir = "/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Mage/Redux/Mage_Nov2015/"
#infiles = ("N2/s0033-2d-210/s0033-2d-210-015sum.fits", "N2/s0033-2d-211/s0033-2d-211-015sum.fits", "N3/s0033-2d-297/s0033-2d-297-015sum.fits", "N3/s0033-2d-298/s0033-2d-298-015sum.fits")

# for running in the Dropbox/S0033 folder
homedir = expanduser("~")
imdir = homedir + "/Dropbox/S0033/MagE/2D_Lya/"
infiles = ("s0033-2d-210-015sum.fits", "s0033-2d-211-015sum.fits", "s0033-2d-297-015sum.fits", "s0033-2d-298-015sum.fits")


frame_names = [210, 211, 297, 298]

d = pyds9.DS9('foo1')  # start ds9.  'd' is the way to call ds9
plt.ion()
plt.close("all")

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

fig2, ax2 = plt.subplots(4, 1, figsize=(8,5))

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)

fig5 = plt.figure()
ax5 = fig5.add_subplot(111)


pix_scale = 0.3 # arcsec per pixel for MagE
stupid_colors=('k', 'blue', 'green', 'red')
stupid_colors=("lightgray", "darkgray", "gray", "black")

# Was working on summed file.  But that blurs things out.  Do on indy exposures.
#sum_im, sum_header = fits.getdata(imdir + "s0033_sum.fits", header=True)
#sum_im = sum_im.astype(np.float64)
#d.set_np2arr(sum_im) # sending ndarray im directly to ds9
#adjust_ds9(d)

#sns.palplot(sns.color_palette("Blues"))
##sns.set_style("white")
#sns.set_style("whitegrid", {'axes.grid' : False})

df = {}
popt = {}  # best fit continuum shapes
ftweak     = {}

for jj, thisfile in enumerate(infiles): 
    print imdir + thisfile
    data_in, header = fits.getdata(imdir + thisfile, header=True)
    im = data_in.astype(np.float64)
    d.set("frame "+ str(jj+1))
    d.set_np2arr(im) # sending ndarray im directly to ds9
    adjust_ds9(d)

    sum_in_Lya = []
    theiis     = []

    lyareg =  (1049, 1120, 11, 28)
    contreg = (1120, 1500)
    
    # Deal with the continuum.  Using pandas
    df[jj] = pandas.DataFrame(data={'cont' : np.nanmedian(im[ : , contreg[0]:contreg[1]], axis=1)})
    guess_pars = (5.0, 17., 3., 0.0)   # amp, x0, sigma, continuum
    (popt[jj], df[jj]['contfit']) = jrr.spec.fit_quick_gaussian(df[jj], guess_pars, colwave='index', colf='cont')
    thisx = (df[jj].index - popt[jj][1]) * pix_scale

    df[jj]['sumLya'] = [np.nansum(im[ii : ii+1, lyareg[0] : lyareg[1]]) for ii in range(0, im.shape[0])]  # wizard FDR list comprehension
    # What the above step does is sum the Lya profile for every row of the iamge

    ax1.plot( df[jj].index,   df[jj]['sumLya'] / df[jj]['sumLya'].sum() )
    width_Lya = lyareg[1] - lyareg[0]   # how many pixels is Lya summed over
    df[jj]['contscaled'] = df[jj]['cont'] * width_Lya  # Get the fluxing the same for cont, Lya.  Integrated over same # of pixels
    df[jj]['sumLya_contsub'] =  df[jj]['sumLya'] - df[jj]['contscaled']
    
    ftweak[jj] = df[jj]['sumLya_contsub'].sum() / df[0]['sumLya_contsub'].sum()
    ax5.plot(thisx,    df[jj]['sumLya_contsub'] /ftweak[jj],  color=stupid_colors[jj], label="")
    ax5.scatter(thisx, df[jj]['sumLya_contsub'] /ftweak[jj],  color=stupid_colors[jj], label="frame "+ str(frame_names[jj]))

    ax4.plot(thisx,    df[jj]['contscaled']/ftweak[jj], color=stupid_colors[jj], label="")
    ax4.scatter(thisx, df[jj]['contscaled']/ftweak[jj], color=stupid_colors[jj], label="frame "+ str(frame_names[jj]))

    dslice = 3
    lyarange = range(lyareg[2] +2, lyareg[3] -2 -dslice, dslice)
    slices = [ np.nansum(im[ii : ii+dslice, lyareg[0] : lyareg[1]], axis=0) for ii in lyarange]  # wizard FDR list comprehension
    for ii, thisslice in enumerate(slices) :
        label = "slice pixs"+ str(lyarange[ii]) + " to " + str(lyarange[ii] + dslice)
        ax2[jj].plot(thisslice, label=label, color = plt.cm.gist_rainbow( ii *1.0 / len(slices)))
        ax2[jj].annotate("frame "+ str(frame_names[jj]), xy=(0.05, 0.89), xycoords='axes fraction', fontsize=10)
 
ax1.set_title("Summed image, scaled by total Lya flux")
#fig2.set_title("Lya profile, not scaled")
ax4.set_title("Continuum spatial profile")
ax5.set_title("Summed (over velocity) Lya spatial profile, continuum subtracted")
ax1.set_xlim(1010,1130)
ax4.set_xlim(-12*pix_scale, 12*pix_scale)
ax5.set_xlim(-12*pix_scale, 12*pix_scale)
ax1.set_ylim(-0.01,0.07)
ax4.set_ylim(-10,700)
ax5.set_ylim(-100,7000)
ax1.set_xlabel("x pixel value")
ax2[-1].set_xlabel("x pixel value - " + str(lyareg[0]))
ax4.set_xlabel("offset from center of continuum (\")")
ax5.set_xlabel("offset from center of continuum (\")")
ax1.legend()
ax2[0].legend()
ax4.legend()
ax5.legend()



