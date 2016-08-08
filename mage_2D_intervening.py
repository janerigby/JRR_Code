''' Looking at 2D spatial profile of interveing absorber.  This is a test case, can generalize later
if this works.  jrigby, apr 2016.'''
import os
import jrr
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def make_one_2D_absline_plot(filename, linewave, linewin, contwave, contwin):
    plt.clf()
    ax1 = plt.subplot2grid((1,5), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((1,5), (0, 2), colspan=1, sharey=ax1)
    ax3 = plt.subplot2grid((1,5), (0, 3), colspan=1, sharey=ax1)
    ax4 = plt.subplot2grid((1,5), (0, 4), colspan=1, sharey=ax1)
    label = str(filename) + "_" + str(thiswave) + "_z" + str(zz)
    plt.annotate(label, (0.1,0.96), xycoords="figure fraction")

    xx = jrr.Mage_2D_iswave(filename, contwave)  # Show the continuum
    cont, header = fits.getdata(filename, header=True) 
    ccutout = cont[0:18, xx-contwin : xx+contwin]
    ax1.imshow(ccutout, cmap='gray', interpolation='none', origin='lower')
    cprofile = np.median(ccutout, axis=1)
    ax1.annotate("Continuum", (0.1,0.96), xycoords="axes fraction")
    ax1.autoscale(False)

    xx = jrr.Mage_2D_iswave(filename, thiswave)  # Show the absorption line
    data, header = fits.getdata(filename, header=True) 
    cutout = data[0:18, xx-win : xx+win]
    profile = np.median(cutout, axis=1)
    ax2.imshow(cutout, cmap='gray', interpolation='none', origin='lower')
    ax2.annotate("Line", (0.1,0.96), xycoords="axes fraction")
    ax2.autoscale(False)

    # show the profiles
    ax3.step(profile, np.arange(profile.size))
    ax3.step(cprofile, np.arange(cprofile.size))
    ax3.annotate("Median", (0.1,0.96), xycoords="axes fraction")

    # show the ratio of the profiles
    ax4.step( profile / cprofile , np.arange(cprofile.size), color='k')
    ax4.step( np.zeros_like(cprofile), np.arange(cprofile.size), color='r')
    ax4.step( np.ones_like(cprofile), np.arange(cprofile.size), color='r')
    ax4.annotate("Ratio", (0.1,0.96), xycoords="axes fraction")
    ax4.set_xlim(-0.5,1.5)

    outfile = label + ".png"
    plt.savefig(outfile)
    plt.show()

    
# Let's look at the interveing MgII absorber in S1527
zz = 1.2823
wave = np.array((2796., 2803)) * (1.0+zz)
contwave = np.array((6340.))
win = 10 # window around the wavelength.  *KLUDGE
contwin = win *2
filename = ("s1527-2d-278/s1527-2d-278-009sum.fits", "s1527-2d-278/s1527-2d-278-010sum.fits", "s1527-2d-280/s1527-2d-280-009sum.fits", "s1527-2d-280/s1527-2d-280-010sum.fits")
cwd = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Mage/Redux/Mage_July2008/Mage_pipe/Slit2_bin1x2"
os.chdir(cwd)

for thisfile in filename :
    for thiswave in wave :
        make_one_2D_absline_plot(thisfile, thiswave, win, contwave, contwin)  # do all the work here.
