# This is a code to empirically measure the spectral resolution, aka
# the line spread function, for MagE spectra.  It measures this from
# skylines that were processed in the same way as the science spectra.
# This code finds skylines via a peak detection method, tosses peaks
# with another peak nearby, fits a single Gaussian to each isolated
# peak, and dumps the resolution R === fwhm/lambda to a file.
# For calling sequence, see Mage-redux-tools/README
# RUN THIS FROM IPYTHON WITH %PYLAB, OR YOU HAVE TO CLICK THROUGH EVERY PLOT.
# jrigby, Feb 2016.  

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import jrr   # Jane's library.  Contains MagE-specific functions
import scipy.signal as signal
from  scipy.optimize import curve_fit
from time import sleep
from math import log, sqrt
from os import system
import subprocess
import re
from sys import exit

GOSLOW = False
mage_mode = "reduction"  # Look at copies of the files on Satchmo
#mage_mode = "released"   

def MyGaussian(x, *params):   # Define a Gaussian with a linear continuum under it, to fit each skyline
    y = np.zeros_like(x)
    for i in range(0, len(params), 4):
        bb = params[i]
        aa = params[i+1]
        cc = params[i+2]
        cont =params[i+3]
        y = y + (aa * np.exp((x - bb)**2/(-2*cc**2))  + cont)
    return y

batch2 = jrr.mage.organize_labels('batch2')   # Adding "Friends of Megasaura" sample  (batch 2)
batch3 = jrr.mage.organize_labels('batch3')   # Adding Apr2018 and Aug2018 runs

#temp_filename =  "RCS0327/KnotB/rcs0327-knotB-comb1.txt"  # Prototyped on this file.  Now, run for all.
#temp_filename = "RCS0327/KnotE/rcs0327-knotE-allres-comb1.txt"

run_just_some = True
if run_just_some:     # A useful way to run just one spectrum, instead of a bunch
# Adding "Friends of Megasaura" sample
    labels = batch2
    labels = batch3  # now doing April 2018 and Aug 2018 runs
    specs = jrr.mage.wrap_getlist(mage_mode, which_list='labels', labels=labels, MWdr=False)
    
else :    # Normal mode
    (specs) = jrr.mage.getlist(mage_mode)  # get list of MagE spectrum filenames and redshifts
    # specs holds the filenames and redshifts, for example   specs['filename'], specs['z_stars']
    infile = specs['filename']

for jj in range(0, len(specs)) :                  #flam_stack[jj] will be jj spectrum
    (specdir, linedir) = jrr.mage.getpath(mage_mode)
    label     = specs['short_label'][jj]
    filename  = specs['filename'][jj]
    zz =  specs['z_neb'][jj]

    if re.search("-combw?C?1.txt", filename) :  # Megasaura with continuum
        outfile   =  specdir + re.sub("-combw?C?1.txt", "-specresln.txt", filename)
        print "DEBUG, preparing to write outfile from original filename ", outfile, filename
        # Error checking, don't want to overwrite filename
        if (outfile == filename) or  (re.search("-combw?C?1.txt", filename)<1) :
            exit("DIED, filename not what was expected")  # die if it didn't change the filename
    else :
        outfile   =  specdir + re.sub("-comb.txt", "-specresln.txt", filename)  # new files no fit cont
        print "DEBUG, preparing to write outfile from original filename ", outfile, filename
        # Error checking, don't want to overwrite filename
        if (outfile == filename) or  (re.search("-comb.txt", filename)<1) :
            exit("DIED, filename not what was expected")  # die if it didn't change the filename
            
    outpng    = outfile.replace(".txt", "")
    (sp_temp, resoln, dresoln) = jrr.mage.open_spectrum(filename, zz, mage_mode)
    print("I opened file ", filename) 

    sp = sp_temp[sp_temp['wave_sky'].gt(5000)].copy()   # Don't look for skylines in the blue
    #sp['fnu_sky'].fillna(-99, None, inplace=True)  # correct for some infinite values
    sp.fnu_sky.replace([np.inf, -np.inf], np.nan, inplace=True)  # replace inf w NaN to keep mad happy
    delta = sp.fnu_sky.mad()  * 2 # delta for peakfinder
    maxtab, mintab = jrr.peakdet.peakdet(sp.fnu_sky,delta)  # Find peaks.  
    peak_ind =  [np.int(p[0]) for p in maxtab]
    ym1      =  [p[1] for p in maxtab]
    plt.clf()
    #plt.ylim(0,5E-15)
    plt.xlim(5000,8500)
    plt.plot(sp.wave_sky, sp.fnu_sky)  # Show the sky spectrum, and flag detected peaks
    plt.scatter(sp['wave_sky'].iloc[peak_ind], sp['fnu_sky'].iloc[peak_ind])
    plt.xlabel("Observed wavelength (A)")
    plt.ylabel("fnu of sky spectrum")
    plt.draw()
    sleep(2.)

    goodpeaks = [] # empty list of good peaks
    RR = []        # empty list of resolutions
    #DEBUG print "Found peaks at ", sp.wave_sky[peak_ind] 
    guess = []
    win = 15 # window of pixels around the 
#    print "#lambda(A)  FWHM(A)   R(lambda/fwhm)";
    for index,item in enumerate(peak_ind) :   # step through the peaks
        plt.plot(sp.wave_sky.iloc[item-win:item+win], sp.fnu_sky.iloc[item-win:item+win], color='b')
        plt.title("Fitting peak " + str(item))

        if ((item == 0) or (item == peak_ind[-1])) :
            dist2blue = 1000.
            dist2red  = 1000.
        else : 
            dist2blue = sp.wave_sky.iloc[item] - sp.wave_sky.iloc[peak_ind[index-1]] 
            dist2red  = sp.wave_sky.iloc[peak_ind[index+1]] - sp.wave_sky.iloc[item]
        minsep = 10.0 #A, test
        #print "DEBUGGING dist2blue dist2red ", dist2blue, dist2red
        if((dist2blue < minsep) or (dist2red < minsep)) :
            pass
            #print "DEBUG Skipping peak at ", sp.wave_sky[item], " because other peak nearby ", dist2blue, dist2red
        else :
            try:
                guess = [sp.wave_sky.iloc[item], sp.fnu_sky.iloc[item], 3, delta*5]
                popt, pcov = curve_fit(MyGaussian, sp.wave_sky.iloc[item-win:item+win], sp.fnu_sky.iloc[item-win:item+win], p0=guess)
                fit = MyGaussian(sp.wave_sky.iloc[item-win:item+win], *popt)
                fwhm = 2*sqrt(2.*log(2)) * abs(popt[2])   # convert sig to FWHM
                RR.append(popt[0]/fwhm)  # save the R  for later.
                goodpeaks.append(sp.wave_sky.iloc[item])
                plt.plot(sp.wave_sky.iloc[item-win:item+win], fit, color='r')
                #plt.show()
                plt.scatter(sp.wave_sky.iloc[peak_ind], sp.fnu_sky.iloc[peak_ind])
                plt.xlim(sp.wave_sky.iloc[item-win],sp.wave_sky.iloc[item+win])
                plt.ylim(0, sp.fnu_sky.iloc[item]*1.3)
                #plt.draw()
                if GOSLOW :
                    sleep(0.5)
            except:
                print "skipping peak at ", sp.wave_sky.iloc[item], "because bad fit"
    plt.clf()
    # Now, plot the measured R vs wavelength, and fit a simple func to it.
    plt.scatter(goodpeaks,  RR)
    median = round(np.median(RR))
    mad = round(jrr.util.mad(RR))
    plt.plot([5000,9000], [median, median])
    plt.ylim(1000,7000)
    plt.text(5000,5500, "Median " + str(median) + " +/- " + str(mad))
    plt.text(5000,5700, filename)
    #plt.show()
    #plt.draw()
    plt.savefig(outpng)
    sleep(1.)
    if GOSLOW :
        sleep(2.)
    ascii.write([goodpeaks,  RR], outfile, names=['wave_sky', 'Resoln'])
    with open(outfile, "a") as myfile:
        myfile.write("#Measured from isolated skylines, from sky spectrum which was combined with same weighting as object spectra, ")
        myfile.write("\n#  but using observed wavelength, with no barycentric correction.") 
        myfile.write("\n#Observed wavelengths are in angstroms.  R=== lambda/fwhm")
        myfile.write("\n#Measured median spectral resoln was:\n")
        myfile.write("MEDIAN  " + filename + "   " + str(round(median)))
        myfile.write("  +/-  " + str(mad) + "  (MEDIAN +/- MAD)\n")  
    system("cat " + outfile)

    # More pythonic way to write the measured spectral resolution to the header
    jrr.util.replace_text_in_file("RESOLUTIONGOESHERE", str(median), filename)
    jrr.util.replace_text_in_file("MADGOESHERE", str(mad), filename)
    if re.search("wC1", filename) :         #  Change the *comb1.txt file as well.
        comb1 = re.sub("wC1", "1", filename)
        jrr.util.replace_text_in_file("RESOLUTIONGOESHERE", str(median), comb1)
        jrr.util.replace_text_in_file("MADGOESHERE", str(mad), comb1)
        
