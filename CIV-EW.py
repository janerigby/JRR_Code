''' Measure the C IV equivalent width (EW) in MageE spectra, following simple index approach
    of Heckman et al. 1998, same used by Crowther et al 2006 and Eldridge & Stanway 2012.
    jrigby, july, oct 2015
'''
from __future__ import print_function
from builtins import str
from builtins import range
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import string
import subprocess
import re
from  astropy.io import ascii
import jrr   # Jane's library.  Contains equivalent width calculator and fnu-to-flambda converter

scale = 1E17  # for plotting, multiply by scale to get y axis ~ unity
mage_mode = "reduction"  # Look for files on satchmo, not released versions on Dropbox

# What should the He II range be?
HeII_range = (1630.0, 1650.0)

# Wavelength range over which to measure CIV index.
# This is the BROAD rest-frame wavelength range to measure CIV emission, from Heckman 1998
# This huge spectral range is good for comparing to Eldridge & Stanway 2012, but DUMB for
# spectra with decent spectral resolution.
CIV_range1 = (1507.0, 1553.0);

# I came up with a good NARROW range, from looking at higher-resoln spectra:
# Shapley et al 2003 composite:  range of 1551--1565 A works well.
# Pettini et al. 2002 cB58:      range of 1551--1560 A is even better.
# RCS0327 Knot E:                cB58 range works well here.
###  Tweaking this to 1551.5 to avoid absorption in 0900
CIV_range2  = (1551.5, 1560.0)  

test_photosph_range = (1280.,1380.) 

ranges = (CIV_range1, CIV_range2, HeII_range, test_photosph_range)

outdir = "REDO/"  
outfile = outdir + "CIV_HeII.out"
png = string.replace(outfile, ".out", ".png")


#Functions:


def plot_a_subpanel( feature, rest_wave, rest_flam, rest_cont, short_label) :    # feature is civ or heii
    pad = 50    # plot continuum on either side of the line
    start = feature[0][0] - pad
    end   = feature[0][-1] + pad
    plt.step(rest_wave[start:end], scale*rest_flam[start:end], color='k')  # plot spectrum as steps
    plt.plot(rest_wave[feature], scale*rest_cont[feature], color='g')          # plot feature
    plt.locator_params(nbins=4)  # sane tick separataion on axes
    plt.title(short_label, fontsize='small', verticalalignment='bottom')
    plt.tick_params(axis='x', which='both', bottom='on', top='off', labelbottom='off')
    plt.xlim( rest_wave[start], rest_wave[end])   # don't fuck with my xrange pyplot!

def plot_one_spec(rest_wave, rest_flam, rest_flam_u, rest_cont, rest_cont_u, short_label, counter) :
    Ncol = len(ranges)    
    for indx, thisrange in enumerate(ranges) :
        ind = np.where(  (rest_wave > thisrange[0])  & (rest_wave < thisrange[1])  )   # indices
        plot_number = (Ncol * counter) + (indx + 1)
        print("Plotting subplot", plot_number)
        plt.subplot(Nspectra, Ncol, plot_number)   # Plot a panel, for example civ
        plot_a_subpanel( ind, rest_wave, rest_flam, rest_cont, short_label)
        if (counter == Nspectra - 1) :   # last subplot gets labels
            plt.xlabel(r'wavelength ($\AA$)')
            plt.tick_params(axis='x', which='both', bottom='on', top='off', labelbottom='on')
            scaleformat = '%.0E' %  (scale**-1)
            if(indx == 0) :
                plt.ylabel(r'f$_\lambda$ (' + scaleformat + r'$\ erg\ s^{-1}\ cm^{-2}\ \AA^{-1}$)')
        if (short_label == "S2111-0114" or short_label == "S1429+1202") :  # Mark HeII that has Skyline contamination
            xmin,xmax = plt.ylim()
            ymin,ymax = plt.ylim()
            put_at_y = ymin + (ymax - ymin)*0.12
            plt.text(1635., put_at_y, "Skyline contam",fontsize='small')
    return(0)

def header_EW_outfile(out) :
    out.write("#index ranges were " + str(ranges))
    out.write("\n#Measured rest-frame EWs and uncertainties, in units of Angstroms.  Negative is emission, positive is absorption\n")

def measure_the_EWs(rest_wave, rest_flam, rest_flam_u, rest_cont, rest_cont_u, short_label, zz, out) :
    civ  = np.where(  (rest_wave > CIV_range[0])  & (rest_wave < CIV_range[1])  )   # indices 
    heii = np.where(  (rest_wave > HeII_range[0]) & (rest_wave < HeII_range[1])  )           
    # Calculate dispersion.  Kludgy but works
    dciv  = np.average(np.array([  (rest_wave[civ][1] - rest_wave[civ][0]),   (rest_wave[civ][-1] - rest_wave[civ][-2])  ]))
    dheii = np.average(np.array([  (rest_wave[heii][1] - rest_wave[heii][0]),   (rest_wave[heii][-1] - rest_wave[heii][-2])  ]))
    # Calculate EW
    (EW_CIV, uncert_EWCIV)   = jrr.calc_EW( rest_flam[civ],  rest_flam_u[civ],  rest_cont[civ],  rest_cont_u[civ],  dciv, 0.0)
    (EW_HeII, uncert_EWHeII) = jrr.calc_EW( rest_flam[heii], rest_flam_u[heii], rest_cont[heii], rest_cont_u[heii], dheii, 0.0)
    out.write("%30s    %5.3f      %5.3f      %5.3f      %5.3f" % (short_label, EW_CIV, uncert_EWCIV, EW_HeII, uncert_EWHeII))
    return(0)

def measure_EW_better(rest_wave, rest_flam, rest_flam_u, rest_cont, rest_cont_u, out) :
    for thisrange in ranges :
        ind = np.where(  (rest_wave > thisrange[0])  & (rest_wave < thisrange[1])  )   # indices 
        dispersion = np.average(np.array([  (rest_wave[ind][1] - rest_wave[ind][0]),   (rest_wave[ind][-1] - rest_wave[ind][-2])  ]))
        (EW, EW_uncert) = jrr.calc_EW( rest_flam[ind],  rest_flam_u[ind],  rest_cont[ind],  rest_cont_u[ind],  dispersion, 0.0)
        out.write(" %5.3f  %5.3f   " % (EW, EW_uncert))
    return(0)
    
def simple_automated_CIVHeII_Cont(wave, flam, flam_u) :  #A simple way to get continuum over the CIV, HeII region
    x1 = 1485.
    x2 = 1690.
    dx = 10.
    ind1 = np.where((wave > (x1-dx)) & (wave < (x1+dx)) )  # clean continuum on blue end
    ind2 = np.where((wave > (x2-dx)) & (wave < (x2+dx)) )  # clean continuum on red end
    ind3 = np.where((wave > (x1-dx)) & (wave < (x2+dx)) )  # span over CIV + HeII + continuum
    y1 = np.average(flam[ind1])
    y2 = np.average(flam[ind2])
    slope = (y1-y2)/(x1-x2)
    b = y1 - slope*x1
    fake_cont = slope*wave[ind3] + b
    return(wave[ind3], flam[ind3], flam_u[ind3], fake_cont)

def measure_EWs_for_a_Bpass_model(infile, outfile) :
    out = open(outfile, 'w')
    header_EW_outfile(out)
    out.write( ("#filename is actually time for model " + infile + "\n"))
    print("Measuring CIV and HeII EWs from bigass bpass model", infile)
    s = ascii.read(infile, comment="#")   # input has units of **flam**
    wave = s['col1']  # wave for all times
    plt.close()
    plt.figure(num=1, figsize=(15,5))
    for nn in range(2,40) :
        colname = "col" + str(nn)
        flam = s[colname]  
        flam_u = np.zeros(shape=flam.shape) 
        time = 10**(6+0.1*(nn -2))
        # Make a simple linear continuum fit from clean regions for CIV and HeII
        (wave_x, flam_x, flam_u_x, cont_x) = simple_automated_CIVHeII_Cont(wave, flam, flam*0.)
        plt.plot(wave_x, flam_x/np.average(flam_x), color="black")
        plt.plot(wave_x, cont_x/np.average(flam_x), color="green")
        pretty_time = "%.3E" % (time)
        out.write( (pretty_time + "  "))
        measure_EW_better(wave_x, flam_x, flam_x*0., cont_x, cont_x*0., out)
        out.write("\n")
    out.close()
    return(0)

            
# MAIN: for each MagE spectrum, measure EW(CIV) and EW(HeII), and plot CIV and HeII
override_file = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Mage/Analysis/Stacked_spectra/spectra-to-stack.txt"
(specs, Nspectra) =jrr.Mage_get_filename_redshift_list(mage_mode, override_file)  # get list of MagE spectrum filenames and redshifts
Nspectra +=1 # Add one so that we can ad the MagE composite spectrum
plt.close('all')
f, axs = plt.subplots(Nspectra, 1, figsize=(6,20))
f.subplots_adjust(left=0.12, bottom=0.06, right=0.91, top=0.96, wspace=0.2, hspace=0.2) #from slider adjust subplot param
out = open(outfile, 'w')
header_EW_outfile(out)
out.write("#First set of EWs uses continuum fits.  Second set uses automated slope continuum.\n")
print("Measuring CIV and HeII EWs for the MagE spectra")
for ii in range(0, Nspectra - 1) :  # Measure EW for each MagE spectrum
    infile      = specs['filename'][ii]
    short_label = specs['short_label'][ii]
    zz          = specs['z_neb'][ii]
    if re.search("comb1.txt", infile) :
        pass 
    else :
        (wave, flam, flam_u, obswave, temp_avgsky, cont, cont_u)  = jrr.Mage_open_spectrum(infile, "flam", mage_mode);
        (rest_wave, rest_flam, rest_flam_u) = jrr.convert2restframe(wave, flam, flam_u, zz)
        (junk     , rest_cont, rest_cont_u) = jrr.convert2restframe(wave, cont, cont_u, zz)
        # should now have arrays of rest wavelength, flambda, and continuum, as
        # rest_wave, rest_flam, rest_flam_u, rest_cont, rest_cont_u
        # Measure the EWs using my beautiful continuum fits
        plot_one_spec(  rest_wave, rest_flam, rest_flam_u, rest_cont, rest_cont_u, short_label, ii)  # split plotting from measurements
        out.write( (short_label + "  "))
        measure_EW_better(rest_wave, rest_flam, rest_flam_u, rest_cont, rest_cont_u, out)
        # Measure the EWs using an automated continuum fitting
        (wave_x, flam_x, flam_u_x, cont_x) = simple_automated_CIVHeII_Cont(rest_wave, rest_flam, rest_flam_u)
        measure_EW_better(wave_x, flam_x, flam_u_x, cont_x, cont_x*0., out)
        out.write("\n")
    
# Want to add the Shapley composite, but can't until I write a continuum for it.

# Measure the MagE composite spectrum.
print("Measuring CIV and HeII EWs for the MagE composite spectrum")
short_label = "MagE-stack"
(restwave, X_avg, X_clipavg, X_median, X_sigma, X_jack_std, Ngal)=  jrr.Mage_open_stacked_spectrum(mage_mode)  # in fnu
restflam = jrr.fnu2flam(restwave, X_avg)
restflam_u = jrr.fnu2flam(restwave, X_sigma)
restcont = np.ones_like(restflam)   # Already continuum normalized
restcont_u = 0.1 * restflam_u  # kludge an uncertainty on continuum
plot_one_spec(  restwave, restflam/scale, restflam_u/scale, restcont/scale, restcont_u/scale, short_label, Nspectra-1)
out.write( (short_label + "  "))
measure_EW_better(restwave, restflam/scale, restflam_u/scale, restcont/scale, restcont_u/scale, out)
(wave_x, flam_x, flam_u_x, cont_x) = simple_automated_CIVHeII_Cont(restwave, restflam, restflam_u)
measure_EW_better(wave_x, flam_x, flam_u_x, cont_x, cont_x*0., out) 
out.write("\n")
out.close()
plt.savefig(png)
plt.show()


# Turn off bpass until i get the EW measurement recoded

# Measure CIV and HeII indices from BPASS models
# develpped in bpass-test.py, now moving over
bpassdir = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Mage/Models/BPASS/"
thisfile = "BPASSv2_imf135_100/OUTPUT_CONT/spectra.z020.dat.gz"
dothis = "ls " + bpassdir + "*/*/*gz"   # temp just solar models
modelnames = subprocess.check_output(dothis, shell=True).split()  # grab the model names

for bpassfile in modelnames :
    temp = string.replace(bpassfile, bpassdir, "") # remove the path the path
    outfile = outdir + string.replace( string.replace(temp, "/", "_"), ".dat.gz", ".out")
    bpasspng = string.replace(outfile, ".out", ".png")
    measure_EWs_for_a_Bpass_model(bpassfile, outfile)
    plt.plot(wave_x, flam_x/np.average(flam_x), color="blue")  # overplot the MagE composite spectrum
#    plt.xlim(1450,1700)
    plt.xlim(1280,1700)
    plt.ylim(0,2)
    plt.savefig(bpasspng)
#    plt.show()
 
