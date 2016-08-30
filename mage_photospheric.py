import jrr
import os
import sys
import re
from matplotlib import pyplot as plt
import pandas
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

mage_mode = "reduction"
#mage_mode = "released"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
color1 = 'k'     # color for spectra
color2 = '0.65'  # color for uncertainty spectra
color3 = '0.5'   # color for continuum
color4 = 'b'     # color for 2nd spectrum, for comparison
outdir = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Mage/Analysis/Photospheric/"  # where to save plots

# Photospheric features that John Chisholm and JR identified as strong in the stack.
line_label_p  = ("CIII 1247", "C_III 1296", "CII 1323", "OIV 1343", "SiIII 1417", "SiII 1485", "SV 1501",
                 "CIII 1620", "FeV 1662", "FeIV 1717", "FeIII 1930 complex", "FeIII 1956 complex", "CIII 2297")
line_center_p = np.array((1247.38, 1296.33,      1323.93,   1343.514,  1417.24,    1485.40,        1501.76,
                  1620.40,     1662.32,   1717.90,      1930.39,               1953.33, 2297.58))
line_label_renorm = line_center_p > 1600.   # Normalization issues w S99 models for wave>1600.  So, renormalize within the plot window.

S99fits = ('stack-A', 'S0004-0103', 'S0957+0509', 'S2111-0114', 'S0033+0242', 'S1226+2152', 'cosmiceye', 'S0108+0624', 'S1429+1202', 'Horseshoe', 'S0900+2234', 'S1527+0652', 'rcs0327-knotE', 'rcs0327-knotG', 'rcs0327-knotU')

sorted_by_age = ('S0033+0242', 'rcs0327-knotE', 'rcs0327-knotG', 'S0108+0624', 'S0957+0509', 'Horseshoe', 'rcs0327-knotU', 'S2111-0114', 'S1429+1202', 'S0004-0103','S0900+2234', 'stack-A', 'S1527+0652', 'S1226+2152') # 'S1458-0023', 


def local_s99_compare_manyspectra(list_S99files, line_cen, line_label, win, label, Ncol=1, vel_plot=False, mage_mode="reduction", LL=[], size=(8,8)) :
    ''' Plot one transition, versus (wavelength or velocity?), for many MagE spectra.  Plot the S99 fit too (a second wave, fnu pair) too
    line_cen:     rest-wavelength of line to plot, in Angstroms.
    line_label:   label for that line
    win:          window around that line.  If vel_plot, then units are km/s.  If not vel_plot, units are rest-frame Angstroms
    vel_plot:     (optional, Boolean): If True, x-axis is velocity.  If False, x-axis is wavelength.
    specs:        (optional) Pandas dataframe of MagE filenames & redshifts.'''
    
    Nspectra = len(list_S99files)
    Nrow = int(np.ceil( (len(list_S99files)*1.0) / Ncol))  # Calculate how many rows to generate
    fig = plt.figure(figsize=size)
    plt.suptitle(label, fontsize=18)

    for ii, rootname in enumerate(list_S99files) :
        print rootname,
        df = jrr.mage.open_S99_spectrum(rootname, 0.0)
        ax = fig.add_subplot(Nrow, Ncol, ii+1)
        if(vel_plot) : # Do velocity plot first, since already coded up.  Later, can add versus wavelength
            vel = jrr.spec.convert_restwave_to_velocity(df.wave, line_cen)   # velocity in km/s
            in_window = vel.between(-1*win, win)
            plt.step(vel[in_window], df.data_fnu[in_window], color=color1)
            plt.step(vel[in_window], df.data_fnu_u[in_window], color=color2)
            plt.step(vel[in_window], df.s99_fnu[in_window], color=color4)
            plt.plot( (0., 0.), (0.0,2), color=color3, linewidth=2)  # plot tics at zero velocity
            plt.xlim(-1*win, win)
        else :
            plt.annotate(rootname, (0.5,0.8), xycoords="axes fraction")
            in_window = df.wave.between(line_cen - win, line_cen + win)
            if(line_cen > 1600.0) :
                plt.step(df.wave[in_window], df.data_fnu[in_window]   /np.median(df.data_fnu[in_window]), color=color1)
                plt.step(df.wave[in_window], df.data_fnu_u[in_window]  /np.median(df.data_fnu[in_window]), color=color2)
                plt.step(df.wave[in_window], df.s99_fnu[in_window]/np.median(df.s99_fnu[in_window]), color=color4)
            else :
                plt.step(df.wave[in_window], df.data_fnu[in_window]   /np.median(df.data_fnu[in_window]), color=color1)
                plt.step(df.wave[in_window], df.data_fnu_u[in_window]  /np.median(df.data_fnu[in_window]), color=color2)
                plt.step(df.wave[in_window], df.s99_fnu[in_window]/np.median(df.data_fnu[in_window]), color=color4)
            plt.plot( (line_cen, line_cen), (0.0,2), color=color3, linewidth=2)  # plot tics at zero velocity
            plt.xlim(line_cen - win, line_cen + win)
        plt.ylim(0.0, 1.5)  # May need to change these limits
        jrr.mage.plot_linelist(LL, 0.0, True, False)  # plot the line IDs
        
    if vel_plot :
        plt.xlabel("rest-frame velocity (km/s)")  
    else :
        plt.xlabel(r'rest-frame wavelength($\AA$)')                                
#        fig.(hspace=0)
    return(0)

            
# housekeeping
zz = 0.0  # stacked spectrum is already in rest frame wavelength
xwin = 12. # +- velocity window (km/s) to consider a line
Ncol = 4

# Open stacked spectrum and associated linelist
linelist = line_path + "stacked.linelist"
(LL, z_sys) = jrr.mage.get_linelist(linelist)
(sp) = jrr.mage.open_stacked_spectrum(mage_mode)


# What I need to do next.
# All this below was written in a frenzy, and I didn't make it re-usable.  Also, 1 of the 2 functions above
# is superceded by mage.boxplot_Nspectra, So I should switch to that.  The one difference is that above
# does a local renormalization.  Add that functionality, either as an option in box_plot, (median of local window), or from autocont?
# what to do about the second routine?  It's useful, but badly hardcoded?  use as is and move on?
 
plot_one_page_per_spectrum = False
if plot_one_page_per_spectrum :
    # This should generate same output as first run of loop below.  Check.
    # Grab S99 fit to MagE Stack A
    the_pdf = "S99-photospheric-bygalaxy.pdf"
    pp = PdfPages(the_pdf)  # output
    S99 = jrr.mage.open_S99_spectrum("stack-A", 0.0)
    alt_file = "magestack_byneb_ChisholmstackA_spectrum.txt"  # this is the spectrum that JChisholm fit
    altsp = jrr.mage.open_stacked_spectrum(mage_mode, alt_file)
    print "STATUS: Making an experimental plot of the photospheric absorption lines"
    jrr.mage.boxplot_Nspectra((altsp.wave,S99.wave), (altsp.fnu,S99.fnu,), (altsp.fnu_u, S99.fnu_u), (0.0,0.0), line_label_p, line_center_p,  xwin, Ncol, LL, extra_label="StackA+S99 for sanity-checking", figsize=(16,16), vel_plot=False, plot_xaxis=True, ylims=(0.0,1.5))  # renorm=line_label_renorm,
    pp.savefig()

    for rootname in S99fits :
        df = jrr.mage.open_S99_spectrum(rootname)
        local_boxplot_2spectra(df.wave, df.data_fnu, df.data_fnu_u, df.wave, df.s99_fnu, S99.s99_fnu*-0.01, line_label_p, line_center_p, redshift, redshift, xwin, Ncol, LL, rootname, (16,16), vel_plot=False, plot_xaxis=True)
        pp.savefig()    
    pp.close()

    # Repeat, but sorted by age
    the_pdf = "MageES99/S99-photospheric-bygalaxy-sortbyage.pdf"
    pp = PdfPages(the_pdf)  # output
    for rootname in sorted_by_age :
        df = jrr.mage.open_S99_spectrum(rootname)
        local_boxplot_2spectra(df.wave, df.data_fnu, df.data_fnu_u, df.wave, df.s99_fnu, S99.s99_fnu*-0.01, line_label_p, line_center_p, redshift, redshift, xwin, Ncol, LL, rootname, (16,16), False, True)
        pp.savefig()    
    pp.close()

    
    
plot_one_page_per_line = True
if plot_one_page_per_line :
    the_pdf = "MageES99/S99-photospheric-bylines.pdf"
    pp = PdfPages(the_pdf)  # output
    for ii in range(0, len(line_center_p)) :
        print "\n", line_label_p[ii], line_center_p[ii],
        local_s99_compare_manyspectra(S99fits, line_center_p[ii], line_label_p[ii], xwin, line_label_p[ii], Ncol, size=(16,8), LL=LL)
        pp.savefig()
    pp.close()

if plot_one_page_per_line :
    the_pdf = "MageES99/S99-photospheric-bylines-sortbyage.pdf"
    pp = PdfPages(the_pdf)  # output
    for ii in range(0, len(line_center_p)) :
        print "\n", line_label_p[ii], line_center_p[ii],
        local_s99_compare_manyspectra(sorted_by_age, line_center_p[ii], line_label_p[ii], xwin, line_label_p[ii], Ncol, size=(16,8), LL=LL)
        pp.savefig()
    pp.close()
    

