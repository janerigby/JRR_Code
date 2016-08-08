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
S99dir = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/MagE_atlas/Contrib/S99/"

# Photospheric features that John Chisholm and JR identified as strong in the stack.
line_label_p  = ("CIII 1247", "C_III 1296", "CII 1323", "OIV 1343", "SiIII 1417", "SiII 1485", "SV 1501",
                 "CIII 1620", "FeV 1662", "FeIV 1717", "FeIII 1930 complex", "FeIII 1956 complex", "CIII 2297")
line_center_p = np.array((1247.38, 1296.33,      1323.93,   1343.514,  1417.24,    1485.40,        1501.76,
                  1620.40,     1662.32,   1717.90,      1930.39,               1953.33, 2297.58))

S99fits = ('stack-A', 'S0004-0103', 'S0957+0509', 'S2111-0114', 'S0033+0242', 'S1226+2152', 'cosmiceye', 'S0108+0624', 'S1429+1202', 'Horseshoe', 'S0900+2234', 'S1527+0652', 'rcs0327-knotE', 'rcs0327-knotG', 'rcs0327-knotU')

sorted_by_age = ('S0033+0242', 'rcs0327-knotE', 'rcs0327-knotG', 'S0108+0624', 'S0957+0509', 'Horseshoe', 'rcs0327-knotU', 'S2111-0114', 'S1429+1202', 'S0004-0103','S0900+2234', 'stack-A', 'S1527+0652', 'S1226+2152') # 'S1458-0023', 

def local_boxplot_2spectra(wave1, fnu1, dfnu1, wave2, fnu2, dfnu2, line_label, line_center, redshift1, redshift2, win, Ncol, LL, spec_label="",figsize=(8,8), vel_plot=True, plot_xaxis=True) :
    ''' Same as mage_boxplot_spectra, but overplot two different spectra, for comparison.'''
    Nrow = int(np.ceil( (len(line_center)*1.0) / Ncol))  # Calculate how many rows to generate
    print "Now plotting", spec_label
    fig = plt.figure(figsize=figsize)
    restwave1 = wave1 / (1.0 + redshift1)
    restwave2 = wave2 / (1.0 + redshift2)
    plt.suptitle(spec_label, fontsize=18)

    for ii, dum in enumerate(line_label) :
        print "    Plotting ", line_label[ii], " at ", line_center[ii]
        ax = fig.add_subplot(Nrow, Ncol, ii+1)
        plt.annotate( line_label[ii], (0.3,0.9), xycoords="axes fraction")
        if(vel_plot) :
            vel1 = jrr.spec.convert_restwave_to_velocity(restwave1, line_center[ii])   # velocity in km/s
            vel2 = jrr.spec.convert_restwave_to_velocity(restwave2, line_center[ii])   # velocity in km/s
            in_window1 = vel1.between(-1*win, win)
            in_window2 = vel2.between(-1*win, win)
            plt.step(vel1[in_window1], fnu1[in_window1], color=color1)   # assumes input is continuum-normalized
            plt.step(vel1[in_window1], dfnu1[in_window1], color=color2)  # plot uncertainty
            plt.step(vel2[in_window2], fnu2[in_window2], color=color4)   # assumes input is continuum-normalized
            plt.step(vel2[in_window2], dfnu2[in_window2], color=color2)  # plot uncertainty
            #plt.plot( (-1*win, win), (1.0,1.0), color=color3)        # plot unity continuum. 
            plt.plot( (0., 0.), (0.0,2), color=color2, linewidth=2)  # plot tics at zero velocity
            plt.ylim(0.0, 1.5)  # May need to change these limits
            plt.xlim(-1*win, win)
        else :
            in_window1 = restwave1.between((line_center[ii] - win), (line_center[ii] + win))
            in_window2 = restwave2.between((line_center[ii] - win), (line_center[ii] + win))
            if(line_center[ii] > 1600.0) :
                # if >1600A, S99 fit is suspect, so normalize both to their median
                plt.step(restwave1[in_window1], fnu1[in_window1]/np.median(fnu1[in_window1]),  color=color1)
                plt.step(restwave1[in_window1], dfnu1[in_window1]/np.median(fnu1[in_window1]),  color=color2)
                plt.step(restwave2[in_window2], fnu2[in_window2]/np.median(fnu2[in_window2]),  color=color4)
            else:  # <1600A so trust the relative fluxing, but scale them both by the same amt to fit onto the plot
                plt.step(restwave1[in_window1], fnu1[in_window1]/np.median(fnu1[in_window1]),  color=color1)
                plt.step(restwave1[in_window1], dfnu1[in_window1]/np.median(fnu1[in_window1]),  color=color2)
                plt.step(restwave2[in_window2], fnu2[in_window2]/np.median(fnu1[in_window1]),  color=color4)                
            #plt.plot((line_center[ii] - win, line_center[ii] + win), (1.0,1.0), color=color3)        # plot unity continuum. 
            plt.plot( (line_center[ii], line_center[ii]), (0.0,2), color=color3, linewidth=2)  # plot tics at zero velocity
            plt.xlim(line_center[ii] - win, line_center[ii] + win)
            plt.ylim(0.0, 1.4)  # May need to change these limits
            
        jrr.mage.plot_linelist(LL, redshift1, True, vel_plot, line_center[ii])  # plot the line IDs
     # plot_linelist(L_all, z_systemic=0.0, restframe=False, velplot=False, line_center=0.0) :  
        if ii == len(line_label) -1 :
            if vel_plot : 
                plt.xlabel("rest-frame velocity (km/s)")  # if last subplot, make xlabel
            else :
                plt.xlabel(r'rest-frame wavelength($\AA$)')
        if (not plot_xaxis) and (ii < len(line_label)-1) :
            ax.axes.xaxis.set_ticklabels([])  # if not last subplot, suppress  numbers on x axis
            # But warning, this will disable x= on interactive matplotlib.  comment out above line to measure numbers interactively on graph
        if not plot_xaxis :
            fig.subplots_adjust(hspace=0)

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
redshift = 0.00  # stacked spectrum is already in rest frame wavelength
xwin = 12. # +- velocity window (km/s) to consider a line
Ncol = 4

# Open stacked spectrum and associated linelist
linelist = line_path + "stacked.linelist"
(LL, z_sys) = jrr.mage.get_linelist(linelist)
(sp) = jrr.mage.open_stacked_spectrum(mage_mode)

plot_one_page_per_spectrum = False
if plot_one_page_per_spectrum :
    # This should generate same output as first run of loop below.  Check.
    # Grab S99 fit to MagE Stack A
    the_pdf = "MageES99/S99-photospheric-bygalaxy.pdf"
    pp = PdfPages(the_pdf)  # output
    S99file = S99dir + "stack-A-sb99-fit.txt"
    S99 = jrr.mage.open_S99_spectrum("stack-A")
    alt_file = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Mage/Analysis/Stacked_spectra/Stack_by_zneb/mage_stack_ChisholmstackA_spectrum.txt"  # this is the spectrum that JChisholm fit
    altsp = jrr.mage.open_stacked_spectrum(mage_mode, alt_file)
    print "STATUS: Making an experimental plot of the photospheric absorption lines"
    local_boxplot_2spectra(altsp.restwave, altsp.X_avg, altsp.X_sigma, S99.wave, S99.s99fit, S99.s99fit*-0.01, line_label_p, line_center_p, redshift, redshift, xwin, Ncol, LL, "StackA+S99 for sanity-checking", (16,16), False, True)
    pp.savefig()    

    for rootname in S99fits :
        df = jrr.mage.open_S99_spectrum(rootname)
        local_boxplot_2spectra(df.wave, df.data_fnu, df.data_fnu_u, df.wave, df.s99_fnu, S99.s99_fnu*-0.01, line_label_p, line_center_p, redshift, redshift, xwin, Ncol, LL, rootname, (16,16), False, True)
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
    

