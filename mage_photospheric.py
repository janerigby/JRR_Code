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

S99fits = ('stack-A', 'stack-bystars', 'stack-highz', 'stack-lowz', 'stack-young', 'stack-mid', 'stack-old', 'S0004-0103', 'S0957+0509', 'S2111-0114', 'S0033+0242', 'S1226+2152', 'cosmiceye', 'S0108+0624', 'S1429+1202', 'Horseshoe', 'S0900+2234', 'S1527+0652', 'rcs0327-knotE', 'rcs0327-knotG', 'rcs0327-knotU')

sorted_by_age = ('S0033+0242', 'rcs0327-knotE', 'rcs0327-knotG', 'S0108+0624', 'S0957+0509', 'Horseshoe', 'rcs0327-knotU', 'S2111-0114', 'S1429+1202', 'S0004-0103','S0900+2234', 'stack-A', 'S1527+0652', 'S1226+2152') # 'S1458-0023', 


def local_s99_compare_manyspectra(list_S99files, line_cen, line_label, win, label, Ncol=1, vel_plot=False, mage_mode="reduction", size=(8,8)) :
    ''' Plot one transition for many MagE galaxies on a page, one page per transition.  Compare spectra and S99 fit
    Wow, this runs REALLY slow. '''
    Nrow = int(np.ceil( (len(list_S99files)*1.0) / Ncol))  # Calculate how many rows to generate
    fig = plt.figure(figsize=size)
    plt.suptitle(label, fontsize=18)

    for ii, rootname in enumerate(list_S99files) :
        print rootname,
        (df, LL) = jrr.mage.open_S99_spectrum(rootname, 0.0)
        ax = fig.add_subplot(Nrow, Ncol, ii+1)
        if(vel_plot) : # Do velocity plot first, since already coded up.  Later, can add versus wavelength
            vel = jrr.spec.convert_restwave_to_velocity(df.rest_wave, line_cen)   # velocity in km/s
            in_window = vel.between(-1*win, win)
            plt.step(vel[in_window], df.rest_fnu_data[in_window]/df.rest_fnu_data_autocont[in_window], color=color1)
            plt.step(vel[in_window], df.rest_fnu_data_u[in_window]/df.rest_fnu_data_autocont[in_window], color=color2)
            plt.step(vel[in_window], df.rest_fnu_s99[in_window]/df.rest_fnu_s99_autocont[in_window], color=color4)
            plt.plot( (0., 0.), (0.0,2), color=color3, linewidth=2)  # plot tics at zero velocity
            plt.xlim(-1*win, win)
        else :
            plt.annotate(rootname, (0.5,0.85), xycoords="axes fraction")
            in_window = df.rest_wave.between(line_cen - win, line_cen + win)
            plt.step(df.rest_wave[in_window], df.rest_fnu_data[in_window]/df.rest_fnu_data_autocont[in_window], color=color1)
            plt.step(df.rest_wave[in_window], df.rest_fnu_data_u[in_window]/df.rest_fnu_data_autocont[in_window], color=color2)
            plt.step(df.rest_wave[in_window], df.rest_fnu_s99[in_window]/df.rest_fnu_s99_autocont[in_window], color=color4)
            plt.plot( (line_cen, line_cen), (0.0,2), color=color3, linewidth=2)  # plot tics at zero velocity
            plt.xlim(line_cen - win, line_cen + win)
        plt.ylim(0.0, 1.5)  # May need to change these limits
        jrr.mage.plot_linelist(LL, 0.0, True, False)  # plot the line IDs
    if vel_plot :
        plt.xlabel("rest-frame velocity (km/s)")  
    else :
        plt.xlabel(r'rest-frame wavelength ($\rm \AA$)')                                
#        fig.(hspace=0)
    return(0)

            
# housekeeping
zz = 0.0  # stacked spectrum is already in rest frame wavelength
xwin = 12. # +- velocity window (km/s) to consider a line
Ncol = 4
GALAXY_PER_PAGE = False   # False so I can debug the LINE_PER_PAGE part
LINE_PER_PAGE = True

if GALAXY_PER_PAGE :
    print "STATUS:  Plotting the photospheric lines, using one page per galaxy."
    the_pdf = "S99-photospheric-bygalaxy.pdf"
    pp = PdfPages(the_pdf)  # output
    for rootname in S99fits :
        print rootname, "Plotting photospheric lines"
        (df, LL) = jrr.mage.open_S99_spectrum(rootname, 0.0)
        jrr.plot.boxplot_Nspectra((df.rest_wave, df.rest_wave), (df.rest_fnu_data/df.rest_fnu_data_autocont, df.rest_fnu_s99/df.rest_fnu_s99_autocont), (df.rest_fnu_data_u/df.rest_fnu_data_autocont, ()), (0.0,0.0), line_label_p, line_center_p, xwin, Ncol, LL, figsize=(16,8),  vel_plot=False, plot_xaxis=True)
        plt.suptitle(rootname, fontsize=18)
        pp.savefig()    
    pp.close()

if GALAXY_PER_PAGE :
    print "STATUS:  Plotting the photospheric lines, using one page per galaxy, but now sorting by light-weighted age."
    the_pdf = "S99-photospheric-bygalaxy-sortbyage.pdf"
    pp = PdfPages(the_pdf)  # output
    for rootname in sorted_by_age :
        print rootname, "Plotting photospheric lines"
        (df, LL) = jrr.mage.open_S99_spectrum(rootname, 0.0)
        jrr.plot.boxplot_Nspectra((df.rest_wave, df.rest_wave), (df.rest_fnu_data/df.rest_fnu_data_autocont, df.rest_fnu_s99/df.rest_fnu_s99_autocont), (df.rest_fnu_data_u/df.rest_fnu_data_autocont, ()), (0.0,0.0), line_label_p, line_center_p, xwin, Ncol, LL, figsize=(16,8),  vel_plot=False, plot_xaxis=True)
        plt.suptitle(rootname, fontsize=18)
        pp.savefig()    
    pp.close()
    
#if LINE_PER_PAGE :
if False:
    print "STATUS:  Plotting the photospheric lines, one page per transition."
    the_pdf = "S99-photospheric-bylines.pdf"
    pp = PdfPages(the_pdf)  # output
#    for ii in range(0, 2) :  # for debugging
    for ii in range(0, len(line_center_p)) :
        print "\n", line_label_p[ii], line_center_p[ii], "***",
        local_s99_compare_manyspectra(S99fits, line_center_p[ii], line_label_p[ii], xwin, line_label_p[ii], Ncol, size=(16,8))
        pp.savefig()
    pp.close()

if LINE_PER_PAGE :
    print "STATUS:  Plotting the photospheric lines, one page per transition, now sorting by light-weighted age."
    the_pdf = "S99-photospheric-bylines-sortbyage.pdf"
    pp = PdfPages(the_pdf)  # output
    for ii in range(0, len(line_center_p)) :
        print "\n", line_label_p[ii], line_center_p[ii],"***",
        local_s99_compare_manyspectra(sorted_by_age, line_center_p[ii], line_label_p[ii], xwin, line_label_p[ii], Ncol, size=(16,8))
        pp.savefig()
    pp.close()
    
