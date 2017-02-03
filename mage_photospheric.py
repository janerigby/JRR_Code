''' Script to make the photospheric absorption plots that I can't eaisly make from within mage_winds.py
jrigby, dec 2016.  Kinda kludgy.  Run from /Mage/Analysis/Wind_plots/Photospheric'''

import jrr
import os
import sys
import re
from matplotlib import pyplot as plt
import pandas
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

#mage_mode = "reduction"
mage_mode = "released"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
color1 = 'k'     # color for spectra
color2 = '0.65'  # color for uncertainty spectra
color3 = '0.5'   # color for continuum
color4 = 'blue'     # color for 2nd spectrum, for comparison

# Photospheric features that John Chisholm and JR identified as strong in the stack.
line_label_p  = ("CIII 1247", "C_III 1296", "CII 1323", "OIV 1343", "SiIII 1417", "SiII 1485", "SV 1501",
                 "CIII 1620", "FeV 1662", "FeIV 1717", "FeIII 1930 complex", "FeIII 1956 complex", "CIII 2297")
line_center_p = np.array((1247.38, 1296.33,      1323.93,   1343.514,  1417.24,    1485.40,        1501.76,
                  1620.40,     1662.32,   1717.90,      1930.39,               1953.33, 2297.58))

S99fits = ('S0004-0103', 'S0033+0242', 'S0108+0624', 'rcs0327-E', 'rcs0327-G', 'rcs0327-U',
'S0900+2234','S0957+0509', 'Horseshoe', 'S1226+2152','S1429+1202', 'S1527+0652', 'S2111-0114', 'Cosmic~Eye','Stack-A', 'chuck')

sortedbyage = ('S2111-0114', 'rcs0327-E' , 'S1458-0023', 'Cosmic~Eye', 'Horseshoe', 'S0108+0624', 'rcs0327-G', 'S1527+0652', 'rcs0327-U', 'S0004-0103', 'S0957+0509', 'chuck', 'S0033+0242', 'Stack-A', 'S1429+1202', 'S0900+2234', 'S1226+2152')
# sorted by new ages, 9 dec 2016

S99_sortbyage = ( 'S2111-0114', 'rcs0327-E', 'Cosmic~Eye',  'Horseshoe', 'rcs0327-G', 'S1527+0652', 'rcs0327-U', 'S0004-0103','S0957+0509', 'chuck', 'S0033+0242', 'S0108+0624', 'Stack-A', 'S1429+1202', 'S0900+2234', 'S1226+2152')
    
# Picking out the higher SNR subset from those w S99 fits
too_noisy_for_photospheric = ('S0957+0509','Horseshoe', 'S2111-0114')
gallist_justhighSNR = ('S0004-0103', 'S0033+0242', 'S0108+0624', 'rcs0327-E', 'rcs0327-G', 'rcs0327-U', 'S0900+2234','S1226+2152','S1429+1202', 'S1527+0652', 'Cosmic~Eye','Stack-A', 'chuck')

sortedbyage_justhighSNR = ( 'rcs0327-E', 'Cosmic~Eye',  'S0108+0624', 'rcs0327-G', 'S1527+0652', 'rcs0327-U', 'S0004-0103', 'chuck', 'S0033+0242', 'Stack-A', 'S1429+1202','S0900+2234','S1226+2152')

print "STATUS:  Loading all the spectra (and their S99 fits) into big honking dataframes."
(df, resoln, dresoln, big_LL, big_zz_sys, specs) = jrr.mage.open_many_spectra(mage_mode, which_list='wcont')  # open honking
#(df, resoln, dresoln, big_LL, big_zz_sys, specs) = jrr.mage.open_many_spectra(mage_mode, which_list='labels', labels=('rcs0327-E','rcs0327-B'))  # DEBUGGING TESTING
df['chuck']  = jrr.mage.read_chuck_UVspec(addS99=True, autofitcont=True)
(df['Stack-A'], big_LL['Stack-A']) = jrr.mage.open_stacked_spectrum(mage_mode, which_stack="Stack-A", addS99=True)
big_LL['chuck'] = big_LL['Stack-A']    # Kludge
big_zz_sys['chuck'] = 0.0  ;  big_zz_sys['Stack-A'] = 0.0
gallist_to_process = [val for val in specs['short_label']]
gallist_to_process.append("chuck")
gallist_to_process.append("Stack-A")
print "DEBUGGING", gallist_to_process


 
def local_s99_compare_manyspectra(labels, line_cen, line_label, win, label, Ncol=1, vel_plot=False, mage_mode="reduction", size=(8,8)) :
    ''' Plot one transition for many MagE galaxies on a page, one page per transition.  Compare spectra and S99 fit'''
    Nrow = int(np.ceil( (len(labels)*1.0) / Ncol))  # Calculate how many rows to generate
    fig = plt.figure(figsize=size)
    plt.suptitle(label, fontsize=18)
    for ii, label in enumerate(labels) :
        print "DEBUG", label, ii
        sp = df[label]
        ax = fig.add_subplot(Nrow, Ncol, ii+1)
        if(vel_plot) :   # x axis is velocity (km/s)
            vel = jrr.spec.convert_restwave_to_velocity(sp.rest_wave, line_cen)   # velocity in km/s
            in_win = vel.between(-1*win, win)
            plt.step(vel[in_win], sp.rest_fnu[in_win]          / sp.rest_fnu_autocont[in_win], color=color1)
            plt.step(vel[in_win], sp.rest_fnu_u[in_win]        / sp.rest_fnu_autocont[in_win], color=color2)
            if jrr.mage.getfullname_S99_spectrum(label) :
                plt.step(vel[in_win], sp.rest_fnu_s99model[in_win] / sp.rest_fnu_autocont[in_win], color=color4)
            plt.plot( (0., 0.), (0.0,2), color=color3, linewidth=2)  # plot tics at zero velocity
            plt.xlim(-1*win, win)
        else :           # x axis is wavelength
            if label == "chuck"     : pretty_label = "Steidel et al. 2016"
            elif label == "Stack-A" : pretty_label = r'$\lambda_{pivot}$-norm. stack'
            else                    : pretty_label = label
            plt.annotate(pretty_label, (0.52,0.85), xycoords="axes fraction")
            in_win = sp.rest_wave.between(line_cen - win, line_cen + win)
            plt.step(sp.rest_wave[in_win], sp.rest_fnu[in_win]          / sp.rest_fnu_autocont[in_win], color=color1)
            plt.step(sp.rest_wave[in_win], sp.rest_fnu_u[in_win]        / sp.rest_fnu_autocont[in_win], color=color2)
            if jrr.mage.getfullname_S99_spectrum(label) :
                plt.step(sp.rest_wave[in_win], sp.rest_fnu_s99model[in_win] / sp.rest_fnu_autocont[in_win], color=color4)
            plt.plot( (line_cen, line_cen), (0.0,2), color=color3, linewidth=2)  # plot tics at zero velocity
            plt.xlim(line_cen - win, line_cen + win)
            plt.locator_params(axis='x', nbins=4)
            plt.locator_params(axis='y', nbins=3)
        plt.ylim(0.0, 1.5)  # May need to change these limits
        jrr.mage.plot_linelist(big_LL[label], big_zz_sys[label], True, False)  # plot the line IDs
        if ii == len(labels)-1 or ii == len(labels)-2 :
            if vel_plot :  plt.xlabel("rest-frame velocity (km/s)")  
            else :         plt.xlabel(r'rest-frame wavelength ($\rm \AA$)')                                
    fig.subplots_adjust(hspace=0)
    return(0)
            
# housekeeping
zz = 0.0  # stacked spectrum is already in rest frame wavelength
xwin = 12. # +- velocity window (km/s) to consider a line
Ncol = 4


def line_per_page(whichgals, the_pdf):
    pp = PdfPages(the_pdf)  # output
    for ii in range(0, len(line_center_p)) :
        print "\n", line_label_p[ii], line_center_p[ii], "***",
        local_s99_compare_manyspectra(whichgals, line_center_p[ii], line_label_p[ii], xwin, line_label_p[ii], Ncol, size=(16,8))
        pp.savefig()
    pp.close()
    plt.clf()
    return(0)

def plot_all_winds(whichgals, the_pdf, size) :
    lab_c     =    ('Al II 1670', 'C II 1334', 'Al III 1854', 'Si IV 1393', 'C IV 1548', 'N V 1238')
    cen_c = np.array(( 1670.7874,   1334.5323,   1854.72,      1393.76,      1548.19,    1238.82))            
    pp = PdfPages(the_pdf)
    for ii, thisone in enumerate(cen_c) :
        local_s99_compare_manyspectra(whichgals,  cen_c[ii], lab_c[ii], 20, lab_c[ii], Ncol=2, vel_plot=False, mage_mode="reduction", size=(10,13))
        pp.savefig()
    pp.close()
    plt.clf()
    
###############################
# Actually run things
print "STATUS:  Plotting the photospheric lines, one page per transition."
line_per_page(gallist_to_process,  "S99-photospheric-bylines.pdf")
line_per_page(gallist_justhighSNR, "S99-photospheric-bylines_highSNR.pdf")

print "STATUS:  Plotting the photospheric lines, one page per transition, now sorting by light-weighted age."
line_per_page(sortedbyage, "S99-photospheric-bylines-sortbyage.pdf")
line_per_page(sortedbyage_justhighSNR, "S99-photospheric-bylines-sortbyage_highSNR.pdf")

# Make 1 page of CIV plots, with S99:
plot_all_winds(S99fits[0:16], 'windlines_all_wSteidelStack.pdf', size=(10,13))
plot_all_winds(S99fits[0:14], 'windlines_all.pdf', size=(10,12))
plot_all_winds(S99_sortbyage[0:16], 'windlines_sortbyage_wSteidelStack.pdf', size=(10,13))

###############################
