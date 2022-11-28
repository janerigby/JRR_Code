from __future__ import print_function
import jrr
import numpy as np
import pandas
import matplotlib 
import  matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator, ScalarFormatter
import re

def calc_smoothSNR(sp, resoln, smooth_window, colwave='wave', colf='fnu', colfu='fnu_u', coldisp='disp') :
    # For a mage spectrum, calculate the SNR *per resoln element*, and smooth it.
    sp['SNRperres'] = sp[colf]/sp[colfu] * np.sqrt(sp[colwave] / resoln / sp[coldisp])
    sp['SNRperres'].interpolate(inplace=True) # remove NaNs, they'll screw up the smoothing
    sp['smoothSNRperres'] = sp['SNRperres'].rolling(window=smooth_window, center=True).median()
    return(0)


mage_mode = 'reduction'
matplotlib.rcParams.update({'font.size': 14})
smooth_window = 201 #301
figsize = (8,4)
#pp = PdfPages("snr_1Feb2018_all.pdf")
pp = PdfPages("snr_Feb2019_all.pdf")

plot_indy_SNR = True
if plot_indy_SNR :
    lab1 = [ 'rcs0327-E', 'rcs0327-G', 'rcs0327-U', 'rcs0327-B', 'rcs0327-counterarc',]
    lab2 = [ 'S1527+0652', 'S1527+0652-fnt']
    lab3 = [ 'S0004-0103', 'S0900+2234', 'S1226+2152']
    lab4 = [ 'Cosmic~Eye', 'S0108+0624',  'S1429+1202']
    lab5 = [ 'S0033+0242', 'S2111-0114', 'S2243-0935' , 'S0957+0509']
    lab6 = ['S1458-0023',  'Horseshoe',  'S1050+0017']
    oldlab7 = [ 'planckarc_h9', 'planckarc_pos1', 'planckarc_f', 'planckarc_h4', 'planckarc_h1a',  'planckarc_h3', 'planckarc_h5', 'planckarc_h1', 'planckarc_h2']#, 'planckarc_slit4a', 'planckarc_slit4bc']
    lab7 = jrr.mage.replace_all_oldstarburstnames_list(oldlab7) # Upgrade to M-X names
    
    lab8 = ['SPT2325', 'SPT0310_slitA', 'SPT0310_slitB', 'PSZ0441_slitA', 'PSZ0441_slitB', 'SPT0142', 'SPT0356', 'S1226+2152image3']
    labs = (lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)
    for labels in labs:
        (df, resoln, dresoln, LL, zz_sys, speclist) = jrr.mage.open_many_spectra(mage_mode, which_list="labels", labels=labels, verbose=True, zchoice='stars', addS99=False, MWdr=True)
        #fig = plt.figure(figsize=figsize)
        fig, ax  = plt.subplots(figsize=figsize)
        for label in labels:
            sp = df[label]
            print(label, resoln[label])
            calc_smoothSNR(sp, resoln[label], smooth_window)
            labelnew = jrr.mage.prettylabel_from_shortlabel(label)
            plt.plot(sp['wave'], sp['smoothSNRperres'], label=labelnew)
        plt.legend(fontsize=10, frameon=True, labelspacing=0, loc='upper left')
        #plt.ylim(0,47) Used for Rigby+2018
        plt.ylim(0,30) 
        plt.xlim(3200, 8250)
        plt.xlabel(r'observed wavelength ($\mathrm{\AA}$)')
        plt.ylabel('SNR per resoln. element')
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        #plt.grid()
        pp.savefig(bbox_inches='tight', pad_inches=0.1)

plot_stack_SNR = True
if plot_stack_SNR :
    fig, ax  = plt.subplots(figsize=figsize)
    (sp, LL) = jrr.mage.open_stacked_spectrum(mage_mode)
    RR = 3300
    calc_smoothSNR(sp, RR, 101, colwave='rest_wave')
    sp['pixperfwhm'] = sp['rest_wave'] / RR / 0.1
    plt.plot(sp['rest_wave'], sp['smoothSNRperres'], color='b')
    #plt.xlim(1400,1500)
    lo = np.array((1440., 1460., 1680., 1760.))
    hi = np.array((1450., 1470., 1700., 1800.))
    print("SNR_method1, SNR_method2, wavelo, wavehi")
    for ii, low in enumerate(lo) :
        subset =  sp.loc[sp['rest_wave'].between(lo[ii],hi[ii])]
        over_resolved = subset['pixperfwhm'].median()
        SNR_med_pp = (subset['fnu']/subset['fnu_u']).median()
        #print SNR_med_pp, SNR_med_pp * np.sqrt( over_resolved ), over_resolved
        sp['temp'] = (sp['fnu']/sp['fnu_u']).rolling(window=int(np.floor(over_resolved)), center=True).apply(jrr.util.add_in_quad)
        SNR_method1 = np.round(SNR_med_pp * np.sqrt( over_resolved ), 0)
        SNR_method2 = np.round((sp['temp'].loc[sp['rest_wave'].between(lo[ii],hi[ii])]).median(), 0)
        print("", SNR_method1, SNR_method2, lo[ii], hi[ii])
    plt.title("stacked spectrum")
    plt.xlabel(r'rest-frame wavelength ($\mathrm{\AA}$)')
    plt.ylabel('SNR per resoln. element')
    pp.savefig(bbox_inches='tight', pad_inches=0.1)
pp.close()
plt.close("all")

#plot_indy_spectra = True
#if plot_indy_spectra:   # Now, make plots like those made by Analysis/Spectra_thumbnails/spectra.gnu, but in python

figsize = (20,3)
labels_RAorder = ['rcs0327-E', 'rcs0327-U', 'rcs0327-B', 'rcs0327-G', 'rcs0327-counterarc', 'S1527+0652', 'S1527+0652-fnt', 'S0004-0103',  'S0033+0242', 'S0108+0624', 'S0900+2234', 'S0957+0509', 'S1050+0017', 'Horseshoe',   'S1226+2152', 'S1429+1202', 'S1458-0023', 'S2111-0114',  'Cosmic~Eye', 'S2243-0935']
labels_to_plot = labels_RAorder+lab7+lab8
(df, resoln, dresoln, LL, zz_sys, speclist) = jrr.mage.open_many_spectra(mage_mode, which_list="labels", labels=labels_to_plot, verbose=True, zchoice='neb', addS99=False, MWdr=True)

scalefactor = 1E28
#pp = PdfPages("spectra-snapshots-Feb2018_all.pdf")
pp = PdfPages("spectra-snapshots-Feb2019_all.pdf")
#for label in labels_RAorder :
for label in (labels_to_plot) :
    sp = df[label].loc[~df[label]['badmask']]   # one galaxy, dont plot bad points
    #sp = df[label]
    fig, ax  = plt.subplots(figsize=figsize)
    #print "Plotting", label
    labelnew = jrr.mage.prettylabel_from_shortlabel(label)
    plt.plot(sp['wave'], sp['fnu']   * scalefactor,   label=labelnew, color='k', linewidth=0.5)
    plt.plot(sp['wave'], sp['fnu_u'] * scalefactor, label='_nolegend_', color='lightblue', linewidth=1)
    if re.match("planckarc_pos1", label) or label == "planckarc" or label == 'sunburst_M-0' :    plt.ylim(-0.5E-28*scalefactor,  5E-27*scalefactor)
    elif re.match("planckarc_slit", label)  :    plt.ylim(-0.5E-28*scalefactor,  1E-27*scalefactor)
    else :   plt.ylim(-0.5E-28*scalefactor,  2.5E-28*scalefactor)
    x1 = max(3200, 912*(1.+zz_sys[label]))  # plot down to 3200, unless its blueward of lyman limit
    x2=8200.
    plt.xlim(x1, x2)
    plt.ylabel(r'$f_{\nu}$ ($10^{-28}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
    ax.set_xlabel(r"observed wavelength ($\mathrm{\AA}$)")
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax2 = ax.twiny()
    up_x1 = x1/(1. + zz_sys[label]);    up_x2 =  x2/(1. +  zz_sys[label])
    # try this? ax2.set_xlim([ax.min(), ax.max())
    ax2.set_xlim( up_x1, up_x2)
    print("DEBUGGING", label, zz_sys[label], " : ",  x1, x2, up_x1, up_x2) 
    ax2.set_xlabel(r"rest-frame wavelength ($\mathrm{\AA}$)")
    ax2.xaxis.set_major_locator(MultipleLocator(200))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax2.get_xaxis().set_tick_params(which='both', direction='in')
    plt.annotate(labelnew, (0.45,0.88), xycoords="axes fraction", fontsize=16)
    pp.savefig(bbox_inches='tight', pad_inches=0.1)
pp.close()
plt.close("all")
