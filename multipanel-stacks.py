''' Same as multipanel-spectrum-ids.py, but plot multiple spectra.  Using for sub-stacks (highZ, loZ),
young, middle-aged, old.
jrigby, Sept 2016. '''

import jrr
import matplotlib.pyplot as plt
import numpy as np
mage_mode = "reduction"  
mage_mode = "released"
methods = ('bystars', 'byneb')  # method of determining systemic redshift

# For sanity checking, make plots of the weighted avg(std), but also the median and sqrt(jackknife variance)
colfnu  = ('X_avg', 'X_median')
colfnuu = ('X_sigma', 'X_jack_std')
avg_or_med = ('wtdavg', 'median') 
for ii in range(0, len(colfnu)):
    # Compare high and low metallicity stacks
    for method in methods :
        stacks = ['magestack_'+method+'_highZ_spectrum.txt', 'magestack_'+method+'_lowZ_spectrum.txt']
        sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0], colfnu=colfnu[ii], colfnuu=colfnuu[ii])
        sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1], colfnu=colfnu[ii], colfnuu=colfnuu[ii])
        the_dfs = [sp1, sp2]
        the_zzs = [0.0, 0.0]  # Stacked spectra are already in rest frame.
        the_pdf =  'PDF_Out2/multipanel_stacked_'+avg_or_med[ii]+'_'+method+'_byZ.pdf'
        (spec_path, line_path) = jrr.mage.getpath(mage_mode)
        (LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, colortab=('black', 'mediumblue'), title="Compare highZ (black) and lowZ (blue) stacks "+method)
    plt.clf()

    # Compare young, middle-age, and old stacks
    label_age = "Compare young (blue), middle-age (black), and old (red) stacks"
    colortab_age=('blue', 'black', 'red')
    for method in methods:
        stacks = ['magestack_'+method+'_younglt8Myr_spectrum.txt', 'magestack_'+method+'_midage8to16Myr_spectrum.txt', 'magestack_'+method+'_oldgt16Myr_spectrum.txt']
        sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0], colfnu=colfnu[ii], colfnuu=colfnuu[ii])
        sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1], colfnu=colfnu[ii], colfnuu=colfnuu[ii])
        sp3 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[2], colfnu=colfnu[ii], colfnuu=colfnuu[ii])
        the_dfs = [sp1, sp2, sp3]
        the_zzs = [0.0, 0.0, 0.0]  # Stacked spectra are already in rest frame.
        the_pdf =  'PDF_Out2/multipanel_stacked_'+avg_or_med[ii]+'_'+method+'_byage.pdf'
        (spec_path, line_path) = jrr.mage.getpath(mage_mode)
        (LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title=label_age + " "+method, colortab=colortab_age, waverange=(1000,3000))
    plt.clf()

    # Plot the standard stacks
    for method in methods:
        stacks = ['magestack_'+method+'_standard_spectrum.txt']
        sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0], colfnu=colfnu[ii], colfnuu=colfnuu[ii])
        the_dfs = [sp1]
        the_zzs = [0.0]  # Stacked spectra are already in rest frame.
        the_pdf =  'PDF_Out2/multipanel_stacked_'+avg_or_med[ii]+'_'+method+'_standard.pdf'
        (spec_path, line_path) = jrr.mage.getpath(mage_mode)
        (LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title="Standard stack"+method, waverange=(1000,3000))
        # make plain stack to drop into paper
        # figure out new outfile for this, so it doesn't overwrite
        the_pdf =  'PDF_Out2/multipanel_stacked_'+avg_or_med[ii]+'_'+method+'_standard_nolabels.pdf'
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, waverange=(1000,3000), annotate=())
    plt.clf()




