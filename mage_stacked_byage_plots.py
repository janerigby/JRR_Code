'''Plotting interesting features in sub-stacks (highZ, loZ),young, middle-aged, old.
jrigby, Oct 2016. '''

import jrr
import matplotlib.pyplot as plt
import numpy as np
mage_mode = "released"  
methods = ('bystars', 'byneb')  # method of determining systemic redshift

line_label_a  = ["Lya", "He II 1640", "C III] 1906", "C II] 2326", "[O II] 2470", "Fe II 2365"]#, "Fe II 2396"]
line_center_a = [1215.6701,  1640.417, 1906.683, 2326.113, 2471.027, 2365.552]#, 2396.1497]

line_label_b  = ["C III 1334", "Si II 1526", "Fe II 1608", "Fe II 2382"]
line_center_b = [1334.5323,    1526.7066,     1608.4511,    2382.765]

#label_age = "Compare young (blue), middle-age (black), and old (red) stacks"
label_age = ""
colortab_age=('blue', 'black', 'red')
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
ymax_a = [20., 1.5, 2., 1.5, 1.5, 1.5]
ymax_b = [1.4, 1.4, 1.4, 1.4]
kinds = ['X_avg', 'X_median']

for kind in kinds :
    for method in methods:
        stacks = ['magestack_'+method+'_younglt8Myr_spectrum.txt', 'magestack_'+method+'_midage8to16Myr_spectrum.txt', 'magestack_'+method+'_oldgt16Myr_spectrum.txt']
        sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
        sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
        sp3 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[2])

        # Plot the emission lines for young, middle, old, to show they get stronger w youth
        thewaves = (sp1.rest_wave, sp2.rest_wave, sp3.rest_wave)
        thefnus  = (sp1[kind]/sp1.fnu_autocont, sp2[kind]/sp2.fnu_autocont, sp3[kind]/sp3.fnu_autocont)
        thedfnus = (sp1.X_sigma/sp1.fnu_autocont, sp2.X_sigma/sp2.fnu_autocont, sp3.X_sigma/sp3.fnu_autocont)
        thezs = (0., 0., 0.)    
        the_pdf =  'Boxplots_for_stacks/boxplot_'+method+'_'+kind+'_byage_emission.pdf'
        jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_a, line_center_a, win=14., Ncol=1, LL=(), extra_label=label_age, figsize=(4,12), vel_plot=False, ymax=ymax_a, colortab=colortab_age)
        plt.savefig(the_pdf, orientation='portrait', bbox_inches='tight', pad_inches=0.1)
        plt.close()
        the_pdf = "Boxplots_for_stacks/boxplot_"+method+'_'+kind+'_byage_abs.pdf'
        jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_b, line_center_b, win=2000., Ncol=1, LL=(), extra_label=label_age, figsize=(4,12), vel_plot=True, plot_xaxis=False,  ymax=ymax_b, colortab=colortab_age)
        plt.savefig(the_pdf, orientation='portrait', bbox_inches='tight', pad_inches=0.1)
        plt.close()

    #repeat, for metallicity
    for method in methods:
        stacks = ['magestack_'+method+'_lowZ_spectrum.txt', 'magestack_'+method+'_highZ_spectrum.txt']
        sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
        sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
        thewaves = (sp1.rest_wave, sp2.rest_wave)
        thefnus  = (sp1[kind]/sp1.fnu_autocont, sp2[kind]/sp2.fnu_autocont)
        thedfnus = (sp1.X_sigma/sp1.fnu_autocont, sp2.X_sigma/sp2.fnu_autocont)
        thezs = (0., 0.)    
        the_pdf =  'Boxplots_for_stacks/boxplot_'+method+'_'+kind+'_byZ_emission.pdf'
        jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_a, line_center_a, win=14., Ncol=1, LL=(), extra_label=label_age, figsize=(4,12), vel_plot=False, ymax=ymax_a, colortab=('blue', 'black'))
        plt.savefig(the_pdf, orientation='portrait', bbox_inches='tight', pad_inches=0.1)
        plt.close()
        the_pdf = 'Boxplots_for_stacks/boxplot_'+method+'_'+kind+'_byZ_abs.pdf'
        jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_b, line_center_b, win=2000., Ncol=1, LL=(), extra_label=label_age, figsize=(4,12), vel_plot=True, plot_xaxis=False,  ymax=ymax_b, colortab=colortab_age)
        plt.savefig(the_pdf, orientation='portrait', bbox_inches='tight', pad_inches=0.1)
        plt.close()
