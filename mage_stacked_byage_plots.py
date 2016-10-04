'''Plotting interesting features in sub-stacks (highZ, loZ),young, middle-aged, old.
jrigby, Oct 2016. '''

import jrr
import matplotlib.pyplot as plt
import numpy as np
mage_mode = "reduction"  
methods = ('bystars', 'byneb')  # method of determining systemic redshift

line_label_a  = ["Lya", "C III] 1906", "He II 1640", "C II] 2323", "[O II] 2470", "Fe II 2365", "Fe II 2396"]
line_center_a = [1215.6701, 1906.683, 1640.417, 2324.214, 2471.027, 2365.552, 2396.1497]

label_age = "Compare young (blue), middle-age (black), and old (red) stacks"
colortab_age=('blue', 'black', 'red')
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.

for method in methods:
    stacks = ['magestack_'+method+'_younglt8Myr_spectrum.txt', 'magestack_'+method+'_midage8to16Myr_spectrum.txt', 'magestack_'+method+'_oldgt16Myr_spectrum.txt']
    sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
    sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
    sp3 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[2])

    # Plot the emission lines for young, middle, old, to show they get stronger w youth
    thewaves = (sp1.rest_wave, sp2.rest_wave, sp3.rest_wave)
    thefnus  = (sp1.X_avg/sp1.fnu_autocont, sp2.X_avg/sp2.fnu_autocont, sp3.X_avg/sp3.fnu_autocont)
    thedfnus = (sp1.X_sigma/sp1.fnu_autocont, sp2.X_sigma/sp2.fnu_autocont, sp3.X_sigma/sp3.fnu_autocont)
    thezs = (0., 0., 0.)
    the_pdf =  'Boxplots_for_stacks/boxplot_'+method+'_byage_emission.pdf'
    jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_a, line_center_a, win=10., Ncol=1, LL=LL, extra_label=label_age, figsize=(8,12), vel_plot=False, ylims=(0,3), colortab=colortab_age)
    plt.savefig(the_pdf, orientation='portrait', bbox_inches='tight', pad_inches=0.1)
    plt.close()
