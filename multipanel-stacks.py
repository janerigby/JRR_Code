''' Same as multipanel-spectrum-ids.py, but plot multiple spectra.  Using for sub-stacks (highZ, loZ),
young, middle-aged, old.
jrigby, Sept 2016. '''

import jrr
import matplotlib.pyplot as plt
import numpy as np
mage_mode = "reduction"  
mage_mode = "released"
methods = ('bystars', 'byneb')  # method of determining systemic redshift

# Compare high and low metallicity stacks
for method in methods :
    stacks = ['magestack_'+method+'_highZ_spectrum.txt', 'magestack_'+method+'_lowZ_spectrum.txt']
    sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
    sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
    the_dfs = [sp1, sp2]
    the_zzs = [0.0, 0.0]  # Stacked spectra are already in rest frame.
    the_pdf =  'PDF_Out2/multipanel_stacked_'+method+'_byZ.pdf'
    (spec_path, line_path) = jrr.mage.getpath(mage_mode)
    (LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
    jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, colortab=('black', 'mediumblue'), title="Compare highZ (black) and lowZ (blue) stacks "+method)
plt.clf()

# Compare young, middle-age, and old stacks
for method in methods:
    stacks = ['magestack_'+method+'_younglt8Myr_spectrum.txt', 'magestack_'+method+'_midage8to16Myr_spectrum.txt', 'magestack_'+method+'_oldgt16Myr_spectrum.txt']
    sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
    sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
    sp3 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[2])
    the_dfs = [sp1, sp2, sp3]
    the_zzs = [0.0, 0.0, 0.0]  # Stacked spectra are already in rest frame.
    the_pdf =  'PDF_Out2/multipanel_stacked_'+method+'_byage.pdf'
    (spec_path, line_path) = jrr.mage.getpath(mage_mode)
    (LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
    jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title="Compare young (blue), middle-age (black), and old (red) stacks "+method, colortab=('blue', 'black', 'red'), waverange=(1000,3000))
plt.clf()

# Plot the standard stacks
for method in methods:
    stacks = ['magestack_'+method+'_standard_spectrum.txt']
    sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
    the_dfs = [sp1]
    the_zzs = [0.0]  # Stacked spectra are already in rest frame.
    the_pdf =  'PDF_Out2/multipanel_stacked_'+method+'_standard.pdf'
    (spec_path, line_path) = jrr.mage.getpath(mage_mode)
    (LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
    jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title="Standard stack"+method, waverange=(1000,3000))
    # make plain stack to drop into paper
    # figure out new outfile for this, so it doesn't overwrite
    the_pdf =  'PDF_Out2/multipanel_stacked_'+method+'_standard_nolabels.pdf'
    jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, waverange=(1000,3000), annotate=())
plt.clf()
