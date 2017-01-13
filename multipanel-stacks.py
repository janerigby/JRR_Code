''' Same as multipanel-spectrum-ids.py, but plot multiple spectra.  Using for sub-stacks (highZ, loZ),
young, middle-aged, old.
jrigby, Sept 2016. '''

import jrr
import matplotlib.pyplot as plt
import numpy as np
mage_mode = "reduction"  
#mage_mode = "released"
methods = ('bystars', 'byneb')  # method of determining systemic redshift

# For sanity checking, make plots of the weighted avg(std), but also the median and sqrt(jackknife variance)
colfnu  = ('X_avg', 'X_median')
colfnuu = ('X_sigma', 'X_jack_std')
avg_or_med = ('wtdavg', 'median')

# New Dec 2016:  Compare Steidel et al 2016 to the MagE stack.
for ii in range(0, len(colfnu)):
    for method in methods: 
        stacks = ['magestack_'+method+'_standard_spectrum.txt']
        sp1, dummyLL   = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0], colfnu=colfnu[ii], colfnuu=colfnuu[ii])
        chuck = jrr.mage.read_chuck_UVspec(autofitcont=True)
        shap  = jrr.mage.read_shapley_composite() 
        the_dfs = [sp1, chuck, shap]
        the_zzs = [0.0, 0.0]
        colortab = ('black', 'blue', 'purple')
        the_pdf =  'PDF_Out2/multipanel_stacked_'+avg_or_med[ii]+'_'+method+'_standard_comparetoSteidel2016.pdf'
        (spec_path, line_path) = jrr.mage.getpath(mage_mode)
        (LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, colwave='rest_wave', colfnu='rest_fnu', colfnu_u='rest_fnu_u', colcont='rest_fnu_autocont', title="Standard stack"+method, waverange=(1000,3000), colortab=colortab)
    plt.clf()

# Plot the standard stacks
for ii in range(0, len(colfnu)):
    for method in methods:
        stacks = ['magestack_'+method+'_standard_spectrum.txt']
        sp1, dummyLL = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0], colfnu=colfnu[ii], colfnuu=colfnuu[ii])
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
