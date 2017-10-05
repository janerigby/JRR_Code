''' Same as multipanel-spectrum-ids.py, but plot multiple spectra.
jrigby, Sept 2016.
Run from /Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Mage/Analysis/Plot-all'''

import jrr
import matplotlib.pyplot as plt
import numpy as np
mage_mode = "reduction"  
methods = ('bystars', 'byneb')  # method of determining systemic redshift

# For sanity checking, make plots of the weighted avg(std), but also the median and sqrt(jackknife variance)
colfnu  = ('fweightavg', 'fmedian')
colfnuu = ('fweightavg_u', 'fjack_std')
avg_or_med = ('wtdavg', 'median')

# Dec 2016:  Compare Steidel et al 2016 to the MagE stack.
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
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, colwave='rest_wave', colfnu='rest_fnu', colfnu_u='rest_fnu_u', colcont='rest_fnu_autocont', title="Standard stack"+method, waverange=(1000,2999), colortab=colortab)
        plt.clf()
plt.close("all")

# May 2017, compare our COS stack to the MagE stack.
sp1, dummyLL   = jrr.mage.open_stacked_spectrum(mage_mode, colfnu=colfnu[0], colfnuu=colfnuu[0])
cos  = jrr.mage.read_our_COS_stack(resoln="matched_mage")
cos['rest_fnu_autocont'] = 1.0 
the_dfs = [sp1, cos]
the_zzs = [0.0, 0.0]
colortab = ('black', 'blue')
jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile="MagE_vs_COSR3300.pdf", plot_cont=True, norm_by_cont=False, apply_bad=False, colwave='rest_wave', colfnu='rest_fnu', colfnu_u='rest_fnu_u', colcont='unity', title="MagE_vs_COS", waverange=(1000,2999), colortab=colortab)
plt.close("all")

cos  = jrr.mage.read_our_COS_stack(resoln="full")
the_dfs = [sp1, cos]
the_zzs = [0.0, 0.0]
colortab = ('black', 'blue')
jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile="MagE_vs_COSR1.4E4.pdf", plot_cont=True, norm_by_cont=False, apply_bad=False, colwave='rest_wave', colfnu='rest_fnu', colfnu_u='rest_fnu_u', colcont='unity', title="MagE_vs_COS", waverange=(1000,2999), colortab=colortab)
plt.close("all")


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
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title="Standard stack"+method, waverange=(1000,2999))
        # make plain stack to drop into paper
        # figure out new outfile for this, so it doesn't overwrite
        the_pdf =  'PDF_Out2/multipanel_stacked_'+avg_or_med[ii]+'_'+method+'_standard_nolabels.pdf'
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, colcont='unity', plot_cont=True, norm_by_cont=False, apply_bad=True, waverange=(1000,2999), plotx2=False, ylim=(0,1.3),  annotate=())
        plt.clf()
plt.close("all")


# Repeat for StackA
for ii in range(0, len(colfnu)):
    for method in methods:
        stacks = ['magestack_'+method+'_ChisholmstackA_spectrum.txt']
        sp1, dummyLL = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0], colfnu=colfnu[ii], colfnuu=colfnuu[ii])
        the_dfs = [sp1]
        the_zzs = [0.0]  # Stacked spectra are already in rest frame.
        the_pdf =  'PDF_Out2/multipanel_stacked_'+avg_or_med[ii]+'_'+method+'_stackA.pdf'
        (spec_path, line_path) = jrr.mage.getpath(mage_mode)
        (LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title="Stack A"+method, waverange=(1000,2999))
        the_pdf =  'PDF_Out2/multipanel_stacked_'+avg_or_med[ii]+'_'+method+'_stackA_nolabels.pdf'
        jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, waverange=(1000.,2999.), plotx2=False, annotate=())
        plt.clf()
    plt.clf()
plt.close("all")
