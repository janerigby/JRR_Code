''' Rewriting the MagE stacking tool to use jrr.spec.stack_spectra().  The history here is that I wrote 
the MagE stacker first, then generalized the code for reuse as jrr.spec.stack_spectra().  Now I am going
back and rewriting the MagE stacking tool so that it uses this function, for two reasons: 
1) I can cite the MagE stack paper for the methodology, and 
2) doing so lets me triple-check the stacking.
jrigby, July 2017.

Methodology:

- Create the new wavelength array.  
- Read in the individual MagE spectra
- Mask out skylines
- Mask out intervening absorbers
- Apply Galactic Extinction (added per referee's suggestion, correcting omission on our part)
- De-redshift the MagE spectra by their systemic redshifts
- Optionally, remove internal extinction
- Normalize the spectra by the preferred method
- Feed to jrr.spec.stack_spectra() a dataframe of rest-frame spectra.  Stack them.
- Write output.'''


import jrr
from   astropy.stats import sigma_clip
import pandas
import matplotlib.pyplot as plt
import numpy as np
import time
import re
import pdb
from os.path import expanduser

def make_linear_wavelength_array(wavelo, wavehi, disp) :
    return ( np.arange(wavelo, wavehi, disp) )

def mask_skylines(sp, colwave='wave', colcont='fnu_cont') :
    # set uncertainties high near skylines [O I] 5577\AA\ and [O I]~6300\AA,
    skyline = (5577., 6300.)
    skywidth = 10.0  # flag spectrum +- skywidth of the skyline
    sp.loc[sp[colwave].between(skyline[0]-skywidth, skyline[0]+skywidth), 'fnu_u'] = 1.0 # huge uncert near skylines
    sp.loc[sp[colwave].between(skyline[1]-skywidth, skyline[1]+skywidth), 'fnu_u'] = 1.0 # huge uncert near skylines
    sp.loc[sp[colcont].eq(9999), 'fnu_u'] = 1.0  # Set huge uncertainties where continuum undefined
    return(0)

def mask_intervening(sp, LL, colint='interve') :
    vmask = 200. # +-200km/s
    sp[colint] = False
    LL['vmask'] = vmask
    jrr.spec.flag_near_lines(sp, LL, linetype=('INTERVE',), colmask=colint)
    sp.loc[sp[colint], 'fnu_u'] = 1. # Set huge uncertainties at positions of known intervening absorbers
    print "DEBUGGING, I masked N pixels for intervening absorbers", sp[colint].sum() 
    return(0)

def make_header_for_stack(rootname, zchoice, norm_method_text) :
    long_rootname = "magestack_by" + zchoice + "_" + rootname
    head  = '# Stack(s) of the MagE spectral atlas of lensed galaxies, named   ' + long_rootname + '\n'
    head += '# Generated on ' + time.strftime("%d/%m/%Y") + '(dd/mm/yyr)\n'
    head += '# Generated by mage_stack_v3.py\n'
    head += '# Generated from ' + str(labels) + "\n"
    head += '# ' + norm_method_text + "\n"
    return(head)

def prep_spectra_for_stacking(df, zchoice, EBV=[], deredden=False, dereddenMW=False, colcont="fnu_cont") : #Do this once per zchoice
    for short_label in df.keys() :  # step through each spectrum in labels
        #print "DEBUGGING, working on", short_label
        df[short_label]['zeros'] = 0.0
        mask_skylines(df[short_label])
        mask_intervening(df[short_label], LL[short_label])
        if   zchoice == 'stars' : zz = specs.ix[short_label]['z_syst']
        elif zchoice == 'neb'   : zz = specs.ix[short_label]['z_neb']
        if dereddenMW :    jrr.mage.deredden_MW_extinction(df[short_label], specs['EBV_MW'][short_label])  # DEPRECATED. INPUT SPECTRA NOW ALREADY DEREDDENED
        jrr.mage.convert_spectrum_to_restframe(df[short_label], zz)
        if deredden and short_label in EBV.index :
            jrr.mage.deredden_internal_extinction(df[short_label], EBV.ix[short_label]['E(B-V)'], colcont)
    return(0)


def stack_mage_spectra(df, LL, specs, labels, rootname, norm_region, norm_func, norm_method_text, mage_mode, zchoice, outwave, deredden=False, EBV=[], colcont='rest_fnu_cont', colcontu='rest_fnu_cont_u') :
    ndf = {}
    # df is a big honking dictionary of dataframes containing the spectra.  LL is a dict of dataframes containing the linelist
    for short_label in labels :  # step through each spectrum in labels
        ndf[short_label] = pandas.DataFrame(data=pandas.Series(outwave), columns=('rest_wave',))  # New data frame, rebinned wavelength
        # Normalize the spectrum and error spectrum
        (temp_norm_fnu, temp_sig) = norm_func(df[short_label]['rest_wave'], df[short_label]['rest_fnu'], df[short_label]['rest_fnu_u'], df[short_label][colcont], df[short_label][colcontu], norm_region)
        ndf[short_label]['fnorm']    = jrr.spec.rebin_spec_new(df[short_label]['rest_wave'], temp_norm_fnu,  outwave)  # normalized fnu, rebinned
        ndf[short_label]['fnorm_u']  = jrr.spec.rebin_spec_new(df[short_label]['rest_wave'], temp_sig,  outwave)       # uncertainty on above
        ndf[short_label]['interve']  = jrr.spec.rebin_spec_new(df[short_label]['rest_wave'], df[short_label]['interve'],  outwave).astype('bool') #intervening
    (stacked, nf, nf_u) = jrr.spec.stack_spectra(ndf, colwave='rest_wave', colf='fnorm', colfu='fnorm_u', colmask='interve', output_wave_array=outwave)
    long_rootname = "magestack_by" + zchoice + "_" + rootname
    outfile = long_rootname + "_spectrum.txt"
    stacked.drop(['fsum', 'fsum_u', 'fmedianxN'], axis=1, inplace=True)  # These columns arent meaningful
    stacked.to_csv('temp', index=False, sep='\t')
    head = make_header_for_stack(rootname, zchoice, norm_method_text)
    jrr.util.put_header_on_file('temp', head, outfile)
    ax = stacked.plot(x='rest_wave', y='fweightavg', color='black')
    stacked.plot(x='rest_wave', y='fjack_std', color='grey', ax=ax)
    plt.ylim(0, stacked['fweightavg'].median() * 3)
    plt.draw()
    return(stacked)

##### SETUP #######
mage_mode = 'reduction'
#mage_mode = "released" # while on plane. ** TEMP
stacklo =  800. #A      # Create a rest-frame wavelength array to stack into
stackhi = 3000. #A
disp = 0.1 # Angstroms  # observed-frame wavelength binning is ~0.3A pper pix for RCS0327.  So, want ~0.1A in rest-frame

#### Which spectra to stack  #######
labels = ['rcs0327-E', 'S0004-0103', 'S0108+0624',  'S0033+0242', 'S0900+2234',  'S0957+0509', 'S1050+0017',  'Horseshoe', 'S1226+2152', 'S1429+1202', 'S1458-0023', 'S1527+0652', 'S2111-0114', 'Cosmic~Eye']
ws99labels = ['rcs0327-E', 'S0004-0103', 'S0108+0624',  'S0033+0242', 'S0900+2234',  'S0957+0509', 'Horseshoe', 'S1226+2152', 'S1429+1202', 'S1458-0023', 'S1527+0652', 'S2111-0114', 'Cosmic~Eye']
labels_censoredB = ['S0033+0242', 'S0900+2234',  'S1050+0017',  'Horseshoe', 'S1429+1202', 'S1458-0023', 'S2111-0114']  # dropped s1226 and s1527 from the stack, bc we're considering them individually.  Dropped Cosmic eye bc of DLA there.
labels_censoredC = ['S0033+0242', 'S0900+2234',  'S1050+0017',  'Horseshoe', 'S1429+1202', 'S1458-0023', 'S2111-0114', 'S1527+0652']
deredden_labels = labels
deredden_labels.remove('S1050+0017')
testlabels = ['rcs0327-E', 'S0004-0103', 'S0108+0624']

# Normalization regions (Angstroms)
norm_region_dum = (1000.0, 1001.0)
norm_regionA = jrr.mage.Chisholm_norm_regionA()
norm_regionB = (1040.0, 1045.0)  # Region where John Chisholm says to normalize


#### Make a bunch of stacks

outwave1 = make_linear_wavelength_array(stacklo, stackhi, disp)  # Not ideal, but same as before.
systemic_methods = ['stars', 'neb']

dereddenMW=False   # Correct for MW reddening in prep_spectra_for_stacking()?  Should be FALSE, instead using input spectra that have been corrected
MWdr = True        # Get spectra that have already been corrected for MW reddening?  SHOULD BE TRUE

for zchoice in systemic_methods :
    (df, resoln, dresoln, LL, zz_sys, specs) = jrr.mage.open_many_spectra(mage_mode, which_list='wcont', zchoice='stars', MWdr=MWdr, addS99=True) # set systemic redshift to stars for all.  Use z_neb when want neb
    prep_spectra_for_stacking(df, zchoice, deredden=False, dereddenMW=dereddenMW)
    # The standard stack.  Normalize the values and shape of each spectrum by the spline continuum, and stack that
    # Load the list of MagE spectrum filenames and redshifts, just for the desired spectra in labels
    rootname = "standard"
    norm_method_text = "Normalized by Janes hand-fit spline continuua, so both value and shape are normalized."
    stacked = stack_mage_spectra(df, LL, specs, labels, rootname, norm_region_dum, jrr.spec.byspline_norm_func, norm_method_text, mage_mode, zchoice, outwave1, colcont='rest_fnu_cont')

    rootname = "divbyS99"  # for each input spectrum, divide by the s99 continuum to normalize.
    norm_method_text = "Normalized by John Chisholms S99 fits to each spectrum, so both value and shape are normalized."
    stacked = stack_mage_spectra(df, LL, specs, ws99labels, rootname, norm_region_dum, jrr.spec.byspline_norm_func, norm_method_text, mage_mode, zchoice, outwave1, colcont='rest_fnu_s99model', colcontu='zeros')

    # Stack A for John Chisholm: normalize flux but not shape of continuum.  May have trouble w spectral tilt at red and blue ends.
    # May be safe near the norm_region
    # Use same labels as the standard stack
    rootname = "ChisholmstackA"
    norm_regionA = jrr.mage.Chisholm_norm_regionA()
    norm_method_textA = "Flux normalized to median in spectral region " + str(norm_regionA) + " but spectral shape not flattened."
    stacked = stack_mage_spectra(df, LL, specs, labels, rootname, norm_regionA, jrr.spec.norm_by_median, norm_method_textA, mage_mode, zchoice, outwave1)

    # Stack B for John Chisholm: normalize flux but not shape of continuum.  May have trouble w spectral tilt at red and blue ends.
# May be safe near the norm_region
    # labels_censored are targets at high enough z that have flux at 1010A
    rootname = "ChisholmstackB"
    norm_method_textB = "Flux normalized to median in spectral region " + str(norm_regionB)  + "but spectral shape not flattened."
    stacked = stack_mage_spectra(df, LL, specs, labels_censoredB, rootname, norm_regionB, jrr.spec.norm_by_median, norm_method_textB, mage_mode, zchoice, outwave1)

    # Stack C for John Chisholm: normalize flux but not shape of continuum.  May have trouble w spectral tilt at red and blue ends.
    # May be safe near the norm_region
    # labels_censored are targets at high enough z that have flux at 1010A
    # dropped s1226 from the stack, bc we're considering it individually.  Dropped Cosmic eye bc of DLA there.
    rootname = "ChisholmstackC"
    stacked = stack_mage_spectra(df, LL, specs, labels_censoredC, rootname, norm_regionB, jrr.spec.norm_by_median, norm_method_textB, mage_mode, zchoice, outwave1)

    # Note:  I am passing a function to make_a_stack, which is the function that says how to
    # normalize each input spectrum and uncertainty spectrum.

    
for zchoice in systemic_methods:
    # Stack the derreddened spectra.  Same as Stack A, but deredden first.
    #   Get E(B-V) values that J Chisholm measured
    S99_sumfile = expanduser('~') + "/Dropbox/MagE_atlas/Contrib/S99/New_right_flam/sb99_overview2.txt"
    S99 = pandas.read_table(S99_sumfile, delim_whitespace=True, comment="#")
    S99.set_index('name', inplace=True, drop=False)
    rootname = "dereddened_StackA"
    norm_method_textA = "Flux normalized to median in spectral region " + str(norm_regionA) + " but spectral shape not flattened."
    norm_method_textA += "\nBefore stacking, dereddened by E(B-V) values measured by Chisholm, as of 9 Dec 2016, in sb99_overview2.txt ."
    # Reload everything fresh, because need to deredden these.
    (df, resoln, dresoln, LL, zz_sys, specs) = jrr.mage.open_many_spectra(mage_mode, which_list='wcont', zchoice='stars', MWdr=MWdr) # set systemic redshift to stars for all.  Use z_neb when want neb
    prep_spectra_for_stacking(df, zchoice, S99, deredden=True, dereddenMW=dereddenMW)
    stacked = stack_mage_spectra(df, LL, specs, deredden_labels, rootname, norm_regionA, jrr.spec.norm_by_median, norm_method_textA, mage_mode, zchoice, outwave1)
