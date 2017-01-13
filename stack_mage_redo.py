''' Stack the MagE spectra together.
    jrigby, oct 2015.  Revised 3/2016, 7/2016, 8/2016, 9/2016
'''
import jrr
from   astropy.stats import sigma_clip
import extinction
import pandas
import matplotlib.pyplot as plt
import numpy as np
import time
import re
import pdb
from os.path import expanduser

debug = True

def byspline_norm_func(wave, rest_fnu, rest_fnu_u, rest_cont, rest_cont_u, norm_region) :
    # Normalization method by the spline fit continuum
    temp_norm_fnu = rest_fnu / rest_cont
    temp_norm_sig = jrr.util.sigma_adivb(rest_fnu, rest_fnu_u,   rest_cont, rest_cont_u) # propogate uncertainty in continuum fit.
    return(temp_norm_fnu, temp_norm_sig)

# July 2016, I want to make many different stacks.  Therefore, I'm going to rewrite this code,
# with the stacking done in function make_a_stack, so that I can call it multiple times.
# I also split the normalization out as a function, norm_func, so that I can modify it.
# Inputs:
#   labels      # labels of spectra that I want to stack.
#   rootname    # what to call this stack.  products will include this rootname
#   norm_region # rest-frame wavelength region to normalize by
#   norm_func   # function that says how to nornmalize each spectrum.
#   norm_method_text # string that describes how the individual spectra were normalized.  For header.
#   mage_mode   # same as other mage functions.  Where to look for spectra
#   zchoice     # How to set systemic redshift.  Choices are "stars", or "neb"
def make_a_stack(labels, rootname, norm_region, norm_func, norm_method_text, mage_mode, zchoice, deredden=False, EBV=[]) :
    plt.close('all')
    plt.ion()
    plt.figure(figsize=(20,5))
#    specs = jrr.mage.getlist_labels(mage_mode, labels)
    specs = jrr.mage.wrap_getlist(mage_mode, which_list="labels", labels=labels)
    Nspectra = len(specs)
    print "DEBUG, Nspectra is", Nspectra
    stacklo =  800. #A      # Create a rest-frame wavelength array to stack into
    stackhi = 3000. #A
    disp = 0.1 # Angstroms  # observed-frame wavelength binning is ~0.3A pper pix for RCS0327.  So, want ~0.1A in rest-frame
    nbins = int((stackhi - stacklo)/disp)
    wave_stack    = np.linspace(stacklo, stackhi, num=nbins)
    nfnu_stack    = np.zeros(shape=(Nspectra, nbins))   # create array that will hold all the spectra
    nfnu_u_stack  = np.zeros(shape=(Nspectra, nbins))
    jackknife     = np.zeros(shape=(Nspectra, nbins))

    print "Filename    label      N_pixels     rest-frame wavelength range (A)"
    for ii in range(0, Nspectra) :                  #nfnu_stack[ii] will be ii spectrum
        label     = specs['short_label'][ii]
        filename  = specs['filename'][ii]

        if(zchoice == "stars") :     # Update 8/2016, enabling either method to to set systemic redshift, selectable as zchoice
            zz =  specs['z_syst'][ii] # New, using john chisholm's S99 fits to photospheric absorption lines where possible.
        elif(zchoice == "neb") :
            zz =  specs['z_neb'][ii]  # Old, what I used prior to 21 july 2016.
        else : raise ValueError('Error, I do not recognize input zchoice (choice to set systemic redshift) as stars or neb')

        (sp, resoln, dresoln)  = jrr.mage.open_spectrum(filename, zz, mage_mode)
        # set uncertainties high near skylines [O I] 5577\AA\ and [O I]~6300\AA,
        skyline = (5577., 6300.)
        skywidth = 10.0  # flag spectrum +- skywidth of the skyline
        sp.fnu_u[sp['wave'].between(skyline[0]-skywidth, skyline[0]+skywidth)] = 1.0 # huge uncert near skylines
        sp.fnu_u[sp['wave'].between(skyline[1]-skywidth, skyline[1]+skywidth)] = 1.0 # huge uncert near skylines
        sp.fnu_u[sp['fnu_cont'].eq(9999)] = 1.0  # Set huge uncertainties where continuum undefined

        (rest_wave, rest_fnu, rest_fnu_u) = jrr.spec.convert2restframe(sp.wave, sp.fnu,  sp.fnu_u,  zz, 'fnu')
        (junk   , rest_cont, rest_cont_u) = jrr.spec.convert2restframe(sp.wave, sp.fnu_cont, sp.fnu_cont_u, zz, 'fnu')
        # should now have arrays of rest wavelength, fnubda, and continuum, as
        # rest_wave, rest_fnu, rest_fnu_u, rest_cont, rest_cont_u
        print filename, label, len(rest_wave),
        print "  %.2f  %.2f" % (  rest_wave[0], rest_wave[-1:])

        if deredden and len(EBV) :
            this_ebv =  S99.loc[label,  'E(B-V)']
            Rv = 4.05
            Av = -1 * Rv * this_ebv  # stupidity w .Series and .as_matrix() is bc extinction package barfs on pandas. pandas->np->pandas
            rest_fnu    = pandas.Series(extinction.apply(extinction.calzetti00(rest_wave.as_matrix(), Av, Rv), rest_fnu.as_matrix()))
            rest_fnu_u  = pandas.Series(extinction.apply(extinction.calzetti00(rest_wave.as_matrix(), Av, Rv), rest_fnu_u.as_matrix()))
            rest_cont   = pandas.Series(extinction.apply(extinction.calzetti00(rest_wave.as_matrix(), Av, Rv), rest_cont.as_matrix()))
            rest_cont_u = pandas.Series(extinction.apply(extinction.calzetti00(rest_wave.as_matrix(), Av, Rv), rest_cont_u.as_matrix()))
        
        # Normalize the spectrum and error spectrum
        (temp_norm_fnu, temp_sig) = norm_func(rest_wave, rest_fnu, rest_fnu_u, rest_cont, rest_cont_u, norm_region)
        nfnu_stack[ii]   = jrr.spec.rebin_spec_new(rest_wave, temp_norm_fnu,   wave_stack) # normalized fnu, rebinned
        nfnu_u_stack[ii] = jrr.spec.rebin_spec_new(rest_wave, temp_sig, wave_stack)  # uncertainty on above
        plt.step(wave_stack, nfnu_stack[ii])

    # fnu_u_stack is 1sigma uncertainty spectrum.  weight is 1/sigma**2
    nfnu_u_stack += 0.00001            # add a very small number to avoid 1/zero errror
    weight_stack = nfnu_u_stack**-2
    sig2clip = 2  # clip at N sigma,  in astropy.stats.sigma_clip
    # Use a masked average to deal with the whole wavelength range.  numpy.ma is masking
    mask1 = np.isnan(nfnu_stack)    # interp1d writes NaNs in non-overlap region.  Mask these
    mask2 = np.isnan(nfnu_u_stack)  # ditto for uncertainty
    crazy_high = 1000.
    mask3 = np.greater(nfnu_stack, crazy_high) + np.less(nfnu_stack, -1*crazy_high)  # flag crazy flux values.
    mask = mask1 + mask2 + mask3
    print "DEBUGGING masks", mask1.sum(), mask2.sum(), mask3.sum(), mask.sum(), mask.shape
    masked_spectrum   = np.ma.array(nfnu_stack, mask=mask)
    nfnu_clip  = sigma_clip(masked_spectrum, sig=sig2clip, iters=None, axis=0)   ## Sigma clipping
    X_avg,     sumweight1   = np.ma.average(masked_spectrum, axis=0, weights=weight_stack, returned=True) # weighted avg of continuum-normalized spectra
    X_clipavg, sumweight2   = np.ma.average(nfnu_clip,      axis=0, weights=weight_stack, returned=True) # weighted avg of cont-normalized spectra, w sig clip
    X_Ngal = np.sum((1 - mask), axis=0)        # Number of galaxies that went into the stack at that wavelength
    X_sigma     = sumweight1**-0.5
    X_clipsigma = sumweight2**-0.5  # spiky, dunno why.  Presumably legacy of sigma_clip?  Don't use it
    X_median = np.ma.median(masked_spectrum, axis=0)

    plt.step(wave_stack, X_avg, color="black", linewidth=3)
    plt.ylim( -1, 3)
    plt.xlim(1000, 1200)
    plt.draw()
    plt.show()
    #pdb.set_trace()
    if GOSLOW :   plt.pause(4)

    # Jack-knife: drop one spectrum and make weighted average, repeat.  Measure the stdev.
    jack_var = np.zeros_like(X_avg)
    for ii in range(0, Nspectra) :
        print "Jackknife, dropping ", specs['short_label'][ii], " from the stack"
        dropit = np.zeros(shape=nfnu_stack.shape, dtype=bool)
        dropit[ii] = True  # Mask out the ii-th spectrum
        jack_mask = mask + dropit
        masked_spectrum   = np.ma.array(nfnu_stack, mask=jack_mask)
        jackknife[ii], weight = np.ma.average(masked_spectrum, axis=0, weights=weight_stack, returned=True)  # all the work is done here.
        jack_var += (jackknife[ii] - X_avg)**2
    jack_var *= ((Nspectra -1.0)/float(Nspectra))
    X_jack_std = np.sqrt(jack_var)
    # Jackknife variance estimation, from wikipedia: var = (n-1)/n * sum(i=1 to n) of (xi  - x.)**2

    #Output the stacked spectrum
    long_rootname = "magestack_by" + zchoice + "_" + rootname
    head  = '# Stack(s) of the MagE spectral atlas of lensed galaxies, named   ' + long_rootname + '\n'
    head += '# Generated on ' + time.strftime("%d/%m/%Y") + '(dd/mm/yyr)\n'
    head += '# Generated from ' + str(labels) + "\n"
    head += '# ' + norm_method_text + "\n"
    head += '# Columns are: \n'
    head += '# wave    :   Rest-frame vacuum wavelength, in Angstroms.  Covers ' + str(stacklo) + " to " + str(stackhi) + '\n'
    head += '# X_avg   :   Weighted avg of continuum-normalized fnu spectra.\n'
    head += '# X_clipavg: same as X_avg, but weighted avg done with sigma clipping at sigma=' + str(sig2clip) + '\n'
    head += '# X_median:  median of continuum-normalized fnu spectra.  No weighting.\n'
    head += '# X_sigma :  uncertainty on X_avg and X_clipavg, from propagating individual error spectra\n'
    head += '# X_jack_std : sqrt of variance estimate from jackknife test. \n'
    head += '# Ngal    :  Number of galaxies that went into the stack at that wavelength.\n'
    head += 'restwave    X_avg    X_clipavg  X_median  X_sigma   X_jack_std   Ngal'
    outfile = long_rootname + "_spectrum.txt"
    np.savetxt(outfile, np.transpose([wave_stack, X_avg, X_clipavg, X_median, X_sigma, X_jack_std,  X_Ngal]), "%.3f  %.2E  %.2E  %.2E  %.2E %.2E  %d", header=head, comments="")
    maskout_head = "#This mask shows which galaxies were stacked at each wavelength. 0=gal used.  1=gal masked out.\n"
    outfile = long_rootname + "_maskused.txt"
    jack_head = "wave   " + re.sub("\n", "", str(specs.short_label.values))
    jack_head = re.sub("'", "", jack_head)
    jack_head = re.sub("\[", "", jack_head)
    jack_head = re.sub("\]", "", jack_head)
    maskout_head += jack_head
    format_string = "%.3f " + "%.0f " * mask.shape[0]
    np.savetxt(outfile, np.transpose(np.vstack((wave_stack, mask))), fmt=format_string, header=maskout_head, comments="")  # so much pain to figure out this syntax
    # The above np magic allows us to concatenate wave_stack (which is 1D, and mask, which is 2D, and the 2nd dimension isn't defined. painful"

    # output the jackknife stacks
    outfile = long_rootname + "_jackknife.txt"
    head2 =  "# Weighted-averages of MagE spectra, each missing 1 spectrum, for jackknife tests\n#wave(A)   weighted_means_w_jackknife\n"
    head2 += "# Each jackknife is missing the following spectrum:\n" + jack_head
    np.savetxt(outfile, np.transpose(np.vstack((wave_stack, jackknife))), fmt="%.5E", header=head2, comments="")

    # plot the stacked spectrum
    plt.clf()
    plt.step(wave_stack, X_avg, color="black")
    plt.step(wave_stack, X_median, color="blue")
    plt.step(wave_stack, X_sigma, color="orange")
    plt.step(wave_stack, X_jack_std, color="red")
    plt.ylabel(r'$f_\nu$')
    #plt.xlim(stacklo, stackhi)
    plt.xlim(800,3000)
    plt.ylim(-1,2)
    figname = long_rootname + "_quicklook.pdf"
    plt.savefig(figname)
    plt.show()
    if GOSLOW :    plt.pause(4)

    PLOT_NGAL = True
    if PLOT_NGAL :
        plt.clf()
        plt.xlabel(r'wavelength ($\AA$)')
        plt.ylabel("number of galaxies that went into the stack")
        plt.plot(wave_stack, X_Ngal, color="red")
        plt.ylim(0,16)
        figname = long_rootname + "_Ngals_in_stack.pdf"
        plt.savefig(figname)
        plt.show()
    return(0)


# Preparation
GOSLOW = True
mage_mode = "reduction"  # Look for files on satchmo, not released versions on Dropbox

# The standard stack.  Normalize the values and shape of each spectrum by the spline continuum, and stack that
# Load the list of MagE spectrum filenames and redshifts, just for the desired spectra in labels
labels = ['rcs0327-E', 'S0004-0103', 'S0108+0624',  'S0033+0242', 'S0900+2234',  'S0957+0509', 'S1050+0017',  'Horseshoe', 'S1226+2152', 'S1429+1202', 'S1458-0023', 'S1527+0652', 'S2111-0114', 'Cosmic~Eye']
rootname = "standard"
norm_method_text = "Normalized by Janes hand-fit spline continuua, so both value and shape are normalized."
norm_region_dum = (1000.0, 1001.0) # dummy value, in this case not used by byspline_norm_func
#make_a_stack(labels, rootname, norm_region_dum, byspline_norm_func,  norm_method_text, mage_mode, "stars")
#make_a_stack(labels, rootname, norm_region_dum, byspline_norm_func,  norm_method_text, mage_mode, "neb")

# Stack A for John Chisholm: normalize flux but not shape of continuum.  May have trouble w spectral tilt at red and blue ends.
# May be safe near the norm_region
# Use same labels as the standard stack
rootname = "ChisholmstackA"
norm_regionA = jrr.mage.Chisholm_norm_regionA()
norm_method_textA = "Flux normalized to median in spectral region " + str(norm_regionA) + " but spectral shape not flattened."
#make_a_stack(labels, rootname, norm_regionA, jrr.mage.norm_by_median, norm_method_textA, mage_mode, 'stars')
#make_a_stack(labels, rootname, norm_regionA, jrr.mage.norm_by_median, norm_method_textA, mage_mode, 'neb')

# Stack B for John Chisholm: normalize flux but not shape of continuum.  May have trouble w spectral tilt at red and blue ends.
# May be safe near the norm_region
# labels_censored are targets at high enough z that have flux at 1010A
labels_censored = ['S0033+0242', 'S0900+2234',  'S1050+0017',  'Horseshoe', 'S1429+1202', 'S1458-0023', 'S2111-0114']
# dropped s1226 and s1527 from the stack, bc we're considering them individually.  Dropped Cosmic eye bc of DLA there.
rootname = "ChisholmstackB"
norm_regionB = (1040.0, 1045.0)  # Region where John Chisholm says to normalize
norm_method_textB = "Flux normalized to median in spectral region " + str(norm_regionB)  + "but spectral shape not flattened."
#make_a_stack(labels_censored, rootname, norm_regionB, jrr.mage.norm_by_median, norm_method_textB, mage_mode, 'stars')
#make_a_stack(labels_censored, rootname, norm_regionB, jrr.mage.norm_by_median, norm_method_textB, mage_mode, 'neb')

# Stack C for John Chisholm: normalize flux but not shape of continuum.  May have trouble w spectral tilt at red and blue ends.
# May be safe near the norm_region
# labels_censored are targets at high enough z that have flux at 1010A
labels_censored = ['S0033+0242', 'S0900+2234',  'S1050+0017',  'Horseshoe', 'S1429+1202', 'S1458-0023', 'S2111-0114', 'S1527+0652']
# dropped s1226 from the stack, bc we're considering it individually.  Dropped Cosmic eye bc of DLA there.
rootname = "ChisholmstackC"
#make_a_stack(labels_censored, rootname, norm_regionB, jrr.mage.norm_by_median, norm_method_textB, mage_mode, 'stars')
#make_a_stack(labels_censored, rootname, norm_regionB, jrr.mage.norm_by_median, norm_method_textB, mage_mode, 'neb')

# Note:  I am passing a function to make_a_stack, which is the function that says how to
# normalize each input spectrum and uncertainty spectrum.
#(temp_nfnu, temp_sig) = norm_func(rest_wave, rest_fnu, rest_fnu_u, rest_cont, rest_cont_u, norm_region)


# Stack the derreddened spectra.  Same as Stack A, but deredden first.
#   Get E(B-V) values that J Chisholm measured
S99_sumfile = expanduser('~') + "/Dropbox/MagE_atlas/Contrib/S99/New_right_flam/sb99_overview2.txt"
S99 = pandas.read_table(S99_sumfile, delim_whitespace=True, comment="#")
S99.set_index('name', inplace=True, drop=False)
deredden_labels = labels
deredden_labels.remove('S1050+0017')
rootname = "dereddened_StackA"
norm_method_textA = "Flux normalized to median in spectral region " + str(norm_regionA) + " but spectral shape not flattened."
norm_method_textA += "\nBefore stacking, dereddened by E(B-V) values measured by Chisholm, as of 9 Dec 2016, in sb99_overview2.txt ."
make_a_stack(deredden_labels, rootname, norm_regionA, jrr.mage.norm_by_median, norm_method_textA, mage_mode, 'stars', deredden=True, EBV=S99)
make_a_stack(deredden_labels, rootname, norm_regionA, jrr.mage.norm_by_median, norm_method_textA, mage_mode, 'neb',   deredden=True, EBV=S99)
