import pandas
import numpy as np
from astropy.stats import sigma_clip
import extinction
import jrr
import time
import re
mage_mode = "reduction"
import matplotlib.pyplot as plt

# Example of using extinction package
Av = -1.0  # negative to deredden; positive to redden
Rv = 4.05
wave = np.arange(4000.,8000.,500.)
ff   = np.ones_like(wave)
result = extinction.apply(extinction.calzetti00(wave, Av, Rv), ff)
#print wave, result

norm_regionA = (1267.0, 1276.0)  # Region where John Chisholm says to normalize
def norm_by_median(wave, rest_fnu, rest_fnu_u, rest_cont, rest_cont_u, norm_region) :
    normalization = np.median(rest_fnu[wave.between(*norm_region)])
    #print "normalization was", normalization, type(normalization)
    return(rest_fnu / normalization,  rest_fnu_u / normalization)
    
# De-redden the mage spectra by the E(B-V) that J Chisholm measured.  
S99_sumfile = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/MagE_atlas/Contrib/S99/sb99-parameters.txt"
S99 = pandas.read_table(S99_sumfile, delim_whitespace=True, comment="#")
S99 = S99[S99['stackit'].eq(1)]
S99.set_index('name', inplace=True, drop=False)
# Set index by name, so that i can use same as short_label, to get EBV['short_label'], and de-redden
labels = S99['name'].values
(specs) = jrr.mage.getlist_labels(mage_mode, labels)
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
rootname = "dereddened"
# filter by S99.stackit

# prepare the stack
Nspectra = len(specs)
stacklo =  800. #A      # Create a rest-frame wavelength array to stack into
stackhi = 3000. #A
disp = 0.1 # Angstroms  # observed-frame wavelength binning is ~0.3A pper pix for RCS0327.  So, want ~0.1A in rest-frame
nbins = int((stackhi - stacklo)/disp)
wave_stack    = np.linspace(stacklo, stackhi, num=nbins)
nfnu_stack    = np.zeros(shape=(Nspectra, nbins))   # create array that will hold all the spectra
nfnu_u_stack  = np.zeros(shape=(Nspectra, nbins))  
jackknife     = np.zeros(shape=(Nspectra, nbins))

plt.ion()
  
for ii in range(0, len(specs)) :           
    label     = specs['short_label'][ii]
    filename  = specs['filename'][ii]
    zz =  np.float64(specs['z_syst'][ii])
    (sp, resoln, dresoln)  = jrr.mage.open_spectrum(filename, zz, mage_mode)
    linelist = jrr.mage.get_linelist_name(filename, line_path)   # convenience function
    (LL, zz1) = jrr.mage.get_linelist(linelist)  # Load the linelist into Pandas data frame
    EBV = S99['EBV'][label]
    norm_method_text = "Flux normalized to median in spectral region " + str(norm_regionA) + " but spectral shape not flattened."
    norm_method_text += "\nSpectrum shifted to rest-frame wavelenth, and E(B-V) dereddening applied."
    norm_method_text += "\nDereddened by E(B-V)=" + str(EBV) + "Rv=" + str(Rv)

    Av = -1 * Rv * EBV
    sp['x']  = sp.wave[~sp.badmask] / (1.+zz)
    sp['y']  = sp.fnu[~sp.badmask]  / (1.+zz)
    sp['dy'] = sp.fnu_u[~sp.badmask]/ (1.+zz)
    (sp.y, sp.dy) = norm_by_median(sp.x, sp.y, sp.dy, None, None, norm_regionA)
    sp['dered']   = extinction.apply(extinction.calzetti00(sp.x.as_matrix(), Av, Rv), sp.y.as_matrix())
    sp['dered_u'] = extinction.apply(extinction.calzetti00(sp.x.as_matrix(), Av, Rv), sp.y.as_matrix())
    (sp.dered, sp.dered_u) = norm_by_median(sp.x, sp.dered, sp.dered_u, None, None, norm_regionA)
    plt.figure(1) # before
    plt.plot(sp.x, sp.y, linewidth=0.1)
    plt.ylim(-2,4)
    plt.xlim(1000,2000)
    plt.figure(2) # after
    plt.plot(sp.x, sp.dered, linewidth=0.1)
    plt.ylim(-2,4)
    plt.xlim(1000,2000)
    plt.show()
    nfnu_stack[ii]   = jrr.spec.rebin_spec_new(sp.x, sp.dered,   wave_stack) # normalized fnu, rebinned
    nfnu_u_stack[ii] = jrr.spec.rebin_spec_new(sp.x, sp.dered_u, wave_stack)  # uncertainty on above

nfnu_u_stack += 0.00001            # add a very small number to avoid 1/zero errror
weight_stack = nfnu_u_stack**-2
sig2clip = 2  # clip at N sigma,  in astropy.stats.sigma_clip
mask1 = np.isnan(nfnu_stack)    # interp1d writes NaNs in non-overlap region.  Mask these
mask2 = np.isnan(nfnu_u_stack)  # ditto for uncertainty 
mask = mask1 #+ mask2 + mask3 
masked_spectrum   = np.ma.array(nfnu_stack, mask=mask)
nfnu_clip  = sigma_clip(masked_spectrum, sig=sig2clip, iters=None, axis=0)   ## Sigma clipping
X_avg,     sumweight1   = np.ma.average(masked_spectrum, axis=0, weights=weight_stack, returned=True) # weighted avg of continuum-normalized spectra
X_clipavg, sumweight2   = np.ma.average(nfnu_clip,      axis=0, weights=weight_stack, returned=True) # weighted avg of cont-normalized spectra, w sig clip
X_Ngal = np.sum((1 - mask), axis=0)        # Number of galaxies that went into the stack at that wavelength
X_sigma     = sumweight1**-0.5
X_clipsigma = sumweight2**-0.5  # spiky, dunno why.  Presumably legacy of sigma_clip?  Don't use it
X_median = np.ma.median(masked_spectrum, axis=0)


# Jack-knife: drop one spectrum and make weighted average, repeat.  Measure the stdev.
jack_var = np.zeros_like(X_avg)
for ii in range(0, Nspectra) :
    print "Jackknife, dropping ", specs['short_label'][ii], " from the stack"
    dropit = np.zeros(shape=nfnu_stack.shape, dtype=bool)  
    dropit[ii] = True  # Mask out the ii-th spectrum
    jack_mask = mask1 + dropit
    masked_spectrum   = np.ma.array(nfnu_stack, mask=jack_mask)
    jackknife[ii], weight = np.ma.average(masked_spectrum, axis=0, weights=weight_stack, returned=True)  # all the work is done here.
    jack_var += (jackknife[ii] - X_avg)**2
jack_var *= ((Nspectra -1.0)/float(Nspectra))
X_jack_std = np.sqrt(jack_var)
# Jackknife variance estimation, from wikipedia: var = (n-1)/n * sum(i=1 to n) of (xi  - x.)**2

#Output the stacked spectrum
head  = '# Stack(s) of the MagE spectral atlas of lensed galaxies.  \n'
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
head += 'restwave    X_avg    X_clipavg  X_median  X_sigma   X_jack_std   Ngal\n'
outfile = "mage_stack_" + rootname + "_spectrum.txt"
np.savetxt(outfile, np.transpose([wave_stack, X_avg, X_clipavg, X_median, X_sigma, X_jack_std,  X_Ngal]), "%.3f  %.2E  %.2E  %.2E  %.2E %.2E  %d", header=head)

maskout_head = "#This mask shows which galaxies were stacked at each wavelength. 0=gal used.  1=gal masked out.\n"
outfile = "mage_stack_" + rootname + "_maskused.txt"
jack_head = "wave   " + re.sub("\n", "", str(specs.short_label.values))
jack_head = re.sub("'", "", jack_head)
jack_head = re.sub("\[", "", jack_head)
jack_head = re.sub("\]", "", jack_head)
maskout_head += jack_head
format_string = "%.3f " + "%.0f " * mask.shape[0]
np.savetxt(outfile, np.transpose(np.vstack((wave_stack, mask))), fmt=format_string, header=maskout_head)  # so much pain to figure out this syntax
# The above np magic allows us to concatenate wave_stack (which is 1D, and mask, which is 2D, and the 2nd dimension isn't defined. painful"

# output the jackknife stacks
outfile = "mage_stack_" + rootname + "_jackknife.txt"
head2 =  "# Weighted-averages of MagE spectra, each missing 1 spectrum, for jackknife tests\n#wave(A)   weighted_means_w_jackknife\n"
head2 += "# Each jackknife is missing the following spectrum:\n" + jack_head
np.savetxt(outfile, np.transpose(np.vstack((wave_stack, jackknife))), fmt="%.5E", header=head2)

# plot the stacked spectrum
plt.figure(3)
plt.step(wave_stack, X_avg, color="black")
plt.step(wave_stack, X_median, color="blue")
plt.step(wave_stack, X_sigma, color="orange")
plt.step(wave_stack, X_jack_std, color="red")
plt.ylabel(r'$f_\nu$')
#plt.xlim(stacklo, stackhi)
plt.xlim(800,3000)
plt.ylim(-1,2)
figname = "mage_stack_" + rootname + "_quicklook.pdf"
plt.savefig(figname)
plt.show()

PLOT_NGAL = True
if PLOT_NGAL :    
    plt.figure(4)
    plt.xlabel(r'wavelength ($\AA$)')
    plt.ylabel("number of galaxies that went into the stack")
    plt.plot(wave_stack, X_Ngal, color="red")
    plt.ylim(0,16)
    figname = "mage_stack_" + rootname + "_Ngals_in_stack.pdf"
    plt.savefig(figname)
    plt.show()
