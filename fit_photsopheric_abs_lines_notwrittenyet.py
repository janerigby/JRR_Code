import jrr
import numpy as np
import pandas
from matplotlib import pyplot as plt
mage_mode = "reduction"
#mage_mode = "released"
import astropy.convolution

labels = ['rcs0327-E'] #['S1226+2152'] #'rcs0327-B', 'S0957+0509', 'S2111-0114']
(specs) = jrr.mage.getlist_labels(mage_mode, labels)
#(specs) = jrr.mage.getlist_wcont(mage_mode)
print "Grabbed this many spectra:", len(specs)
(spec_path, line_path) = jrr.mage.getpath(mage_mode)

for ii in range(0, len(specs)) :                  #nfnu_stack[ii] will be ii spectrum
    label     = specs['short_label'][ii]
    filename  = specs['filename'][ii]
    zz =  specs['z_neb'][ii]
    (sp, resoln, dresoln)  = jrr.mage.open_spectrum(filename, zz, mage_mode)
    linelist = jrr.mage.get_linelist_name(filename, line_path)   # convenience function
    (LL, zz1) = jrr.mage.get_linelist(linelist)  # Load the linelist into Pandas data frame

    # Automatically fit continuum.  results written to sp.fnu_autocont, sp.flam_autocont.
    jrr.mage.auto_fit_cont(sp, LL, zz)  # Fit the continuum.
    # Will want to run this on S99, and on the MagE data, for the photospehric lines, in the same way.

    # Plot the results.
    plt.figure(figsize=(20,5))
    plt.step(sp.wave[sp.badmask.eq(False)], sp.fnu[sp.badmask.eq(False)], color='b')
#    plt.step(sp.wave, sp.fnu_u, color='y')
    plt.step(sp.wave, sp.fnu_cont, color='y')
    plt.plot(sp.wave, sp.fnu_autocont, color='k')
    plt.ylim(0, 2E-28)
    #plt.xlim(5000,7000)
    plt.title(label + "  z=" + str(zz))
    plt.scatter((5577., 6300.), (0, 0))
    plt.show()

