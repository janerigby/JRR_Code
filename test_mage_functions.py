# Test the MagE routines in jrr.py
import jrr
import matplotlib.pyplot as plt
from time import sleep

mage_mode = "reduction"  # Look at copies of the files on Satchmo
#mage_mode = "released"   # Look at the released version of MagE spectra, on Dropbox

print "Here are some examples of grabbing filenames and redshifts for MagE spectra."
print "Results are stored in a Pandas data frame"

print "Example 1)  Grab info for all spectra"
(specs) = jrr.mage.getlist(mage_mode)  # get list of MagE spectrum filenames and redshifts

print "Example 2) Grab all that have continuum fits"
(specs_wcont) = jrr.mage.getlist_wcont(mage_mode)
print "This pandas data frame contains filenames & redshifts for all the spectra w continuum fits."
print "Number of spectra: ", len(specs_wcont)
print "Column names: ", specs_wcont.keys()
print "Other summary info: "
print specs_wcont.info()

print "Example 3) Grab info for a list of spectra (short_label), in this case, all spectra of RCS0327"
labels = ['rcs0327-B', 'rcs0327-E', 'rcs0327-Ehires', 'rcs0327-Elores', 'rcs0327-G', 'rcs0327-U', 'rcs0327-BDEFim1', 'rcs0327-counterarc']
(specs) = jrr.mage.getlist_labels(mage_mode, labels)

(specs) = jrr.mage.getlist_wcont(mage_mode)   # Get all w continuua

(specs) = jrr.mage.getlist(mage_mode)

print "It's then easy to iterate through this list, open spectra, and plot"
 
for jj in range(0, len(specs)) :                  #flam_stack[jj] will be jj spectrum
    label     = specs['short_label'][jj]
    filename  = specs['filename'][jj]
    zz =  specs['z_neb'][jj]
    (specdir, linedir) = jrr.mage.getpath(mage_mode)
    (sp, resoln, dresoln)  = jrr.mage.open_spectrum(filename, zz, mage_mode)
    print("I opened file ", filename) 
    plt.plot(sp.wave, sp.fnu)
    plt.xlim(4000,8000)
    plt.ylim(-1E-28,3E-28)
print "All done"
plt.title("All the MagE spectra")
plt.show()

(sp, LL) = jrr.mage.open_stacked_spectrum(mage_mode)
plt.plot(sp.restwave, sp.X_avg)
plt.title("stacked spectrum")
plt.show()

