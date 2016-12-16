import jrr
import matplotlib.pyplot as plt

mage_mode = "reduction"  # Look at copies of the files on Satchmo

(specs_wcont) = jrr.mage.getlist_wcont(mage_mode)
print "This pandas data frame contains filenames & redshifts for all the spectra w continuum fits."
print "Number of spectra: ", len(specs_wcont)
print "Column names: ", specs_wcont.keys()
print "Other summary info: "
print specs_wcont.info()

for label in specs_wcont['short_label'] :
    (sp, resoln, dresoln, LL, zz_syst) = jrr.mage.wrap_open_spectrum(label, mage_mode) 
    print sp['disp'].mean(), sp['disp'].median(), sp['disp'].min(), sp['disp'].max(), "IQR:", jrr.util.IQR(sp['disp'])
