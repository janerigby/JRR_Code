''' Want to read Ayan's equivalent width measurements, and plot as a function of
stellar age and metallicity as fit by John Chisholm.  Rigby, Oct 2016'''

from matplotlib import pyplot as plt
import jrr
import pandas
import numpy as np
from os.path import expanduser

pandas.set_option('display.width', 2000)

# Read Ayan's file
homedir = expanduser("~")
flux_file = homedir + "/Dropbox/MagE_atlas/Contrib/EWs/allspec_fitted_emission_linelist.txt"
dfL = pandas.read_table(flux_file, delim_whitespace=True, comment="#")

line_lab1 = "CIII]1906"
line_lab2 = "CIII]1908"

subset1 = dfL[dfL['line_lab'].eq(line_lab1)]
subset2 = dfL[dfL['line_lab'].eq(line_lab2)]
ax = subset1.plot(x='EWr_fit', y='EWr_sum', kind='scatter', xerr='EWr_fit_u', yerr='EWr_sum_u', color='blue') # native pandas plotting.  pandas.DataFrame.plot()
subset2.plot(x='EWr_fit', y='EWr_sum', kind='scatter', xerr='EWr_fit_u', yerr='EWr_sum_u', xlim=(-2,2), ylim=(-2,2), ax=ax, color='red') 
plt.plot( (-2,2), (-2,2), color='green')
plt.savefig('EWrfit_EWrsum.pdf')
plt.show()
# Why don't the EWr_fit and EWr_sum agree?  Need to ask Ayan
plt.close()
