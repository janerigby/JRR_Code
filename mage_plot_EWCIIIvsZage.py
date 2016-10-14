''' Want to read Ayan's equivalent width measurements, and plot as a function of
stellar age and metallicity as fit by John Chisholm.  Rigby, Oct 2016'''

from matplotlib import pyplot as plt
import jrr
import pandas
import numpy as np
from os.path import expanduser

pandas.set_option('display.width', 2000)
line_lab1 = "CIII]1906"
line_lab2 = "CIII]1908"
line_lab3 = "Ly-alpha"
# Need to re-write this so that it's a loop over several lines of interest.
# Then, do the CIII1907+1909 by hand (as done below), since it's a sum.  Special case.



# Read Ayan's file
homedir = expanduser("~")
flux_file = homedir + "/Dropbox/MagE_atlas/Contrib/EWs/allspec_fitted_emission_linelist.txt"
dfL = pandas.read_table(flux_file, delim_whitespace=True, comment="#")

# Read John's file with ages, metallicities
ageZ_file = homedir + "/Dropbox/MagE_atlas/Contrib/S99/sb99-parameters-sigfigs.txt"
dftZ = pandas.read_table(ageZ_file, delim_whitespace=True, comment="#")
dftZ.set_index('name', inplace=True)

subset1 = dfL[dfL['line_lab'].eq(line_lab1)]
subset1.set_index('label', inplace=True)
subset2 = dfL[dfL['line_lab'].eq(line_lab2)]
subset2.set_index('label', inplace=True)

subset3 = dfL[dfL['line_lab'].eq(line_lab3)]
subset3.set_index('label', inplace=True)
dfLya = pandas.concat([dftZ, subset3], axis=1)


sum_EW   = subset1['EWr_fit'] + subset2['EWr_fit']
sum_EW_u = np.sqrt(subset1['EWr_fit_u']**2 + subset2['EWr_fit_u']**2)
EWsum = pandas.concat([sum_EW, sum_EW_u], axis=1).reset_index()
EWsum.set_index('label', inplace=True)

plot_debug = False
if plot_debug :
    ax = subset1.plot(x='EWr_fit', y='EWr_sum', kind='scatter', xerr='EWr_fit_u', yerr='EWr_sum_u', color='blue') # native pandas plotting.  pandas.DataFrame.plot()
    #subset2.plot(x='EWr_fit', y='EWr_sum', kind='scatter', xerr='EWr_fit_u', yerr='EWr_sum_u', xlim=(-2,2), ylim=(-2,2), ax=ax, color='red') 
    plt.plot( (-2,2), (-2,2), color='green')
    plt.savefig('EWrfit_EWrsum.pdf')
    plt.show()
    # Why don't the EWr_fit and EWr_sum agree?  Have asked Ayan.  Maybe b/c _sum is affected by blend?

# Let's move on.  Need to sum EWr for the two lines, and plot it, versus age, Z from stellar fits
df1 = pandas.concat([dftZ, subset1], axis=1)
df2 = pandas.concat([dftZ, subset2], axis=1)
df3 = pandas.concat([dftZ, EWsum], axis=1)
ax = df1.plot(x='age', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='blue')
df2.plot(x='age', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='red', ax=ax)
df3.plot(x='age', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='green', ax=ax)
plt.savefig('CIII_vs_age.pdf')

plt.figure()
ax2 = df1.plot(x='metallicity', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='blue')
df2.plot(x='metallicity', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='red', ax=ax2)
df3.plot(x='metallicity', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='green', ax=ax2)
plt.savefig('CIII_vs_Z.pdf')

# Do the same, for Lyman apha
ax = dfLya.plot(x='age', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='blue', ylim=(-100,0))
plt.savefig('Lya_vs_age.pdf')
plt.figure()
ax2 = dfLya.plot(x='metallicity', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='blue',ylim=(-100,0))
plt.savefig('Lya_vs_Z.pdf')
plt.show()
