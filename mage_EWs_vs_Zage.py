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

# Read John's file with ages, metallicities
ageZ_file = homedir + "/Dropbox/MagE_atlas/Contrib/S99/sb99-parameters-sigfigs.txt"
dftZ = pandas.read_table(ageZ_file, delim_whitespace=True, comment="#")
dftZ.set_index('name', inplace=True)

#First, do the sum of CIII 1907, 1909 as a special case.
line_lab = ("CIII]1906", "CIII]1908")
subset1 = dfL[dfL['line_lab'].eq(line_lab[0])]
subset1.set_index('label', inplace=True)
subset2 = dfL[dfL['line_lab'].eq(line_lab[1])]
subset2.set_index('label', inplace=True)
sum_EW   = subset1['EWr_fit'] + subset2['EWr_fit']
sum_EW_u = np.sqrt(subset1['EWr_fit_u']**2 + subset2['EWr_fit_u']**2)
EWsum = pandas.concat([sum_EW, sum_EW_u], axis=1).reset_index()
EWsum.set_index('label', inplace=True)
df1 = pandas.concat([dftZ, subset1], axis=1)
df2 = pandas.concat([dftZ, subset2], axis=1)
df3 = pandas.concat([dftZ, EWsum], axis=1)
#ax = df1.plot(x='age', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='blue')
#df2.plot(x='age', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='red', ax=ax)
df3.plot(x='age', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='green',  title='CIII sum')
plt.savefig('CIIIsum_vs_age.pdf')
#
#ax2 = df1.plot(x='metallicity', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='blue')
#df2.plot(x='metallicity', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='red', ax=ax2)
df3.plot(x='metallicity', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='green', title='CIII sum')
plt.savefig('CIIIsum_vs_Z.pdf')

plot_debug = False
if plot_debug :
    ax = subset1.plot(x='EWr_fit', y='EWr_sum', kind='scatter', xerr='EWr_fit_u', yerr='EWr_sum_u', color='blue') # native pandas plotting.  pandas.DataFrame.plot()
    #subset2.plot(x='EWr_fit', y='EWr_sum', kind='scatter', xerr='EWr_fit_u', yerr='EWr_sum_u', xlim=(-2,2), ylim=(-2,2), ax=ax, color='red') 
    plt.plot( (-2,2), (-2,2), color='green')
    plt.savefig('EWrfit_EWrsum.pdf')
    plt.show()
    # Why don't the EWr_fit and EWr_sum agree?  Have asked Ayan.  Maybe b/c _sum is affected by blend?

# Now, do a bunch of other lines as a loop
line_labs = ("CIII]1906", "CIII]1908", "Ly-alpha", "HeII1640", "CII2325c", "OII2470mid")
for line_lab in line_labs :
    subset = dfL[dfL['line_lab'].eq(line_lab)]
    subset.set_index('label', inplace=True)
    df = pandas.concat([dftZ, subset], axis=1)
    ax = df.plot(x='age', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='blue', title=line_lab)
    if line_lab == "Ly-alpha" : plt.ylim(-100,1)
    else : plt.ylim(-3,0.5)
    thepdf = line_lab + "_vs_age.pdf"
    plt.savefig(thepdf)
    plt.figure()
    ax2 = df.plot(x='metallicity', y='EWr_fit', yerr='EWr_fit_u', kind='scatter', color='blue', title=line_lab)
    if line_lab == "Ly-alpha" : plt.ylim(-100,1)
    else : plt.ylim(-3,0.5)
    thepdf = line_lab + "_vs_Z.pdf"
    plt.savefig(thepdf)

