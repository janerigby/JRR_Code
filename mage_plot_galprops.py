import pandas
import matplotlib.pyplot as plt
import numpy as np
import jrr

# Compare the derived nebular metallicities and stellar metallicities
# Cosmic Eye gets plotted twice, because it has 2 different nebular Zs


def add_other_Zneb_CosmicEye(df):
    starkZ = 8.6;  starkZu = 0.25 # starkZu is a guess
    altZneb =  10**(starkZ - 8.69)
    altZnebsig = 10**(starkZ + starkZu - 8.69) - 10**(starkZ  - 8.69)
    plt.scatter(altZneb, df.loc['Cosmic~Eye']['Z_S99'], color='b')
    plt.errorbar(altZneb, df.loc['Cosmic~Eye']['Z_S99'], xerr=altZnebsig, ls='none', color='k')
    print "DEBUGGING", altZneb, altZnebsig
    return(0)

# Get results from mage_calc_galprops.py
galprops =  pandas.read_table("mage_galproperties.txt", delim_whitespace=True, comment="#")
galprops.set_index('label', inplace=True, drop=False)
plt.ion()

hasZneb   = galprops.loc[galprops['sigZ']>0]
uplimZneb = galprops.loc[galprops['sigZ']==-99]

plt.scatter(hasZneb['Zneb'], hasZneb['Z_S99'], color='b')
plt.errorbar(hasZneb['Zneb'], hasZneb['Z_S99'], xerr=hasZneb['sigZ'], ls='none', color='k')
plt.scatter(uplimZneb['Zneb'], uplimZneb['Z_S99'], color='b')
plt.errorbar(uplimZneb['Zneb'], uplimZneb['Z_S99'], xerr=0.03, ls='none', xuplims=np.ones_like(uplimZneb['Z_S99']), color='k')
jrr.plot.annotate_from_dataframe(galprops, xcol='Zneb', ycol='Z_S99', text='label')
add_other_Zneb_CosmicEye(hasZneb)

plt.xlabel('Nebular Z (fraction of solar)')
plt.ylabel('S99 Z (fraction of solar)')
plt.xlim(0,1)
plt.ylim(0,1)
plt.plot((0,1), (0,1))
plt.show()

