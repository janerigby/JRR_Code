import pandas
import matplotlib.pyplot as plt
import numpy as np
import jrr


# Get results from mage_calc_galprops.py
galprops =  pandas.read_table("mage_galproperties.txt", delim_whitespace=True, comment="#")
plt.ion()

hasZneb   = galprops.loc[galprops['sigZ']>0]
uplimZneb = galprops.loc[galprops['sigZ']==-99]

plt.scatter(hasZneb['Zneb'], hasZneb['Z_S99'], color='blue')
plt.errorbar(hasZneb['Zneb'], hasZneb['Z_S99'], xerr=hasZneb['sigZ'], ls='none')
plt.scatter(uplimZneb['Zneb'], uplimZneb['Z_S99'], color='red')
plt.errorbar(uplimZneb['Zneb'], uplimZneb['Z_S99'], xerr=0.05, ls='none', xuplims=np.ones_like(uplimZneb['Z_S99']))
jrr.plot.annotate_from_dataframe(galprops, xcol='Zneb', ycol='Z_S99', text='label')

plt.xlabel('Nebular Z (fraction of solar)')
plt.ylabel('S99 Z (fraction of solar)')
plt.xlim(0,1)
plt.ylim(0,1)
plt.plot((0,1), (0,1))
plt.show()

