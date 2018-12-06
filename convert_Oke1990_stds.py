''' The standard star files on /lib/onedstds/oke1990/ are not really compatible with the iraf
task standard. So, I am converting them.  jrigby, Dec 2018.'''

import jrr
import numpy as np
import pandas
import  matplotlib.pyplot  as plt

plt.ion()
starfile = "/Users/jrrigby1/Ureka-x/Ureka/variants/common/iraf/gmisc/lib/onedstds/oke1990/gd108.dat"
outfile  = "gd108_try.dat"
colnames = ('wave', 'mag', 'bp')
df = pandas.read_table(starfile, names=colnames, delim_whitespace=True, comment="#")
df['flux'] = 10**(-0.4 * df['mag'])


startwave = 3200.
stopwave  = 9000.
bandpass = 100.
df2 = pandas.DataFrame(data=pandas.Series(np.arange(startwave, stopwave, bandpass)), columns=('wave',))
df2['newf'] = jrr.spec.rebin_spec_new(df['wave'], df['flux'], df2['wave'], fill='extrapolate')
df2['bp'] = bandpass

ax = df.plot(x='wave', y='flux', kind='scatter', color='blue')
df2.plot(x='wave', y='newf', ax=ax, kind='scatter', color='red')
plt.xlim(3100,8500)
plt.ylim(0,1E-5)
plt.show()

df2.to_csv(outfile, sep='\t', index=False, header=False)
