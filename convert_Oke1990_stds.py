''' The standard star files on /lib/onedstds/oke1990/ are not really compatible with the iraf
task standard. So, I am converting them to 100A bandpasses.  jrigby, Dec 2018.'''

import jrr
import glob
from os.path import basename
import numpy as np
import pandas
import  matplotlib.pyplot  as plt


def mag2flux(mag) :
    return (10**(-0.4 * mag))

def flux2mag(flux) :
    return (-2.5 * np.log10(flux))

#plt.ion()
indir  = '/iraf/iraf/noao/lib/onedstds/oke1990/'
outdir = '/iraf/iraf/noao/lib/onedstds/oke1990jrr/'
#starfile = "/Users/jrrigby1/Ureka-x/Ureka/variants/common/iraf/gmisc/lib/onedstds/oke1990/gd108.dat"
#outfile  = "gd108_try.dat"

startwave = 3200.
stopwave  = 9000.
bandpass = 100.

filenames = [ basename(x) for x in glob.glob(indir + '*.dat')]
for filename in filenames :
    colnames = ('wave', 'mag', 'bp')
    df = pandas.read_table(indir + filename, names=colnames, delim_whitespace=True, comment="#")
    df['flux'] = mag2flux(df['mag'])
    df2 = pandas.DataFrame(data=pandas.Series(np.arange(startwave, stopwave, bandpass)), columns=('wave',))
    df2['flux'] = jrr.spec.rebin_spec_new(df['wave'], df['flux'], df2['wave'], fill='extrapolate')
    df2['mag'] = flux2mag(df2['flux'])
    df2['bp'] = bandpass
    ax = df.plot(x='wave', y='mag', kind='scatter', color='blue')
    df2.plot(x='wave', y='mag', ax=ax, kind='scatter', color='red')
    plt.xlim(3100,8500)
    plt.show()
    df2.drop(['flux'], axis=1, inplace=True)
    df2.to_csv(outdir + filename, sep='\t', index=False, header=False)
