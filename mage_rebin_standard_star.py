import pandas
import numpy as np
import jrr
thefile = "/Users/jrrigby1/SCIENCE/Lensed-LBGs/Mage/Redux/Mage_Feb2020/1D/Final-Products/gd153.dat"
columns = ('wave', 'flux', 'width')
df = pandas.read_csv(thefile, delim_whitespace=True, skiprows=1, names=columns)

new_wave = np.arange(3000., 10000., 50)
new_flux = jrr.spec.rebin_spec_new(df['wave'], df['flux'], new_wave)
df2 = pandas.DataFrame({'wave' : new_wave,  'flux' : new_flux})
df2['width'] = 50.0
df2.to_csv("gd153_binned.dat", sep="\t", index=False)
