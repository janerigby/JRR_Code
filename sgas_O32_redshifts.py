'''Just one last hack for the day:  the redshifts where we can both OII 3727 and OIII 5007 in a single IFU spectrum are  z=1.60--2.77 for the G140[M or H]/F100LP grating/filter combo, and z=3.45--5.32 for the G235[M or H]/F170LP.  And z=0.61--9.5 for the prism/clear combo.'''
 
import pandas
import numpy as np
infile = 'jwst_grating_wavelengths.txt'
df = pandas.read_csv(infile, delim_whitespace=True, names=('grating', 'waveon', 'waveoff'), comment='#')
np.array((0.3727, 0.5008))
#waves = np.array((0.3727, 0.5008))  # O2 and O3
waves = np.array((0.4861, 0.6563)) # Hbeta to Halha
waves = np.array((0.375, 0.6750)) # Hbeta to Halha
waves = np.array((0.3725, 0.9080)) # OII to SII 

df['z1'] = df.waveon/waves[0] - 1
df['z2'] = df.waveoff/waves[1] - 1
