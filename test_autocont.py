import jrr
import os
import sys
import re
from matplotlib import pyplot as plt
import pandas
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
#plt.ion()
mage_mode = "reduction"
#mage_mode = "released"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)

 # temporary, determine nice boxcar in rest wavelength 
labels = ['rcs0327-E', 'S0004-0103', 'S0108+0624',  'S0033+0242', 'S0900+2234',  'S0957+0509', 'S1050+0017',  'Horseshoe', 'S1226+2152', 'S1429+1202', 'S1458-0023', 'S1527+0652', 'S2111-0114', 'Cosmic~Eye']

labels = ['S1226+2152']  # Open a MagE spectrum
labels = ['SPT0310',]
for label in labels :
    (st, resoln, dresoln, LL, z_systemic)=jrr.mage.wrap_open_spectrum(label, mage_mode)

# Open its S99 fit
linelist = line_path + "stacked.linelist"
(S99, S99LL) = jrr.mage.open_S99_spectrum(labels[0], 0.0)

plt.plot(S99.wave, S99.fnu, color='b')
plt.plot(S99.wave, S99.data_fnu, color='g')
norm_regionA = (1267.0, 1276.0)
normalization = np.median(st.rest_fnu[st.rest_wave.between(*norm_regionA)])
plt.plot(st.rest_wave, st.rest_fnu/normalization, color='k')
plt.plot(st.rest_wave, st.rest_fnu_autocont/normalization, color='pink')
plt.plot(S99.wave, S99.fnu_autocont, color='y')
plt.xlim(800,2200)
plt.ylim(0,2)
plt.show()
