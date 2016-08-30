import jrr
import os
import sys
import re
from matplotlib import pyplot as plt
import pandas
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
plt.ion()
mage_mode = "reduction"
#mage_mode = "released"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
 
# Testing, does the autocont for a single spectrum look ok?  Yes
(st, resoln, dresoln, LL, z_systemic)=jrr.mage.wrap_open_spectrum('S1226+2152', mage_mode)
plt.plot(st.wave, st.fnu, color='b')
plt.plot(st.wave, st.fnu_autocont, color='y')
plt.scatter(st.wave, st.linemask, 1)
plt.xlim(4000,7000)
plt.ylim(0,1E-28)
plt.clf()

# again, for an S99 spectrum.
linelist = line_path + "stacked.linelist"
(S99, S99LL) = jrr.mage.open_S99_spectrum("S1226+2152", 0.0)
plt.plot(S99.wave, S99.fnu, color='b')
plt.plot(S99.wave, S99.fnu_autocont, color='y')
plt.scatter(S99.wave, S99.linemask, 1)
#plt.xlim(4000,7000)
plt.ylim(0,2)


 #(st) = jrr.mage.open_stacked_spectrum(mage_mode)
