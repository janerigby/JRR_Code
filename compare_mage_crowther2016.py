from numpy import array, concatenate
import jrr
import os
import sys
import re
from numpy import zeros_like
from matplotlib import pyplot as plt
mage_mode = "reduction"
#mage_mode = "released"

(sp, LL) = jrr.mage.open_stacked_spectrum(mage_mode)
plt.plot(sp.restwave, sp.X_avg)
plt.plot(sp.restwave, sp.X_sigma)
plt.xlim(1150, 1700)
plt.ylim(0, 1.4)
jrr.mage.plot_linelist(LL) 
plt.show()
