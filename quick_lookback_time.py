from astropy.cosmology import Planck18 as cosmo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 18})

plt.ion()
plt.clf()
zz = np.arange(0, 15, 0.1)
lookback =  cosmo.lookback_time(zz)
plt.plot(lookback, zz, lw=3)

pts2highlight = (1.0, 5.0)
plt.scatter(cosmo.lookback_time(pts2highlight), pts2highlight)

plt.show()
plt.ylabel('redshift', fontsize=20)
plt.xlabel('lookback time (Gyr)', fontsize=20)
plt.grid()
plt.ylim(0,15)
plt.xlim(0,13.8)
