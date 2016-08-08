from astropy.io import fits
data, header = fits.getdata("icwx05020_drz.fits", header=True)
import matplotlib.pyplot as plt
import numpy as np
# Make a quick cutout of the spectrum
x1 = 647
x2 = 784
y1 = 553
y2 = 559
cut3 = data[y1:y2,x1:x2]
plt.imshow(cut3, origin='lower')

bkg = data[y1+15:y2+15, x1:x2]
plt.clf()
plt.imshow(cut3 - bkg)
oneD = np.sum(cut3 - bkg, axis=0)   # subtract a rough background 
plt.plot(oneD*10)
plt.show()
