from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import pandas
# used http://www.astropy.org/astropy-tutorials/UVES.html as a guide

infile = "KBSS-LM1.uv.fnu.fits"
uncert_file = "KBSS-LM1.uv.fnu.sig.fits"

sp = fits.open(infile)
sp.info()
header = sp[0].header
wcs = WCS(header)
wcs2 = wcs.dropaxis(1)  # Kill the WCS's dummy 2nd dimension
index = np.arange(header['NAXIS1'])
temp =  (np.array(wcs2.wcs_pix2world(index, 0))).T
wavelength = (10**temp)[:,0]
fnu = sp[0].data

# Get the uncertainty
sp2 = fits.open(uncert_file)
fnu_u = sp2[0].data

# Make a pandas data frame
foo = np.array((wavelength, fnu, fnu_u))
df = pandas.DataFrame(foo.T, columns=("wave", "fnu", "fnu_u"))
df.head()
