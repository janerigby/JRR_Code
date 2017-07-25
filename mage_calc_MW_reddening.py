import jrr
import numpy as np
import query_argonaut
import astropy.io.ascii
from astropy.coordinates import SkyCoord
from astropy import units as u

mage_mode = 'reduction'

coords_file = "/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Mage/Papers/Mage-sample/SUBMITTED/tab_sample_v2.tex"
table = astropy.io.ascii.read(coords_file, format='aastex')
df = table.to_pandas()
df['DEC (J2000)']  = df['DEC (J2000)'].str.replace('$', '')
df['Source name']  = df['Source name'].str.replace('$', '')
df.fillna(value=0.0, inplace=True)  # drop the nans
print df.head()

c = SkyCoord(ra=df['RA (J2000)'] , dec=df['DEC (J2000)'], unit=(u.hourangle, u.deg))  # Convert from segidecimal to decimal
result = query_argonaut.query(c.ra.value.tolist(), c.dec.value.tolist(), coordsys='equ', mode='sfd')  # query only takes lists as input
df['EBV_Green2015'] = (np.array(list(result.values()))).flatten()
print df[['Source name', 'EBV_Green2015']]


planckarc = SkyCoord(ra="15:50:04.57", dec="-78:10:59.5", unit=(u.hourangle, u.deg))   # Coords for S1 in Table 1 of Dahle et al. 2016
result_planck = query_argonaut.query(planckarc.ra.value, planckarc.dec.value, coordsys='equ', mode='sfd')  # query only takes lists as input
print "E(B-V) for the Planck arc is", result_planck


# Then by-hand emacs-ing to get these EBV values in to spectra-filenames-redshifts.txt

