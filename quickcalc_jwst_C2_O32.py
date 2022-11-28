# Quick inputs to the ETC, for doing IFU spectroscopy of lensed galaxies to find high O3/O2 clumps.  Noodling around before Cycle 2 proposals, jrigby, 
import extinction
from astropy.cosmology import Planck15 as cosmo
from astropy import units
from math import pi
import numpy as np


##### Change these parameters
Av = 1.0
Rv  = 4.05
zz = 3.0
SFR = 1.  # inuts are solar masses per year
magnification = 10. 

#####################
# SFR (Msol/yr) = 7.9E-42 * L(Halpha) (erg/s)  # Eqn 2 of Kennicutt 1998
SFR2LHa = 7.9E-42  

#  Compute the extinction
print("Flux diference for this reddening at O2, Hbeta, O3:")
flux = 1.0

# Get me some of that sweet sweet cosmology
flux_to_lum = 4 * pi * (cosmo.luminosity_distance(zz).cgs.value)**2     # Apply luminosity distance
print("Scale factor to convert flux to luminosity at zz", zz, flux_to_lum)  # in cgs units

flux_HB = SFR / SFR2LHa / 2.8 / flux_to_lum   # SFR to LHa, LHa to LHbeta,  LHbeta to f(Hbeta)
print("Here goes nothing, unextincted Hbeta line flux is", flux_HB)

flux_all_lines = np.array((86.2,  121.,  156.,   282., 849., 449.)) / 156.  # dereddened line fluxes for S1723, scaled to Hbeta

#                             O2       O2         Hbeta    O3       O3         Halpha our friends
wave_all_lines = np.array([3727.092, 3729.875, 4862.683, 4960.295, 5008.240, 6564.61])
flux_reddening_all_lines = extinction.apply(extinction.calzetti00(wave_all_lines, Av, Rv), flux)

predicted_fluxes = flux_HB * flux_all_lines * flux_reddening_all_lines  * magnification
print(predicted_fluxes)  # Predicted reddennnnnned fluxes for the emission lines
print( (1.+zz)*wave_all_lines / 1E4) # in microns you're welcome ETC

