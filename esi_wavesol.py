#! /usr/bin/env python
'''Quickly estimate wavelength of feature on ESI spectrum.
   Run on command line as:  >./esi_wavesol.py  order Yval [redshift]
   for example, ./esi_wavesol.py 7 2010.2 [3.14]
   Returns approx observed wavelength, and if redshift entered, approx rest wavelength'''
import numpy as np
import sys

#print "DEBUG, arguments received were: ", sys.argv, len(sys.argv)
if(len(sys.argv) == 3 or len(sys.argv) == 4):
    order =  str(sys.argv[1])  # input the ESI spectral order
    pixel = np.float64(sys.argv[2])   # input the Y value of a pixel in that order
    if len(sys.argv)== 4 :
        zz = np.float64(sys.argv[3])  # Optional, redshift
else:
    print "ERROR!  Two arguments required:  [spectral order]  and [Y pixel value]\nOptional third argument is redshift\n"
    sys.exit(1)

#Dictionary of ESI wavelength coefficients, from http://www2.keck.hawaii.edu/inst/esi/QuickRef.html#Wavelengths
ed = {}
ed['15'] = np.array((4077.46, 0.154482, -1.140e-6, -3.106e-10))
ed['14'] = np.array((4366.24, 0.165128, -2.0205e-6, 5.71e-10))
ed['13'] = np.array((4699.50, 0.179043, -1.912e-6, -8.44e-11))
ed['12'] = np.array((5088.55, 0.194456, -2.140e-6, 4.00e-11))
ed['11'] = np.array((5549.09, 0.212052, -2.365e-6, -1.23e-10))
ed['10'] = np.array((6101.46, 0.233675, -2.593e-6, -1.105e-10))
ed['9']  = np.array((6776.99, 0.259847, -2.826e-6, -1.90e-10))
ed['8']  = np.array((7621.60, 0.29266, -3.203e-6, -2.77e-10))
ed['7']  = np.array((8707.59, 0.334496, -3.6815e-6, -2.58e-10))
ed['6']  = np.array((10156,   0.39, -4.25e-6, 0))

def approx_wave_esi(order, Yval) :
    ''' arguments are order (6 through 15) and Y pixel coordinate 
    returns approx wavelength in Angstroms'''
    coeff = ed[order]
    wave = coeff[0] + coeff[1]*(Yval-2048) + coeff[2]*(Yval-2048)**2 + coeff[3]*(Yval-2048)**3
    return(wave)

wave =  approx_wave_esi(order, pixel)
print "observed wave: ", wave, "  "
if len(sys.argv)== 4 :
    rest_wave = wave / (1.0+zz)
    print "rest wave:     ", rest_wave


