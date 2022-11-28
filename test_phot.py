import jrr
from astropy.io import fits
import numpy as np
import pyregion
import matplotlib.pyplot as plt
import timeit

indir = '/Users/jrrigby1/SCIENCE/Lensed-LBGs/Planck_Arc/JR_narrowband/'
infile = 'OII_contsub_convolved/OII_contsub_convolved.fits'
regionfile = 'box_for_median_imcoords.reg'

print("testing method 1")
phot1 = jrr.phot.pyregion_photometry(indir + infile, indir + regionfile)

print("testing method 2")
phot2 = jrr.phot.pyregion_photometry_ma(indir + infile, indir + regionfile)

