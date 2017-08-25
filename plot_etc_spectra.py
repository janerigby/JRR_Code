''' JWST reports the per-pixel SNR.  What we often want is the SNR of the integrated line flux.  This script 
measures the latter quantity.  jrigby, day before the ERS deadline'''

from astropy.io.fits import getdata
from astropy.io.fits import getdata
import pandas
import numpy as np
import glob
import re
from astropy.io import ascii, fits
from astropy.table import Table
import matplotlib.pyplot as plt
import jrr

thefiles =  glob.glob("*/*/lineplot/lineplot_sn.fits")
df = {}
for ii, thisfile in enumerate(thefiles):
    df[thisfile] = Table.read(thisfile).to_pandas()
    df[thisfile]['totalsnr'] = 0.0
    window = 31
#   window = 61
    df[thisfile]['totalsnr'] =  df[thisfile]['sn'].rolling(window=window, center=True).apply(jrr.util.add_in_quad)

#print "Comparing SGAS targets, medium versus high resolution grating"
#pair1= ("ETC_SGAS/calc37/lineplot/lineplot_sn.fits", "ETC_SGAS/calc58/lineplot/lineplot_sn.fits")
#pair1a= ("ETC_SGAS/calc58/lineplot/lineplot_sn.fits", "ETC_SGAS/calc58short/lineplot/lineplot_sn.fits")
#pair2= ("ETC_SGAS/calc35/lineplot/lineplot_sn.fits", "ETC_SGAS/calc59/lineplot/lineplot_sn.fits")
#pair3= ("ETC_SGAS/calc34/lineplot/lineplot_sn.fits", "ETC_SGAS/calc60/lineplot/lineplot_sn.fits")

print("Double-checking latest integration times for SGAS")
p1  = "ETC_SGASnew/calc37/lineplot/lineplot_sn.fits"
p2  = "ETC_SGASnew/calc35/lineplot/lineplot_sn.fits"
p3  = "ETC_SGASnew/calc34/lineplot/lineplot_sn.fits"
p4  = "ETC_SGASnew/calc19/lineplot/lineplot_sn.fits"

print "Comparing SPT targets, medium versus high resolution grating"
p5 = 'ETC_SPT/s2147_calc1/lineplot/lineplot_sn.fits'
p6 = 'ETC_SPT/s2134_calc1/lineplot/lineplot_sn.fits'

print "Calculating for Travis"
p7 =  'ETC_Travis/Travis1/lineplot/lineplot_sn.fits'
p8 =  'ETC_Travis/Travis25/lineplot/lineplot_sn.fits'
p9 =  'ETC_Travis/Travis27/lineplot/lineplot_sn.fits'
p10 =  'ETC_Travis/Travis22/lineplot/lineplot_sn.fits'
p11=  'ETC_Travis/Travis24/lineplot/lineplot_sn.fits'

title = ("SGAS1723 shorter Pa + PB", "S1226 Ha, HB", "S1723 Ha, HB", "S1226 MIRI Pa", "SPT2147", "SPT2134", "Travis1", "Travis25", "Travis27", "TravisPAH", "Travis [Ne III]")

for ii, thisone in enumerate((p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11)) :
    ax = df[thisone].plot(x='WAVELENGTH', y='sn', label=thisone, color='blue')
    df[thisone].plot(x='WAVELENGTH', y='totalsnr', label="totalsnr 0", ax=ax, color='black')
    df[thisone].plot.scatter(x='WAVELENGTH', y='sn', label=thisone, ax=ax, color='blue')
    plt.title(title[ii])
    plt.show()


# If doing pairs
'''for ii, pair in enumerate((pair1, pair1a, pair2, pair3, pair4, pair5, pair6, pair7, pair8)) :
    ax = df[pair[0]].plot(x='WAVELENGTH', y='sn', label=pair[0], color='blue')
    df[pair[0]].plot(x='WAVELENGTH', y='totalsnr', label="totalsnr 0", ax=ax, color='black')
    df[pair[0]].plot.scatter(x='WAVELENGTH', y='sn', label=pair[0], ax=ax, color='blue')
    df[pair[1]].plot(x='WAVELENGTH', y='sn', ax=ax, label=pair[1], color='red')
    df[pair[1]].plot(x='WAVELENGTH', y='totalsnr', ax=ax, label="totalsnr 1", color='orange')
    df[pair[1]].plot.scatter(x='WAVELENGTH', y='sn', ax=ax, label=pair[1], color='red')
    plt.title(title[ii])
    plt.show()
'''
