from __future__ import print_function
from builtins import range
from astropy.io import fits
import  matplotlib.pyplot as plt
import pandas
import numpy as np
import jrr
from astropy.table import Table
import astropy
from os.path import expanduser

homedir = expanduser("~")
workdir = homedir + "/WORK/S2222_QSO_intervening/all_1D_spectra/"
#spectra = ('merged_spectra/manual_combspec_sgas2222_qsoA_N.fits', 'merged_spectra/manual_combspec_sgas2222_qsoB_N.fits', "sgas2222_15BQ27m02_21_F.fits")
spectra = ('merged_spectra/manual_combspec_sgas2222_qsoA_F.fits', 'merged_spectra/manual_combspec_sgas2222_qsoB_N.fits', "sgas2222_15BQ27m02_21_F.fits")

label = ['QSO_A', 'QSO_B', 'QSO_C']
colors= ('k', 'b', 'g')

def auto_fit_cont(sp, LL, zz, vmask=500, boxcar=1001, flag_lines=True) : # Automatically fit the continuum
    ''' Automatically fits a smooth continuum to a spectrum.
     Inputs:  sp,  a Pandas data frame containing the spectra, opened by mage.open_spectrum or similar
              LL,  a Pandas data frame containing the linelist, opened by mage.get_linelist(linelist) or similar
              zz,  the systemic redshift. Used to set window threshold around lines to mask
              vmask (Optional), velocity around which to mask lines, km/s. 
              boxcar (Optional), the size of the boxcar smoothing window, in pixels
              flag_lines (Optional) flag, should the fit mask out known spectral features? Can turn off for debugging
    '''
    # First, mask out big skylines. Done by mage.flag_skylines, which is called by mage.open_spectrum
    # Second, flag regions with crazy high uncertainty.  done in flag_huge_uncert, called by mage.open_spectrum
    # Third, mask out regions near lines.  Flagged in sp.linemask
    if flag_lines : flag_near_lines(sp, LL, zz, vmask)
    # Populate flam_autocont with flam, unless pixel is bad or has a spectral feature, in which case it stays nan.
    temp_flam = sp.flam.copy(deep=True)
    temp_flam.loc[(sp.badmask | sp.linemask).astype(np.bool)] = np.nan
    # Last, smooth with a boxcar
    smooth1 = astropy.convolution.convolve(np.array(temp_flam), np.ones((boxcar,))/boxcar, boundary='fill', fill_value=np.nan)
    sp['flam_autocont'] = pandas.Series(smooth1)  # Write the smooth continuum back to data frame
    return(0)

    
fig = plt.figure(figsize=(8,10))
zz = 2.2987981  # redshift of intervening absorber.  Measured below as np.median(zzlist)
name_line = ("Mg II 2796,2803", "Fe II 2600")
wav_line = (2800., 2600.)
win = 30.0 #Angstroms, plotting window
Ncol = 1  #columns of subplot
Nrow = len(wav_line)  # rows of subplot

# Read in the spectra to a dataframe df_all
df = {}  # data frame to use later
for ii, thisfile in enumerate(spectra) :
    # This deals with the weird inflg=2 x_readspec formatting.  Wavelength in hdu2 extension
    hdulist = fits.open(workdir + thisfile)
    flam  = hdulist[0].data
    flam_u = hdulist[1].data
    wave = hdulist[2].data
    rest_wave = wave / (1.+zz)  # this is rest-frame of the intervening absorber, not the QSO
    badmask = np.zeros_like(flam).astype(np.bool)
    linemask = np.zeros_like(flam).astype(np.bool)

    temptab = Table( (wave,rest_wave, flam,flam_u, badmask, linemask), names=('wave', 'rest_wave','flam','flam_u', 'badmask', 'linemask'))
    df[label[ii]] = temptab.to_pandas()  #dictionary of pandas data frames.
    #df.sort_values(by='wave', inplace=True)
    #print df[label[ii]].head()
    auto_fit_cont(df[label[ii]], None, 0.0, vmask=500, boxcar=101, flag_lines=False)
df_all = pandas.concat(df)

# Plot the results
for jj, thisline in enumerate(wav_line) :  # Loop over spectral features
    ax = fig.add_subplot(Nrow, Ncol, 1+jj)
    for ii, thisfile in enumerate(spectra) :
        plt.step(df[label[ii]].rest_wave, df[label[ii]].flam/df[label[ii]].flam_autocont, color=colors[ii])
        plt.step(df[label[ii]].rest_wave, df[label[ii]].flam_u/df[label[ii]].flam_autocont, color=colors[ii])
    plt.xlim( thisline - win,  thisline + win)
    plt.plot((0.,1.E4), (1.,1.), color='k')
    plt.annotate(name_line[jj], xy=(thisline,1.37), color="black", fontsize=18)
    plt.ylim(0.0,1.5)
    plt.xlabel(r'rest wavelength ($\AA$)', fontsize=18)

plt.savefig("/Users/jrrigby1/WORK/S2222_QSO_intervening/s2222_intervening.pdf")
plt.show()

print("#obj   line    cuton  cutoff  EW_r(A)       EW_r_u(A)     redshift")
# Measure EWs for Mg II 2803 and FeII 2600.  Also, get rest wavelength
cuton  = (2801.6, 2598.5)
cutoff = (2807.7, 2603.3)
centers = (2803.531, 2600.1729)
names  = ("MgII2803", "FeII2600")
zzlist = []
EWar   = np.zeros(shape=(len(spectra), len(cuton)))
EWar_u = np.zeros(shape=(len(spectra), len(cuton)))
for jj in range(0, len(cuton)) :
    for ii, thisfile in enumerate(spectra) :
        cutout = df[label[ii]].rest_wave.between(cuton[jj], cutoff[jj])
        dispersion = df[label[ii]][cutout].rest_wave.diff().median()
        (EWar[ii][jj], EWar_u[ii][jj]) = jrr.spec.calc_EW(df[label[ii]][cutout].flam, df[label[ii]][cutout].flam_u, df[label[ii]][cutout].flam_autocont, df[label[ii]][cutout].flam_autocont*0.0, dispersion, zz)
        
        # Measure the redshift of line, from the mode (lowest point) of the absorption line
        cutout2 = df[label[ii]].wave.between(cuton[jj]*(1+zz), cutoff[jj]*(1+zz))
        idmin = df[label[ii]][cutout2].flam.idxmin()
        line_zz = df[label[ii]].iloc[idmin].wave / centers[jj] - 1.0
        zzlist.append(line_zz)
        print(label[ii], names[jj], cuton[jj], cutoff[jj], EWar[ii][jj], EWar_u[ii][jj], line_zz)

print("mean and std of lines:", np.mean(EWar, axis=0), np.std(EWar, axis=0))        
print("Median redshift and mad for fitted lines: ",  np.median(zzlist), jrr.util.mad(zzlist))

# Need to concatenate the data frames, so that I can deal with all 3 spectra at once.


