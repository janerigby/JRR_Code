from os.path import expanduser, basename
import glob
import pandas
import numpy as np
import astropy.convolution
import jrr
from matplotlib import pyplot as plt

def check_ESI_fluxing(thisdir) :
    myfiles = ['s1723_arc_a_esi.txt', 's1723_arc_b_esi.txt', 's1723_side_a_esi.txt', 's1723_side_b_esi.txt', 's1723_center_a_esi.txt', 's1723_center_b_esi.txt', 's1723_counter_a_esi.txt', 's1723_counter_b_esi.txt']
    for thisfile in myfiles :
        print "DEBUG", thisdir, thisfile
        df = pandas.read_table(thisdir + thisfile, delim_whitespace=True, comment="#")
        boxcar = 201
        df['smooth'] = df['flam'].rolling(window=boxcar,center=False).median()
        #smooth = astropy.convolution.convolve(df['flam'].as_matrix(), np.ones((boxcar,))/boxcar, boundary='extend', fill_value=np.nan) # boxcar smooth
        scaleby = df[df['obswave'].between(7000,7200)]['flam'].median()
        print "scaleby", scaleby
        df['smoothscaled'] = df['smooth'] / scaleby
        plt.plot(df['obswave'], df['smoothscaled'], label=thisfile)
    plt.ylim(-3,10) ; plt.xlim(4000, 1.E4)
    return(myfiles)

#### Directories and filenames
home = expanduser("~")
wdir = home + '/Dropbox/SGAS-shared/Grism_S1723/'  
file_ESI  = home + '/Dropbox/MagE_atlas/Contrib/ESI_spectra/2016Aug/s1723_ESI_JRR_sum.txt'
file_MMT  = wdir + 'MMT_BlueChannel/spec_s1723_14b_final_F.txt'  # Updated Jan 23 2017
file_GNIRS = wdir + 'GNIRS/s1723_gnirs_residual_subtracted_spectrum_jrr.csv'
file_G102 = wdir + 'WFC3_fit_1Dspec/FULL_G102_coadded.dat' 
file_G141 = wdir + 'WFC3_fit_1Dspec/FULL_G141_coadded.dat'

#### Read ESI spectra
# Units are:  wave: barycentric corrected vacuum Angstroms;  flambda in erg/cm2/s/angstrom
sp_ESI = pandas.read_table(file_ESI, delim_whitespace=True, comment="#")
sp_ESI.rename(columns={'fsum_jrr' : 'flam'}, inplace=True)
sp_ESI.rename(columns={'fsum_u' : 'flam_u'}, inplace=True)
sp_ESI['flam_u'] = pandas.to_numeric(sp_ESI['flam_u'], errors='coerce') # convert to float64

#### Read WFC3 spectra (w continuua, both grisms)
names = ('wave', 'flam', 'flam_u', 'cont', 'flam_contsub')  # assumed flam. **Check w Michael
sp_G102 = pandas.read_table(file_G102, delim_whitespace=True, comment="#", names=names)
sp_G141 = pandas.read_table(file_G141, delim_whitespace=True, comment="#", names=names)

#### Read MMT Blue Channel spectrum. # wave in vacuum Ang, flam in 1E-17 erg/s/cm^2/A
names=('wave', 'flam', 'flam_u')  
sp_MMT = pandas.read_table(file_MMT, delim_whitespace=True, comment="#", names=names)
sp_MMT['flam'] *= 1.0E-17   # flam was in units of 1E-17
sp_MMT['flam_u'] *= 1.0E-17
#sp_MMT.replace([np.inf, -np.inf], np.nan, inplace=True)

#### Read GNIRS spectrum
sp_GNIRS = pandas.read_csv(file_GNIRS, comment="#")  # This is in counts.  NOT FLUXED

# Calculate some scaling factors, and scale the spectra
# Scale everything to G102
wavcut = np.array((8000., 9500., 4700., 5150.))
f_G102A = sp_G102.loc[sp_G102['wave'].between(wavcut[0], wavcut[1])]['flam'].median()
f_ESIA  = sp_ESI.loc[sp_ESI['wave'].between(wavcut[0], wavcut[1])]['flam'].median()
sp_ESI['flam_cor']   = sp_ESI['flam']   * f_G102A / f_ESIA 
sp_ESI['flam_u_cor'] = sp_ESI['flam_u'] * f_G102A / f_ESIA 

f_ESIB =  sp_ESI.loc[sp_ESI['wave'].between(wavcut[2], wavcut[3])]['flam_cor'].median()  # Scale MMT from already-scaled ESI
f_MMTB =  sp_MMT.loc[sp_MMT['wave'].between(wavcut[2], wavcut[3])]['flam'].median()
sp_MMT['flam_cor']   = sp_MMT['flam']   * f_ESIB / f_MMTB
sp_MMT['flam_u_cor'] = sp_MMT['flam_u'] * f_ESIB / f_MMTB

# Report how I am scaling the fluxes
print "Scaled ESI spectrum to match G102, from median flam in range", wavcut[0], "to", wavcut[1]
print "   by factor",  f_G102A / f_ESIA
print "Scaled MMT BC spectrum to match scaled ESI, from median flam in range", wavcut[2], 'to', wavcut[3]
print "   by factor", f_ESIB / f_MMTB
# ESI flux was high b/c I was summing multiple exposures. See MagE_atlas/Contrib/ESI_spectra/2016Aug/readme.txt

# Plot the scaled spectra
plt.clf()
plt.plot(sp_ESI['wave'],  sp_ESI['flam_cor'],    color='blue', label='Keck ESI')
plt.plot(sp_ESI['wave'],  sp_ESI['flam_u_cor'],  color='lightblue', label='_nolegend_', linewidth=0.2)
plt.plot(sp_MMT['wave'],  sp_MMT['flam_cor'],    color='green', Label='MMT Blue Channel')
plt.plot(sp_MMT['wave'],  sp_MMT['flam_u_cor'],  color='lightgreen', label='_nolegend_', linewidth=0.2)
plt.plot(sp_G102['wave'], sp_G102['flam'],   color='red', label='WFC3 G102')
plt.plot(sp_G102['wave'], sp_G102['flam_u'], color='red',label='_nolegend_', linewidth=0.2)
plt.plot(sp_G141['wave'], sp_G141['flam'],   color='purple', label='WFC3 G141')
plt.plot(sp_G141['wave'], sp_G141['flam_u'], color='purple', label='_nolegend_', linewidth=0.2)
plt.legend()
plt.ylim(0,1E-16)
plt.xlim(3200,16000)
plt.xlabel("Vacuum wavelength (Angstroms)")
plt.ylabel(r'$f_{\lambda}$ ($erg$ $s^{-1}$ $cm^{-2}$ $\AA^{-1}$)')
plt.show()


### Check the fluxing of the ESI sum
plt.clf()
test = check_ESI_fluxing(home + '/Dropbox/MagE_atlas/Contrib/ESI_spectra/2016Aug/')
boxcar = 21
sp_ESI['smooth'] = sp_ESI['flam_cor'].rolling(window=boxcar,center=False).median()
#sp_ESI['smooth'] = astropy.convolution.convolve(sp_ESI['flam_cor'].as_matrix(), np.ones((boxcar,))/boxcar, boundary='extend', fill_value=np.nan) # boxcar smooth
scaleby = sp_ESI[sp_ESI['wave'].between(7000,7200)]['smooth'].median()
print "scaleby", scaleby
plt.plot(sp_ESI['wave'], sp_ESI['smooth']/scaleby,  color='black', label='JRR SUM', linewidth=1.5)
plt.legend()
plt.show()
print "***CAUTION: something is deeply wrong with some of these input spectra shortward of 5000A.  Have written to ask Ayan"
