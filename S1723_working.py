from os.path import expanduser, basename
import subprocess 
import glob
import pandas
import numpy as np
import astropy.convolution
import jrr
from matplotlib import pyplot as plt


def smooth_the_noise(sp, win=21, colwave='wave', colf='flam_cor', colfu='flam_u_cor', outcol='flam_u_col_smooth') :
    sp[outcol] =  sp[colfu].rolling(window=win, center=True).median()
    return(0)

def flag_noise_peaks(sp, delta=0.15, colfu='flam_u', contmask='contmask') :                         
    maxtab, mintab = jrr.peakdet.peakdet(sp[colfu], delta)  # Find peaks.
    peak_ind =  [np.int(p[0]) for p in maxtab] # The maxima
    sp[contmask].iloc[peak_ind] = True
    return(peak_ind)

def wrap_fit_continuum(sp, LL, zz, boxcar, colwave='wave', colf='flam_cor', colfu='flam_u_cor', colcont='flamcor_autocont', new_way=False, label="") :
    (smooth1, smooth2) =  jrr.spec.fit_autocont(sp, LL, zz, boxcar=boxcar, colf=colf,  colcont=colcont, new_way=False)
    plt.plot(sp[colwave],  sp[colf],    color='green', label=label)
    plt.plot(sp[colwave],  sp[colfu],   color='lightgreen')
    plt.plot(sp[colwave],  sp['contmask']*sp[colf].median(),   color='yellow', label='masked')
    plt.plot(sp[colwave],  sp[colcont], color='k', label='Auto continuum fit')
    plt.ylim(sp[colf].median() * -0.1, sp[colf].median() * 5)
    plt.legend()
    return(smooth1, smooth2)

def check_ESI_fluxing(thisdir) :
    myfiles = ['s1723_arc_a_esi.txt', 's1723_arc_b_esi.txt', 's1723_side_a_esi.txt', 's1723_side_b_esi.txt', 's1723_center_a_esi.txt', 's1723_center_b_esi.txt']#, 's1723_counter_a_esi.txt', 's1723_counter_b_esi.txt']
    for thisfile in myfiles :
        #print "DEBUG", thisdir, thisfile
        df = pandas.read_table(thisdir + thisfile, delim_whitespace=True, comment="#")
        boxcar = 11
        df['smooth'] = df['flam'].rolling(window=boxcar,center=True).median()
        #smooth = astropy.convolution.convolve(df['flam'].as_matrix(), np.ones((boxcar,))/boxcar, boundary='extend', fill_value=np.nan) # boxcar smooth
        scaleby = df[df['obswave'].between(7000,7200)]['flam'].median()
        print "scaleby", scaleby
        df['smoothscaled'] = df['smooth'] / scaleby
        plt.plot(df['obswave'], df['smoothscaled'], label=thisfile)
    plt.ylim(-3,10) ; plt.xlim(4000, 1.E4)
    return(myfiles)

def load_linelists(linelistdir, zz) :
    LL_uv  = pandas.read_table(linelistdir + "rest_UV_emission_linelist.txt",      delim_whitespace=True, comment="#")
    LL_opt = pandas.read_table(linelistdir + "rest_optical_emission_linelist.txt", delim_whitespace=True, comment="#")
    (spec_path, line_path) = jrr.mage.getpath('released')
    (LL_temp, zz_notused) =  jrr.mage.get_linelist(line_path + 's1723.linelist')
    LL_uvabs = LL_temp.loc[LL_temp['type'].eq("ISM")]  # This already has redshift from .linelist.  may be out of synch****
    LL_uv['zz']  = zz  ;   LL_opt['zz'] = zz   # Load the redshift
    LL_uv.sort_values( by='restwav', inplace=True)
    LL_opt.sort_values(by='restwav', inplace=True)
    LL_uvabs.sort_values(by='restwav', inplace=True)
    LL = pandas.concat([LL_uv, LL_opt, LL_uvabs], ignore_index=True)  # Merge the two linelists
    LL.sort_values(by='restwav', inplace=True)  # Sort by wavelength
    LL.reset_index(drop=True, inplace=True)
    LL['obswav'] = LL['restwav'] * (1.0 + LL['zz'])
    LL['fake_wav'] = 0
    LL['fake_v'] = 0
    LL['vmask'] = 300.  # Dummy for now
    LL.drop_duplicates(subset=('restwav', 'lab1'), inplace=True)  # drop duplicate entries, w same rest wavelength, same ion label    
    return(LL)

##########  END OF FUNCTIONS ###############


#######  Housekeeping  ##############
plt.close("all")
figsize=(12,4)

#### Directories and filenames
home = expanduser("~")
wdir = home + '/Dropbox/Grism_S1723/'  
file_ESI  = home + '/Dropbox/MagE_atlas/Contrib/ESI_spectra/2016Aug/s1723_arc_ESI_JRR_sum.txt'  # currently just arc_a, arc_b
#file_ESI  = home + '/Dropbox/MagE_atlas/Contrib/ESI_spectra/2016Aug/s1723_arc_a_esi.txt'  #TEMP KLUDGE CHECK
file_MMT  = wdir + 'MMT_BlueChannel/spec_s1723_14b_final_F.txt'  # Updated Jan 23 2017
file_GNIRS = wdir + 'GNIRS/s1723_gnirs_residual_subtracted_spectrum_jrr.csv'
file_G102 = wdir + 'WFC3_fit_1Dspec/FULL_G102_coadded.dat' 
file_G141 = wdir + 'WFC3_fit_1Dspec/FULL_G141_coadded.dat'

#### Load the linelists, and set the redshift
zHa   = 1.329279 ;  zHa_u = 0.000085  # From GNIRS, see gnirs_writeup.txt
LL = load_linelists(wdir+'Linelists/', zHa)

#### Read ESI spectra
# Units are:  wave: barycentric corrected vacuum Angstroms;  flambda in erg/cm2/s/angstrom
sp_ESI = pandas.read_table(file_ESI, delim_whitespace=True, comment="#")
sp_ESI.rename(columns={'fsum_jrr' : 'flam'}, inplace=True)  # TEMP COMMENTED OUT ****
sp_ESI.rename(columns={'fsum_u' : 'flam_u'}, inplace=True) # TEMP COMMENTED OUT ****
#sp_ESI['wave'] = sp_ESI['obswave'] 
sp_ESI['flam_u'] = pandas.to_numeric(sp_ESI['flam_u'], errors='coerce') # convert to float64
sp_ESI['contmask'] = False
sp_ESI.loc[sp_ESI['wave'].between(7582.,7750.), 'contmask'] = True   # Mask telluric A-band
sp_ESI.loc[sp_ESI['wave'].between(6300.,6303.), 'contmask'] = True   # Mask sky line
sp_ESI.loc[sp_ESI['wave'].between(5047.,5057.), 'contmask'] = True   # Mask sky line
sp_ESI.loc[sp_ESI['wave'].between(6300.,6303.), 'contmask'] = True   # Mask sky line


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
sp_MMT['contmask'] = False
sp_MMT.loc[sp_MMT['wave'].between(4999.,5017.), 'contmask'] = True   # Mask these noisy regions from continuum fitting.
sp_MMT.loc[sp_MMT['wave'].gt(5200.), 'contmask'] = True   # Mask these noisy regions from continuum fitting.
peak_ind_MMT = flag_noise_peaks(sp_MMT, colfu='flam_u', delta=1E-17)


#### Read GNIRS spectrum
sp_GNIRS = pandas.read_csv(file_GNIRS, comment="#")  # This is in counts.  NOT FLUXED

# Calculate scaling factors, and scale the spectra
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

#### Flux the GNIRS spectrum
fig = plt.figure(4, figsize=figsize)
first_guess_cont = sp_GNIRS.loc[sp_GNIRS['wave'].between(1.5E4, 1.56E4)]['mean'].median()
guesspars = (100., 6564.61 *(1+zHa), 10., first_guess_cont)
(popt, fit) = jrr.spec.fit_quick_gaussian(sp_GNIRS.interpolate(), guesspars, colwave='wave', colf='mean')
quick_flux_Ha_GNIRS = jrr.spec.sum_of_gaussian(popt)  # In counts
flux_Ha_G141 = 200.0E-17  # from SDSSJ1723+3411_G141.fit
sp_GNIRS['flam'] = sp_GNIRS['mean'] * flux_Ha_G141 / quick_flux_Ha_GNIRS
sp_GNIRS['flam_u'] = sp_GNIRS['errinmean'] * flux_Ha_G141 / quick_flux_Ha_GNIRS
plt.plot(sp_GNIRS['wave'], sp_GNIRS['flam'], color='k', label="GNIRS flux")
plt.plot(sp_GNIRS['wave'], sp_GNIRS['flam_u'], color='grey', label="GNIRS uncert")
plt.legend()
sp_GNIRS_cutout = sp_GNIRS.loc[sp_GNIRS['wave'].between(1.5E4,1.55E4)]

# Report how I am scaling the fluxes
print "Scaled ESI spectrum to match G102, from median flam in range", wavcut[0], "to", wavcut[1]
print "   by factor",  f_G102A / f_ESIA
print "Scaled MMT BC spectrum to match scaled ESI, from median flam in range", wavcut[2], 'to', wavcut[3]
print "   by factor", f_ESIB / f_MMTB
# ESI flux was high b/c I was summing multiple exposures. See MagE_atlas/Contrib/ESI_spectra/2016Aug/readme.txt

print "Scaled GNIRS spectrum by ratio of Halpha flux, by factor", flux_Ha_G141 / quick_flux_Ha_GNIRS

# De-redshift the spectra  (should loop this...)
jrr.spec.convert2restframe_df(sp_MMT, zHa, units='flam', colwave='wave', colf='flam_cor', colf_u='flam_u_cor')
jrr.spec.convert2restframe_df(sp_ESI, zHa, units='flam', colwave='wave', colf='flam_cor', colf_u='flam_u_cor')
#jrr.spec.convert2restframe_df(sp_G102, zHa, units='flam', colwave='wave', colf='flam', colf_u='flam_u_cor')
#jrr.spec.convert2restframe_df(sp_G141, zHa, units='flam', colwave='wave', colf='flam', colf_u='flam_u_cor')

# Plot the scaled spectra
fig = plt.figure(0, figsize=figsize)
plt.plot(sp_ESI['wave'],  sp_ESI['flam_cor'],    color='green', label='Keck ESI')
plt.plot(sp_ESI['wave'],  sp_ESI['flam_u_cor'],  color='lightgreen', label='_nolegend_', linewidth=0.2)
plt.plot(sp_MMT['wave'],  sp_MMT['flam_cor'],    color='blue', Label='MMT Blue Channel')
plt.plot(sp_MMT['wave'],  sp_MMT['flam_u_cor'],  color='lightblue', label='_nolegend_', linewidth=0.2)
plt.plot(sp_GNIRS_cutout['wave'], sp_GNIRS_cutout['flam'], color='orange', label="GNIRS")
plt.plot(sp_GNIRS_cutout['wave'], sp_GNIRS_cutout['flam_u'], color='lightgrey', label="_nolegend_")
plt.plot(sp_G102['wave'], sp_G102['flam'],   color='red', label='WFC3 G102')
plt.plot(sp_G102['wave'], sp_G102['flam_u'], color='red',label='_nolegend_', linewidth=0.2)
plt.plot(sp_G141['wave'], sp_G141['flam'],   color='purple', label='WFC3 G141')
plt.plot(sp_G141['wave'], sp_G141['flam_u'], color='purple', label='_nolegend_', linewidth=0.2)
plt.legend()
plt.ylim(0,1E-16)
plt.xlim(3200,16000)
plt.xlabel("Vacuum wavelength (Angstroms)")
plt.ylabel(r'$f_{\lambda}$ ($erg$ $s^{-1}$ $cm^{-2}$ $\AA^{-1}$)')
plt.title("Have scaled continuua.")
#plt.show()


### Check the fluxing of the ESI 
fig = plt.figure(1, figsize=figsize)
test = check_ESI_fluxing(home + '/Dropbox/MagE_atlas/Contrib/ESI_spectra/2016Aug/')
boxcar = 5
sp_ESI['smooth'] = sp_ESI['flam_cor'].rolling(window=boxcar,center=False).median()
#sp_ESI['smooth'] = astropy.convolution.convolve(sp_ESI['flam_cor'].as_matrix(), np.ones((boxcar,))/boxcar, boundary='extend', fill_value=np.nan) # boxcar smooth
scaleby = sp_ESI[sp_ESI['wave'].between(7000,7200)]['smooth'].median()
print "scaleby", scaleby
plt.plot(sp_ESI['wave'], sp_ESI['smooth']/scaleby,  color='black', label='JRR arc sum', linewidth=1.5)
plt.legend()
plt.title("Diagnosing problems with ESI fluxing, ESI SUM")
print "***CAUTION: something is deeply wrong with some of these input spectra shortward of 5000A.  Have written to ask Ayan"

### Fit nice continuum for MMT.  The boxcar # is arbitrary; can make it more physical later
fig = plt.figure(2, figsize=figsize)
plt.title("Fit continuum for MMT bluechannel")
(smooth1, smooth2) = wrap_fit_continuum(sp_MMT, LL, zHa, boxcar=151, colf='flam_cor',  colcont='flamcor_autocont', label="MMT Blue Channel")

### First stab at a cont fit for ESI.  Need to diagnose problems w ESI extraction, then come back to this
peak_ind = flag_noise_peaks(sp_ESI, colfu='flam_u_cor', delta=0.3E-17)
fig = plt.figure(3, figsize=figsize)
plt.title("First attempt to fit continuum to ESI.  Not done yet")
(smooth1, smooth2) = wrap_fit_continuum(sp_ESI, LL, zHa, boxcar=251, colf='flam_cor',  colcont='flamcor_autocont', label="ESI with problems at blue end")

plt.show()  # Show all plots at once, each in a separate window

print "Fitting line fluxes.  MMT first"
sp_MMT.to_csv("s1723_MMT_wcont.txt", sep='\t')
# Ayan first runs a translator to get it into format EW_fitter.py likes
subprocess.call(home + "/Python/AYAN_Code/convert_spectra_format.py --inpath ~/Dropbox/Grism_S1723/JRR_Working/ --infile s1723_MMT_wcont.txt --flamconst 1. --flamcol flam_cor --flamucol flam_u_cor --wavecol wave --flamcontcol flamcor_autocont --z 1.32952 --zu 4e-4")
#
subprocess.call(home + "/Python/AYAN_Code/EW_fitter.py --short s1723_MMT_wcont_new-format --spec_list_file ~/Dropbox/Grism_S1723/JRR_Working/other-spectra-filenames-redshifts.txt --silent --path ~/Dropbox/Grism_S1723/JRR_Working/ --linelistpath ~/Dropbox/Grism_S1723/JRR_Working/ --fout s1723_MMT_emission_measured.out --useflamcont flam_cont --savepdf --hide")

#execfile( home + "/Python/AYAN_Code/convert_spectra_format.py    --inpath ./  --infile s1723_MMT_wcont.txt --flamconst# 1.0  --flamcol flam_cor --flamucol flam_u_cor  --wavecol wave   --flamcontcol flamcor_autocont --z 1.329279 --zu 0.000085")

#run ~/Python/AYAN_Code/EW_fitter.py --short s1723_MMT_wcont_new-format  --spec_list_file ./other-spectra-filenames-redshifts.txt --silent --hide --savepdf --fout s1723_MMT_measuredlines.out --useflamcont  flam_cont --linelistpath ./


