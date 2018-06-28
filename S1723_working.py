from os.path import expanduser, basename
import subprocess 
import glob
import jrr  # Moved some functions to jrr.grism
import pandas
import numpy as np
#import query_argonaut
import astropy.convolution
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator

''' This script processes the spatially-integrated spectra of SGAS J1723, so that we can get emission line
fluxes.  jrigby, 2017.  Now updating June 2018 to use new grizli-reduced grism spectra.'''


def annotate_header(header_ESI, header_MMT, header_GNIRS) :
    for header in (header_ESI, header_MMT, header_GNIRS) :
        header += ("#\n# Below is extra processing by S1723_working.py\n")
        header_ESI +=("# Corrected for MW reddening of E(B-V)=" + str(EBV) + "\n")
    return(header_ESI, header_MMT, header_GNIRS)


def MMT_barycentric_correction(sp) :  # Applying barycentric correction to wavelengths
    thistime =  Time('2014-05-05T09:21:00', format='isot', scale='utc')
    thisradec = SkyCoord("17:23:37.23", "+34:11:59.07", unit=(u.hourangle, u.deg), frame='icrs')
    mmt= EarthLocation.of_site('mmt') # Needs internet connection 
    barycor_vel = jrr.barycen.compute_barycentric_correction(thistime, thisradec, location=mmt)
    print "FYI, the barycentric correction factor for s1723 was", barycor_vel
    jrr.barycen.apply_barycentric_correction(sp, barycor_vel, colwav='oldwave', colwavnew='wave') 
    return(0)

def flag_noise_peaks(sp, delta=0.15, neighborpix=2, colfu='flam_u', contmask='contmask') :                         
    maxtab, mintab = jrr.peakdet.peakdet(sp[colfu], delta)  # Find peaks.
    peak_ind =  [np.int(p[0]) for p in maxtab] # The maxima
    for thispeak in peak_ind :
        sp.iloc[thispeak - neighborpix : thispeak + neighborpix][contmask] = True
    return(peak_ind)

def sum_spectraline_wflatcont(sp, linewave1, linewave2, contwave1, contwave2, colwave='wave', colflam='flam') :
    # This measures line flux via direct summation.  Used to apply relative flux scaling to the ESI and MMT spectra.
    # sp is a dataframe
    # linewave1, linewave2 are the begining and end wavelengths of the emission line to sum flux via direct summation
    # contwave1, contwave2 are the begining and end wavelengths of the region to get continuum, via median
    sp['disp'] = sp[colwave].diff()
    contlevel  = sp.loc[sp[colwave].between(contwave1, contwave2)][colflam].median()
    part1 = sp.loc[sp[colwave].between(linewave1, linewave2)][colflam].sum() * sp.loc[sp[colwave].between(linewave1, linewave2)]['disp'].median()
    part2 = (linewave2 -  linewave1) * contlevel
    lineflux = part1 - part2
    return(lineflux)

def adjust_plot() :
    plt.legend()
    x1 = 3200. ; x2 = 16900.
    plt.ylim(0,1E-16)
    plt.xlim(x1, x2)
    plt.xlabel(r"Vacuum wavelength ($\AA$)")
    plt.ylabel(r'$f_{\lambda}$ ($erg$ $s^{-1}$ $cm^{-2}$ $\AA^{-1}$)')
    #LL_notabs = LL.loc[LL['type'] != 'ISM']
    #jrr.mage.plot_linelist(LL_notabs, z_systemic=np.nan, restframe=False, velplot=False, line_center=np.nan) 
    ax2 = ax.twiny()
    ax2.set_xlim( x1/(1. + zHa), x2/(1. + zHa) )
    ax2.set_xlabel(r"rest-frame vacuum wavelength ($\AA$)")
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    fig.subplots_adjust(hspace=0)
    return(0)

##########  END OF FUNCTIONS ###############


#######  Housekeeping  ##############
plt.close("all")
figsize=(16,4)

#### Directories and filenames
home = expanduser("~")
wdir = home + '/Dropbox/Grism_S1723/'  
file_ESI  = home + '/Dropbox/MagE_atlas/Contrib/ESI_spectra/2016Aug/s1723_wholearc_ESI_JRR_sum.txt'   #  just arc_a, arc_b
file_MMT  = wdir + 'MMT_BlueChannel/spec_s1723_14b_final_F.txt'  # Updated Jan 23 2017
file_GNIRS = wdir + 'GNIRS/s1723_gnirs_residual_subtracted_spectrum_jrr.csv'
file_G102 = wdir + 'Grizli_redux/1Dsum/Wcont/sgas1723_1Dsum_bothroll_G102_wcontMWdr.txt' # New grisli redux, both rolls
file_G141 = wdir + 'Grizli_redux/1Dsum/Wcont/sgas1723_1Dsum_bothroll_G141_wcontMWdr.txt' # New grisli redux, both rolls

#### Load the linelists, and set the redshift
EBV = jrr.grism.get_MWreddening_S1723()         # Get the Milky Way reddening value
zHa   = 1.329279 ;  zHa_u = 0.000085  # From GNIRS, see gnirs_writeup.txt
zESI_prelim = 1.3333365               # Why is this discrepant??
LL = jrr.grism.load_linelists(wdir+'Linelists/', zHa)   #Need to check, z may be discrepant **

#### Read ESI spectra
# Units are:  wave: barycentric corrected vacuum Angstroms;  flambda in erg/cm2/s/angstrom
sp_ESI = pandas.read_table(file_ESI, delim_whitespace=True, comment="#")
header_ESI = jrr.util.read_header_from_file(file_ESI, comment="#")  # save the header
sp_ESI.rename(columns={'fsum_jrr' : 'flam'}, inplace=True)  
sp_ESI.rename(columns={'fsum_u' : 'flam_u'}, inplace=True) 
sp_ESI['flam_u'] = pandas.to_numeric(sp_ESI['flam_u'], errors='coerce') # convert to float64
jrr.spec.deredden_MW_extinction(sp_ESI, EBV, colwave='wave', colf='flam', colfu='flam_u')  # Correct for Milky Way reddening
sp_ESI['contmask'] = False
sp_ESI.loc[sp_ESI['wave'].between(7576.,7739.), 'contmask'] = True   # Mask telluric A-band
sp_ESI.loc[sp_ESI['wave'].between(5583.,5592.), 'contmask'] = True   # Mask sky line
sp_ESI.loc[sp_ESI['wave'].between(5611.,5623.), 'contmask'] = True   # Mask sky line
sp_ESI.loc[sp_ESI['wave'].between(6308.,6313.), 'contmask'] = True   # Mask sky line
sp_ESI.loc[sp_ESI['flam'] == 0.00, 'contmask'] = True   # Mask bad column, which Ayan replaced w zero

#### Read WFC3 spectra (w continuua, both grisms)  # MW dereddening was applied at continuum fitting stage.
sp_G102 = pandas.read_csv(file_G102, comment="#")
sp_G141 = pandas.read_csv(file_G141, comment="#")
grism_info_G102 = jrr.grism.get_grism_info("G102") 
grism_info_G141 = jrr.grism.get_grism_info("G141") 


#### Read MMT Blue Channel spectrum. # wave in vacuum Ang, flam in 1E-17 erg/s/cm^2/A
names=('oldwave', 'flam', 'flam_u')  
sp_MMT = pandas.read_table(file_MMT, delim_whitespace=True, comment="#", names=names)
header_MMT =("# MMT Blue Channel spectrum of SGAS J1723.\n")  # Nothing worth grabbing from original header
sp_MMT['flam']   *= 1.0E-17   # flam was in units of 1E-17
sp_MMT['flam_u'] *= 1.0E-17
#sp_MMT.replace([np.inf, -np.inf], np.nan, inplace=True)
MMT_barycentric_correction(sp_MMT) 
jrr.spec.deredden_MW_extinction(sp_MMT, EBV, colwave='wave', colf='flam', colfu='flam_u')  # Correct for Milky Way reddening
sp_MMT['contmask'] = False
sp_MMT.loc[sp_MMT['wave'].between(4999.,5017.), 'contmask'] = True   # Mask these noisy regions from continuum fitting.
sp_MMT.loc[sp_MMT['wave'].gt(5200.), 'contmask'] = True   # Mask these noisy regions from continuum fitting.
peak_ind_MMT = flag_noise_peaks(sp_MMT, colfu='flam_u', delta=1E-17)

#### Read GNIRS spectrum
sp_GNIRS = pandas.read_csv(file_GNIRS, comment="#")  # This is in counts.  NOT FLUXED
header_GNIRS = jrr.util.read_header_from_file(file_GNIRS, comment="#")  # save the header

# Update the headers
(header_ESI, header_MMT,header_GNIRS) = annotate_header(header_ESI, header_MMT, header_GNIRS)

# Scale the ESI spectrum to the HST flux scale, using the flux in the [O II] doublet.  Direct summation of flux, linear local continuua
flux_OII_G102 = (70.991 + 97.251) * 1E-17  # from 1Dsum/sgas1723_1Dsum_bothroll_G102_wcontMWdr_meth2.fitdf
rawflux_OII_ESI = sum_spectraline_wflatcont(sp_ESI, 8684., 8712., 8500., 8900.)

# Scale the MMT spectrum to the ESI flux, using the flux in the summed [C~III] doublet
rawflux_CIII_ESI = sum_spectraline_wflatcont(sp_ESI, 4443., 4456., 4400., 4490)
rawflux_CIII_MMT = sum_spectraline_wflatcont(sp_MMT, 4436., 4450., 4400., 4490)
# ESI flux in file_ESI should be 2x too high, b/c I summed 2 exposures.
sp_ESI['flam_cor']   = sp_ESI['flam']   *  flux_OII_G102/rawflux_OII_ESI 
sp_ESI['flam_u_cor'] = sp_ESI['flam_u'] *  flux_OII_G102/rawflux_OII_ESI
sp_MMT['flam_cor']   = sp_MMT['flam']   *  flux_OII_G102/rawflux_OII_ESI * rawflux_CIII_ESI/rawflux_CIII_MMT
sp_MMT['flam_u_cor'] = sp_MMT['flam_u'] *  flux_OII_G102/rawflux_OII_ESI * rawflux_CIII_ESI/rawflux_CIII_MMT

#### Flux the GNIRS spectrum
flux_Ha_G141 = 408.8E-17  #  kludge, hardcoded, from sgas1723_1Dsum_bothroll_G141_wcontMWdr_meth2.fitdf.  Updated 6/28/2018
first_guess_cont = sp_GNIRS.loc[sp_GNIRS['wave'].between(1.5E4, 1.56E4)]['mean'].median()
guesspars = (100., 6564.61 *(1+zHa), 10.)  # Fix the continuum, don't let curve_fit change it; it's drifting too high
(popt_GNIRS, fit_GNIRS) = jrr.spec.fit_gaussian_fixedcont(sp_GNIRS.interpolate(), guesspars, contlevel=first_guess_cont, colwave='wave', colf='mean')
quick_flux_Ha_GNIRS = jrr.spec.sum_of_gaussian(popt_GNIRS) # In counts
GNIRS_scaleflux = flux_Ha_G141 / quick_flux_Ha_GNIRS
sp_GNIRS['flam'] = sp_GNIRS['mean']        * GNIRS_scaleflux
sp_GNIRS['flam_u'] = sp_GNIRS['errinmean'] * GNIRS_scaleflux
#fig = plt.figure(4, figsize=figsize)
#plt.plot(sp_GNIRS['wave'], sp_GNIRS['flam'], color='k', label="GNIRS flux")
#plt.plot(sp_GNIRS['wave'], sp_GNIRS['flam_u'], color='grey', label="GNIRS uncert")
#plt.plot(sp_GNIRS['wave'], fit_GNIRS * GNIRS_scaleflux, color='blue', label='Ha fit')
#plt.legend()
sp_GNIRS_cutout = sp_GNIRS.loc[sp_GNIRS['wave'].between(1.5E4,1.55E4)]

# Report how scaled the fluxes
howscaled_ESI = "# Scaled flux in ESI so that sum of [O II] 3727+3729 in ESI matches flux of those lines in G102\n"
howscaled_MMT = "# Scaled flux in MMT so that sum of CIII in MMT matches sum in ESI\n"
#howscaled_ESI = "# flam_cor is scaled ESI flambda, to match WFC3 G102, from median flam in range " + str(wavcut[0]) + " to " + str(wavcut[1]) + "\n#    Flux was scaled by factor  " + str(f_G102A / f_ESIA) + "\n"
#howscaled_MMT = "# Scaled MMT BC spectrum to match scaled ESI, from median flam in range " + str(wavcut[2]) + ' to ' + str(wavcut[3]) + "\n#   Flux was scaled by factor " + str(f_ESIB / f_MMTB) + "\n"
howscaled_GNIRS = "# Scaled GNIRS spectrum by ratio of Halpha flux, by factor " + str(GNIRS_scaleflux)
header_ESI   += howscaled_ESI
header_MMT   += howscaled_MMT
header_GNIRS += howscaled_GNIRS 
print "How scaled spectra\n:", howscaled_ESI, "\n", howscaled_MMT, "\n", howscaled_GNIRS, "\n"

### Fit MMT continuum.  The boxcar # is arbitrary, seems to match
fig = plt.figure(2, figsize=figsize)
plt.title("Fit continuum for MMT bluechannel")
(smooth1, smooth2) = jrr.grism.wrap_fit_continuum(sp_MMT, LL, zHa, boxcar=151, colf='flam_cor',  colcont='flamcor_autocont', label="MMT Blue Channel")

### Fit ESI continuum
peak_ind = flag_noise_peaks(sp_ESI, colfu='flam_u_cor', delta=0.05E-17)
fig = plt.figure(3, figsize=figsize)
#plt.scatter(sp_ESI['wave'].iloc[peak_ind], sp_ESI['flam_cor'].iloc[peak_ind]) 
plt.title("Continuum fitting for ESI")
(smooth1, smooth2) = jrr.grism.wrap_fit_continuum(sp_ESI, LL, zHa, boxcar=351, colf='flam_cor',  colcont='flamcor_autocont', label="ESI")
header_ESI += ("# Applied smooth continuum \n")
plt.ylim(0,9E-17)

print "Writing flux-adjusted, continuum-fitted spectra for MMT and ESI."
sp_ESI.to_csv('temp1', sep='\t', na_rep='NaN')
sp_MMT.to_csv("temp2", sep='\t', na_rep='NaN')
jrr.util.put_header_on_file('temp1', header_ESI, "s1723_ESI_wcont.txt")
jrr.util.put_header_on_file('temp2', header_MMT, "s1723_MMT_wcont.txt")

# Plot the scaled spectra
ax = sp_ESI.plot(    x='wave', y='flam_cor', color='green',  label='Keck ESI', figsize=figsize)
sp_GNIRS_cutout.plot(x='wave', y='flam',     color='orange', label="GNIRS",     ax=ax)
sp_MMT.plot( x='wave', y='flam_cor',    color='blue',   label='MMT BC',    ax=ax)
sp_G102.plot(x='wave', y='flam',        color='red',    label='WFC3 G102', ax=ax)
sp_G141.plot(x='wave', y='flam',        color='purple', label='WFC3 G141', ax=ax)
sp_MMT.plot( x='wave', y='flam_u_cor',  color='lightblue', label='_nolegend_', linewidth=0.2, ax=ax)
sp_ESI.plot( x='wave', y='flam_u_cor',  color='lightgreen', label='_nolegend_', linewidth=0.2, ax=ax)
sp_G102.plot(x='wave', y='flam_u',      color='red',label='_nolegend_', linewidth=0.2, ax=ax)
sp_G141.plot(x='wave', y='flam_u',      color='purple', label='_nolegend_', linewidth=0.2, ax=ax)
sp_GNIRS_cutout.plot(x='wave', y='flam_u', color='lightgrey', label="_nolegend_", ax=ax)
plt.title("Flux-scaled spectra")
adjust_plot()

fig = plt.figure(5, figsize=figsize)   # Prettier plot
how2comb = {'flam_cor': 'mean', 'flam_u_cor': jrr.util.convenience1, 'wave': 'mean', 'flamcor_autocont' : 'mean'} #how to bin
bin_MMT =  np.arange(3200, 4800,  5)
bin_ESI =  np.arange(4350, 8900, 10)
cutout_MMT = jrr.spec.bin_boxcar_better(sp_MMT, bin_MMT, how2comb, bincol='wave')
cutout_ESI = jrr.spec.bin_boxcar_better(sp_ESI, bin_ESI, how2comb, bincol='wave')
cutout_G102 = sp_G102.loc[sp_G102['wave'].between(grism_info_G102['x1'], grism_info_G102['x2'])]
cutout_G141 = sp_G141.loc[sp_G141['wave'].between(grism_info_G141['x1'], grism_info_G141['x2'])]
ax = cutout_ESI.plot(    x='wave', y='flam_cor', color='green',  label='Keck ESI', figsize=figsize)
#sp_GNIRS_cutout.plot(x='wave', y='flam',     color='orange', label="GNIRS",     ax=ax)
cutout_MMT.plot( x='wave', y='flam_cor',    color='blue',   label='MMT BC',    ax=ax)
cutout_G102.plot(x='wave', y='flam',        color='red',    label='WFC3 G102', ax=ax)
cutout_G141.plot(x='wave', y='flam',        color='purple', label='WFC3 G141', ax=ax)
cutout_MMT.plot( x='wave', y='flam_u_cor',  color='lightblue', label='_nolegend_', linewidth=0.2, ax=ax)
cutout_ESI.plot( x='wave', y='flam_u_cor',  color='lightgreen', label='_nolegend_', linewidth=1, ax=ax)
cutout_G102.plot(x='wave', y='flam_u',      color='red',label='_nolegend_', linewidth=0.2, ax=ax)
cutout_G141.plot(x='wave', y='flam_u',      color='purple', label='_nolegend_', linewidth=0.2, ax=ax)
plt.plot(sp_ESI['wave'],  sp_ESI['flamcor_autocont'],    color='black', label='_nolegend_')
plt.plot(sp_MMT['wave'],  sp_MMT['flamcor_autocont'],    color='black',label='_nolegend_')
plt.plot(cutout_G102['wave'], cutout_G102['cont'],   color='black', label='_nolegend_')
plt.plot(cutout_G141['wave'], cutout_G141['cont'],   color='black', label='_nolegend_')
sp_GNIRS_cutout.plot(x='wave', y='flam_u', color='lightgrey', label="_nolegend_", ax=ax)
plt.title("Prettier plot.  MMT and ESI have been boxcar smoothed")
adjust_plot()


# Plot just the scaled continuum fits
fig = plt.figure(5, figsize=figsize)
plt.plot(sp_ESI['wave'],  sp_ESI['flamcor_autocont'],    color='green', label='Keck ESI')
plt.plot(sp_MMT['wave'],  sp_MMT['flamcor_autocont'],    color='blue', Label='MMT Blue Channel')
plt.plot(cutout_G102['wave'], cutout_G102['cont'],   color='red', label='WFC3 G102')
plt.plot(cutout_G141['wave'], cutout_G141['cont'],   color='purple', label='WFC3 G141')
plt.title("Plotting the scaled continuua")
adjust_plot()

# Sanity check how I fluxed the spectra
cont_sanity_check1 = sp_ESI.loc[sp_ESI['wave'].between(8000., 8500.)]['flam_cor'].median() / sp_G102.loc[sp_G102['wave'].between(8000., 8500.)]['flam'].median()
cont_sanity_check2 = sp_MMT.loc[sp_MMT['wave'].between(4200., 4300.)]['flam_cor'].median()  / sp_ESI.loc[sp_ESI['wave'].between(4200., 4300.)]['flam_cor'].median()
print "Sanity checks:  ESI/G102:", cont_sanity_check1, "and MMT/ESI:", cont_sanity_check2


plt.show()  # Show all plots at once, each in a separate window



# Ayan first runs a translator to get it into format EW_fitter.py likes
#subprocess.call(home + "/Python/AYAN_Code/convert_spectra_format.py --inpath ~/Dropbox/Grism_S1723/JRR_Working/ --infile s1723_MMT_wcont.txt --flamconst 1. --flamcol flam_cor --flamucol flam_u_cor --wavecol wave --flamcontcol flamcor_autocont --z 1.32952 --zu 4e-4")
#
#subprocess.call(home + "/Python/AYAN_Code/EW_fitter.py --short s1723_MMT_wcont_new-format --spec_list_file ~/Dropbox/Grism_S1723/JRR_Working/other-spectra-filenames-redshifts.txt --silent --path ~/Dropbox/Grism_S1723/JRR_Working/ --linelistpath ~/Dropbox/Grism_S1723/JRR_Working/ --fout s1723_MMT_emission_measured.out --useflamcont flam_cont --savepdf --hide")

##execfile( home + "/Python/AYAN_Code/convert_spectra_format.py    --inpath ./  --infile s1723_MMT_wcont.txt --flamconst# 1.0  --flamcol flam_cor --flamucol flam_u_cor  --wavecol wave   --flamcontcol flamcor_autocont --z 1.329279 --zu 0.000085")

##run ~/Python/AYAN_Code/EW_fitter.py --short s1723_MMT_wcont_new-format  --spec_list_file ./other-spectra-filenames-redshifts.txt --silent --hide --savepdf --fout s1723_MMT_measuredlines.out --useflamcont  flam_cont --linelistpath ./

# Have scaled fluxes by emission lines.  Let's sanity check that continuum isn't too far out of agreement



