from os.path import expanduser, basename, exists
from os import remove
import warnings
import subprocess
import glob
from re import sub
import jrr  # Moved some functions to jrr.grism
import pandas
import numpy as np
#import query_argonaut
import astropy.convolution
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.ticker import AutoMinorLocator
matplotlib.rcParams.update({'font.size': 14})

''' This script processes the spatially-integrated spectra of SGAS J1723, so that we can get 
emission line fluxes, and make figures for the paper.  jrigby, 2017.  
Updated June 2018 to use Michael's new grizli redux of grism spectra.'''

# RUN THIS FROM THE FOLLOWING DIR:  /Dropbox/Grism_S1723/JRR_Working2/

def annotate_header(header_ESI, header_MMT, header_GNIRS) :
    for header in (header_ESI, header_MMT, header_GNIRS) :
        header += ("#\n# Below is extra processing by S1723_working.py\n")
        header += ("# Corrected for MW reddening of E(B-V)=" + str(EBV) + "\n")
    return(header_ESI, header_MMT, header_GNIRS)

def MMT_barycentric_correction(sp, header_MMT) :  # Applying barycentric correction to wavelengths
    thistime =  Time('2014-05-05T09:21:00', format='isot', scale='utc')
    thisradec = SkyCoord("17:23:37.23", "+34:11:59.07", unit=(u.hourangle, u.deg), frame='icrs')
    mmt= EarthLocation.of_site('mmt') # Needs internet connection 
    barycor_vel = jrr.barycen.compute_barycentric_correction(thistime, thisradec, location=mmt)
    header_MMT += ("# The barycentric correction factor for s1723 was" + str(barycor_vel))
    jrr.barycen.apply_barycentric_correction(sp, barycor_vel, colwav='oldwave', colwavnew='wave') 
    return(0)

def flag_noise_peaks(sp, delta=0.15, neighborpix=2, colfu='flam_u', contmask='contmask') :                         
    maxtab, mintab = jrr.peakdet.peakdet(sp[colfu], delta)  # Find peaks.
    peak_ind =  [np.int(p[0]) for p in maxtab] # The maxima
    for thispeak in peak_ind :
        #sp.iloc[thispeak - neighborpix : thispeak + neighborpix][contmask] = True
        sp.iloc[thispeak - neighborpix : thispeak + neighborpix, sp.columns.get_loc(contmask)] = True # Longer, avoids SettingWithCopyWarning
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

def adjust_plot(x1 = 3200., x2 = 16900., y1=0., y2=1.25E-16) :
    plt.legend()
    plt.ylim(y1, y2)
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
    #fig.set_tight_layout(True)
    fig.subplots_adjust(hspace=0.1)
    plt.tight_layout(True)
    return(0)

##########  END OF FUNCTIONS ###############
plt.ion()

#######  Housekeeping  ##############
plt.close("all")
figsize=(16,6)
nolegend="_" # Matplotlib versions handle this differently.  Used to be "_nolegend_"
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")  # nolegend generates cluttering warnings.
plot_MMTESIcont=False

#### Directories and filenames
home = expanduser("~")
wdir = home + '/Dropbox/Grism_S1723/'  
file_ESI  = home + '/Dropbox/MagE_atlas/Contrib/ESI_spectra/2016Aug/s1723_wholearc_ESI_JRR_sum.txt'   #  just arc_a, arc_b
file_MMT  = wdir + 'MMT_BlueChannel/spec_s1723_14b_final_F.txt'  # Updated Jan 23 2017
file_GNIRS = wdir + 'GNIRS/s1723_gnirs_residual_subtracted_spectrum_jrr.csv'

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

#### Read WFC3 spectra (w continuua, both grisms, both rolls, & sum of rolls).  MW dereddening was applied at continuum fitting stage.
grism_info = {}  # Dict of dicts that store the grism info
grism_info['G102'] = jrr.grism.get_grism_info("G102") 
grism_info['G141'] = jrr.grism.get_grism_info("G141")
grism_info['G102']['x1'] = 8500. # Clipping this harsher
sp_grism     = {}   # Load grism spectra into a dict of dataframes
cutout_grism = {}   # The same, for wavelength-trimmed versions
grismfiles = [x for x in glob.glob(wdir + 'Grizli_redux/1Dsum/Wcont/*wcontMWdr.txt')] 
for grismfile in grismfiles :
    parts = jrr.grism.parse_filename(basename(grismfile))
    label = str(parts['roll']) + '_' + parts['grating']
    sp_grism[label] = pandas.read_csv(grismfile, comment="#")
    sp_grism[label].set_index('wave', inplace=True, drop=False)
    if 'both' in parts['roll'] :  jrr.grism.half_the_flux( sp_grism[label] )
    # Retrieve the best-fit generated by grism_fit_spectra_v4.py
    modelfile = sub(".txt", "_meth2.model", sub("Wcont/", "", sub("Grizli_redux", "WFC3_fit_1Dspec", grismfile)))
    sp_temp = pandas.read_csv(modelfile, comment="#")  
    sp_temp.set_index('wave', inplace=True, drop=False)
    sp_grism[label]['bestfit'] = sp_temp['bestfit']  # This should work?
    sp_grism[label]['flamcontsub'] = sp_grism[label]['flam'] - sp_grism[label]['cont']
    cutout_grism[label] = sp_grism[label].loc[sp_grism[label]['wave'].between(grism_info[parts['grating']]['x1'], grism_info[parts['grating']]['x2'])]
    
#### Read MMT Blue Channel spectrum. # wave in vacuum Ang, flam in 1E-17 erg/s/cm^2/A
names=('oldwave', 'flam', 'flam_u')  
sp_MMT = pandas.read_table(file_MMT, delim_whitespace=True, comment="#", names=names)
header_MMT =("# MMT Blue Channel spectrum of SGAS J1723.\n")  # Nothing worth grabbing from original header
sp_MMT['flam']   *= 1.0E-17   # flam was in units of 1E-17
sp_MMT['flam_u'] *= 1.0E-17
#sp_MMT.replace([np.inf, -np.inf], np.nan, inplace=True)
MMT_barycentric_correction(sp_MMT, header_MMT) 
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

# Scale the flux of the ESI spectrum, to match HST G102, using [O II] doublet flux.  Direct summation of flux, linear local continuua.
flux_OII_G102 = (70.991 + 97.251) * 1E-17  # from 1Dsum/sgas1723_1Dsum_bothroll_G102_wcontMWdr_meth2.fitdf.  Updated 17July2018
rawflux_OII_ESI = sum_spectraline_wflatcont(sp_ESI, 8684., 8712., 8500., 8900.)
sp_ESI['flam_cor']   = sp_ESI['flam']   *  flux_OII_G102/rawflux_OII_ESI 
sp_ESI['flam_u_cor'] = sp_ESI['flam_u'] *  flux_OII_G102/rawflux_OII_ESI
# ESI flux in file_ESI should be 2x too high, b/c I summed 2 exposures.
howscaled_ESI = "# Scaled ESI flux so that sum of [O II] 3727+3729 fluxes matches flux of those lines in G102.  Factor was "
howscaled_ESI += (str(np.round((flux_OII_G102/rawflux_OII_ESI), 4)) + "\n")

# Scale the MMT spectrum to the ESI flux, using the flux in the summed [C~III] doublet
rawflux_CIII_ESI = sum_spectraline_wflatcont(sp_ESI, 4443., 4456., 4400., 4490)
rawflux_CIII_MMT = sum_spectraline_wflatcont(sp_MMT, 4436., 4450., 4400., 4490)
sp_MMT['flam_cor']   = sp_MMT['flam']   *  flux_OII_G102/rawflux_OII_ESI * rawflux_CIII_ESI/rawflux_CIII_MMT
sp_MMT['flam_u_cor'] = sp_MMT['flam_u'] *  flux_OII_G102/rawflux_OII_ESI * rawflux_CIII_ESI/rawflux_CIII_MMT
howscaled_MMT = "# Scaled MMT flux so that sum of CIII in MMT matches sum of those lines in ESI, after ESI scaled to G102. Factor was "
howscaled_MMT +=  (str(np.round((flux_OII_G102/rawflux_OII_ESI * rawflux_CIII_ESI/rawflux_CIII_MMT), 4)) + "\n")
# The MMT fluxing gets nudged up by 34%.  That's sensible given a longslit.

#### Flux the GNIRS spectrum
flux_Ha_G141 = 408.79E-17  #  kludge, hardcoded, from sgas1723_1Dsum_bothroll_G141_wcontMWdr_meth2.fitdf.  Updated 17July2018
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
cutout_GNIRS = sp_GNIRS.loc[sp_GNIRS['wave'].between(1.5E4,1.55E4)]

# Report how scaled the fluxes
howscaled_GNIRS = "# Scaled GNIRS spectrum by ratio of Halpha flux, by factor " + str(GNIRS_scaleflux)
header_ESI   += howscaled_ESI
header_MMT   += howscaled_MMT
header_GNIRS += howscaled_GNIRS 
print "How scaled spectra:\n", howscaled_ESI, howscaled_MMT,  howscaled_GNIRS


# Sanity checks on the MMT, ESI, G102 flux scaling
cont_sanity_check1 = sp_ESI.loc[sp_ESI['wave'].between(9100.,1E4)]['flam_cor'].median() / sp_grism['bothroll_G102'].loc[sp_grism['bothroll_G102']['wave'].between(9100.,1E4)]['flam'].median()
cont_sanity_check2 = sp_MMT.loc[sp_MMT['wave'].between(4200., 4300.)]['flam_cor'].median()  / sp_ESI.loc[sp_ESI['wave'].between(4200., 4300.)]['flam_cor'].median()
print "Flux sanity checks:  ESI/G102:", cont_sanity_check1, "and MMT/ESI:", cont_sanity_check2

fig, ax = plt.subplots(figsize=figsize)
order_keys = ['bothroll_G102', 'roll139_G102', 'roll308_G102', 'bothroll_G141', 'roll139_G141', 'roll308_G141']  # plot in this order
colors     = ['k', 'dimgrey',                     'purple']*2 
labels     = ['Average of rolls',   r'Roll $139\degree$',   r'Roll $308\degree$', nolegend, nolegend, nolegend]
for ii, grismkey in enumerate(order_keys) :
    plt.plot(cutout_grism[grismkey]['wave'], cutout_grism[grismkey]['flamcontsub'], label=labels[ii], color=colors[ii], lw=1.5) #, linestyle='steps-mid')
    plt.plot(cutout_grism[grismkey]['wave'], cutout_grism[grismkey]['flamcontsub'], label=labels[ii], color=colors[ii], lw=1.5, marker='o')
    #plt.plot(cutout_grism[grismkey]['wave'], cutout_grism[grismkey]['bestfit'],     '--', label=labels[ii], color=colors[ii], lw=1.)
print "Checking continuum shape in grism vs roll.  Contamination?"
adjust_plot(x1=8000, y2=0.8E-16)
plt.legend()


### Fit MMT continuum.  The boxcar # is arbitrary, seems to match
if plot_MMTESIcont :
    fig = plt.figure(2, figsize=figsize)
    plt.title("Fit continuum for MMT bluechannel")
(smooth1, smooth2) = jrr.grism.wrap_fit_continuum(sp_MMT, LL, zHa, boxcar=151, colf='flam_cor',  colcont='flamcor_autocont', label="MMT Blue Channel", makeplot=plot_MMTESIcont)

### Fit ESI continuum
peak_ind = flag_noise_peaks(sp_ESI, colfu='flam_u_cor', delta=0.05E-17)
if plot_MMTESIcont :
    fig = plt.figure(3, figsize=figsize)
    plt.title("Fit continuum for Keck ESI")
(smooth1, smooth2) = jrr.grism.wrap_fit_continuum(sp_ESI, LL, zHa, boxcar=351, colf='flam_cor',  colcont='flamcor_autocont', label="ESI", makeplot=plot_MMTESIcont)
header_ESI += ("# Applied smooth continuum \n")
if plot_MMTESIcont :  plt.ylim(0,9E-17)

print "Writing flux-adjusted, continuum-fitted spectra for MMT and ESI."
sp_ESI.to_csv('temp1', sep='\t', na_rep='NaN')
sp_MMT.to_csv("temp2", sep='\t', na_rep='NaN')
jrr.util.put_header_on_file('temp1', header_ESI, "s1723_ESI_wcont.txt")
jrr.util.put_header_on_file('temp2', header_MMT, "s1723_MMT_wcont.txt")

# Plot the scaled spectra
fig, ax  = plt.subplots(figsize=figsize)
sp_ESI.plot(    x='wave', y='flam_cor', color='green',  label='Keck ESI', ax=ax)
cutout_GNIRS.plot(x='wave', y='flam',     color='orange', label="GNIRS",     ax=ax)
sp_MMT.plot( x='wave', y='flam_cor',    color='blue',     label='MMT BC',    ax=ax)
sp_grism['bothroll_G102'].plot(x='wave', y='flam',        color='red',    label='WFC3 G102', ax=ax)
sp_grism['bothroll_G141'].plot(x='wave', y='flam',        color='purple', label='WFC3 G141', ax=ax)
sp_MMT.plot( x='wave', y='flam_u_cor',  color='lightblue', label=nolegend, linewidth=0.2, ax=ax)
sp_ESI.plot( x='wave', y='flam_u_cor',  color='lightgreen', label=nolegend, linewidth=0.2, ax=ax)
sp_grism['bothroll_G102'].plot(x='wave', y='flam_u',      color='red',label=nolegend, linewidth=0.2, ax=ax)
sp_grism['bothroll_G141'].plot(x='wave', y='flam_u',      color='purple', label=nolegend, linewidth=0.2, ax=ax)
cutout_GNIRS.plot(x='wave', y='flam_u', color='lightgrey', label=nolegend, ax=ax)
plt.title("Flux-scaled spectra")
adjust_plot()

# Roll this up into one composite spectrum.  For every wavelength, pick the best spectrum.
n1 = sp_MMT[['wave', 'flam_cor', 'flam_u_cor']].loc[sp_MMT['wave'].between(3200., 4750.)]
n1.set_index('wave', inplace=True, drop=False)
n2 = sp_ESI[['wave', 'flam_cor', 'flam_u_cor']].loc[sp_ESI['wave'].between(4750., 8500.)]
n2.set_index('wave', inplace=True, drop=False)
n3 = sp_grism['bothroll_G102'].loc[sp_grism['bothroll_G102']['wave'].between(8500., 11000)][['wave', 'flam', 'flam_u']]
n3.rename(columns={'flam' : 'flam_cor',  'flam_u' : 'flam_u_cor'}, inplace=True)
n4 = sp_grism['bothroll_G141'].loc[sp_grism['bothroll_G141']['wave'].between(11000., 16250.)][['wave', 'flam', 'flam_u']]
n4.rename(columns={'flam' : 'flam_cor',  'flam_u' : 'flam_u_cor'}, inplace=True)
composite = pandas.concat((n1, n2, n3, n4))
composite.to_csv("temp3", na_rep='NaN', index=False)
comphead = "# Composite convenience spectrum for SGAS 1723.  Ignores overlap, picks one spectrum for each wavelength, as follows:\n"
comphead +="# MMT Blue Channel 3200--4750A; Keck ESI 4750--8500A; HST WFC3-IR G102 grism 8500--11000A; HST WFC3-IR G141 grism 11000--16250A\n"
comphead +="# units:  wave is wavelength in vacuum angstroms; f_lambda is erg/s/cm^2/s\n";
jrr.util.put_header_on_file('temp3', comphead, "S1723_composite_spectrum.csv")

fig, ax = plt.subplots(figsize=figsize)   # Prettier plot
how2comb = {'flam_cor': 'mean', 'flam_u_cor': jrr.util.convenience1, 'wave': 'mean', 'flamcor_autocont' : 'mean'} #how to bin
bin_MMT =  np.arange(3200, 4800,  5)
#bin_ESI =  np.arange(4350, 8900, 10)
bin_ESI =  np.arange(4350, 10000, 10)
binned_MMT = jrr.spec.bin_boxcar_better(sp_MMT, bin_MMT, how2comb, bincol='wave')
binned_ESI = jrr.spec.bin_boxcar_better(sp_ESI, bin_ESI, how2comb, bincol='wave')
binned_ESI_clean = binned_ESI.loc[binned_ESI['flam_u_cor'].lt(1E-14)] # Avoid the stuff w big errorbars, for plotting
binned_MMT.plot(  x='wave', y='flam_cor',         color='blue',   label='MMT BC', figsize=figsize, lw=2, ax=ax)
binned_ESI_clean.plot( x='wave', y='flam_cor',         color='green',  label='Keck ESI', ax=ax, lw=2)
cutout_grism['bothroll_G102'].plot(x='wave', y='flam',             color='red',    label='WFC3 G102',ax=ax, lw=2)
cutout_grism['bothroll_G141'].plot(x='wave', y='flam',             color='purple', label='WFC3 G141',ax=ax, lw=2)
binned_MMT.plot(       x='wave', y='flam_u_cor',       color='blue',   label=nolegend, ax=ax, lw=1)
binned_ESI_clean.plot( x='wave', y='flam_u_cor',       color='green',  label=nolegend, ax=ax, lw=1)
cutout_grism['bothroll_G102'].plot(      x='wave', y='flam_u',           color='red',    label=nolegend, ax=ax, lw=1)
cutout_grism['bothroll_G141'].plot(      x='wave', y='flam_u',           color='purple', label=nolegend, ax=ax, lw=1)
binned_MMT.plot(       x='wave', y='flamcor_autocont', color='black',  label=nolegend, ax=ax)
binned_ESI_clean.plot( x='wave', y='flamcor_autocont', color='black',  label=nolegend, ax=ax)
cutout_grism['bothroll_G102'].plot(      x='wave',  y='cont',            color='black',  label=nolegend, ax=ax)
cutout_grism['bothroll_G141'].plot(     x ='wave',  y='cont',            color='black',  label=nolegend, ax=ax)
#plt.annotate("Prettier plot.  MMT and ESI have been boxcar smoothed", (0.3,0.9), xycoords='axes fraction')
adjust_plot()





outfilename1 = 'hahb_rats_G102G141.txt'
outfilename2 = 'hahb_rats_G141only.txt'
plt.show()  # Show all plots at once, each in a separate window
print("\n\nMeasure some basic line ratios.")
subdirs = ('1D_complete_images_A2_A3/', '1Dbyclumps/', '1Dsum/')
G102fits = []
for subdir in subdirs:
    fitdir = home + '/Dropbox/Grism_S1723/WFC3_fit_1Dspec/' + subdir
    G102fits += [ x for x in glob.glob(fitdir + '*G102*2.fitdf')]
print "Ha/Hbeta from G141, G102", subdir
jrr.grism.measure_linerats_usebothgrisms(G102fits, outfilename1, line1='Halpha_G141', line2='Hbeta_G102', verbose=True)
print "Ha/Hbeta from G141 only", subdir
jrr.grism.measure_linerats_usebothgrisms(G102fits, outfilename2, line1='Halpha_G141', line2='Hbeta_G141', verbose=True)

    #print "4363/Hbeta ratios"
    #jrr.grism.measure_linerats_fromfiles(G102fits, fitdir, '[O~III]', 'Hbeta', verbose=True)
    #print "4363/HBeta again, this time using method 1 fits"
    #jrr.grism.measure_linerats_fromfiles(G102fits_method1, fitdir, '[O~III]', 'Hbeta', verbose=True)

# Need to grab both a G102 and G141 2.fitdf, merge them (w linenames + -G102, G141 to be unique), and then
# do the same measure_linerats() query across the 2 grisms.  Doublecheck Ayan's flux ratios, and then try to
# reconcile with the PDFs of the line fits.  What's going on with Ha/Hbeta?



# Ayan first runs a translator to get it into format EW_fitter.py likes
#subprocess.call(home + "/Python/AYAN_Code/convert_spectra_format.py --inpath ~/Dropbox/Grism_S1723/JRR_Working/ --infile s1723_MMT_wcont.txt --flamconst 1. --flamcol flam_cor --flamucol flam_u_cor --wavecol wave --flamcontcol flamcor_autocont --z 1.32952 --zu 4e-4")
#
#subprocess.call(home + "/Python/AYAN_Code/EW_fitter.py --short s1723_MMT_wcont_new-format --spec_list_file ~/Dropbox/Grism_S1723/JRR_Working/other-spectra-filenames-redshifts.txt --silent --path ~/Dropbox/Grism_S1723/JRR_Working/ --linelistpath ~/Dropbox/Grism_S1723/JRR_Working/ --fout s1723_MMT_emission_measured.out --useflamcont flam_cont --savepdf --hide")

##execfile( home + "/Python/AYAN_Code/convert_spectra_format.py    --inpath ./  --infile s1723_MMT_wcont.txt --flamconst# 1.0  --flamcol flam_cor --flamucol flam_u_cor  --wavecol wave   --flamcontcol flamcor_autocont --z 1.329279 --zu 0.000085")

##run ~/Python/AYAN_Code/EW_fitter.py --short s1723_MMT_wcont_new-format  --spec_list_file ./other-spectra-filenames-redshifts.txt --silent --hide --savepdf --fout s1723_MMT_measuredlines.out --useflamcont  flam_cont --linelistpath ./

# Have scaled fluxes by emission lines.  Let's sanity check that continuum isn't too far out of agreement



