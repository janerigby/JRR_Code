'''Making Jane's sensitivity plots comparing JWST to previous observatories.  Now,
with post-commissioning (for Cycle 2) sensitivities from pandeia 2.0, courtesy Klaus.  
I think these are the sensitivities they used for the Pocket Guide.  I need to clean
up the Keck and Spitzer comparisons.
jrigby, Dec 2022'''

from __future__ import print_function
from builtins import range
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter
from matplotlib.backends.backend_pdf import PdfPages
from os.path import expanduser, exists, basename
import glob
import re
import pandas
import jrr
from   astropy.io import fits
from astropy.table import Table

def load_photometry() :
    # Adding the comparison to other observatories. Switching from gnuplot to python
    phot_codes = {'NIRSpec' : 0, 'NIRCam': 1, 'MIRI' : 2, 'HST' : 4, 'WISE' : 5, 'Spitzer' : 6, 'Gemini' :10 , 'Herschel' : 11, 'SOFIA' : 12}
    #13=alma(cycle0), 14=alma(finished), 15=VLA, 16=EVLA
    names = ('name', 'wave', 'limfnu', 'code')                
    phot_df = pandas.read_csv(pandir + '../jwst-phot.dat', comment="#", delim_whitespace=True, usecols=[0,1,2,3], names=names)
    phot_df['wave'] =  pandas.to_numeric(phot_df['wave'])
    return(phot_df, phot_codes)

def load_spectroscopy() :  # input units are W/m^2
    # Adding the comparison to other observatories. Switching from gnuplot to python
    spec_codes = { 'nirspec' : 0, 'NIRCAM' : 1, 'MIRI' : 2, 'HST' : 4, 'WISE' : 5, 'Spitzer' : 6, 'Keck' : 9, 'Gemini' : 10, 'VLT' : 11, 'SOFIA' : 12, 'TMT' : 13, 'ELT' : 14}
    names = ('name', 'wave', 'limflux_wm2', 'R', 'code', 'comment')
    spec_df = pandas.read_csv(pandir + '../jwst-spec.dat', comment='#', delim_whitespace=True, names=names)
    spec_df['limflux_cgs'] = spec_df['limflux_wm2'] * 1E3
    return(spec_df, spec_codes)

def load_lowres_spec() : # input are Jy per spectral resoln element
    names = ('telescope', 'instrument', 'wave', 'limflux_Jy', 'foo', 'R', 'code', 'comment')
    df = pandas.read_csv(pandir + '../jwst-lowres-spec.dat', comment='#', delim_whitespace=True, names=names)
    df['limflux_mJy'] = df.limflux_Jy *1000
    jrr.spec.pssc_to_pssl_df(df, PSSLcol='PSSL', PSSCcol='limflux_mJy', wavecol='wave', Rcol='R') # should do same as below
    df['PSSL2']   = (df.limflux_mJy) * 3.0E-15 / (df.R  * df.wave) * 1E3  #The last 1E3 converts to cgs    
    return(df)
    
def add_annotations():
    plt.annotate("NIRISS", xy=(1,2E-17), color='red', xycoords='data', fontsize=fs2)
    plt.annotate("MIRI", xy=(14,5E-17), color='k', xycoords='data', fontsize=fs2)
    plt.annotate("NIRSpec MOS", xy=(0.6,0.6E-18), color='orange', xycoords='data', fontsize=fs2)
    plt.annotate("& IFU", xy=(1.6, 0.6E-18), color='green', xycoords='data', fontsize=fs2)
    plt.annotate("NIRCam grism", xy=(5,0.3E-17), color='blue', xycoords='data', fontsize=fs2)

def add_annotations_lowres():
    plt.annotate("NIRISS", xy=(1,2E-7), color='red', xycoords='data', fontsize=fs2)
    plt.annotate("MIRI", xy=(12,3E-5), color='k', xycoords='data', fontsize=fs2)
    plt.annotate("NIRSpec MOS", xy=(0.6,3E-8), color='orange', xycoords='data', fontsize=fs2)
    plt.annotate("& IFU", xy=(1.5, 3E-8), color='green', xycoords='data', fontsize=fs2)
    plt.annotate("NIRCam grism", xy=(5,1E-5), color='blue', xycoords='data', fontsize=fs2)
    
def pretty_plot() :
    plt.xlabel("wavelength (micron)", fontsize=fs3)
    plt.ylabel('limiting flux density (Jy)', fontsize=fs3)
    #plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.annotate('Jane.Rigby@nasa.gov, pandeia v2.0', xy=(0.64,0.16), xycoords='figure fraction')
    jrr.plot.force_axisticks_linear((ax.xaxis,))
    plt.xticks(fontsize=fs3)
    plt.yticks(fontsize=fs3)
    plt.tight_layout()


def load_nirspec_postlaunch() : # Load nirspec's post-commissioning curves, from their PASP paper
    # files from Box, via  Nimisha Kumari of STScI
    ndir = '/Users/jrrigby1/MISSIONS/JWST/Sens/Postlaunch/Nirspec_post-commissioning_Sensitivity_Curves_IFS_MOS/'
    IFS_files = [basename(x) for x in glob.glob(ndir + 'IFS-*Jy.src')]
    MOS_files = [basename(x) for x in glob.glob(ndir + 'MOS-*Jy.src')]
    #print(IFS_files, MOS_files)
    df_ifs = {} ; df_mos = {}
    for thisfile in IFS_files :  df_ifs[thisfile] = pandas.read_json(ndir + thisfile)
    for thisfile in MOS_files :  df_mos[thisfile] = pandas.read_json(ndir + thisfile)
    return(df_ifs, df_mos)
    
def load_nircam_grism_postlaunch() : # Load nircam grism post-commissioning curves, from their PASP paper.
    ndir = '/Users/jrrigby1/MISSIONS/JWST/Sens/NIRCam_WFSS_post-commissioning_sens/'
    grism_files = [basename(x) for x in glob.glob(ndir + 'F*.dat')]
    #columns = ('wavelength', 'cont_flux_uJy', 'cont_flux_uJy_old', 'line_flux_Em18', 'line_flux_Em18_old')
    df_nircam = {}
    for thisfile in grism_files:  df_nircam[thisfile] = pandas.read_csv(ndir + thisfile, delim_whitespace=True, comment='#', skiprows=1)
    return(df_nircam)

def grab_R_for_nirspec(grating) : #grating as in 'prism', 'g140m', 'g235h'
    ndir = '/Users/jrrigby1/MISSIONS/JWST/Sens/Postlaunch/NIRSpec_R/'
    Rfilename = 'jwst_nirspec_' + grating.lower() + '_disp.fits'
    RR, Rheader = fits.getdata(ndir + Rfilename, header=True)
    return(RR, Rheader)

def grab_niriss_wfss():
    ndir = '/Users/jrrigby1/MISSIONS/JWST/Sens/Postlaunch/'
    filenames = ['niriss-wfss-line-sens.txt', 'niriss-wfss-cont-sens.txt']
    df_line = pandas.read_csv(ndir + filenames[0], delim_whitespace=True, comment='#')
    df_line['sens_cgs'] = df_line['Sensitivity'] * 1E-20 * 1E3   # sensitivity is in units of (10^-20 W/m^2)
    df_cont = pandas.read_csv(ndir + filenames[1], delim_whitespace=True, comment='#')
    # need to test df_cont
    return(df_line, df_cont)

def spec_compare_to_other_obs(df_lowres_spec):  # compare to other observatories
    otherspec_df, spec_codes = load_spectroscopy()
    shortcodes = ('SOFIA', 'Keck', 'VLT', 'Gemini',  'Spitzer')
    for instr in shortcodes : # spec_codes :
        subset =  otherspec_df.loc[otherspec_df['code'] == spec_codes[instr]]
        plt.plot(subset['wave'], subset['limflux_cgs'], label=instr, marker='o', linestyle='dashed')
    subset = df_lowres_spec.loc[df_lowres_spec.telescope == 'Spitzer']      # Add low-resoln spitzer
    plt.plot(subset.wave, subset.PSSL, color='purple', label='Spitzer lowres', marker='o', linestyle='dashed')
    plt.legend(loc='upper center')
    return(0)


###############################################################
#pandir = '/Users/jrrigby1/MISSIONS/JWST/Sens/pandeia_sensitivities_1.3/'
#pandir = '/Users/jrrigby1/MISSIONS/JWST/Sens/pandeia_sensitivities_1.5/'
pandir = '/Users/jrrigby1/MISSIONS/JWST/Sens/Pandeia_v2.0/'  # POST-COMMISSIONING
fs1 = 14 ; fs2 = 16 ; fs3 = 22 # fontsizes
###############################################################

# Which curves from the IDTs to plot?
phot_compare_IDT = {'nircam': False, 'niriss': False}
spec_compare_IDT = {'nircam': False, 'niriss': False, 'nirspec': False}


df_sens = jrr.instruments.load_Klaus_sensitivities(pandir)        # Klaus's sensitivities from Pandeia.  I *think* these are in mJy
df_disp = jrr.instruments.load_Klaus_dispersions(pandir=None)     # Klaus's spectral resolutions from Pandeia
print("I think I imported all of Klaus's files:", list(df_sens.keys()))
plt.ion()
pp = PdfPages("jwst_sensitivity_jrigby_v2.pdf")  # output
plt.close("all")


########################################################################################################################
#   PHOTOMETRY PLOT 
fig, ax = plt.subplots(figsize=(8.5,6))
modes  = ['nircam_sw', 'nircam_lw', 'miri_imaging', 'niriss_imaging']
labels = ['JWST NIRCam', '_nolegend_', 'JWST MIRI', 'JWST NIRISS']
colors = ['blue', 'blue', 'k', 'red']
# Plot imaging
for ii, mode in enumerate(modes) :
    thismode = mode + '_sensitivity'
    subset = df_sens[thismode].loc[df_sens[thismode]['configs'].astype(str).str.contains('w') &  \
                                       ~df_sens[thismode]['configs'].astype(str).str.contains('w2')].sort_values(by='wavelengths')  # Skip W2
    # above is only W filters, but not W2
    print(subset.head())
    plt.plot(subset['wavelengths'], subset['lim_fluxes']*1E-3, label=labels[ii], color=colors[ii], marker='o', linestyle='-', markersize=7, alpha=1.0)

# Harvesting Jane's sensitivities for previous missions
(phot_df, phot_codes) = load_photometry()
df_lowres_spec = load_lowres_spec() # load the prelaunch lowres spec file
#instrs = ['NIRCam', 'MIRI', 'HST']
instrs = ['HST', 'Spitzer', 'Gemini']
for instr in instrs:
    subset =  phot_df.loc[phot_df['code'] == phot_codes[instr]]
    print(subset.head())
    plt.plot(subset['wave'], subset['limfnu'], label=instr, marker='o', linestyle='-')
plt.xlim(0.5, 28)
plt.ylim(2E-9, 3E-4)
plt.title("Photometric sensitivity, pt src, SNR=10 in " + r'$10^4$' + "s", fontsize=fs2) 
plt.text(1.1, 2.5E-6, "Gemini", color='green', fontsize=fs3)
plt.text(0.8, 8E-8, "Hubble", color='#1f77b4', fontsize=fs3)
plt.text(4.5, 1.5E-8, "NIRISS on-orbit", color='red', fontsize=fs3)
plt.text(4.5, 7E-9, "NIRCam on-orbit", color='blue', fontsize=fs3)
plt.text(8., 1E-7, "MIRI on-orbit", color='black', fontsize=fs3)
plt.text(8.8, 4.5E-5, "Spitzer", color='orange', fontsize=fs3)
pretty_plot()

# Compare to post-commissioning sensitivity estimates from NIRCam and NIRISS IDTs
df_postlaunch_phot = pandas.read_csv("../Postlaunch/jwst_postlaunch.dat", comment='#', delim_whitespace=True)
df_postlaunch_nircam = df_postlaunch_phot.loc[df_postlaunch_phot.instrument.str.contains('nircam')].sort_values(by='wavelength')
df_postlaunch_niriss = df_postlaunch_phot.loc[df_postlaunch_phot.instrument.str.contains('niriss')].sort_values(by='wavelength')
subset_M = df_postlaunch_nircam.loc[df_postlaunch_nircam.filtname.str.contains('M')] # just nircam M
subset_W = df_postlaunch_nircam.loc[df_postlaunch_phot.filtname.str.contains('W')]     # just nircam W
if phot_compare_IDT['nircam'] :
    plt.plot(subset_W.wavelength, subset_W.limiting_sensitivity_Jy, marker='*', color='blue', alpha=0.5, linestyle='dashed')
    #plt.plot(subset_M.wavelength, subset_M.limiting_sensitivity_Jy, marker='*', color='aqua')  # medium band NIRCam
if phot_compare_IDT['niriss'] :
    plt.plot(df_postlaunch_niriss.wavelength, df_postlaunch_niriss.limiting_sensitivity_Jy, marker='*', color='red', linestyle='dashed', alpha=0.5)
pp.savefig()
plt.show()


########################################################################################################################
# SPECTROSCOPY PLOT IN LINE FLUX
#
# Plot limiting sensitivity to an unresolved line (in erg/s/cm^2).  What a spectroscopist wants.
fig2, ax2 = plt.subplots(figsize=(8.5,6))
cuton  = { 'f277w': 2.46, 'f322w2': 2.46, 'f356w': 3.15, 'f444w': 3.94} # clean up plotting of NIRCam grism mode
cutoff = { 'f277w': 3.1,  'f322w2': 3.95, 'f356w': 3.95, 'f444w': 4.93} # clean up plotting of NIRCam grism mode

df_specsens = {}   # make a smaller dict of sensitivities for the spectroscopic modes only
#specmodes = ['nirspec_fs','nirspec_ifu', 'nirspec_msa',  'miri_mrs']# 'nircam_wfgrism', 'miri_lrs',  'niriss_wfss','niriss_soss'
specmodes = ['nirspec_ifu', 'nirspec_msa',  'miri_mrs', 'nircam_wfgrism',  'niriss_wfss']
just_plot_grisms = ('f444w', 'f322w2')

# Loop through the spectroscopic modes, computing PSSL from PSSC
for specmode in specmodes:
    thismode = specmode + '_sensitivity'    
    df_specsens[thismode] = df_sens[thismode].copy(deep=True)
    df_specsens[thismode]['mode'] = mode
    df_specsens[thismode]['R'] = df_specsens[thismode]['wavelengths'] * 0  # initialize
    df_specsens[thismode]['PSSL'] = df_specsens[thismode]['wavelengths'] * 0  # initialize
    for index, row in df_specsens[thismode].iterrows() :     # Go row by row, for each grating, order, etc
        if 'miri' in thismode :     disp_key = 'jwst_miri_' + row.aperture + '-' + row.disperser
        elif 'nircam' in thismode:  disp_key = 'jwst_nircam'
        elif 'nirspec' in thismode: disp_key = 'jwst_' + re.split('_', thismode)[0] + '_' + row.disperser 
        elif 'niriss' in thismode:
                disp_key = 'jwst_niriss_' + row.configs['disperser'] + '-ord1'  # should have an order attached.  Try first
                df_disp[disp_key].rename(columns={'Wavelength': 'WAVELENGTH'}, inplace=True)  # regularize this 
        rebinnedR = jrr.spec.rebin_spec_new(df_disp[disp_key]['WAVELENGTH'], df_disp[disp_key]['R'], row.wavelengths)
        df_specsens[thismode].loc[index, 'R'] = rebinnedR
        # Now, need to convert to limiting line flux for a spectrally unresolved line.
        # Eqn 2.7 of the Spitzer IRS manual:  PSSL = 3.0E-15 PSSC / (R lambda), where PSSC is in units of mJy, and lambda is in microns.  PSSL is in W/m^2
        df_specsens[thismode].loc[index, 'PSSL'] = row.lim_fluxes * 3.0E-15 / ( row.R  * row.wavelengths) * 1E3  #The last 1E3 converts to cgs

# I'm computing PSSL in cgs units.
# For MIRI, Pandeia has a column line_limits: Limiting integrated line flux in W/m^2.   So I don't need to recalculate it.
# But, it disagrees w PSSL.  Have written to Klaus, Maca
        
# Now plot it.  I've separated computing PSSL (above) from plotting it (here)
for specmode in specmodes:
    thismode = specmode + '_sensitivity'    
    for index, row in df_specsens[thismode].iterrows() :     # Go row by row, for each grating, order, etc
        if 'miri' in thismode :
            plt.plot(row.wavelengths, row.PSSL, color='k', alpha=1)
            #plt.plot(row.wavelengths, row.line_limits*1E3, color='pink', alpha=1)
        elif specmode == 'nircam_wfgrism' :
            threshold = 2
            indices = np.where(row.PSSL < threshold * row.PSSL.min())
            if row.configs['filter'] in just_plot_grisms:  # only plot predicted grism sens for wide filters!
                plt.plot(row.wavelengths[indices], row.PSSL[indices], color='blue', alpha=1)
                #plt.annotate(row.configs['filter'], (np.median(row.wavelengths), row.PSSL.min()), fontsize=10, color='blue')
        elif 'nirspec' in thismode:
            if 'ifu' in thismode:       plt.plot(row.wavelengths, row.PSSL, color='green', alpha=1)  
            elif 'msa' in thismode:     plt.plot(row.wavelengths, row.PSSL, color='orange', alpha=1)
        elif 'niriss' in thismode:
            threshold = 2.5
            indices = np.where(row.PSSL < threshold * row.PSSL.min())
            plt.plot(row.wavelengths[indices], row.PSSL[indices] * (3./2.), color='red', alpha=1)
            # the factor of 3/2 is because Steph says they assumed line covered 3 pix not 2
        else:
            plt.plot(row.wavelengths, row.PSSL)

plt.ylabel(r"line flux det. at SNR=10 in $10^4$s (erg s$^-1$ cm$^2$)", fontsize=fs2)
plt.title("Spectroscopic line sensitivity, pt src, SNR=10 in " + r'$10^4$' + "s", fontsize=fs2)

bigspec_df = pandas.concat(df_specsens)

if spec_compare_IDT['nirspec']:
    (df_ifs, df_mos) = load_nirspec_postlaunch()  # post-commissioning nirspec sensitivities
    threshold = 10 # Don't plot the edges of the filter curve, just where it has sensitivity
    for df1 in (df_ifs, df_mos) :
        for key in df1.keys():
            df = df1[key]
            # Parse key to get grating
            grating = re.split('-|_', key)[2]
            RR, Rheader = grab_R_for_nirspec(grating)
            df['R'] = jrr.spec.rebin_spec_new(RR.WAVELENGTH, RR.R, df.wavelength * 1E6)  # Interpolate R onto the sensitivity curve
            df['limflux_mJy'] = df.data *1E3
            df['wave_um'] = df.wavelength *1E6
            df['PSSL2'] = 3.0E-15 * (df.data * 1E3) / (df.R * df.wavelength *1E6) *1E3  # in cgs
            jrr.spec.pssc_to_pssl_df(df, PSSLcol='PSSL', PSSCcol='limflux_mJy', wavecol='wave_um', Rcol='R') # should do same as above.  Just checking
            color = 'green' if ('IFS' in key) else 'orange'
            subset = df.loc[df.PSSL < threshold * df.PSSL.min()]
            plt.plot(subset.wavelength *1E6, subset.PSSL, color=color, linestyle='dashed', alpha=0.5) 

if spec_compare_IDT['nircam']:
    df_nircamgrism = load_nircam_grism_postlaunch() #post-commissioning nircam sensitivities
    threshold = 14 # Don't plot the edges of the filter cuve, just where it has sensitivity
    for key in df_nircamgrism.keys():
        df = df_nircamgrism[key]
        subset = df.loc[ (df['line_flux_E-18']  < threshold * df['line_flux_E-18'].min())]
        plt.plot(subset.wavelength, subset['line_flux_E-18'] * 1E-18, color='lightblue', linestyle='dashed')  # data is in cgs, unit of 10^-18 erg/s/cm^2.
        #plt.annotate(key, (np.median(subset.wavelength), subset['line_flux_E-18'].min() *1E-18), fontsize=10, color='blue', alpha=0.5)

# post-commissioning niriss wfss sensitivity
if spec_compare_IDT['niriss']:
    df_niriss_wfss_line, df_niriss_wfss_cont  = grab_niriss_wfss()
    for thisfilter in df_niriss_wfss_line.Filter.unique():
        subset = df_niriss_wfss_line.loc[df_niriss_wfss_line.Filter == thisfilter]
        plt.plot(subset.Wavelength, subset.sens_cgs, color='red', linestyle='dashed', alpha=0.5)

spec_compare_to_other_obs(df_lowres_spec)

pretty_plot()
plt.ylabel(r'limiting line flux (erg s$^{-1}$ cm$^{-2}$ )', fontsize=fs3)

plt.xlim(0.5,30)
plt.ylim(5E-19,2E-15)
ax2.xaxis.set_major_formatter(ScalarFormatter())
add_annotations()
plt.annotate('will update to Keck/MOSFIRE', (2, 5E-17))
plt.tight_layout()
pp.savefig()
plt.show()



###################################################################################################################
### Spectroscopy plot with continuum sensitivity
# What are the units here? From above, looks like mJy.  Double check!
#
fig, ax = plt.subplots(figsize=(8.5,6))
for specmode in specmodes:
    thismode = specmode + '_sensitivity'    
    for index, row in df_specsens[thismode].iterrows() :     # Go row by row, for each grating, order, etc
        if 'miri' in thismode :
            plt.plot(row.wavelengths, row.lim_fluxes *1E-3, color='k')
        elif specmode == 'nircam_wfgrism' :
            threshold = 1.4 
            indices = np.where(row.lim_fluxes < threshold * row.lim_fluxes.min())
            if row.configs['filter'] in just_plot_grisms:  # only plot predicted grism sens for wide filters!
                plt.plot(row.wavelengths[indices], row.lim_fluxes[indices] *1E-3, color='blue', alpha=0.5)
        elif 'nirspec' in thismode:
            if 'ifu' in thismode:       plt.plot(row.wavelengths, row.lim_fluxes *1E-3, color='green')  
            elif 'msa' in thismode:     plt.plot(row.wavelengths, row.lim_fluxes *1E-3, color='orange')
        elif 'niriss' in thismode:
            threshold = 2
            indices = np.where(row.lim_fluxes < threshold * row.lim_fluxes.min())
            plt.plot(row.wavelengths[indices], row.lim_fluxes[indices] *1E-3, color='red')
        else:
            plt.plot(row.wavelengths, row.lim_fluxes)

# Above are prelaunch predictions.  Now plot the actuals
if spec_compare_IDT['nircam'] :      
    for key in df_nircamgrism.keys():
        df = df_nircamgrism[key]
        subset = df.loc[ (df['cont_flux_uJy']  < threshold * df['cont_flux_uJy'].min())]
        plt.plot(subset.wavelength, subset['cont_flux_uJy'] * 1E-6, color='blue', linestyle='dashed', alpha=0.5)  # convert from microJy to Jy

# post-commissioning nirspec sensitivities
if spec_compare_IDT['nirspec'] :
    threshold = 10 # Don't plot the edges of the filter curve, just where it has sensitivity
    for df1 in (df_ifs, df_mos) :
        for key in df1.keys():
            df = df1[key]
            # Parse key to get grating
            grating = re.split('-|_', key)[2]
            color = 'green' if ('IFS' in key) else 'orange'
            subset = df.loc[df.data < threshold * df.data.min()]
            plt.plot(subset.wavelength *1E6, subset.data, color=color, linestyle='dashed', alpha=0.5) 

# post-commissioning NIRISS sensitivity 
if spec_compare_IDT['niriss']:
    for thisfilter in df_niriss_wfss_cont.Filter.unique():
        subset = df_niriss_wfss_cont.loc[df_niriss_wfss_cont.Filter == thisfilter]
        plt.plot(subset.Wavelength, subset.Sensitivity * 1E-6, color='k', linestyle='dashed', alpha=0.5)

    
# plot comparisons:
subset = df_lowres_spec.loc[df_lowres_spec.telescope == 'HST']
plt.scatter(subset.wave, subset.limflux_Jy, color='red', label='HST')             
#
subset = df_lowres_spec.loc[df_lowres_spec.telescope == 'Spitzer']
plt.scatter(subset.wave, subset.limflux_Jy, color='purple', label='Spitzer')             
#
plt.title("Spectroscopic sensitivity, pt src, per resoln element, SNR=10 in " + r'$10^4$' + "s", fontsize=fs2)
plt.legend()

add_annotations_lowres()
pretty_plot()
pp.savefig()
plt.show()



######### Temp plot
fig, ax = plt.subplots(figsize=(8.5,6))
for specmode in specmodes:
    thismode = specmode + '_sensitivity'    
    for index, row in df_specsens[thismode].iterrows() :     # Go row by row, for each grating, order, etc
        if 'miri' in thismode :
            plt.plot(row.wavelengths, row.PSSL / (row.line_limits*1E3), color='k')
plt.xlabel("wavelength (micron)")
plt.ylabel("PSSL /  (line_limits * 1E3)")
pp.savefig()
plt.show()
            
pp.close()



pandeia_ver = "v2.0"
outfile = 'pandeia_lineflux_sensitivities_' + pandeia_ver + '.readme'  
header_text =  '# PSSL (erg/s/cm^2) is the line flux detectable for a spectrally unresolved line, point source, at SNR=10, in 1E4s, from Pandeia ' + pandeia_ver + '\n'
header_text += '# Adapting Eqn 2.7 of the Spitzer IRS manual:  PSSL (erg/s/cm^2) = 3.0E-15 PSSC / (R lambda) * 1E3,\n'
header_text += '# where PSSC is in units of mJy, lambda is in microns, and PSSL is in erg/s/cm^2\n'
header_text += '# In Python pandas, read this file as:  df = pandas.read_csv(\'' + outfile + '\', index_col=[0,1], comment=\'#\')\n'
header_text += '# Jane.Rigby@nasa.gov, 20 Aug. 2020.\n#\n'
jrr.util.put_header_on_file('/dev/null', header_text, outfile)
bigspec_df.to_pickle(re.sub('readme', 'pcl', outfile))
# read this as  df_sens = pandas.read_pickle('pandeia_lineflux_sensitivities_v1.5.0.pcl')
 

# TODO:  
# DONE Add actual NIRCam and NIRISS sensitivities from Marcia's paper and the SciPerf paper
# DONE Import actual NIRSpec sensitivities from Torsten
# Improve comparison to Keck. Have asked John O'M
# Ask MIRI for the same
# Ask NIRISS for same
# Separate plots for pre-launch, post-launch, so easier to powerpoint.  Dotted lines for predicts?


