from __future__ import print_function
from builtins import range
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter
from os.path import expanduser, exists
import re
import pandas
import jrr
from   astropy.io import fits
from astropy.table import Table          

#pandir = '/Users/jrrigby1/MISSIONS/JWST/Sens/pandeia_sensitivities_1.3/'
pandir = '/Users/jrrigby1/MISSIONS/JWST/Sens/pandeia_sensitivities_1.5/'
fs1 = 14 ; fs2 = 18 ; fs3 = 22 # fontsizes

def load_photometry() :
    # Adding the comparison to other observatories. Switching from gnuplot to python
    phot_codes = {'NIRSpec' : 0, 'NIRCam': 1, 'MIRI' : 2, 'HST' : 4, 'WISE' : 5, 'Spitzer' : 6, 'gemini' :10 , 'herschel' : 11, 'sofia' : 12}
    #13=alma(cycle0), 14=alma(finished), 15=VLA, 16=EVLA
    names = ('name', 'wave', 'limfnu', 'code')                
    phot_df = pandas.read_csv(pandir + '../jwst-phot.dat', comment="#", delim_whitespace=True, usecols=[0,1,2,3], names=names)
    phot_df['wave'] =  pandas.to_numeric(phot_df['wave'])
    return(phot_df, phot_codes)

def load_spectroscopy() :
    # Can I grab the continuum sensitivity folder, so it's easy?  Not what users want, though.  Have asked for spectral R
    # so that I can convert from continuum sensitivity to line sensitivity
    return(0)

def pretty_plot() :
    plt.xlabel("wavelength (micron)", fontsize=fs3)
    plt.ylabel('limiting flux density (Jy)', fontsize=fs3)
    #plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()

    
df_sens = jrr.instruments.load_Klaus_sensitivities(pandir)          # Klaus's sensitivities from Pandeia
df_disp = jrr.instruments.load_Klaus_dispersions(pandir=None)  # Klaus's spectral resolutions from Pandeia
print("I think I imported all of Klaus's files:", list(df_sens.keys()))
plt.ion()


plt.close("all")
fig, ax = plt.subplots(figsize=(8.5,6))
plt.xticks(fontsize=fs3)
plt.yticks(fontsize=fs3)
# MAKE A PHOTOMETRY PLOT
# First, plot Klaus's values.  May need to screen for wide filters
modes  = ['nircam_sw', 'nircam_lw', 'miri_imaging']
labels = ['JWST NIRCam', '_nolegend_', 'JWST MIRI']
# Plot imaging
for ii, mode in enumerate(modes) :
    thismode = mode + '_sensitivity'
    subset = df_sens[thismode].loc[df_sens[thismode]['configs'].astype(str).str.contains('w')]  # Only Wide (broadband) filters  
    plt.plot(subset['wavelengths'], subset['lim_fluxes']*1E-3, label=labels[ii], color='red', marker='o', linestyle='-')

# Harvesting Jane's sensitivities
(phot_df, phot_codes) = load_photometry()
#instrs = ['NIRCam', 'MIRI', 'HST']
instrs = ['HST', 'Spitzer', 'gemini']
for instr in instrs:
    subset =  phot_df.loc[phot_df['code'] == phot_codes[instr]]
    print(subset.head())
    plt.plot(subset['wave'], subset['limfnu'], label=instr, marker='o', linestyle='-')

plt.xlim(0.5, 28)
plt.ylim(2E-9, 3E-4)
plt.title("photometric performance, pt source, SNR=10 in " + r'$10^4$' + "s", fontsize=fs2) 
plt.text(1.1, 2.5E-6, "Gemini", color='green', fontsize=fs3)
plt.text(0.8, 8E-8, "Hubble", color='#1f77b4', fontsize=fs3)
plt.text(2.0, 3E-9, "JWST NIRCam", color='red', fontsize=fs3)
plt.text(11., 3E-7, "JWST MIRI", color='red', fontsize=fs3)
plt.text(8.8, 4.5E-5, "Spitzer", color='orange', fontsize=fs3)
pretty_plot()
jrr.plot.force_axisticks_linear((ax.xaxis,))
plt.show()
fig.savefig(pandir + 'phot_plot.pdf')

# Now, plot spectroscopy
# If we can grab the spectral resolution, can plot limiting sensitivity to an unresolved line (in erg/s/cm^2).  What a spectroscopist wants.
fig2, ax2 = plt.subplots()
modes  = ['nirspec_msa', 'nirspec_ifu', 'miri_mrs'] #, 'miri_lrs'] # 'nircam_wfgrism']

# make a smaller dict of sensitivities for the spectroscopic modes only
df_specsens = {}
specmodes = ['nirspec_fs','nirspec_ifu', 'nirspec_msa',  'miri_mrs']# 'nircam_wfgrism', 'miri_lrs',  'niriss_wfss','niriss_soss'
for specmode in specmodes:
    df_specsens[specmode + "_sensitivity"] = df_sens[specmode + "_sensitivity"].copy(deep=True)
for ii, mode in enumerate(specmodes) :
    thismode = mode + '_sensitivity'
    df_specsens[thismode]['mode'] = mode
    df_specsens[thismode]['R'] = df_specsens[thismode]['wavelengths'] * 0  # initialize
    df_specsens[thismode]['PSSL'] = df_specsens[thismode]['wavelengths'] * 0  # initialize
    for index, row in df_specsens[thismode].iterrows() :     # Go row by row, for each grating, order, etc
        #plt.plot(row.wavelengths, row.lim_fluxes*1E-3, color='red')
        if 'miri' in thismode :    disp_key = 'jwst_miri_' + row.aperture + '-' + row.disperser
        elif 'nircam' in thismode: disp_key = 'jwst_nircam'
        else:                      disp_key = 'jwst_' + re.split('_', thismode)[0] + '_' + row.disperser
        print("DEBUGGING", disp_key)
        rebinnedR = jrr.spec.rebin_spec_new(df_disp[disp_key]['WAVELENGTH'], df_disp[disp_key]['R'], row.wavelengths)
        df_specsens[thismode].loc[index, 'R'] = rebinnedR
    # Now, need to convert to limiting line flux for an unresolved line.
    # Eqn 2.7 of the Spitzer IRS manual:  PSSL = 3.0E-15 PSSC / (R lambda), where PSSC is in units of mJy, and lambda is in microns.  PSSL is in W/m^2
        df_specsens[thismode].loc[index, 'PSSL'] = row.lim_fluxes * 3.0E-15 / ( row.R  * row.wavelengths) * 1E3  #The last 1E3 converts to cgs
        plt.plot(row.wavelengths, row.PSSL)
        plt.yscale('log')
        plt.xlabel("wavelength (micron)")
        plt.ylabel("limiting unresolved line flux (erg/s/cm^2)")
bigspec_df = pandas.concat(df_specsens)

pandeia_ver = "v1.5.0"
outfile = 'pandeia_lineflux_sensitivities_' + pandeia_ver + '.readme'  
header_text =  '# PSSL (erg/s/cm^2) is the line flux detectable for a spectrally unresolved line, point source, at SNR=10, in 1E4s, from Pandeia ' + pandeia_ver + '\n'
header_text += '# Adapting Eqn 2.7 of the Spitzer IRS manual:  PSSL (erg/s/cm^2) = 3.0E-15 PSSC / (R lambda) * 1E3,\n'
header_text += '# where PSSC is in units of mJy, lambda is in microns, and PSSL is in erg/s/cm^2\n'
header_text += '# In Python pandas, read this file as:  df = pandas.read_csv(\'' + outfile + '\', index_col=[0,1], comment=\'#\')\n'
header_text += '# Jane.Rigby@nasa.gov, 20 Aug. 2020.\n#\n'
jrr.util.put_header_on_file('/dev/null', header_text, outfile)
bigspec_df.to_pickle(re.sub('readme', 'pcl', outfile))
# read this as    df_sens = pandas.read_pickle('pandeia_lineflux_sensitivities_v1.5.0.pcl')
 
#nirspec_setttings = [x for x in df_disp.keys() if 'nirspec' in x] 
    



