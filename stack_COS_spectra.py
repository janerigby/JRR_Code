''' Want to stack the z~0-0.1 COS spectra w J Chisholm.
    modeled after mage_stack_redo.py (which is hardcoded for mage kinds of spectra), 
    and esi_stack_spectra.py, which is more general but did a straight stack rather than weighted avg.
    jane.rigby@nasa.gov, 4/2017
'''

indir = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/MagE_atlas/Contrib/Chisholm16/raw/"  # Running in this dir.
#indir = "/Users/jrrigby1/Dropbox/MagE_atlas/Contrib/Chisholm16/raw/" # on milk

import jrr
import glob
import re
from os.path import basename
import pandas
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.stats import gaussian_fwhm_to_sigma
import astropy.convolution

def get_redshifts(filename='stack_samp.txt') :
    colnames = ('galname', 'zz', 'has_G130M', 'has_G160M')  #has=1 is True
    zzlist = pandas.read_csv(filename, comment="#", names=colnames)
    return(zzlist)

def get_MW_geocoronal_linelist(infile='MW_geo.txt') :
    colnames = ('linelabel', 'obswav')
    linelist = pandas.read_csv(infile, comment="#", names=colnames)
    linelist['restwav']  = linelist['obswav']  
    linelist['zz'] = 0.0  # Milky Way and Geocoronal lines have zero redshift
    linelist['unity'] = True
    return(linelist)

def get_stacked_linelist(vmask=300., Lyamask=1000.) :  # need a linelist to flag lines before fitting continuum in spec.fit_autocont
    # Mask most lines by +- vmask km/s ;  mask HI (Lya) by +- Lyamask km/s
    (spec_path, line_path) = jrr.mage.getpath('released')
    (LL, foobar) = jrr.mage.get_linelist(line_path + "stacked.linelist")
    LL['vmask'] = vmask
    LL.loc[LL['lab1'] == 'H I', 'vmask'] = Lyamask
    LL2 = LL.loc[LL['type'].ne("PHOTOSPHERE")]  # don't mask photospheric lines
    return(LL2)
     
def get_spectra(filenames) :
    df = {}
    for this in filenames :
        #print "loading file ", this
        colnames = ('obswave', 'flam', 'flam_u')   # Units are 1E-14 erg/s/cm^2/A
        df[this] = pandas.read_table(this, delim_whitespace=True, comment="#", names=colnames)
        prepare_spectra(df)
    return(df)  # return a dictionary of dataframes of spectra

def prepare_spectra(df) :   # Process the spectra dataframes for use later
    for this in df.keys() :
        df[this]['Nfiles'] = 1 # N of exposures that went into this spectrum
        df[this]['badmask']   = False
        df[this]['linemask']  = False
        df[this]['contmask']  = False
        df[this]['stackmask'] = False
        df[this].loc[df[this]['flam_u'] == 0.00, 'badmask'] = True   # flag bad values
        jrr.spec.calc_dispersion(df[this], colwave='obswave', coldisp='obsdisp')
    return(0)
    
def blur_the_spectra(df, intermed_wave, new_wave, R1=2.0E4, R2=3500.) :  # Convolve and rebin the COS spectra to simulate MagE data.
    idf = {}  ; ndf = {}  # intermediate and new data frames
    for this in df.keys() : #for each spectrum dataframe
        # Input is oversampled spectrum with high spectral resolution R1.
        # This function rebins to critically sampled, then convolves w Gaussian, then rebins again to R=R2.
        # Start by rebinning to an intermediate wavelength array. R=2E4 nyquist=2.2.
        idf[this] = pandas.DataFrame(data=pandas.Series(intermed_wave), columns=('obswave',))
        idf[this]['flam']   = jrr.spec.rebin_spec_new(df[this]['obswave'], df[this]['flam'],  intermed_wave)
        idf[this]['flam_u'] = jrr.spec.rebin_spec_new(df[this]['obswave'], df[this]['flam_u'], intermed_wave)
        # Intermediate wavelength array has pixels to keep R=constant. It is critically sampled.  Now, convolve w Gaussian,
        # assuming destination is also a critically-sampled array
        kern_fwhm = np.sqrt((R1/R2)**2 - 1.)
        kern = astropy.convolution.Gaussian1DKernel(gaussian_fwhm_to_sigma * kern_fwhm)
        idf[this]['smooth']   = astropy.convolution.convolve(idf[this]['flam'].as_matrix(),   kern, boundary='fill', fill_value=np.nan)
        idf[this]['smooth_u'] = astropy.convolution.convolve(idf[this]['flam_u'].as_matrix(), kern, boundary='fill', fill_value=np.nan)
        # Should we be doing soemthing smarter here with smooth_u?
        ndf[this] = pandas.DataFrame(data=pandas.Series(new_wave), columns=('obswave',))
        ndf[this]['flam']   = jrr.spec.rebin_spec_new(idf[this]['obswave'], idf[this]['smooth'],   new_wave)
        ndf[this]['flam_u'] = jrr.spec.rebin_spec_new(idf[this]['obswave'], idf[this]['smooth_u'], new_wave)
    prepare_spectra(idf) ; prepare_spectra(ndf)
    return(idf, ndf)

 
def flag_geoMW_lines(df, geoMW_linelist, vmask, Lyamask) :
    geoMW_linelist['vmask'] = vmask
    geoMW_linelist.loc[geoMW_linelist['linelabel'] == 'HI', 'vmask'] = Lyamask
    for this in df.keys() :  # For each spectrum dataframe
        jrr.spec.flag_near_lines(df[this], geoMW_linelist, colv2mask='vmask', colwave='obswave', colmask='geoMWlinemask')
    return(0)
#        line_cen = np.array(geoMW_linelist.obswav)
#        vmask_ar = np.array(geoMW_linelist.vmask)
#        line_lo   = line_cen * (1.0 - vmask_ar/2.997E5)
#        line_hi   = line_cen * (1.0 + vmask_ar/2.997E5)
#        temp_wave = np.array(df[this]['obswave'])
#        temp_mask = np.zeros_like(temp_wave).astype(np.bool)
#        for ii in range(0, len(line_lo)) :    # doing this in observed wavelength
#            temp_mask[np.where( (temp_wave > line_lo[ii]) & (temp_wave < line_hi[ii]))] = True
#        df[this]['geoMWlinemask'] = temp_mask  # Using temp np arrays is much faster than writing repeatedly to the pandas df
    return(0)

def from_galname_get_redshift(galname, zzlist) :
    redshift = zzlist[zzlist['galname'].eq(galname)]['zz']
    return(redshift.values[0])

def from_filename_get_props(this, zzlist) :
    m = re.search('(\S+)-(G\w\w\wM)', this)
    galname = m.group(1)
    grating = m.group(2)
    redshift = from_galname_get_redshift(galname, zzlist)
    return(galname, grating, redshift)

def deredshift_all_spectra(df, zzlist) :
    for this in df.keys() :
        (galname, grating, redshift) =  from_filename_get_props(this, zzlist)
        jrr.spec.convert2restframe_df(df[this], redshift, units='flam', colwave='obswave', colf='flam', colf_u='flam_u')
        jrr.spec.calc_dispersion(df[this], colwave='rest_wave', coldisp='rest_disp')
    return()

def wrapper_fit_continuua(df, smooth_length, debug=False) : 
    print "Fitting Continuua (target  grating redshift)"
    for this in df.keys() : 
        (galname, grating, redshift) =  from_filename_get_props(this, zzlist)
        print "      ", galname, grating, redshift
        df[this].loc[df[this]['badmask'],        'contmask'] = True   # Prepare mask for continuum fitting
        df[this].loc[df[this]['geoMWlinemask'],  'contmask'] = True
        df[this].loc[df[this]['badmask'],        'stackmask'] = True  # Prepare mask for stacking
        df[this].loc[df[this]['geoMWlinemask'],  'stackmask'] = True

        boxcar = jrr.spec.get_boxcar4autocont(df[this], smooth_length)
        jrr.spec.fit_autocont(df[this], LL, 0.0, colv2mask='vmask', boxcar=boxcar, flag_lines=True, colwave='rest_wave', colf='rest_flam', colmask='contmask', colcont='rest_flam_autocont')
        df[this]['rflam_norm']  = df[this]['rest_flam']  / df[this]['rest_flam_autocont']    # Normalize by continuum here.
        df[this]['rflamu_norm'] = df[this]['rest_flam_u']/ df[this]['rest_flam_autocont']
        if debug : # Plot the continuum fits
            fig = plt.figure(figsize=(20,4))
            plt.clf()
            plt.plot(df[this]['rest_wave'], df[this]['rest_flam'], color='black', linewidth=1)
            used_for_cont = df[this].loc[df[this]['contmask'] == False]
            plt.plot(used_for_cont['rest_wave'], used_for_cont['rest_flam'], color='grey', linewidth=1)
            plt.ylim(0,3)
            plt.plot(df[this]['rest_wave'], df[this]['rest_flam_autocont'], color='green', linewidth=1)
            #plt.plot(df[this]['rest_wave'], df[this]['contmask'], color='red', linewidth=1)
            plt.xlabel("rest wavelength")
            print "Plotting", this, galname, grating, redshift
            plt.show()
    return(0)


# Set parameters
#                          Velocity to mask features for continuum fitting, and stacking, +-, in km/s
geoMW_vmask      = 300.  # Mask Geocoronal and Milky Way emission for most lines.  Mask for continuum and stacking 
Lya_geoMW_mask   = 3000. # Larger mask for Milky Way damped Lyman alpha.  Mask for continuum and stacking.
starburst_vmask  = 500.  # mask lines when fitting continuum.  Leave in stack
Lya_starburst    = 2000. # Mask Lya when fitting continuum.  Leave in stack.
smooth_length_origR = 20. # rest-frame Angstroms.  Smoothing length for continuum fitting.
smooth_length_likeMage= 50.
R2E4_wavear  = jrr.spec.make_wavearray_constant_resoln(1145, 1790, 2E4, 2.2)  # desired wavelength array for COS, nyquist sampled at R=2E4
R3500_wavear = jrr.spec.make_wavearray_constant_resoln(1145, 1790, 3500, 2.2) # appropriate output wavelength array for simulated-MagE stack

# NOW DO THINGS
zzlist = get_redshifts()
filenames =  [ basename(x) for x in glob.glob(indir + "*G*M") ]    # Find the files. 
print "Loading spectra"
#df = get_spectra(filenames)                 # load spectra into dict of dataframes
df = get_spectra(filenames[:20])  # DEBUGGING, subset is faster *****
geoMW_linelist =  get_MW_geocoronal_linelist()  # Line list, prepartory to masking geocoronal & MW lines
flag_geoMW_lines(df, geoMW_linelist, geoMW_vmask, Lya_geoMW_mask)  # Flag the geocoronal and MW emission.
deredshift_all_spectra(df, zzlist)
LL = get_stacked_linelist(starburst_vmask, Lya_starburst) 
wrapper_fit_continuua(df, smooth_length_origR, debug=False)

# To make a stack w spectral resoln matched to MagE, repeat this for ndf
(idf, ndf) = blur_the_spectra(df, R2E4_wavear, R3500_wavear, R1=2.0E4, R2=3500.)
flag_geoMW_lines(ndf, geoMW_linelist, geoMW_vmask, Lya_geoMW_mask)  # Flag the geocoronal and MW emission.
deredshift_all_spectra(ndf, zzlist)
wrapper_fit_continuua(ndf, smooth_length_likeMage, debug=False)


# STACKING STARTS HERE

# Make stack, with output array nyquist sampled for COS
stacked =  jrr.spec.stack_spectra(df, colwave='rest_wave', colf='rflam_norm', colfu='rflamu_norm', colmask='stackmask', output_wave_array=R2E4_wavear)
stacked_output = "stacked_COS_spectrum_R2E4.csv"
stacked.to_csv(stacked_output, index=False)

# Make stack, with spectral resolution matched to MagE
stacked2 =  jrr.spec.stack_spectra(ndf, colwave='rest_wave', colf='rflam_norm', colfu='rflamu_norm', colmask='stackmask', output_wave_array=R3500_wavear)
stacked_output = "stacked_COS_spectrum_R3500.csv"
stacked2.to_csv(stacked_output, index=False)

plt.close("all")
ax = stacked.plot(x='rest_wave', y='fweightavg', linewidth=0.5, color='k', drawstyle="steps-post")
stacked.plot(x='rest_wave', y='fweightavg_u', linewidth=1, color='grey', ax=ax, drawstyle="steps-post")
plt.plot(stacked.rest_wave, stacked.Ngal / stacked.Ngal.max(), color='green', linewidth=1)
plt.plot((1000,2000), (1,1), color='grey', linewidth=1)
plt.ylim(0,1.5)
plt.show()

plt.close("all")
ax = stacked2.plot(x='rest_wave', y='fweightavg', linewidth=0.5, color='k', drawstyle="steps-post")
stacked2.plot(x='rest_wave', y='fweightavg_u', linewidth=1, color='grey', ax=ax, drawstyle="steps-post")
plt.plot(stacked2.rest_wave, stacked2.Ngal / stacked2.Ngal.max(), color='green', linewidth=1)
plt.plot((1000,2000), (1,1), color='grey', linewidth=1)
plt.ylim(0,1.5)
plt.show()

# This works!  Columns are a bit ugly, but it works





# For debugging, plot a bunch of the spectra, to see if MW, Geocoronal features were masked out
if False :    
    for this in df.keys() :
        (galname, grating, redshift) =  from_filename_get_props(this, zzlist)
        fig = plt.figure(figsize=(20,4))
        ax1 = fig.add_subplot(111)
        masked = df[this][df[this]['geoMWlinemask']]
        unmasked = df[this][~df[this]['geoMWlinemask']]
        print "Plotting", this
        plt.plot(unmasked.obswave, unmasked.flam, color='black', linewidth=0.5)
        plt.plot(unmasked.obswave, unmasked.flam_u, color='orange', linewidth=0.5)
        plt.scatter(masked.obswave, masked.flam, color='red', s=1)
        plt.scatter(geoMW_linelist.obswave, geoMW_linelist.unity*0.1, color='blue')
        if   grating == "G130M" :  x1=1200; x2=1350
        elif grating == "G160M" :  x1=1500; x2=1560
        else :                     x1=1200; x2=1800
        ax1.set_xlim(x1,x2)
        plt.ylim(0,2)
        plt.show()


