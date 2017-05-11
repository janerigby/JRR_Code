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

def get_the_redshifts(filename='stack_samp.txt') :
    colnames = ('galname', 'zz', 'has_G130M', 'has_G160M')  #has=1 is True
    zzlist = pandas.read_csv(filename, comment="#", names=colnames)
    return(zzlist)

def get_MW_geocoronal_linelist(infile='MW_geo.txt') :
    colnames = ('linelabel', 'obswave')
    linelist = pandas.read_csv(infile, comment="#", names=colnames)
    linelist['unity'] = True
    return(linelist)

def get_stacked_linelist() :  # needed for fit_autocont
    (spec_path, line_path) = jrr.mage.getpath('released')
    (LL, foobar) = jrr.mage.get_linelist(line_path + "stacked.linelist") 
    return(LL)
     
def get_the_spectra(filenames) :
    df = {}
    for this in filenames :
        print "loading file ", this
        colnames = ('obswave', 'flam', 'flam_u')   # Units are 1E-14 erg/s/cm^2/A
        df[this] = pandas.read_table(this, delim_whitespace=True, comment="#", names=colnames)
        df[this]['Nfiles'] = 1 # N of exposures that went into this spectrum
        df[this]['badmask'] = False
        df[this]['linemask'] = False
        df[this]['contmask'] = False
        df[this].loc[df[this]['flam_u'] == 0.00, 'badmask'] = True   # flag bad values
        jrr.spec.calc_dispersion(df[this], colwave='obswave', coldisp='obsdisp')
    return(df)  # return a dictionary of dataframes of spectra

def flag_geoMW_lines(df, geoMW_linelist, vmask, Lyamask) :
    geoMW_linelist['vmask'] = geoMW_linelist['unity'] * vmask
    geoMW_linelist.loc[geoMW_linelist['linelabel'] == 'HI', 'vmask'] = Lyamask
    for this in df.keys() :  # For each spectrum dataframe
        line_cen = np.array(geoMW_linelist.obswave)
        vmask_ar = np.array(geoMW_linelist.vmask)
        line_lo   = line_cen * (1.0 - vmask_ar/2.997E5)
        line_hi   = line_cen * (1.0 + vmask_ar/2.997E5)
        temp_wave = np.array(df[this]['obswave'])
        temp_mask = np.zeros_like(temp_wave).astype(np.bool)
        for ii in range(0, len(line_lo)) :    # doing this in observed wavelength
            temp_mask[np.where( (temp_wave > line_lo[ii]) & (temp_wave < line_hi[ii]))] = True
        df[this]['geoMWlinemask'] = temp_mask  # Using temp np arrays is much faster than writing repeatedly to the pandas df
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

# Velocity to mask, +-, in km/s
geoMW_vmask      = 300.  # Mask Geocoronal and Milky Way emission for most lines  
Lya_geoMW_mask   = 3000. # Larger mask for Milky Way damped Lyman alpha
starburst_vmask  = 300.  # mask lines when fitting continuum

# NOW DO THINGS
zzlist = get_the_redshifts()
#filenames =  [ basename(x) for x in glob.glob(indir + "*G*M") ]    # Find the files. 
filenames =  [ basename(x) for x in glob.glob(indir + "NGC*G*M") ]  # TEMP KLUDGE FOR SPEED
df = get_the_spectra(filenames)                 # load spectra into dict of dataframes
geoMW_linelist =  get_MW_geocoronal_linelist()  # Line list, prepartory to masking geocoronal & MW lines
flag_geoMW_lines(df, geoMW_linelist, geoMW_vmask, Lya_geoMW_mask)  # Flag the geocoronal and MW emission.
deredshift_all_spectra(df, zzlist)

print "Experimenting from here below:"
LL = get_stacked_linelist() 
last = "SBS1415+437-G130M"
for this in df.keys()[0:5] :
    fig = plt.figure(figsize=(20,4))
    (galname, grating, redshift) =  from_filename_get_props(this, zzlist)
    df[this].loc[df[this]['badmask'],        'contmask'] = True   # Prepare mask for continuum fitting
    df[this].loc[df[this]['geoMWlinemask'],  'contmask'] = True
    smooth_length = 25.
    boxcar = jrr.spec.get_boxcar4autocont(df[this], smooth_length)
    # Need to add larger vmask for Lya here.... ****
    (smooth1, smooth2) = jrr.spec.testing_fit_autocont(df[this], LL, 0.0, vmask=starburst_vmask, boxcar=boxcar, flag_lines=True, colwave='rest_wave', colf='rest_flam', colmask='contmask', colcont='rest_flam_autocont')
    # Plot the results:
    plt.clf()
    plt.plot(df[this]['obswave'], df[this]['rest_flam'], color='black', linewidth=1)
    used_for_cont = df[this].loc[df[this]['contmask'] == False]
    plt.plot(used_for_cont['obswave'], used_for_cont['rest_flam'], color='grey', linewidth=1)
    plt.ylim(0,3)
    plt.plot(df[this]['obswave'], df[this]['rest_flam_autocont'], color='green', linewidth=1)
    plt.xlabel("Observed wavelength")
    print "Plotting", this, galname, grating, redshift
    plt.show()


#stacked =  jrr.spec.stack_spectra(df, straight_sum=True, colwav='rest_wave', colf='rest_flam', colfu='rest_flam_u', stacklo=1144., stackhi=1787., disp=1.0)  #disp is a dummy value*&&*****

# To do next:
# Normalize the spectra.  Let's try autofitcont first.
# Then stack the spectra, using jrr.spec.stack_spectra().  buggy, solving now...

stacked =  jrr.spec.stack_spectra(df, straight_sum=True, colwav='rest_wave', colf='rest_flam', colfu='rest_flam_u', stacklo=1200., stackhi=1400., disp=2.0)  # testing.  Not right stacklo, stackhi, dispersion
# John says should bin by ~7 pixels.  COS spectra are over-sampled.  Do in constant R.



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


