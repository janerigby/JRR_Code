''' Want to stack the z~0-0.1 COS spectra w J Chisholm.
    modeled after mage_stack_redo.py (which is hardcoded for mage kinds of spectra), 
    and esi_stack_spectra.py, which is more general but did a straight stack rather than weighted avg.
    jane.rigby@nasa.gov, 4/2017
'''

indir = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/MagE_atlas/Contrib/Chisholm16/raw/"  # Running in this dir.

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

def get_the_spectra(filenames) :
    df = {}
    for thisfile in filenames :
        print "loading file ", thisfile
        colnames = ('obswave', 'flam', 'flam_u')   # Units are 1E-14 erg/s/cm^2/A
        df[thisfile] = pandas.read_table(thisfile, delim_whitespace=True, comment="#", names=colnames)
        df[thisfile]['Nfiles'] = 1 # N of exposures that went into this spectrum
        df[thisfile]['flam_u'][df[thisfile]['flam_u'].eq(0.000)] = 1E6
    return(df)  # return a dictionary of dataframes of spectra

def flag_near_lines(df, geoMW_linelist) :
    vmask = 300.  # velocity to mask, +-, km/s   # checked that it masks geocoronal emission
    Lyamask = 1000.
    for thisfile in df.keys() :  # For each spectrum dataframe
        line_cen = np.array(geoMW_linelist.obswave)
        line_lo   = line_cen * (1.0 - vmask/2.997E5)
        line_hi   = line_cen * (1.0 + vmask/2.997E5)
        temp_wave = np.array(df[thisfile]['obswave'])
        temp_mask = np.zeros_like(temp_wave).astype(np.bool)
        for ii in range(0, len(line_lo)) :    # doing this in observed wavelength
            temp_mask[np.where( (temp_wave > line_lo[ii]) & (temp_wave < line_hi[ii]))] = True
        df[thisfile]['linemask'] = temp_mask  # Using temp np arrays is much faster than writing repeatedly to the pandas df
    return(0)

def from_filename_get_galname_grating(filename) :
    m = re.search('(\S+)-(G\w\w\wM)', thisfile)
    galname = m.group(1)
    grating = m.group(2)
    return(galname, grating)

def from_galname_get_redshift(galname, zzlist) :
    redshift = zzlist[zzlist['galname'].eq(galname)]['zz']
    return(redshift.values[0])

# NOW DO THINGS
zzlist = get_the_redshifts()
filenames =  [ basename(x) for x in glob.glob(indir + "*G130M") ]    # Find the files.   curerntly doing G130M only
df = get_the_spectra(filenames)  # load spectra into dict of dataframes
geoMW_linelist =  get_MW_geocoronal_linelist()  # Will need this to mask the geocoronal lines
flag_near_lines(df, geoMW_linelist)  # Flag the geocoronal and MW emission.


first = df.keys()[0]
#ax = df[first].plot(x='obswave', y='flam')
figsize=(20,4)
for thisfile in df.keys() :
    (galname, grating) =  from_filename_get_galname_grating(thisfile)
    redshift = from_galname_get_redshift(galname, zzlist)
    print thisfile, galname, grating, redshift
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(111)
#    ax2 = ax1.twiny()
    #df[thisfile].plot(x='obswave', y='flam', ax=ax, linewidth=1)
    masked = df[thisfile][df[thisfile]['linemask']]
    unmasked = df[thisfile][~df[thisfile]['linemask']]
    print "Plotting", thisfile
    plt.plot(unmasked.obswave, unmasked.flam, color='black', linewidth=0.5)
    plt.plot(unmasked.obswave, unmasked.flam_u, color='orange', linewidth=0.5)
    plt.scatter(masked.obswave, masked.flam, color='red', s=1)
    plt.scatter(geoMW_linelist.obswave, geoMW_linelist.unity*0.1, color='blue')
    x1=1200; x2=1350
    ax1.set_xlim(x1,x2)
#    ax2.set_xlim(x1/(1.+redshift), x2/(1.+redshift))
    plt.ylim(0,2)
    plt.show()

# Notes to future Jane:  Here's what to do to finish this
# Add extra masking to Lya, try +-1000 km/s to get rid of MW DLAs.
# Check whether that removed the DLA signature.  If not, change mask window
# If yes, move on, convert to rest frame, and stack the spectra.

