import shutil
import jrr
import glob
import re
from os.path import basename
import pandas
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy import constants, units
from astropy.io import fits
# Needed for barycentric correction
from astropy import coordinates
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation


def get_observation_datetimes() :
    infile  = "obslog.txt"
    obslog = pandas.read_table(infile, delim_whitespace=True, comment="#")
    obslog.set_index("object", inplace=True, drop=True)  # change index
    return(obslog)


# Apply the barycentric correction:
mmt = EarthLocation.of_site('mmt')
obslog = get_observation_datetimes()
prefix = "Before_Barycencor_ascii"
filenames = [ basename(x) for x in glob.glob(prefix+"*txt") ]
print "DEBUGGING", filenames
for thisfile in filenames :
    print "Working on", thisfile
    approx_rootname =  re.sub('s', 'sgas', re.split('_', thisfile)[1]) + "_arc"
    # Above willl get UTtime slightly off for redgals, but does not matter for barycentric cor
    thisobs = obslog.loc[approx_rootname]
    my_target = SkyCoord(thisobs['RA'], thisobs['DEC'], unit=(units.hourangle, units.deg), frame='icrs')
    my_time = Time( thisobs['UTDate'] + "T" + thisobs['UTTime_midpt'] , format='isot', scale='utc')
    barycor_vel = jrr.barycen.compute_barycentric_correction(my_time, my_target, location=mmt)
    print "barycor:", barycor_vel
    df = pandas.read_table(thisfile, delim_whitespace=True, comment="#", columns=('obswave','flam', 'dflam'))
    
# Actually run things:

