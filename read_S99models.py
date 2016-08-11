''' This is a script to read in John Chisholm's Starburst99 models, so that I can
fit auto continuum and equivalent widths.'''

from os.path import expanduser
from astropy.table import Table
from astropy.io import fits
import re
import numpy as np
import pandas
import matplotlib.pyplot as plt

homedir = expanduser("~")
thedir = homedir + "/Dropbox/MagE_atlas/Contrib/S99/models/"

def get_t_Z_from_label(label) :
    label = re.sub('t', '',  label)
    (age, Z) = re.split("Z", label)
    return(age, str(Z)) # returns age in Myr and metallicity where 0.02 is solar.

def get_label_from_t_Z(t,Z) :
    return ("t"+t+"Z"+Z)
    
def get_Z_from_filename(filename) :
     metallicity = np.float64(re.sub(".fits", "", re.sub("sb99", "0.", filename)))
     return(str(metallicity))

def convert_age_to_Myr(inage):
    return(((np.float64(inage)/1E6).astype(np.int)).astype(np.str))

        
     
models = ('sb99001.fits', 'sb99004.fits', 'sb99008.fits', 'sb9902.fits', 'sb9904.fits')  # J. Chisholm's models
test = [re.sub("sb99", "0.", x) for x in models]         # list comprehension
metallicities = [re.sub(".fits", "", x) for x in test]   # list comprehension

# Get the list of ages.
mytab = Table.read(thedir + models[0])
ages = convert_age_to_Myr(mytab['AGE'][0,:])

#These are the unique model names, format t(age, Myr)Z(metallicity, where 0.02 is solar)
indices = ["t"+x+"Z"+y for x in ages for y in metallicities]  # all the model names
df = {}  # empty data frame.  Preparing for dictionary addressing.

for thismod in (models[:]) : 
    metallicity =get_Z_from_filename(thismod)
    mytab = Table.read(thedir + thismod)
    (foo, Nages) = mytab['AGE'].shape
    for tt in range(0, Nages) :
        age = convert_age_to_Myr(mytab['AGE'][0,tt])
        label = get_label_from_t_Z(age,metallicity)
        print thismod, metallicity, age, label
        wave = mytab['WAVE'][0,:]
        flam = mytab['FLUX'][0,tt,:]
        # Make a pandas data frame for this time, this metallicity
        df[label] = pandas.DataFrame(data=np.transpose([wave.data, flam.data]), columns=['wave', 'flam'])
        
df_all = pandas.concat(df)

# OK, read into a big dataframe, df_all.  Use as in 3D-HST NB example.  For each, fit autocont, then use
# Ayan's EW fitter.
