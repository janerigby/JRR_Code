from __future__ import print_function
from builtins import range
import numpy as np
from matplotlib import pyplot as plt
import glob
from os.path import basename, expanduser, exists
import pandas 
plt.ion()

# Load Klaus's sensitivities
pandir = '/Users/jrrigby1/MISSIONS/JWST/Sens/pandeia_sensitivities_1.3/'
filenames = [ basename(x) for x in glob.glob(pandir + '*npz')]
pan = {} # empty dictionary of pandeia files (in weird numpy format)
df = {}  # empty dict of dataframes
for filename in filenames :
    (rootname, suffix) = filename.split('.')
    pan[rootname] = np.load(pandir + filename)

    uglyintermed = []
    for ii in range(0, len(list(pan[rootname].items()))) :
        uglyintermed.append( list(pan[rootname].items())[ii][1].flatten())  #oh god it burns
    tempdf = pandas.DataFrame(data=uglyintermed)
    #tempdf = pandas.DataFrame(data=[pan[rootname].items()[0][1].flatten(), pan[rootname].items()[1][1].flatten(), pan[rootname].items()[2][1].flatten(), pan[rootname].items()[3][1].flatten(), pan[rootname].items()[4][1].flatten()])
    df[rootname] = tempdf.transpose()
    df[rootname].columns = list(pan[rootname].keys())

print("I think I imported all of Klaus's files:", list(df.keys()))

# Plot the photometry 
modes = ['nircam_sw', 'nircam_lw', 'miri_imaging']  # replace reptition below with loop
# Plot imaging
for mode in modes :
    thekey = mode + '_sensitivity'
    plt.scatter(df[thekey]['wavelengths'], df[thekey]['lim_fluxes'], label=mode)
plt.xlim(0.6, 30)
plt.ylim(0.0, 0.05)
plt.xlabel("wavelength (micron)")
plt.ylabel('limiting flux density (microJy)')
plt.legend()
plt.show()
