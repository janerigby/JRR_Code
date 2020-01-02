from __future__ import print_function
from builtins import range
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter
import glob
from os.path import basename, expanduser, exists
import pandas 

#pandir = '/Users/jrrigby1/MISSIONS/JWST/Sens/pandeia_sensitivities_1.3/'
pandir = '/Users/jrrigby1/MISSIONS/JWST/Sens/pandeia_sensitivities_1.5/'
fs1 = 14 ; fs2 = 18 ; fs3 = 22 # fontsizes

def load_photometry() :
    # Adding the comparison to other observatories. Switching from gnuplot to python
    phot_codes = {'NIRSpec' : 0, 'NIRCam': 1, 'MIRI' : 2, 'HST' : 4, 'WISE' : 5, 'Spitzer' : 6, 'gemini' :10 , 'herschel' : 11, 'sofia' : 12}
    #13=alma(cycle0), 14=alma(finished), 15=VLA, 16=EVLA
    names = ('name', 'wave', 'limfnu', 'code')                
    phot_df = pandas.read_table(pandir + '../jwst-phot.dat', comment="#", delim_whitespace=True, usecols=[0,1,2,3], names=names)
    phot_df['wave'] =  pandas.to_numeric(phot_df['wave'])
    return(phot_df, phot_codes)

def load_spectroscopy() :
    # Can I grab the continuum sensitivity folder, so it's easy?  Not what users want, though.  Have asked for spectral R
    # so that I can convert from continuum sensitivity to line sensitivity
    return(0)


def force_axisticks_linear(whichaxis) :
    for axis in whichaxis :  #ax.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
       
def pretty_plot() :
    plt.xlabel("wavelength (micron)", fontsize=fs3)
    plt.ylabel('limiting flux density (Jy)', fontsize=fs3)
    #plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()


        
plt.ion()
# Load Klaus's sensitivities
filenames = [ basename(x) for x in glob.glob(pandir + '*npz')]
pan = {} # empty dictionary of pandeia files (in weird numpy format)
df = {}  # empty dict of dataframes
for filename in filenames :
    (rootname, suffix) = filename.split('.')
    pan[rootname] = np.load(pandir + filename, allow_pickle=True)

    uglyintermed = []
    for ii in range(0, len(list(pan[rootname].items()))) :
        uglyintermed.append( list(pan[rootname].items())[ii][1].flatten())  #oh god it burns
    tempdf = pandas.DataFrame(data=uglyintermed)
    #tempdf = pandas.DataFrame(data=[pan[rootname].items()[0][1].flatten(), pan[rootname].items()[1][1].flatten(), pan[rootname].items()[2][1].flatten(), pan[rootname].items()[3][1].flatten(), pan[rootname].items()[4][1].flatten()])
    df[rootname] = tempdf.transpose()
    df[rootname].columns = list(pan[rootname].keys())

print("I think I imported all of Klaus's files:", list(df.keys()))

plt.clf()
fig, ax = plt.subplots(figsize=(8.5,6))
plt.xticks(fontsize=fs3)
plt.yticks(fontsize=fs3)
# MAKE A PHOTOMETRY PLOT
# First, plot Klaus's values.  May need to screen for wide filters
modes  = ['nircam_sw', 'nircam_lw', 'miri_imaging']
labels = ['JWST NIRCam', '_nolegend_', 'JWST MIRI']
# Plot imaging
for ii, mode in enumerate(modes) :
    thekey = mode + '_sensitivity'
    subset = df[thekey].loc[df[thekey]['configs'].astype(str).str.contains('w')]  # Only Wide (broadband) filters  
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
force_axisticks_linear((ax.xaxis,))
plt.show()
fig.savefig(pandir + 'phot_plot.pdf')

# Now, plot spectroscopy 
fig2, ax2 = plt.subplots()
modes  = ['nirspec_msa', 'miri_mrs'] #, 'miri_lrs'] # 'nircam_wfgrism']
for ii, mode in enumerate(modes) :
    thekey = mode + '_sensitivity'
    for index, row in df[thekey].iterrows() :     # Ah, shit, need go row by row.
        plt.plot(row.wavelengths, row.lim_fluxes*1E-3, color='red')
plt.xlim(0.5, 30)
plt.ylim(3E-9, 1E-2)
pretty_plot()
force_axisticks_linear((ax2.xaxis,))
#


