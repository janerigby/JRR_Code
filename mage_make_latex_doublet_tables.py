''' This is an unfinished script to make a pretty latex table of doublets found in the Megasaura spectra, 
for the intervening absorber paper.  jrigby, Feb 2019'''

from os.path import expanduser
import pandas
import numpy as np

def prettify_df(df) :
    df.replace("Cosmic\textasciitildeEye", "Cosmic Eye", inplace=True)
    df.replace("Cosmic~Eye", "Cosmic Eye", inplace=True)
    df.replace("rcs0327-E", "RCS-GA 032727$-$132609", inplace=True)
    df.replace("SPT2325", "SPT2325-4111", inplace=True)
    #df.replace("planckarc", "PSZ1-ARC G311.6602-18.4624", inplace=True)
    df.replace("planckarc", "Sunburst Arc", inplace=True)
    df.replace("PSZ0441_slitA", "PSZ0441-0946", inplace=True)
    df.replace("MgII", "Mg II", inplace=True)
    df.replace("CIV", "C IV", inplace=True)
    df.replace("SiIV", "Si IV", inplace=True)
    df['flag'].replace('y', '-', inplace=True)  # Change all the "y" flags (yes a detection) to a dash
    del df['in_forest']
    del df['2doublets']
    return(0) # acts on df

homedir = expanduser("~") + "/"
doublet_dir = "Dropbox/MagE_atlas/Contrib/Intervening/Doublet_search/"
infile1 = homedir + doublet_dir + "Results_16Feb2018/found_doublets_SNR4_JRRedit.txt"  # results for batch1, batch2
infile2 = homedir + doublet_dir + "Results_21Feb2019/found_doublets_batch3_SNR4_JRRedit.txt"  # Results for batch3

dfa = pandas.read_table(infile1, comment="#", delim_whitespace=True)
dfb = pandas.read_table(infile2, comment="#", delim_whitespace=True)
df = pandas.concat([dfa, dfb])

df['EWr1'] = df['EWobs1'] / (1. + df['zz'])
df['EWr2'] = df['EWobs2'] / (1. + df['zz'])

decimal_places = {'zz' : 5, 'wave1' : 1 , 'wave2' : 1, 'EWr1' : 2, 'EWr2' : 2, 'snr1' : 1, 'snr2' : 1}

df2 = df.round(decimal_places)
del df2['EWobs1']  # remove unnecessary columns
del df2['EWobs2']
new_col_order = ['gal', 'pointing', 'zz', 'doubname', 'wave1', 'wave2', 'EWr1', 'EWr2', 'snr1', 'snr2', 'in_forest', '2doublets', 'flag' ]
df3 = df2[new_col_order]
prettify_df(df3)

# Should change the y values in the flag column to blanks

def_detected = df3.loc[df3['flag'].str.contains('-')]
det_blended  = df3.loc[df3['flag'].eq('b')]
iffy         = df3.loc[df3['flag'].str.contains('p')]

# Check sizes of these subsets
thesizes = np.array((def_detected.shape[0], det_blended.shape[0], iffy.shape[0]))
if  np.sum(thesizes) == df3.shape[0] : print "Good, the sum of the subsets is as big as the big df."
else : print "ERROR, sizes do not match! The sum of", thesizes, np.sum(thesizes), "should sum to the total of", df3.shape[0]


#print df3.head()
df3.to_latex("doublet_table.tex", index=False)
