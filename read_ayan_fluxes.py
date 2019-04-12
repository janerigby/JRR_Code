from __future__ import print_function
import jrr
import pandas
import numpy as np
from os.path import expanduser

def trust_this_line(df,  signi_thresh=3., SNR_thresh=1., newcol='Realdet'):
    df[newcol] =  (df['EW_signi'] > signi_thresh)  &  (df['f_line'] / df['f_line_u'] > SNR_thresh)
    # Above translated:  We should believe the line if it's 3 sigma significant, and has a sane SNR
    return(0)
    
homedir = expanduser("~")
flux_file = homedir + "/Dropbox/MagE_atlas/Contrib/EWs/emission/allspec_fitted_emission_linelist.txt"
df = pandas.read_table(flux_file, delim_whitespace=True, comment="#", thousands=',')
trust_this_line(df)
signifi = 3.0 # sigma
sig2 = 3
Trust = 'Realdet'

# Do this as a loop, it works
#lines = df['line_lab'].unique()  # unique spectral features
#for thisline in lines :
#    print thisline, (df['line_lab'].eq(thisline) & df['EW_significance'].gt(signifi)).sum()

# How many detections of each line?
bylines = pandas.DataFrame(data=None, index=df['line_lab'].unique())  # unique spectral features
bylines['Ndet'] = np.NaN
for thisline in bylines.index :
    Ndet = ((df['line_lab'].eq(thisline) & df[Trust])).sum()
    bylines['Ndet'][thisline] = Ndet
bylines.sort_values(by='Ndet', ascending=False, inplace=True)
print(bylines.head(60))

# How many 3 sigma detections of each galaxy?
bygals = pandas.DataFrame(data=None, index=df['label'].unique()) # unique gal name
bygals['Ndet'] = np.NaN
for thisgal in bygals.index :
    Ndet = (df['label'].eq(thisgal) & df[Trust]).sum()
    bygals['Ndet'][thisgal] = Ndet
bygals.sort_values(by='Ndet', ascending=False, inplace=True)
print(bygals.head(20))

# Below seems super-dodgy, need to update this
lines_to_check = ('CIII1908', "OIII1666", "OIII1660", 'OIII2320')
for checkline in lines_to_check :
    gooddet = df[df['line_lab'].eq(checkline) & df['EW_signi'].gt(sig2)]
    gooddet.sort_values(by='EW_signi', ascending=False, inplace=True)
    print(gooddet.head(30))

