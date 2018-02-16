import jrr
import pandas
import numpy as np
from os.path import expanduser

def trust_this_line(df, line,  signi_thresh=3., SNR_thresh=1.):
    condition1 = df[df['line_lab'].eq(line)].EW_signi.values[0] > signi_thresh
    condition2 = df[df['line_lab'].eq(line)].f_line.values[0]/df[df['line_lab'].eq(line)].f_line_u.values[0] > SNR_thresh)
    return(condition1 and condition2)
    
homedir = expanduser("~")
flux_file = homedir + "/Dropbox/MagE_atlas/Contrib/EWs/emission/allspec_fitted_emission_linelist.txt"
df = pandas.read_table(flux_file, delim_whitespace=True, comment="#", thousands=',')
signifi = 3.0 # sigma
sig2 = 3

# Do this as a loop, it works
#lines = df['line_lab'].unique()  # unique spectral features
#for thisline in lines :
#    print thisline, (df['line_lab'].eq(thisline) & df['EW_significance'].gt(signifi)).sum()

# How many detections of each line?
bylines = pandas.DataFrame(data=None, index=df['line_lab'].unique())  # unique spectral features
bylines['Ndet'] = np.NaN
for thisline in bylines.index :
    Ndet = (df['line_lab'].eq(thisline) & df['EW_signi'].gt(signifi)).sum()
    bylines['Ndet'][thisline] = Ndet
bylines.sort_values(by='Ndet', ascending=False, inplace=True)
print bylines.head(60)

# How many 3 sigma detections of each galaxy?
bygals = pandas.DataFrame(data=None, index=df['label'].unique()) # unique gal name
bygals['Ndet'] = np.NaN
for thisgal in bygals.index :
    Ndet = (df['label'].eq(thisgal) & df['EW_signi'].gt(signifi)).sum()
    bygals['Ndet'][thisgal] = Ndet
bygals.sort_values(by='Ndet', ascending=False, inplace=True)
print bygals.head(20)

lines_to_check = ("OIII1666", "OIII1660", 'OIII2320')
for checkline in lines_to_check :
    gooddet = df[df['line_lab'].eq(checkline) & df['EW_signi'].gt(sig2)]
    gooddet.sort_values(by='EW_signi', ascending=False, inplace=True)
    print gooddet.head(30)

    trust_this_line(df, thisline)

