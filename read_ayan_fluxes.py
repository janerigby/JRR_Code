import jrr
import pandas
import numpy as np
from os.path import expanduser

homedir = expanduser("~")
flux_file = homedir + "/Dropbox/MagE_atlas/Contrib/EWs/fitted_emission_list_allspec.txt"
df = pandas.read_table(flux_file, delim_whitespace=True, comment="#")

# Do this as a loop, it works
#lines = df['line_lab'].unique()  # unique spectral features
#for thisline in lines :
#    print thisline, (df['line_lab'].eq(thisline) & df['EW_significance'].gt(3.0)).sum()

# Trying to get this to work as a pandas data frame, it breaks
lines = pandas.DataFrame(data=None, index=df['line_lab'].unique())  # unique spectral features
lines['Ndet'] = np.NaN
for thisline in lines.index :
    Ndet = (df['line_lab'].eq(thisline) & df['EW_significance'].gt(3.0)).sum()
    lines['Ndet'][thisline] = Ndet 

lines.sort_values(by='Ndet', ascending=False, inplace=True)
print lines.head(60)

