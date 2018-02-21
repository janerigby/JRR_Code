from os.path import expanduser
import pandas

homedir = expanduser("~")
doublet_file = homedir + "/Dropbox/MagE_atlas/Contrib/Intervening/Doublet_search/found_doublets_SNR4.txt"

df = pandas.read_csv(doublet_file, sep='\t', comment="#", index_col=0)
df['EWr1'] = df['EWobs1'] / (1. + df['zz'])
df['EWr2'] = df['EWobs2'] / (1. + df['zz'])

decimal_places = {'zz' : 5, 'wave1' : 3 , 'wave2' : 3, 'EWr1' : 2, 'EWr2' : 2, 'snr1' : 1, 'snr2' : 1}

df2 = df.round(decimal_places)
del df2['EWobs1']
del df2['EWobs2']
new_col_order = ['flag', 'gal', 'zz', 'doubname', 'wave1', 'wave2', 'EWr1', 'EWr2', 'snr1', 'snr2', 'in_forest']
df3 = df2[new_col_order]

print df3.head()
df3.to_latex("doublet_table.tex")
