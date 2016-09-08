''' Same as multipanel-spectrum-ids.py, but plot multiple spectra.  Using for individual galaxies
jrigby, Sept 2016. '''

import jrr
import matplotlib.pyplot as plt
import numpy as np
mage_mode = "reduction"  
#mage_mode = "released"
methods = ('bystars', 'byneb')  # method of determining systemic redshift

specs = jrr.mage.getlist_wcont(mage_mode, drop_s2243=True)
short_labels = specs['short_label'].unique()
for label in short_labels :
    (sp, resoln, dresoln, LL, z_systemic,boxcar) = jrr.mage.wrap_open_spectrum(label, mage_mode) 
    the_dfs = [sp]
    the_zzs = [z_systemic]  
    the_pdf =  "PDF_Out2/multipanel_" + label + '.pdf'
    jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title=(label+" "+str(the_zzs)), topfid=(1.1,1))#, waverange=(1000,3000))
plt.clf()
