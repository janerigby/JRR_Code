import jrr
from os.path import expanduser
import pandas
import numpy as np
import matplotlib.pyplot as plt

mage_mode = 'released'
label = 'S1527+0652'

restwaves = np.array((1037.6167,  6564.61,   4862.683))
linenames =         ('O VI 1037', 'H alpha', 'H beta')

(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(sp_mage, resoln, dresoln, LL, zz_syst) = jrr.mage.wrap_open_spectrum(label, mage_mode)

# Spectrum that contains Halpha
colnames=('wave', 'flam', 'noise', 'f1', 'n1', 'f2', 'n2', 'f3', 'n3', 'f4', 'n4', 'f5', 'n5')
filters = ('N4', 'N6', 'N7') # NIRSPEC filters for S1527
sp_nir = {}  # Empty dictionary of dataframes
for filt in filters :
    nirspec_specfile = expanduser("~") + '/Dropbox/MagE_atlas/Rest-Opt/s1527_' + filt + '.spec'
    sp_nir[filt] = pandas.read_table(nirspec_specfile, delim_whitespace=True, comment="#", names=colnames)
    jrr.spec.flam2fnu_df(sp_nir[filt], colwave='wave', colf='flam', colf_u='noise') 
    jrr.spec.convert2restframe_df(sp_nir[filt], zz_syst, units='fnu', colwave='wave', colf='fnu', colf_u='fnu_u') 


jrr.plot.boxplot_Nspectra((sp_mage['wave'],), (sp_mage['fnu'],), (sp_mage['fnu_u'],), (zz_syst,), (linenames[0],), (restwaves[0],), win=2000, ymax=(1E-28,), figsize=(6,4), vel_plot=True)

jrr.plot.boxplot_Nspectra((sp_nir['N7']['wave'],), (sp_nir['N7']['fnu'],), (sp_nir['N7']['fnu_u'],), (zz_syst,), (linenames[1],), (restwaves[1],), win=2000, ymax=(1.5E-26,), figsize=(6,4), vel_plot=True)

jrr.plot.boxplot_Nspectra((sp_nir['N6']['wave'],), (sp_nir['N6']['fnu'],), (sp_nir['N6']['fnu_u'],), (zz_syst,), (linenames[2],), (restwaves[2],), win=2000, ymax=(0.7E-26,), figsize=(6,4), vel_plot=True)


plt.show()
