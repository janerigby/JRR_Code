import jrr
import numpy as np
import pandas
import string
from os.path import expanduser, basename
import glob
import re
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
mage_mode = "released"
plt.ion()

batch23 = jrr.mage.organize_labels('batch23')   # Adding Apr2018 and Aug2018 runs
labels = ('PSZ0441', 'SPT0310', 'SPT2325', 'SPT0356', 'SPT0142')  # Just the 5 targets of this proposal
Z      = (0.10,       0.23,       0.27,     0.22,       0.23 )  # John didnt fit these, ro at least, they're not in hot stars table 3
age    = (19.3,      15.7,        9.2,      8.6,        10.1)
age_dict = dict(zip(labels, age))
Z_dict =   dict(zip(labels, Z))


from_meg = ('rcs0327-E', 'S1226+2152')
(sp, resoln, dresoln, LL, zz_sys, speclist)  = jrr.mage.open_many_spectra(mage_mode, which_list='labels', labels=from_meg, addS99=True)
for ii, key in enumerate(sp.keys()) :
    print "redshift, object", zz_sys[key], key


#jrr.mage.plot_linelist(LL[key], z_systemic=zz_sys[key],  restframe=True)


figdir = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/SGAS-shared/PROPOSALS/Hubble/Cycle27/Friends_o_Meg/Figure/"
plotorder = ('PSZ0441',  'SPT0310', 'SPT0142', 'SPT2325', 'SPT0356')   # from bottom to top
filenames =  [ basename(x) for x in glob.glob(figdir + "*-fit.txt") ]
colnames=('wave', 'obsflux', 'errflux', 'S99fit')
df = {}
#for ii, filename in enumerate(filenames):
#    myre = re.compile('(\S+)-sb99-fit.txt')
#    key = myre.match(filename).groups()[0]

fig = plt.figure(figsize=(11,12))
ax = fig.add_subplot(1, 1, 1)
nudgear = (1, 1, 1, 1, 1, 1)
scalef = 1.2
nudge  = 1.1
key = from_meg[-1]
plt.plot(sp[key].rest_wave, (sp[key].rest_fnu / sp[key].rest_fnu_autocont) - nudge, label=key, linestyle='steps-mid', color='k')  # Plot S21226
plt.plot(sp[key].rest_wave, (sp[key].rest_fnu_s99model / sp[key].rest_fnu_autocont) - nudge, label='_nolegend_', color='b')
plt.annotate(key +", age=26 Myr, Z=0.34",  xy=(1530, 0 - nudge), xycoords='data', fontsize=14)
for ii, key in enumerate(plotorder):
    filename = key + '-sb99-fit.txt'
    df[key] = pandas.read_table(figdir + filename, delim_whitespace=True, comment="#", names=colnames)
    df[key]['obsflux'].loc[df[key]['obsflux'].lt(-10)] = np.nan
    df[key]['tempy'] = df[key]['obsflux']  /  df[key]['obsflux'][df[key].wave.between(1530,1536)].median()  # renormalize for plot    
    df[key]['tempS99'] = df[key]['S99fit'] /  df[key]['obsflux'][df[key].wave.between(1530,1536)].median()  # renormalize for plot    
    plt.plot(df[key]['wave'], df[key]['tempy'] + ii*scalef, label=key, linestyle='steps-mid')
    plt.plot(df[key]['wave'], df[key]['tempS99'] + ii*scalef, color='b', label='_nolegend_')
    plotlab = key + ", age=" + str(age_dict[key]) + " Myr, Z=" + str(Z_dict[key]) 
    plt.annotate(plotlab, xy=(1530,ii*scalef+0.5), xycoords='data', fontsize=14)
key = from_meg[0]  # Plot RCS0327-E
nudge = ii +  2.5
plt.plot(sp[key].rest_wave, (sp[key].rest_fnu / sp[key].rest_fnu_autocont) + nudge, label=key, linestyle='steps-mid', color='k')
plt.plot(sp[key].rest_wave, (sp[key].rest_fnu_s99model / sp[key].rest_fnu_autocont) + nudge, label='_nolegend_', color='b')
plt.annotate(key +", age=2.5 Myr, Z=0.24",  xy=(1530, nudge+0.3), xycoords='data', fontsize=14)
plt.ylim(-1.2,8.4)
plt.xlim(1530,1570)
#plt.xlim(1500,1600)
#plt.xlim(1350,1450)
plt.xlabel(r'rest wavelength ($\AA$)',  fontsize=18)
plt.ylabel('scaled flux', fontsize=18)

# Bigger x axis, suppress y axis
matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16)
plt.yticks([])

#handles, labels = ax.get_legend_handles_labels()
#ax.legend(reversed(handles), reversed(labels), title='Line', loc='upper left')
