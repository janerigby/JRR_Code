from __future__ import print_function
from builtins import str
from builtins import zip
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
labels = ('PSZ0441', 'SPT0310', 'SPT2325', 'SPT0356', 'SPT0142', 'planckarc_pos1')  #
Z      = (0.10,       0.23,       0.27,     0.22,       0.23,     0.55 )  
age    = (19.3,       15.7,        9.2,      8.6,        10.1,     2.9)
age_dict = dict(list(zip(labels, age)))
Z_dict =   dict(list(zip(labels, Z)))


from_meg = ('rcs0327-E', 'S1226+2152')
(sp, resoln, dresoln, LL, zz_sys, speclist)  = jrr.mage.open_many_spectra(mage_mode, which_list='labels', labels=from_meg, addS99=True)
for ii, key in enumerate(list(sp.keys())) :
    print("redshift, object", zz_sys[key], key)


#jrr.mage.plot_linelist(LL[key], z_systemic=zz_sys[key],  restframe=True)

home= expanduser("~")
#figdir = home + "/Dropbox/SGAS-shared/PROPOSALS/Hubble/Cycle27/Friends_o_Meg/Figure/"
#figdir = home + "/Texts/Proposals/HST/C28/Friends_of_megasaura"
figdir = home + "/Texts/Proposals/HST/C29/Friends_of_Megasaura/Figure/"
plotorder = ('PSZ0441',  'SPT0310', 'SPT0142', 'SPT2325', 'planckarc_pos1')   # from bottom to top
need_HST = ('PSZ0441', 'SPT0142')  

filenames =  [ basename(x) for x in glob.glob(figdir + "*-fit.txt") ]
colnames=('wave', 'obsflux', 'errflux', 'S99fit')
df = {}
#for ii, filename in enumerate(filenames):
#    myre = re.compile('(\S+)-sb99-fit.txt')
#    key = myre.match(filename).groups()[0]

what_to_plot =   "CIV" # "CII" #  

if what_to_plot == "CII" :
    xy1 = 1328 ; xy2 = 1332
    scalef = 2
    xlim = (1315,1355)
    labelshift = 11
elif what_to_plot == "CIV":
    xy1 = 1530 ; xy2 = 1536  
    scalef = 2
    xlim = (1530,1570)
    labelshift = 0
elif what_to_plot == "SiII" :
    xy1 = 1250 ; xy2 = 1255
    scalef =2
    xlim = (1240,1280)
    labelshift = 5


fig = plt.figure(figsize=(11,12))
ax = fig.add_subplot(1, 1, 1)
nudge  = 1.1
key = from_meg[-1]
plt.step(sp[key].rest_wave, (sp[key].rest_fnu / sp[key].rest_fnu_autocont) - nudge, label=key, where='mid')  # Plot S21226
plt.plot(sp[key].rest_wave, (sp[key].rest_fnu_s99model / sp[key].rest_fnu_autocont) - nudge, label='_nolegend_', color='b')
plt.annotate(key +", age=26 Myr, Z=0.34",  xy=(xy1-labelshift, 0 - nudge), xycoords='data', fontsize=14)
for ii, key in enumerate(plotorder):
    if key in need_HST : lw =2 
    else : lw =1
    if key == "planckarc_pos1" : altlab = "Sunburst Arc"
    else :                     altlab = key
    plotlab = altlab + ", age=" + str(age_dict[key]) + " Myr, Z=" + str(Z_dict[key])
    filename = key + '-sb99-fit.txt'
    df[key] = pandas.read_table(figdir + filename, delim_whitespace=True, comment="#", names=colnames)
    df[key]['obsflux'].loc[df[key]['obsflux'].lt(-10)] = np.nan
    df[key]['tempy'] = df[key]['obsflux']  /  df[key]['obsflux'][df[key].wave.between(xy1,xy2)].median()  # renormalize for plot    
    df[key]['tempS99'] = df[key]['S99fit'] /  df[key]['obsflux'][df[key].wave.between(xy1,xy2)].median()  # renormalize for plot
    if key in need_HST:
        plt.step(df[key]['wave'], df[key]['tempy'] + ii*scalef, label=key, where='mid', lw=1.8, color='k')
        #plt.annotate(plotlab, xy=(xy1-labelshift,ii*scalef+0.), xycoords='data', fontsize=18, fontweight='bold')
        plt.annotate(plotlab, xy=(xy1 + 21,ii*scalef+0.), xycoords='data', fontsize=18, fontweight='bold')
    else :
        plt.step(df[key]['wave'], df[key]['tempy'] + ii*scalef, label=key, where='mid', lw=1.5)
        plt.annotate(plotlab, xy=(xy1-labelshift,ii*scalef+0.5), xycoords='data', fontsize=14)
    plt.plot(df[key]['wave'], df[key]['tempS99'] + ii*scalef, color='b', label='_nolegend_')

key = from_meg[0]  # Plot RCS0327-E
nudge = ii +  6
plt.step(sp[key].rest_wave, (sp[key].rest_fnu / sp[key].rest_fnu_autocont) + nudge, label=key, where='mid')
plt.plot(sp[key].rest_wave, (sp[key].rest_fnu_s99model / sp[key].rest_fnu_autocont) + nudge, label='_nolegend_', color='b')
plt.annotate(key +", age=2.5 Myr, Z=0.24",  xy=(xy1-labelshift, nudge+0.3), xycoords='data', fontsize=14)
plt.ylim(-1.2,12)
plt.xlim(xlim)
#plt.xlim(1350,1450)
plt.xlabel(r'rest wavelength ($\AA$)',  fontsize=20)
plt.ylabel('scaled flux', fontsize=20)
plt.plot((1260.4221,1260.4221), (-20,20), color='grey', lw=2) # Draw a vert line at SiII
plt.plot((1335.7077, 1335.7077), (-20,20), color='grey', lw=2) # Draw a vert line at C II
plt.plot((1548.195, 1548.195), (-20,20), color='grey', lw=2) # Draw a vert line at CIV
plt.plot((1550.770, 1550.770), (-20,20), color='grey', lw=2) # Draw a vert line at CIV
plt.suptitle(what_to_plot, fontsize=24, fontweight='bold')
#plt.annotate(what_to_plot, xy=(0.6,0.01), xycoords='axes fraction', fontsize=24, fontweight='bold',horizontalalignment='right')

# Bigger x axis, suppress y axis
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)
plt.yticks([])
plt.tight_layout()
fig.savefig(what_to_plot + "_C29extended.pdf")


#handles, labels = ax.get_legend_handles_labels()
#ax.legend(reversed(handles), reversed(labels), title='Line', loc='upper left')
