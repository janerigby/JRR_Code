import jrr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas

def pretty_the_plot() :
    plt.xlabel("rest-frame wavelength (A)", fontsize=18)
    plt.ylabel("normalized flux", fontsize=18)
    plt.xlim(1200, 2000)
    plt.ylim(0,1.6)
    plt.rcParams.update({'font.size': 16})
    return(0)

def plot_the_spectrum(ax, df,  windlabels_df, photlabels_df):
    ax.plot(df['wave'], df['flam'], color='k')
    ax.plot(df['wave'], df['S99fit'], color='r')
    ax.plot(df['wave'], df['flam_u'], color='grey')
    jrr.plot.annotate_from_dataframe(windlabels_df, xcol='wave', ycol='dum', text='label1', ha='center', fontsize=16, ax=ax)
    jrr.plot.annotate_from_dataframe(photlabels_df, xcol='wave', ycol='dum', text='label1', ha='center', fontsize=16, ax=ax)
    ax.stem(windlabels_df['wave'], windlabels_df['dum'], linefmt='k--', markerfmt=' ', bottom=1.0, basefmt=' ', use_line_collection=True)
    ax.stem(photlabels_df['wave'], photlabels_df['dum'], linefmt='k--', markerfmt=' ', bottom=0.6, basefmt=' ', use_line_collection=True)
    pretty_the_plot()
    return(0)

plt.ion()
windlabels_df = pandas.read_csv('wind_labels.txt', comment="#", delim_whitespace=True)
photlabels_df = pandas.read_csv('photospheric_labels.txt', comment="#", delim_whitespace=True)

pp1 = PdfPages('S1723.pdf')
fig1, ax1 = plt.subplots(figsize=(10,3.5))
cols = ('wave', 'flam', 'flam_u', 'S99fit')
df_s1723 = pandas.read_csv('S1723-sb99-fit.txt', comment="#", delim_whitespace=True, names=cols)
df_s1723_good =  df_s1723.loc[df_s1723['flam_u'].lt(3)]     
plot_the_spectrum(ax1, df_s1723_good,   windlabels_df, photlabels_df)
ax1.annotate("S1723", xy=(0.85,0.05), xycoords='axes fraction',  fontsize=20)
plt.tight_layout()
pp1.savefig()

pp2 = PdfPages('S1226.pdf')
fig2, ax2 = plt.subplots(figsize=(10,3.5))
s1226file = '/Users/jrrigby1/Dropbox/MagE_atlas/Contrib/S99/Chisholm19/ApJ/Megasaura/S1226+2152-sb99-fit.txt'
df_s1226 = pandas.read_csv(s1226file, comment="#", delim_whitespace=True, names=cols)
df_s1226_good =  df_s1226.loc[df_s1226['flam_u'].lt(100)]     
plot_the_spectrum(ax2, df_s1226_good,   windlabels_df, photlabels_df)
plt.annotate("S1226",  xy=(0.85,0.05), xycoords='axes fraction', fontsize=20)
plt.tight_layout()
pp2.savefig()

# Plot on same axis
pp3 = PdfPages('S1723_S1226.pdf')
fig, (ax3, ax4) = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(12, 6))
plot_the_spectrum(ax3, df_s1723_good,   windlabels_df, photlabels_df)
ax3.annotate("MMT/BC, Keck/ESI for S1723", xy=(0.6, 0.02), xycoords='axes fraction',  fontsize=20)
plot_the_spectrum(ax4, df_s1226_good,   windlabels_df, photlabels_df)
ax4.annotate("Magellan/MagE for S1226", xy=(0.6, 0.05), xycoords='axes fraction', fontsize=20)
plt.tight_layout()
plt.subplots_adjust(hspace=0)
pp3.savefig()
pp1.close() ; pp2.close(); pp3.close()
