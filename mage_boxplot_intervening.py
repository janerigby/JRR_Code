from __future__ import print_function
from builtins import str
import jrr
#import ayan
from re import sub
import numpy as np
import pandas
import string
from os.path import expanduser
from astropy.io.ascii import read
mage_mode = "released"
# Make a boxplot for each isolated Intervening system in the MagE Megasaura spectra
import warnings
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages


warnings.filterwarnings('ignore')  # maybe this is sloppy, but I can't see the output for all the pandas warnings


def make_boxplots_of_doublets(doublet_df, trans_list, pp) :  # Plot a lot of transitions at once, one boxplot per page
    for row in doublet_df.itertuples():
        if "y" in row.flag or "p" in row.flag or '-' in row.flag:
            (sp, resoln, dresoln, LL, zz_syst) = jrr.mage.wrap_open_spectrum(row.shortname, mage_mode, addS99=False)  # load spectrum
            lims = jrr.spec.get_waverange_spectrum(sp)
            transin = trans_list[trans_list.wave.between(lims[0]/(1.+row.zz), lims[1]/(1.+row.zz))]  # Only plot transitions covered by the spectrum
            transin.reset_index(inplace=True, drop=True)
            jrr.plot.boxplot_Nspectra( (sp.wave,), (sp.fnu/sp.fnu_autocont,), (sp.fnu_u/sp.fnu_autocont,), (row.zz,), transin.label, transin.wave, win=2000, Ncol=1, LL=LL, vel_plot=True) 
            extra_label = row.shortname + ",  intervening absorber z="+ str(np.round(row.zz, decimals=3)) + " (" + row.flag + ")"
            print("Working on " + extra_label)
            plt.suptitle(extra_label, fontsize=20)
            pp.savefig()
    return(doublet_df)


def nested_subplots(doublet_df, trans_list2, pp, plot_linelist=True):  # Same as plot_only_detected_doublets, but in grid format on PDF
    nrows = 5 
    ncols = 2 
    Npages = 7
    counter = 0
    # There are 66 plotted doublets, so ideally want to plot in grids, 2 columns x 5 rows, so that's 7 pages
    for pg in range(0, Npages) :
        fig = plt.figure(figsize=(8, 10))
        outer = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1)
        for i in range(nrows*ncols):
            inner = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[i], wspace=0.0, hspace=0.0)
            for j in range(1):
                #ax = plt.Subplot(fig, inner[j])
                #t = ax.text(0.5,0.5, 'pg=%d out=%d in=%d c=%d' % (pg, i, j, counter))
                #t.set_ha('center')
                #ax.set_xticks([])
                #ax.set_yticks([])
                #plt.plot((0,1), (1,3))  # debugging
                plot_one_doublet(counter, doublet_df, trans_list2, plot_linelist=False, makenewfig=True)
                #fig.add_subplot(ax)
                counter += 1
        plt.tight_layout()
        pp.savefig(fig)
    return(0)

def plot_one_doublet(ii, doublet_df, trans_list2, plot_linelist=True, makenewfig=False) :  # make a plot for just one row of the doublet dataframe
    pretty_trans = { 'MgII' : 'Mg II', 'CIV' : 'C IV', 'SiIV' : 'Si IV'}
    row = doublet_df.iloc[ii]
    if "y" in row.flag or "p" in row.flag or '-' in row.flag:
        (sp, resoln, dresoln, LL, zz_syst) = jrr.mage.wrap_open_spectrum(row.shortname, mage_mode, addS99=False)  # load spectrum
        lims = jrr.spec.get_waverange_spectrum(sp)
        transin = trans_list2[trans_list2.wave.between(lims[0]/(1.+row.zz), lims[1]/(1.+row.zz))]  # Only plot transitions covered by the spectrum
        transin.reset_index(inplace=True, drop=True)
        ## Just plot the detections. Added second_doublet argument to jrr.plot.boxplot_Nspectra(), to mark second line of doublet
        doubname = sub(' ', '', row.doubname)
        if    doubname == 'MgII' :
            trans_det = transin.loc[transin['label'] == 'Mg_II_2796']
            doublet2nd = 2803.5310
        elif  doubname == 'CIV':
            trans_det = transin.loc[transin['label'] == 'C_IV_1548']
            doublet2nd = 1550.7700
        elif  doubname == 'SiIV':
            trans_det = transin.loc[transin['label'] == 'Si_IV_1393']
            doublet2nd =  1402.770
        else : raise Exception("ERROR, unrecognized doublet", doubname)
        if plot_linelist : inLL = LL
        else             : inLL=()
        jrr.plot.boxplot_Nspectra( (sp.wave,), (sp.fnu/sp.fnu_autocont,), (sp.fnu_u/sp.fnu_autocont,), (row.zz,), trans_det.label.values, trans_det.wave.values, win=15, Ncol=1, LL=inLL, vel_plot=False, figsize=(9,4), second_doublet = doublet2nd, ymax=(2,), plot_linelabel=False, makenewfig=makenewfig)
        # Clean up the title
        extra_label = jrr.mage.prettylabel_from_shortlabel(row.shortname)  + ',  ' + pretty_trans[doubname]
        extra_label = sub("image3", " image 3", extra_label)
        round_zz = np.format_float_positional(row.zz, precision=3, trim='k', unique=False)
        extra_label += r" at $z_{abs} = $" + round_zz
        if row.flag != '-' : extra_label += " (" + row.flag + ")"
        print("Working on " + extra_label)
        plt.title(extra_label, fontname='Times')
        plt.tight_layout()
        return(0)


# Want to plot same order as in the tables. So, grab the as-published table
batch = 'aspublished'
doublet_df = jrr.mage.load_doublet_df(batch)
make_big_boxplots = True
make_pub_boxplots_1perpage = False
try_nested = False

trans_file  = expanduser("~") + "/Dropbox/MagE_atlas/Linelists/MINE/interven_short.lst"
trans_list = pandas.read_table(trans_file, delim_whitespace=True, comment="#", names=("wave", "lab1", "lab2", "f", "d", "type"))
trans_list['label'] = trans_list.lab1 + "_" + trans_list.lab2.astype('str')

trans_file2 = expanduser("~") + "/Dropbox/MagE_atlas/Linelists/MINE/interven_shorter.lst"
trans_list2 = pandas.read_table(trans_file2, delim_whitespace=True, comment="#", names=("wave", "lab1", "lab2", "f", "d", "type"))
trans_list2['label'] = trans_list2.lab1 + "_" + trans_list2.lab2.astype('str')

if make_big_boxplots :  # Make boxplots with lots of emission lines. For data exploration, not for publication.
    matplotlib.rcParams.update({'font.size': 14})
    the_pdf  = "intervening_doublets_boxplots_" + batch + ".pdf"
    pp  = PdfPages(the_pdf)
    result  =   make_boxplots_of_doublets(doublet_df, trans_list,  pp)
    pp.close()

if make_pub_boxplots_1perpage : # Plot only the detected doublets.  For publication.  1 plot per page
    matplotlib.rcParams.update({'font.size': 20})
    the_pdf2 = "intervening_doublets_justdetected_" + batch + ".pdf"
    pp2 = PdfPages(the_pdf2)
    for ii, row in enumerate(doublet_df.itertuples()):
        plot_one_doublet(ii, doublet_df, trans_list2, plot_linelist=False)
        pp2.savefig()
    pp2.close()
    
if try_nested:  # Try making plots in multiple pages of 2x5 grid
    matplotlib.rcParams.update({'font.size': 20})
    the_pdf3 = "intervening_doublets_justdetected_grid_" + batch + ".pdf"
    pp3 = PdfPages(the_pdf3)
    nested_subplots(doublet_df, trans_list2, pp3, plot_linelist=False)
    pp3.close()
plt.clf()


# TO DO:
# write a wrapper for plot_only_detected_doublets, so it makes a pretty 2xN plot on several pages.  Publish and move on!
# Remove the line IDs?  Most are distracting and not useful...  I think I can do this by just feeding LL=False  Try this
