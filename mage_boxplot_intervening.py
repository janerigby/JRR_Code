import jrr
#import ayan
import numpy as np
import pandas
from os.path import expanduser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
mage_mode = "released"
# Make a boxplot for each isolated Intervening system in the MagE Megasaura spectra

def make_boxplots_of_doublets(doublet_file) :
    doublet_list = pandas.read_table(doublet_file, delim_whitespace=True, comment="#")
    #foo = doublet_list[0:5]
    for row in doublet_list.itertuples():
    #for row in foo.itertuples():  # TEMP, only run a few doublets, for debugging
        if row.flag == "y" :
            (sp, resoln, dresoln, LL, zz_syst) = jrr.mage.wrap_open_spectrum(row.gal, mage_mode, addS99=True)  # load spectrum
            lims = jrr.spec.get_waverange_spectrum(sp)
            # Only plot the transitions that are covered by the spectrum, at intervening absorber redshift row.zz
            transin = trans_list[trans_list.wave.between(lims[0]/(1.+row.zz), lims[1]/(1.+row.zz))] 
            transin.reset_index(inplace=True, drop=True)
            jrr.plot.boxplot_Nspectra( (sp.wave,), (sp.fnu/sp.fnu_autocont,), (sp.fnu_u/sp.fnu_autocont,), (row.zz,), transin.label, transin.wave, win=2000, Ncol=1, LL=LL, vel_plot=True)            
            extra_label = row.gal + " " + row.pointing + " Intervening absorber at z="+ str(np.round(row.zz, decimals=4)) + " (" + row.flag + ")"
            print "Working on " + extra_label
            plt.suptitle(extra_label, fontsize=20)
            pp.savefig()
    return(0)


batches = ('batch1', 'batch23')
doubletfiles = {}
doubletfiles['batch1']  = expanduser("~") + "/Dropbox/MagE_atlas/Contrib/Intervening/Doublet_search/Results_16Feb2018/found_doublets_SNR4_JRRedit.txt"
doubletfiles['batch23'] = expanduser("~") + "/Dropbox/MagE_atlas/Contrib/Intervening/Doublet_search/Results_21Feb2019/found_doublets_batch3_SNR4_JRRedit.txt"


trans_file = expanduser("~") + "/Dropbox/MagE_atlas/Linelists/MINE/interven.lst"
trans_list = pandas.read_table(trans_file, delim_whitespace=True, comment="#", names=("wave", "lab1", "lab2", "f", "d", "type"))
trans_list['label'] = trans_list.lab1 + "_" + trans_list.lab2.astype('str')

for batch in batches :
    the_pdf = "intervening_doublets_boxplots_" + batch + ".pdf"
    pp = PdfPages(the_pdf)
    foo =  make_boxplots_of_doublets(doubletfiles[batch])
    pp.close()
    

'''         plt.clf()
            
            print "debugging, working on this row of doublet list: ", row
            plt.step(sp.wave, sp.fnu/sp.fnu_autocont, color='black')
            plt.step(sp.wave, sp.fnu/sp.fnu_s99model, color='blue')
            plt.step(sp.wave, sp.fnu_u/sp.fnu_autocont, color='grey')
            plt.title(row.doubname + " at z="+ str(np.round(row.zz, decimals=5)))
            plt.xlim(row.wave1 - 30., row.wave2+30.)
            plt.ylim(-1,2)
            plt.vlines( (row.wave1, row.wave2), -3, 3, colors='red')
            plt.show()
            plt.pause(1)
            plt.clf()'''
