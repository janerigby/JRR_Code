''' Same as multipanel-spectrum-ids.py, but plot multiple spectra.  Using for sub-stacks (highZ, loZ),
young, middle-aged, old.
jrigby, Sept 2016. '''

import jrr
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from commands import getoutput
import sys
from re import split, sub, search
mage_mode = "reduction"  
#mage_mode = "released"

def annotate_the_first_page() :
    plt.annotate("f_nu (erg/s/cm^2/Hz)",  xy=(0.05,0.055), color="black", xycoords="figure fraction")
    plt.annotate("Line color coding:", xy=(0.05,0.043), color="black", xycoords="figure fraction")
    plt.annotate("Photosphere", xy=(0.05,0.03), color="blue", xycoords="figure fraction")            
    plt.annotate("Emission",    xy=(0.17,0.03), color="red", xycoords="figure fraction")
    plt.annotate("ISM",         xy=(0.25,0.03), color="green", xycoords="figure fraction")
    plt.annotate("Wind",        xy=(0.32,0.03), color="cyan", xycoords="figure fraction")
    plt.annotate("Fine structure", xy=(0.4,0.03), color="purple", xycoords="figure fraction")
    plt.annotate("Intervening", xy=(0.55,0.03), color="orange", xycoords="figure fraction")


# Compare high and low metallicity stacks
stacks = ['magestack_bystars_highZ_spectrum.txt', 'magestack_bystars_lowZ_spectrum.txt']
sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
the_dfs = [sp1, sp2]
the_zzs = [0.0, 0.0]  # Stacked spectra are already in rest frame.
the_pdf =  "multipanel_stacked_bystars_byZ.pdf"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title="Compare lowZ and highZ stacks", annotate=annotate_the_first_page)
plt.clf()

stacks = ['magestack_byneb_highZ_spectrum.txt', 'magestack_byneb_lowZ_spectrum.txt']
sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
the_dfs = [sp1, sp2]
the_zzs = [0.0, 0.0]  # Stacked spectra are already in rest frame.
the_pdf =  "multipanel_stacked_byneb_byZ.pdf"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title="Compare lowZ and highZ stacks", annotate=annotate_the_first_page)
plt.clf()

# Compare young, middle-age, and old stacks
stacks = ['magestack_bystars_younglt8Myr_spectrum.txt', 'magestack_bystars_midage8to16Myr_spectrum.txt', 'magestack_bystars_oldgt16Myr_spectrum.txt']
sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
sp3 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[2])
the_dfs = [sp1, sp2, sp3]
the_zzs = [0.0, 0.0, 0.0]  # Stacked spectra are already in rest frame.
the_pdf =  "multipanel_stacked_bystars_byage.pdf"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title="Compare young (blue), middle-age (black), and old (red) stacks", annotate=annotate_the_first_page, colortab=('blue', 'black', 'red'), waverange=(1000,3000))
plt.clf()

# Compare young, middle-age, and old stacks
stacks = ['magestack_byneb_younglt8Myr_spectrum.txt', 'magestack_byneb_midage8to16Myr_spectrum.txt', 'magestack_byneb_oldgt16Myr_spectrum.txt']
sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
sp2 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
sp3 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[2])
the_dfs = [sp1, sp2, sp3]
the_zzs = [0.0, 0.0, 0.0]  # Stacked spectra are already in rest frame.
the_pdf =  "multipanel_stacked_byneb_byage.pdf"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=True, title="Compare young (blue), middle-age (black), and old (red) stacks", annotate=annotate_the_first_page, colortab=('blue', 'black', 'red'), waverange=(1000,3000))
plt.clf()
