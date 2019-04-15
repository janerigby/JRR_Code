from __future__ import print_function
# CALL THIS from command line, as:
#   ipython ~/Python/JRR_Code/mage_calc_galprops.py > mage_galproperties.txt   #clunky I/O, I know

from os.path import expanduser
import pandas
import numpy as np
import jrr
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
from astropy import units
from math import pi
mage_mode = "released"
homedir = expanduser("~")

# Load the S99 summary file
S99file = homedir + "/Dropbox/MagE_atlas/Contrib/S99/MW_Corrected/sb99_overview2.txt"
names = ("name", "EBV", "sigmaEBV", "ZS99", "tS99")
S99 = pandas.read_table(S99file, delim_whitespace=True, comment="#", names=names)
S99.set_index('name', inplace=True)

# Load the nebular metallicities from Table 2 of MagE sample paper
nebfile = homedir + "/Dropbox/MagE_atlas/Rest-Opt/megasaura_metallicities_samplepaper.txt"
Zneb =  pandas.read_table(nebfile, delim_whitespace=True, comment="#")
Zneb.set_index('label', inplace=True)
Zneb['Zneb_solar'] = 10**(Zneb['12+logOH'] - 8.69)
Zneb['Zneb_sig']   = 10**(Zneb['12+logOH'] + Zneb['sigZ'] - 8.69) - 10**(Zneb['12+logOH']  - 8.69)
Zneb.loc[Zneb['sigZ'] == -99, 'Zneb_sig'] = np.int(-99)
Zneb.loc[Zneb['sigZ'] == -9, 'Zneb_sig'] = np.int(-9)

# Load the Megasaura spectra
(sp, resoln, dresoln, LL, zz_sys, speclist) = jrr.mage.open_many_spectra(mage_mode, which_list="wcont", zchoice='neb', addS99=True, MWdr=True, silent=True, verbose=False)

wave_targ = 1500.0  # where to measure rest-frame flux density, Angstroms
wave_win = 10.      # window +-1500, to measure rest-frame flux density.  units are rest-frame Angstroms

print("label   EBV_S99  sigma_EBV  Zneb  sigZ  Z_S99   age_S99  fnu1500  fnu1500_dered  SFR_raw  SFR_dered")
print("#(SFRs and fnus not yet corrected for magnification)")
for label in speclist['short_label'] :
    fnu1500 = sp[label]['rest_fnu'].loc[sp[label]['rest_wave'].between(wave_targ - wave_win, wave_targ + wave_win)].median()
    if label in S99.index :  # if there's a EBV from S99 fit
        jrr.spec.deredden_internal_extinction(sp[label], S99.loc[label]['EBV'])
        fnu1500_dered = sp[label]['rest_fnu_dered'].loc[sp[label]['rest_wave'].between(wave_targ - wave_win, wave_targ + wave_win)].median()
        zz = speclist.loc[label]['z_syst']
        flux_to_lum = 4 * pi * (cosmo.luminosity_distance(zz).cgs.value)**2     # Apply luminosity distance
        K98_convert = 1.4E-28  # Eqn 1 of Kennicutt et al. 1998, units are erg/s/Hz
        SFR_raw   =  fnu1500 * flux_to_lum * K98_convert        # ** Magnification should be divided here
        SFR_dered =  fnu1500_dered * flux_to_lum * K98_convert  # ** Magnification should be divided here
        print(label, S99.loc[label]['EBV'],  S99.loc[label]['sigmaEBV'], Zneb.loc[label]['Zneb_solar'],  Zneb.loc[label]['Zneb_sig'], end=' ')
        print(S99.loc[label]['ZS99'],  S99.loc[label]['tS99'], fnu1500, fnu1500_dered, SFR_raw, SFR_dered)


# Make a new data frame w results, and dump it to a file.
galprops =  pandas.read_table("mage_galproperties.txt", delim_whitespace=True, comment="#")
