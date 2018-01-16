import pandas
import jrr
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
from astropy import units
from math import pi
mage_mode = "reduction"


# Load the S99 summary file
S99file = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/MagE_atlas/Contrib/S99/MW_Corrected/sb99_overview2.txt"
names = ("name", "EBV", "sigmaEBV", "ZS99", "tS99")
S99 = pandas.read_table(S99file, delim_whitespace=True, comment="#", names=names)
S99.set_index('name', inplace=True)

# Load the Megasaura spectra
(sp, resoln, dresoln, LL, zz_sys, speclist) = jrr.mage.open_many_spectra(mage_mode, which_list="wcont", zchoice='neb', addS99=True, MWdr=True)

wave_targ = 1500.0  # where to measure rest-frame flux density, Angstroms
wave_win = 10.      # window +-1500, to measure rest-frame flux density.  units are rest-frame Angstroms

print "label   EBV_S99  Z_S99   ZS99  fnu1500  fnu1500_dered  SFR_raw  SFR_dered (SFRs and fnus not yet corrected for magnification; missing nebular Z)"
for label in speclist['short_label'] :
    fnu1500 = sp[label]['rest_fnu'].loc[sp[label]['rest_wave'].between(wave_targ - wave_win, wave_targ + wave_win)].median()
    if label in S99.index :  # if there's a EBV from S99 fit
         jrr.mage.deredden_internal_extinction(sp[label], S99.ix[label]['EBV'], 'rest_fnu_cont')
         fnu1500_dered = sp[label]['rest_fnu_dered'].loc[sp[label]['rest_wave'].between(wave_targ - wave_win, wave_targ + wave_win)].median()
         #*** Check whether Calzetti is implmented properly!  Stars or gas!
         zz = speclist.ix[label]['z_syst']
         flux_to_lum = 4 * pi * (cosmo.luminosity_distance(zz).cgs.value)**2     # Apply luminosity distance
         K98_convert = 1.4E-28  # Eqn 1 of Kennicutt et al. 1998, units are erg/s/Hz
         SFR_raw   =  fnu1500 * flux_to_lum * K98_convert        # ** Multiplication should be divided here
         SFR_dered =  fnu1500_dered * flux_to_lum * K98_convert  # ** Multiplication should be divided here
         print label, S99.ix[label]['EBV'],  S99.ix[label]['ZS99'],  S99.ix[label]['tS99'], fnu1500, fnu1500_dered, SFR_raw, SFR_dered
         # Still missing nebular metallicity -- grab it from Table 2 of Sample paper
# Double-check whether Extinction-correction is right.... **
