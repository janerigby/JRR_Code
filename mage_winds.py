from numpy import array, concatenate
import jrr
import os
import sys
import re
from numpy import zeros_like
from matplotlib import pyplot as plt
import pandas
import numpy as np
#mage_mode = "reduction"
mage_mode = "released"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)


def plot_winds_neutral_stellar(prefix, thewaves, thefnus, thedfnus, thezs, vwin, Ncol, label="", LL=[], z_sys=0.0, colortab=False) :
    ''' thewaves, thefnus, thedfnus, thezs are TUPLES of arrays of wavelength, fnu, sigma, and redshift.  If only plotting one, use thewaves=(wave_array,) '''
    jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_a, line_center_a, vwin, Ncol, LL, extra_label=label, colortab=colortab)
    plt.savefig(prefix + "a.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close()
    jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_b, line_center_b, vwin, Ncol, LL, extra_label=label, colortab=colortab)
    plt.savefig(prefix + "b.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close()
    jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_c, line_center_c, vwin, Ncol, LL, extra_label=label, colortab=colortab)
    plt.savefig(prefix + "c.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close()

        
# Define the lines to plot
line_label_a  = ('Lya', 'C II 1334',  'Si II 1260', 'Si II 1526', 'Al II 1670')
  #, 'O I 1302')  # Nino and Marc R. say these trace neutral gas.  
  #Added SiII1260, AlII, AlIII at Rongmon's suggestion
line_center_a = array((1215.6701, 1334.5323,  1260.4221, 1526.7066, 1670.7874))
  ##, 1302.1685))  Temporarily commented out OI for C23 proposal

line_label_b  =  ('MgII2796', 'FeII2344', 'FeII2383') #, 'FeII2600' )
line_center_b = array((2796.352,  2344.214,  2382.765)) #,  2600.1729,))
line_label_c  = ("N V 1238", "Si IV 1393", "C IV 1548", "Al III 1854")
line_center_c = array((1238.82, 1393.76, 1548.19, 1854.72))


line_label_extra  = ("OVI 1031", "OVI 1037", "CIII] 1907", "C II] 2326", "OIII] 1660", "O II] 2470", "SiII] 2335", "HeII1640", "SiIII 1882", "SiIII 1892", "LyB", "SiIV")
line_center_extra =array((1031.9261, 1037.6167, 1906.68, 2326.00, 1660.81, 2470.97, 2335.123, 1640.42,  1883.00, 1892.03, 1025.7223, 1402.770))

line_label_all  = line_label_a  + line_label_b  + line_label_c + line_label_extra
line_center_all = concatenate((line_center_a, line_center_b, line_center_c, line_center_extra))
    
# housekeeping
redshift = 0.0  # stacked spectrum is already in rest frame wavelength
vwin = 4000. # +- velocity window (km/s) to consider a line
Ncol = 1

print "STATUS:  Plotting wind lines for MagE stack"
linelist = line_path + "stacked.linelist"
(LL, z_sys) = jrr.mage.get_linelist(linelist)
(sp) = jrr.mage.open_stacked_spectrum(mage_mode)
plot_winds_neutral_stellar("MagEstack/", (sp.wave,), (sp.X_avg,), (sp.X_sigma,), (redshift,), vwin, Ncol, "MagE stack", LL, z_sys)
plt.clf()
print "STATUS:  Repating, for median MagE rather than X_avg"
plot_winds_neutral_stellar("MagEmedian/", (sp.wave,), (sp.X_median,), (sp.X_jack_std,), (redshift,), vwin, Ncol, "MagE median", LL, z_sys)
plt.clf

print "STATUS:  Plotting wind lines for S99 fit to MagE Stack A"
(S99, zz) = jrr.mage.open_S99_spectrum("stack-A", 0.0)  
plot_winds_neutral_stellar("S99fit/", (S99.rest_wave,), (S99.rest_fnu_s99,), (S99.rest_fnu_s99_u,), (0.0,0.0), vwin, Ncol, "S99 fit", LL, z_sys)
plt.clf

# Renormalizing the spectra with autocont, for better plotting...
print "STATUS: Overplotting MagE stack and S99 fit to Mage Stack, for wind lines"
alt_file = "magestack_byneb_ChisholmstackA_spectrum.txt"  # this is the spectrum that JChisholm fit
altsp = jrr.mage.open_stacked_spectrum(mage_mode, alt_file)
plot_winds_neutral_stellar("MageES99/", (altsp.wave, S99.rest_wave), (altsp.X_avg/altsp.fnu_autocont, S99.rest_fnu_s99/S99.rest_fnu_s99_autocont), (altsp.X_sigma/altsp.fnu_autocont, S99.rest_fnu_s99*-0.01), (0.0, 0.0), vwin, Ncol, "Stack and S99 fit", LL, z_sys)

# Same, but for median spectrum.  Not sure how well the autocont will line up
altsp = jrr.mage.open_stacked_spectrum(mage_mode, alt_file,  colfnu='X_median', colfnuu='X_jack_std')
# Above should put the median in _fnu, and do autocont on it.  Let's see...
plot_winds_neutral_stellar("MageES99_median/", (altsp.wave, S99.rest_wave), (altsp.X_median/altsp.fnu_autocont, S99.rest_fnu_s99/S99.rest_fnu_s99_autocont), (altsp.X_jack_std/altsp.fnu_autocont, S99.rest_fnu_s99*-0.01), (0.0, 0.0), vwin, Ncol, "Stack and S99 fit", LL, z_sys)
                           

### Obsolete, a bad alley to walk down.
## Same, but by age
#st1 = jrr.mage.open_stacked_spectrum(mage_mode, "magestack_bystars_younglt8Myr_spectrum.txt")
#st2 = jrr.mage.open_stacked_spectrum(mage_mode, "magestack_bystars_midage8to16Myr_spectrum.txt")
#st3 = jrr.mage.open_stacked_spectrum(mage_mode, "magestack_bystars_oldgt16Myr_spectrum.txt")
#plot_winds_neutral_stellar("MagEstack_byage/", (st1.wave, st2.wave, st3.wave), (st1.X_avg/st1.fnu_autocont,st2.X_avg/st2.fnu_autocont,st3.X_avg/st3.fnu_autocont), (st1.X_sigma/st1.fnu_autocont, st2.X_sigma/st2.fnu_autocont,st3.X_sigma/st3.fnu_autocont), (0.0, 0.0, 0.0), vwin, Ncol, "stacked by age: young (blue), middle (black), old (red)", LL, 0.0, colortab=('blue', 'black', 'red'))
#
## same, but by metallicity
#st1 = jrr.mage.open_stacked_spectrum(mage_mode, "magestack_bystars_lowZ_spectrum.txt")
#st2 = jrr.mage.open_stacked_spectrum(mage_mode, "magestack_bystars_highZ_spectrum.txt")
#plot_winds_neutral_stellar("MagEstack_byZ/", (st1.wave, st2.wave), (st1.X_avg/st1.fnu_autocont,st2.X_avg/st2.fnu_autocont), (st1.X_sigma/st1.fnu_autocont, st2.X_sigma/st2.fnu_autocont), (0.0, 0.0, 0.0), vwin, Ncol, "stacked by Z: lowZ (blue), highZ (black)", LL, 0.0, colortab=('blue', 'black'))


print "STATUS:  Making same figure as Heckman et al. Figure 1 but for our sample"
line_label_Heck  = ('S II 1260', 'C II 1334', 'Si III 1206', 'Si IV 1393', 'N II 1084')  
line_center_Heck = array((1260.4221, 1334.5323,   1206.500, 1393.76, 1084.5659))  
jrr.plot.velocity_overplot(sp.wave, sp.X_avg, line_label_Heck, line_center_Heck, redshift, -1600, 500, (8,5))
plt.savefig("MagEstack/like-heckman2015fig1.pdf", bbox_inches='tight', pad_inches=0.1)
#sys.exit()
print "STATUS: again, but removing the transitions blueward of Lya."
line_label_Heck  = ('S II 1260', 'C II 1334',  'Si IV 1393')  
line_center_Heck = array((1260.4221, 1334.5323, 1393.76))  
jrr.plot.velocity_overplot(sp.wave, sp.X_avg, line_label_Heck, line_center_Heck, redshift, -1600, 500, (8,5))
plt.savefig("MagEstack/like-heckman2015fig1-onlytrans_redwardlya.pdf", bbox_inches='tight', pad_inches=0.1)


# Conclusion:  it takes about 2500 km/s for the continuum to recover blueward of the SiIV,
# C IV, and Al III wind lines.  This must be related to the max or term velocity of the wind.
# Is this a new measurement?  What does it mean?
# How quickly do the ISM lines recover?  (get from FeII 2600, MgII)
# Can Separating ISM from wind kinematics.
# And then go read stis papers; is the high-vel tail interesting?


#print "STATUS:  Plotting  the Leitherer et al. 2011 composite for comparison"
#(leith) = jrr.mage.open_Leitherer_2011_stack() 
#leith['Luncert'] = zeros_like(leith.avg_flux)
#plot_winds_neutral_stellar("Leitherer2011/", (leith.wave,), (leith.avg_flux,), (leith.Luncert,), redshift, vwin, Ncol, "Leitherer et al. 2011")

#print "STATUS:  making velocity plots of the stacked Crowther et al. 2016 spectrum"
#sp = jrr.mage.open_Crowther2016_spectrum()
# Stopped here, got stuck -- not continuum normalized.  ** resume from here.
#plot_winds_neutral_stellar("Crowther2016/", (sp.wave,), (sp.fnu,), (sp.X_jack_std,), (redshift,), vwin, Ncol, "MagE median", LL, z_sys)
#plt.clf


print "STATUS:  Making velocity plots of wind lines  for all the MagE spectra.  May see something in those w high SNR"
(specs) = jrr.mage.getlist_wcont(mage_mode) 
for ii in range(0, len(specs)) :
    zz =  specs['z_neb'][ii]
    label = specs['short_label'][ii]
    filename  = specs['filename'][ii]
    linelist = line_path + re.sub(".txt", ".linelist", (re.split("/", filename)[-1]))
    (LL, z_sys) = jrr.mage.get_linelist(linelist)
    (sp, resoln, dresoln) = jrr.mage.open_spectrum(filename, zz, mage_mode)
    fnu_norm   = sp.fnu / sp.fnu_cont
    fnu_norm_u = jrr.util.sigma_adivb(sp.fnu, sp.fnu_u, sp.fnu_cont, sp.fnu_cont_u)
    prefix = "All_Mage/" + label
    plot_winds_neutral_stellar(prefix, (sp.wave,), (fnu_norm,), (fnu_norm_u,), (zz,), vwin, Ncol, label, LL, z_sys)    

plot_onepagers = True  # This step is SLOW.  Turn off while debugging above.
if plot_onepagers :
    print "STATUS: Making velocity plots of all the MagE spectra on one page, for a bunch of lines"
    for ii in range(0, len(line_label_all)) :
        print "   Plotting ", line_label_all[ii]
        jrr.mage.plot_1line_manyspectra(line_center_all[ii], line_label_all[ii], 4000., True, mage_mode) 
        foo = line_label_all[ii]
        outfile = "Each_line_all_spectra/" + foo.replace(" ", "") + ".pdf"
        plt.savefig(outfile, orientation='portrait', bbox_inches='tight', pad_inches=0.1)
        plt.close()
    #    jrr.mage.plot_1line_manyspectra(1334.5323, 'C II 1334', 4000.)   #velocity plot


print "STATUS: Just plot CIV for a few spectra"
labels = ['rcs0327-E', 'S0900+2234', 'S1527+0652', 'S1429+1202', 'S1226+2152', 'Cosmic~Eye', 'S1458-0023']
specs = jrr.mage.getlist_labels(mage_mode, labels)

line_label   = "C IV 1548"
line_center  = array((1548.19))
jrr.mage.plot_1line_manyspectra(line_center, line_label, 4000., True, mage_mode, specs)
plt.savefig("CIV_excerpts.pdf", orientation='portrait', bbox_inches='tight', pad_inches=0.1) 

print "STATUS:  Plotting C III] for Cycle 24 HST proposal"
line_label = ("C III] 1907")
line_center  = array((1906.68))
jrr.mage.plot_1line_manyspectra(line_center, line_label, 10., False, mage_mode, specs)
plt.savefig("CIII_C24.pdf", bbox_inches='tight', pad_inches=0.1)
plt.close()

jrr.mage.plot_1line_manyspectra(line_center, line_label, 10., False, mage_mode)
plt.savefig("CIII_all.pdf", bbox_inches='tight', pad_inches=0.1)
plt.close()

# Plot OVI for John Chisholm, just those w decent SNR.
(line_label, line_center) = ("OVI 1031", 1031.9261)
labels = ['S1226+2152', 'S1527+0652', 'S1527+0652-fnt', 'S1429+1202', 'S1458-0023', 'S2111-0114',  'S0033+0242']
specs =  jrr.mage.getlist_labels(mage_mode, labels)
jrr.mage.plot_1line_manyspectra(line_center, line_label, 20., False, mage_mode, specs, size=(8,16), fontsize=20)
plt.tight_layout()
plt.savefig("OVI_just_decentSNR.pdf", orientation='portrait', bbox_inches=None, pad_inches=0.) 
 
