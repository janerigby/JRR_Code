from numpy import array, concatenate
import jrr
import sys
import re
from numpy import zeros_like
from matplotlib import pyplot as plt
import pandas
import numpy as np
#mage_mode = "reduction"
mage_mode = "released"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)

def plot_winds_neutral_stellar(prefix, thewaves, thefnus, thedfnus, thezs, label="", LL=[], z_sys=0.0, colortab=False) :
    ''' thewaves, thefnus, thedfnus, thezs are TUPLES of arrays of wavelength, fnu, sigma, and redshift.  If only plotting one, use thewaves=(wave_array,) '''
    ymax= [1.5]
    vwin = 4000. # +- velocity window (km/s) to consider a line
    Ncol = 1
    jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_a, line_center_a, vwin, Ncol, LL, extra_label=label, colortab=colortab, ymax=ymax, verbose=False)
    plt.savefig(prefix + "a.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.clf()
    jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_b, line_center_b, vwin, Ncol, LL, extra_label=label, colortab=colortab, ymax=ymax, verbose=False)
    plt.savefig(prefix + "b.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.clf()
    jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_c, line_center_c, vwin, Ncol, LL, extra_label=label, colortab=colortab, ymax=ymax, verbose=False)
    plt.savefig(prefix + "c.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.clf()
    plt.close("all")
    return(0)

def plot_photospheric_lines(prefix, thewaves, thefnus, thedfnus, thezs, label="", LL=[], z_sys=0.0, colortab=False) :
    ''' thewaves, thefnus, thedfnus, thezs are TUPLES of arrays of wavelength, fnu, sigma, and redshift.  If only plotting one, use thewaves=(wave_array,) '''
    ymax= [1.5]
    win = 13
    Ncol = 4
    jrr.plot.boxplot_Nspectra(thewaves, thefnus, thedfnus, thezs, line_label_p, line_center_p, win, Ncol, LL, extra_label=label, figsize=(16,8), vel_plot=False, colortab=colortab, ymax=ymax, verbose=False)
    plt.suptitle(prefix)
    plt.savefig(prefix + "photospheric.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.clf()
    plt.close("all")
    return(0)
    
# Define the lines to plot
line_label_a  = ('Lya', 'Si II 1260', 'O I 1302',  'C II 1334', 'Si II 1526', 'Al II 1670') # Nino and Marc R. say these trace neutral gas.
  #Added SiII1260, AlII, AlIII at Rongmon's suggestion
line_center_a = array((1215.6701, 1260.4221, 1302.1685, 1334.5323,  1526.7066, 1670.7874))
line_label_b  =  ('MgII2796', 'FeII2344', 'FeII2383')
line_center_b = array((2796.352,  2344.214,  2382.765))
line_label_c  = ("N V 1238", "Si IV 1393", "C IV 1548", "Al III 1854")
line_center_c = array((1238.82, 1393.76, 1548.19, 1854.72))
line_label_extra  = ("OVI 1031", "OVI 1037", "CIII] 1907", "C II] 2326", "OIII] 1660", "O II] 2470", "SiII] 2335", "HeII1640", "SiIII 1882", "SiIII 1892", "LyB", "SiIV1402", "Fe II 2600", "S II 1259")
line_center_extra =array((1031.9261, 1037.6167, 1906.68, 2326.00, 1660.81, 2470.97, 2335.123, 1640.42,  1883.00, 1892.03, 1025.7223, 1402.770, 2600.1729, 1259.519))
line_label_all  = line_label_a  + line_label_b  + line_label_c + line_label_extra
line_center_all = concatenate((line_center_a, line_center_b, line_center_c, line_center_extra))

# Photospheric features that John Chisholm and JR identified as strong in the stack.
line_label_p  = ("CIII 1247", "C_III 1296", "CII 1323", "OIV 1343", "SiIII 1417", "SiII 1485", "SV 1501",
                 "CIII 1620", "FeV 1662", "FeIV 1717", "FeIII 1930 complex", "FeIII 1956 complex", "CIII 2297")
line_center_p = np.array((1247.38, 1296.33,      1323.93,   1343.514,  1417.24,    1485.40,        1501.76,
                  1620.40,     1662.32,   1717.90,      1930.39,               1953.33, 2297.58))
    
# housekeeping
zz = 0.0  # stacked spectrum is already in rest frame wavelength
pdir = "Photospheric"

def plot_wind_stack(sanitycheck_S99=False) :
    stack_choices = ("standard", "Stack-A")  
    outdir = ("StdStack", "StackA")
    linelist = line_path + "stacked.linelist"
    (LL, zz) = jrr.mage.get_linelist(linelist)
    for ii, which_stack in enumerate(stack_choices) :
        print "STATUS:  Plotting wind lines for MagE stack ", outdir[ii]
        sp = jrr.mage.open_stacked_spectrum(mage_mode, which_stack=which_stack)
        plot_winds_neutral_stellar(outdir[ii]+"/"+outdir[ii]+"-wtdavg-", (sp.wave,), (sp.X_avg,), (sp.X_sigma,), (zz,), "MagE wtdavg "+outdir[ii], LL, zz)       
        plot_winds_neutral_stellar(outdir[ii]+"/"+outdir[ii]+"-median-", (sp.wave,), (sp.X_median,), (sp.X_jack_std,), (zz,), "MagE median "+ outdir[ii], LL, zz)
        plot_photospheric_lines(pdir+"/"+outdir[ii]+"-wtdavg-", (sp.wave,), (sp.X_avg,), (sp.X_sigma,),       (zz,), "", LL, zz)
        plot_photospheric_lines(pdir+"/"+outdir[ii]+"-median-", (sp.wave,), (sp.X_median,), (sp.X_jack_std,), (zz,), "", LL, zz)
        if which_stack == "Stack-A"  :  # S99 was fit to Stack-A, not to the standard stack
            plot_winds_neutral_stellar(outdir[ii]+"/"+outdir[ii]+"-wtdavgwS99-", (sp.wave,sp.wave), (sp.X_avg,    sp.fnu_s99model), (sp.X_sigma,    sp.wave*0), (zz,zz), "MagE wtdavg StackA", LL, zz)
            plot_photospheric_lines(pdir+"/"+outdir[ii]+"-wtdavgwS99-",          (sp.wave,sp.wave), (sp.X_avg,    sp.fnu_s99model), (sp.X_sigma,    sp.wave*0), (zz,zz), "", LL, zz)
            plot_winds_neutral_stellar(outdir[ii]+"/"+outdir[ii]+"-medianwS99-", (sp.wave,sp.wave), (sp.X_median, sp.fnu_s99model), (sp.X_jack_std, sp.wave*0), (zz,zz), "MagE median StackA", LL, zz)
            plot_photospheric_lines(pdir+"/"+outdir[ii]+"-medianwS99-",          (sp.wave,sp.wave), (sp.X_median, sp.fnu_s99model), (sp.X_jack_std, sp.wave*0), (zz,zz), "", LL, zz)
            print "STATUS:  Making plot.echelle_spectrum multipage plots for Stack-A, w S99 fits"
            sp2 = sp.copy(deep=True)
            sp['temp_fnu']   = sp['fnu_norm']           # Plot echelle spectrum w the S99 fit  
            sp['temp_fnu_u'] = sp['fnu_norm_u']
            sp['temp_cont']  = sp['fnu_cont']
            sp2['temp_fnu']   = sp2['fnu_s99model']
            sp2['temp_fnu_u'] = sp2['fnu_s99model']*0
            sp2['temp_cont']  = sp2['fnu_cont']
            print "STATUS:  Plotting echelle format spectra, w S99 fits"
            jrr.plot.echelle_spectrum((sp[~sp['badmask']],sp2[~sp2['badmask']]), (0.0,0.0), outfile=outdir[ii]+"/"+outdir[ii]+"-multipanel_wtdavg_stack_wS99_norm.pdf",   title="", norm_by_cont=True, plot_cont=True, colwave='rest_wave', colfnu='temp_fnu', colfnu_u='temp_fnuu', colcont='temp_cont', verbose=False)
            jrr.plot.echelle_spectrum((sp[~sp['badmask']],sp2[~sp2['badmask']]), (0.0,0.0), outfile=outdir[ii]+"/"+outdir[ii]+"-multipanel_wtdavg_stack_wS99_nonorm.pdf", title="", norm_by_cont=False, plot_cont=True, colwave='rest_wave', colfnu='temp_fnu', colfnu_u='temp_fnuu', colcont='temp_cont', verbose=False)
            sp['temp_fnu']  = sp['X_median']  
            sp['temp_fnuu'] = sp['X_jack_std']
            jrr.plot.echelle_spectrum((sp[~sp['badmask']],sp2[~sp2['badmask']]), (0.0,0.0), outfile=outdir[ii]+"/"+outdir[ii]+"-multipanel_median_stack_wS99_norm.pdf", title="", norm_by_cont=True, plot_cont=True, colwave='rest_wave', colfnu='temp_fnu', colfnu_u='temp_fnuu', colcont='temp_cont', verbose=False)
            jrr.plot.echelle_spectrum((sp[~sp['badmask']],sp2[~sp2['badmask']]), (0.0,0.0), outfile=outdir[ii]+"/"+outdir[ii]+"-multipanel_median_stack_wS99_nonorm.pdf", title="", norm_by_cont=False, plot_cont=True, colwave='rest_wave', colfnu='temp_fnu', colfnu_u='temp_fnuu', colcont='temp_cont', verbose=False)
            
            if sanitycheck_S99 :
                (S99, ignore_linelist) = jrr.mage.redo_open_S99_spectrum("Stack-A", denorm=True)
                plot_winds_neutral_stellar(outdir[ii]+"/"+outdir[ii]+"-s99-", (S99.rest_wave,), (S99.rest_fnu_s99,), (S99.rest_fnu_s99_u,), (zz,), "S99 fit", LL, zz, colortab=("blue"))
                plot_photospheric_lines(pdir+"/"+outdir[ii]+"-s99-",          (S99.rest_wave,), (S99.rest_fnu_s99,), (S99.rest_fnu_s99_u,), (zz,), "", LL, zz, colortab=("blue"))
                plot_winds_neutral_stellar(outdir[ii]+"/"+outdir[ii]+"-s99sanity-",  (S99.rest_wave,S99.rest_wave), (S99.rest_flam_data_norm, S99.rest_flam_s99_norm), (S99.rest_flam_data_u_norm,S99.rest_flam_data_u_norm), (zz,zz), "Sanity-check S99 norm", LL, zz)
                plot_photospheric_lines(pdir+"/"+outdir[ii]+"-s99sanity-",  (S99.rest_wave,S99.rest_wave), (S99.rest_flam_data_norm, S99.rest_flam_s99_norm), (S99.rest_flam_data_u_norm,S99.rest_flam_data_u_norm), (zz,zz), "", LL, zz)
                S99['temp_fnu']    = S99['rest_fnu_S99']
                S99['temp_fnuu']   = S99['rest_fnu_s99'] * 0.0
                S99['temp_cont']   =  jrr.spec.rebin_spec_new(sp['rest_wave'], sp['rest_fnu_autocont'], S99['rest_wave'])
                sp['temp_fnu']  = sp['X_avg']  
                sp['temp_fnuu'] = sp['X_sigma']
                sp['temp_cont'] = sp['rest_fnu_autocont']
                jrr.plot.echelle_spectrum((sp[~sp['badmask']],S99), (0.0,0.0), outfile=outdir[ii]+"/"+outdir[ii]+"-multipanel_wtdavg_stack_wS99_nonorm.pdf", title="", norm_by_cont=False, plot_cont=True, colwave='rest_wave', colfnu='temp_fnu', colfnu_u='temp_fnuu', colcont='temp_cont', verbose=False)
                sp['temp_fnu']  = sp['X_median']  
                sp['temp_fnuu'] = sp['X_jack_std']
                jrr.plot.echelle_spectrum((sp[~sp['badmask']],S99), (0.0,0.0), outfile=outdir[ii]+"/"+outdir[ii]+"-multipanel_median_stack_wS99_nonorm.pdf", title="", norm_by_cont=False, plot_cont=True, colwave='rest_wave', colfnu='temp_fnu', colfnu_u='temp_fnuu', colcont='temp_cont', verbose=False)

                
        print "STATUS:  Making same figure as Heckman et al. Figure 1 but for our sample"
        line_label_Heck  = ('S II 1260', 'C II 1334', 'Si III 1206', 'Si IV 1393', 'N II 1084')  
        line_center_Heck = array((1260.4221, 1334.5323,   1206.500, 1393.76, 1084.5659))  
        jrr.plot.velocity_overplot(sp.wave, sp.X_avg, line_label_Heck, line_center_Heck, zz, -1600, 500, (8,5))
        plt.savefig(outdir[ii]+"/"+outdir[ii]+"-likeheckman2015fig1.pdf", bbox_inches='tight', pad_inches=0.1)
        plt.clf()
        print "STATUS: again, but removing the transitions blueward of Lya."
        line_label_Heck  = ('S II 1260', 'C II 1334',  'Si IV 1393')  
        line_center_Heck = array((1260.4221, 1334.5323, 1393.76))  
        jrr.plot.velocity_overplot(sp.wave, sp.X_avg, line_label_Heck, line_center_Heck, zz, -1600, 500, (8,5))
        plt.savefig(outdir[ii]+"/"+outdir[ii]+"-like-heckman2015fig1-onlytrans_redwardlya.pdf", bbox_inches='tight', pad_inches=0.1)
        plt.clf()
    plt.close("all")
    return(0)
    
# Conclusion:  it takes about 2500 km/s for the continuum to recover blueward of the SiIV,
# C IV, and Al III wind lines.  This must be related to the max or term velocity of the wind.
# Is this a new measurement?  What does it mean?
# How quickly do the ISM lines recover?  (get from FeII 2600, MgII)
# Can Separating ISM from wind kinematics.
# And then go read stis papers; is the high-vel tail interesting?


def plot_wind_indy(sanitycheck_S99=False) :
    outdir_ech = "../Plot-all/PDF_Out2_S99/"
    print "STATUS:  Making velocity plots of wind lines  for all the MagE spectra.  May see something in those w high SNR"
    (specs) = jrr.mage.getlist_wcont(mage_mode) 
    for ii in range(0, len(specs)) :
        label = specs['short_label'][ii]
        (sp, resoln, dresoln, LL, z_sys, boxcar) =  jrr.mage.wrap_open_spectrum(label, mage_mode, addS99=True)
        sp['fnu_norm']   = sp.fnu / sp.fnu_cont
        sp['fnu_norm_u'] = jrr.util.sigma_adivb(sp.fnu, sp.fnu_u, sp.fnu_cont, sp.fnu_cont_u)    
        prefix = "All_Mage/" + label
        print "   plotting ", prefix
        plot_winds_neutral_stellar(prefix, (sp.wave,), (sp.fnu_norm,), (sp.fnu_norm_u,), (z_sys,), label, LL, z_sys)
        plot_photospheric_lines(pdir+"/"+label+"-s99-", (sp.wave,), (sp.fnu_norm,), (sp.fnu_norm_u,), (z_sys,), "", LL, z_sys)
        plt.clf()
        if jrr.mage.getfullname_S99_spectrum(label) :  # If a S99 fit file exists
            # This is the preferred method:  gets the S99 fit from mage.wrap_open_spectrum
            plot_winds_neutral_stellar(prefix+"wS99",      (sp.wave,sp.wave), (sp.fnu_norm, sp.fnu_s99model/sp.fnu_cont), (sp.fnu_norm_u, sp.fnu_s99model*-1), (z_sys, z_sys), label, LL, z_sys)
            plot_photospheric_lines(pdir+"/"+label+"wS99", (sp.wave,sp.wave), (sp.fnu_norm, sp.fnu_s99model/sp.fnu_cont), (sp.fnu_norm_u, sp.fnu_s99model*-1), (z_sys, z_sys), "", LL, z_sys)
            sp2 = sp.copy(deep=True)
            sp['temp_fnu']   = sp['fnu']           # Plot echelle spectrum w the S99 fit  
            sp['temp_fnu_u'] = sp['fnu_u']
            sp['temp_cont']  = sp['fnu_cont']
            sp2['temp_fnu']   = sp2['fnu_s99model']
            sp2['temp_fnu_u'] = sp2['fnu_s99model']*-1
            sp2['temp_cont']  = sp2['fnu_cont']
            jrr.plot.echelle_spectrum((sp[~sp['badmask']],sp2[~sp2['badmask']]), (z_sys, z_sys), LL=LL, outfile=outdir_ech+"multipanel_"+label+"_wS99_norm.pdf", title="", norm_by_cont=True, plot_cont=True, colwave='wave', colfnu='temp_fnu', colfnu_u='temp_fnu_u', colcont='temp_cont', verbose=False)
            jrr.plot.echelle_spectrum((sp[~sp['badmask']],sp2[~sp2['badmask']]), (z_sys, z_sys), LL=LL, outfile=outdir_ech+"multipanel_"+label+"_wS99_nonorm.pdf", title="", norm_by_cont=False, plot_cont=True, colwave='wave', colfnu='temp_fnu', colfnu_u='temp_fnu_u', colcont='temp_cont', verbose=False)
            if sanitycheck_S99 :
                # This is the deprecated, less clean method, used for sanity-checking
                (S99, ignore_linelist) = jrr.mage.redo_open_S99_spectrum(label)
                S99['fnu_s99_normasdata'] = S99['fnu_s99'] / jrr.spec.rebin_spec_new(sp['wave'], sp['fnu_cont'], S99['wave'])
                S99['fnu_data_normasdata'] = S99['fnu_data'] / jrr.spec.rebin_spec_new(sp['wave'], sp['fnu_cont'], S99['wave'])
                S99['fnu_u_normasdata'] = S99['fnu_data_u'] / jrr.spec.rebin_spec_new(sp['wave'], sp['fnu_cont'], S99['wave'])
                plot_winds_neutral_stellar(prefix+"testS99", (S99.wave,S99.wave), (S99.fnu_data_normasdata,S99.fnu_s99_normasdata), (S99.fnu_u_normasdata,S99.wave*0), (z_sys,z_sys), label, LL, z_sys)
                plot_photospheric_lines(pdir+"/"+label+"testS99", (S99.wave,S99.wave), (S99.fnu_data_normasdata,S99.fnu_s99_normasdata), (S99.fnu_u_normasdata,S99.wave*0), (z_sys,z_sys), "", LL, z_sys)
                S99['temp_fnu']   = S99['fnu_s99_normasdata']
                S99['temp_fnu_u'] = S99['fnu_s99_normasdata'] * -1
                S99['temp_cont']  =  jrr.spec.rebin_spec_new(sp['wave'], sp['fnu_cont'], S99['wave'])
                jrr.plot.echelle_spectrum((sp[~sp['badmask']],S99), (z_sys, z_sys), LL=LL, outfile=outdir_ech+"multipanel_"+label+"_wS99_norm.pdf", title="", norm_by_cont=True, plot_cont=True, colwave='wave', colfnu='temp_fnu', colfnu_u='temp_fnu_u', colcont='temp_cont', verbose=False)
                jrr.plot.echelle_spectrum((sp[~sp['badmask']],S99), (z_sys, z_sys), LL=LL, outfile=outdir_ech+"multipanel_"+label+"_wS99_nonorm.pdf", title="", norm_by_cont=False, plot_cont=True, colwave='wave', colfnu='temp_fnu', colfnu_u='temp_fnu_u', colcont='temp_cont', verbose=False)
    plt.close("all")
    return(0)

            
def plot_onepagers() : # This step is SLOW.  Turn off while debugging rest of script
    print "STATUS: Making velocity plots of all the MagE spectra on one page, for a bunch of lines"
    for ii in range(0, len(line_label_all)) :
#    for ii in range(0, 1) :  # For Lya only
        print "   Plotting ", line_label_all[ii]
        jrr.mage.plot_1line_manyspectra(line_center_all[ii], line_label_all[ii], 4000., True, mage_mode, fontsize=10) 
        foo = line_label_all[ii]
        outfile = "Each_line_all_spectra/" + foo.replace(" ", "") + ".pdf"
        #outfile = "Each_line_all_spectra/" + foo.replace(" ", "") + "_bigger-yrange.pdf"  # For Lya only
        plt.savefig(outfile, orientation='portrait', bbox_inches='tight', pad_inches=0.1)
        plt.clf()
    plt.close("all")
    return(0)

def plot_some_CIV() :        
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
    return(0)

def plot_OVI_forJC() :
    # Plot OVI for John Chisholm, just those w decent SNR.
    (line_label, line_center) = ("OVI 1031", 1031.9261)
    labels = ['S1226+2152', 'S1527+0652', 'S1527+0652-fnt', 'S1429+1202', 'S1458-0023', 'S2111-0114',  'S0033+0242']
    specs =  jrr.mage.getlist_labels(mage_mode, labels)
    jrr.mage.plot_1line_manyspectra(line_center, line_label, 20., False, mage_mode, specs, size=(8,16), fontsize=20)
    plt.tight_layout()
    plt.savefig("OVI_just_decentSNR.pdf", orientation='portrait', bbox_inches=None, pad_inches=0.) 
    return(0)


    
#######################################################
# What do I want to run today?  Running all of them is slow; I did one at a time
#plot_wind_stack()    # outdir = ("StdStack", "StackA")
plot_wind_indy()      # All_Mage/,  ../Plot-all/PDF_Out2_S99/
#plot_onepagers()    # Each_line_all_spectra/
#plot_some_CIV()     # .
#plot_OVI_forJC()     # .
#######################################################

#print "STATUS:  Plotting  the Leitherer et al. 2011 composite for comparison"
#(leith) = jrr.mage.open_Leitherer_2011_stack() 
#leith['Luncert'] = zeros_like(leith.avg_flux)
#plot_winds_neutral_stellar("Leitherer2011/", (leith.wave,), (leith.avg_flux,), (leith.Luncert,), redshift, "Leitherer et al. 2011")

#print "STATUS:  making velocity plots of the stacked Crowther et al. 2016 spectrum"
#sp = jrr.mage.open_Crowther2016_spectrum()
# Stopped here, got stuck -- not continuum normalized.  ** resume from here.
#plot_winds_neutral_stellar("Crowther2016/", (sp.wave,), (sp.fnu,), (sp.X_jack_std,), (redshift,), "MagE median", LL, z_sys)
#plt.clf
