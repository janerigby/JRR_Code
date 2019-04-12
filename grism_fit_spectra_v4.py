from __future__ import print_function
# This is a script to fit the 1D G102 and G141 WFC/HST grism spectra for S1723 and S2340.
# It used to be in IDL (fit_g141.pro and fit_G102.pro), which was used to fit NIRSPEC
# spectra of RCS0327 knots, as published in Wuyts et al 2014.
# Now ported to Python, using LMFIT as the fitting package.
# ** This code assumes the continuum has already been removed.**
# Run from each of these 2 dirs:  Dropbox/Grism_S2340/WFC3_fit_1Dspec,  Dropbox/Grism_S1723/WFC3_fit_1Dspec
# Run twice in each dir, for G102, and G141.  So, a total of 4 infiles to process both galaxies.
#; jrigby, feb 2012, march 2012, oct 2016, may 2018

from builtins import str
import numpy as np
import pandas
import jrr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import lmfit 
from os.path import expanduser, basename, abspath, exists
from os import makedirs
import glob
import re
import time

matplotlib.rcParams.update({'font.size': 14})

def get_filenames(dir, which_grism) :
    myfiles = [ basename(x) for x in glob.glob(dir + "*" + which_grism + "*.txt") ]
    return(myfiles)

def error_unknown_grism(which_grism) :
    raise Exception('ERROR: Cannot understand which_grism', which_grism)
    return(0)

def get_redshift(which_gal) :
    if which_gal == 'S1723' :
        zz = 1.331366
        #datadir = expanduser("~") + '/Dropbox/Grism_S1723/Grizli_redux/1Dsum/Wcont/'
    elif which_gal == 'S2340' :
        zz = 1.4275 
        #datadir = expanduser("~") + '/Dropbox/Grism_S2340/Grizli_redux/1Dsum/Wcont/'
    else : raise Exception("Unrecognized galaxy:", which_gal)
    return(zz)

def get_line_wavelengths(which_grism) :
    if which_grism == "G141" :
        restwaves   = np.array((4862.683,  4960.295,  5008.240,  5877.59,  6310.,            6549.85,  6564.61,  6585.28, 6725.))
        yaxcoords   =         [0.25,       -1,       0.96,        0.1,     0.1,               -1,       0.8,      -1,       0.2] 
        plotnames   =         [r'H$\beta$', '',    '[O III]',    'He I',  '[O I]+[S III]',   '',      r'H$\alpha$', '',       '[S II]']
        linenames   =         ['Hbeta',   '[O~III]', '[O~III]',  'He~I',  '[O~I]+[S~III]',  '[N~II]', 'Halpha', '[N~II]', '[S~II]']
        #                       0             1            2         3          4             5           6           7         8

    elif which_grism == "G102" :
        restwaves = np.array((3727.092,   3729.875,   3869.86,  3890.151,  3968.593,   3971.195,   4025.,        4102.892, 4341.684, 4364.436, 4472.7,  4687.02,   4741.4489, 4862.683))
        linenames =         ['[O~II]',    '[O~II]',  '[Ne~III]', 'H8+HeI', '[Ne~III]', 'Heps',     'HeI+HeII',  'Hdelt',  'Hgam',  '[O~III]', 'He~I', 'He~II',    '[Ar~IV]', 'Hbeta']
        #                       0             1            2         3          4         5           6           7          8        9        10        11         12           13
        plotnames =         ['[O II]',    '',        '[Ne III],H8,He I', '', r'[Ne III], H$\epsilon$', '',  'HeI+HeII',  r'H$\delta$',  '', r'H$\gamma$, [O III]', 'He I', 'He II',    '[Ar IV]', r'H$\beta$']
        yaxcoords =         [  0.95,       -1,        0.35,              0.,    0.2,                   0.,   -0.13,       0.2,           0., 0.35,                     0.1,     0.1,      0.2,          0.7]

    else : error_unknown_grism(which_grism)
    label_df = pandas.DataFrame({'restwave' : restwaves, 'yaxcoord' : yaxcoords, 'plotname' : plotnames})  # For convenience, make dataframe of labels to plot
    return(restwaves, linenames, label_df)

# SHould have done the below functions with just a single "params" var that gets unpcked, as in v=params.valuesdict().  See https://lmfit.github.io/lmfit-py/parameters.html
# Also, probably should have done all this in astropy, which has similar functionality. Moving on.

def fitfunc_G141(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8) : # A custom G141 fitting function
    # Parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, its clearer to have params be FLUXES.
    # This version fits a global redshift zz, and allows no wavelength miscalibration between lines.
    fluxes          = np.array((f0,        f1,        f2,        f3,         f4,          f5,          f6,       f7,  f8))
    (restwaves, linenames, label_df) = get_line_wavelengths("G141")
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G141'))

def fitfunc_G141_waveoff(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, d0, d1, d2, d3, d4, d5, d6, d7, d8) : # A custom G141 fitting function
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, its clearer to have params be FLUXES.
    # This version lets each line have a dwave offset, to deal w relative wavelength calibration problems in grism spectra
    fluxes          = np.array((f0,        f1,        f2,        f3,      f4,      f5,      f6,       f7, f8))
    dwaves       = np.array((d0, d1, d2, d3, d4, d5, d6, d7, d8))
    (restwaves, linenames, label_df) = get_line_wavelengths("G141")
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G141', dwaves=dwaves))

def fitfunc_G102(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13) : # A custom G102 fitting function
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, it is much easier to have params be FLUXES.  Easier to report, too
    fluxes    = np.array((f0,         f1,        f2,        f3,       f4,       f5,        f6,       f7,       f8,       f9,       f10,      f11,      f12, f13))
    (restwaves, linenames, label_df) = get_line_wavelengths("G102")
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G102'))

def fitfunc_G102_waveoff(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13) : 
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, it is much easier to have params be FLUXES.  Easier to report, too
    fluxes    = np.array((f0,         f1,        f2,        f3,       f4,       f5,        f6,       f7,       f8,       f9,       f10,      f11,      f12, f13))
    dwaves       = np.array((d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13))
    (restwaves, linenames, label_df) = get_line_wavelengths("G102")
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G102', dwaves=dwaves))

def fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism, dwaves=()) :  # This part is the same for either grism
    if len(dwaves) == 0 :  dwaves = np.zeros_like(restwaves)
    grism_info = jrr.grism.get_grism_info(which_grism) 
    sigmas      = restwaves * (1+zz) / ( grism_info['R'] *2.35482/morph_broad)
    y = np.zeros_like(wave)
    for ii, restwave in enumerate(restwaves) :
        aa = fluxes[ii] / sigmas[ii] / np.sqrt(2.*np.pi)  #  Convert from gaussian flux to gaussian height
        y += jrr.spec.onegaus(wave, aa, (restwaves[ii]*(1+zz) + dwaves[ii]), sigmas[ii], 0.0)
    return(y)

def pick_fitting_function(which_grism, waveoff=False) :
    if   which_grism == 'G141':
        if waveoff: return(fitfunc_G141_waveoff)
        else :      return(fitfunc_G141)
    elif which_grism == 'G102':
        if waveoff : return(fitfunc_G102_waveoff)
        else:        return(fitfunc_G102)
    else : error_unknown_grism(which_grism)

def prep_params(which_grism, waveoff=False, scale=1.) : # Create param container for fitfunc_G141 or fitfunc_G102
    if   which_grism == 'G141': 
        guesses =  np.array((             250.,  500., 1500., 25.,  20., 16.,  800.,  50.,  60.))*scale
        if waveoff :
            parnames = ('zz', 'morph_broad', 'f0',  'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'd0', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8') #fitfunc_G141_waveoff
            guesses  = np.concatenate((guesses, np.zeros(len(guesses))))   # initialize the dwavelengths
        else:
            parnames = ('zz', 'morph_broad', 'f0',  'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8')  # Par names from fitfunc_G141

    elif which_grism == 'G102' :
        guesses = np.array((              144,  200, 100,   20,   50, 30,   20,   50,  120,   10, 10,  10,    5,   250))*scale       
        if waveoff : 
            parnames = ('zz', 'morph_broad', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11', 'f12', 'f13', 'd0', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'd10', 'd11', 'd12', 'd13') #fitfunc_G102_waveoff
            guesses  = np.concatenate((guesses, np.zeros(len(guesses))))  # initialize the dwavelengths
        else : 
            parnames = ('zz', 'morph_broad', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11', 'f12', 'f13')  # Par names from fitfunc_G102
    else : error_unknown_grism(which_grism)
    return(guesses, parnames)

def set_params(mymodel, parnames, guesses, zz, morph_broad=1.) :  # Set initial values of params for fitfunc_G141 or fitfunc_G102
    mypars = mymodel.make_params()
    mypars.add('zz', value=zz, vary=True)
    mypars.add('morph_broad', value=morph_broad, vary=True)
    for ii, parname in enumerate(parnames[2:]) :
        mypars.add(parname, value=guesses[ii], vary=True)
    return(mypars)

def lock_params(mypars, which_grism, which_gal, waveoff=False, sigoff=0.0, limit_wave=False) :
    grism_info = jrr.grism.get_grism_info(which_grism)
    if which_grism == 'G141' and waveoff :
        mypars['d0'].set(expr='d2')  # Hbeta, 4959, 5007 all shift in wavelength together
        mypars['d1'].set(expr='d2')  # ditto
        mypars['d4'].set(expr='d6')  # 6310 shifts w HA
        mypars['d5'].set(expr='d6')  # Halpha and [N II] doublet shift together in wavelength
        mypars['d7'].set(expr='d6')
        mypars['d8'].set(expr='d6')
    if which_grism == 'G102' and waveoff :
        mypars['d1'].set(expr='d0')  # 3727, 3729 shift in wavelength together
        mypars['d3'].set(expr='d2')  #  #Other closely spaced line groups should move together
        mypars['d5'].set(expr='d4')  # 
        mypars['d6'].set(min=-30, max=30) # Don't let the 4025 line wander too far
        mypars['d9'].set(expr='d8')  # 
        mypars['d10'].set(expr='d8')  # 
        mypars['d11'].set(expr='d13')  #
        mypars['d12'].set(expr='d13')  #
        if which_gal == 'S2340'    : mypars['d13'].set(value=0., vary=False)  # keep this from wandering
            
    if waveoff and limit_wave :   # May be necessary to bound the delta wavelengths, but bounds makes lmfit really slow, so try not to use it. 
        dwave_keys = [x for x in list(mypars.keys()) if re.match('d', x)]  # find all the d0... dX wavelength offset pars
        for dwave in dwave_keys :
            mypars[dwave].set(min= -1. * sigoff * grism_info['wave_unc'])
            mypars[dwave].set(max=       sigoff * grism_info['wave_unc'])
            
    if which_grism == 'G141' :
        mypars['f1'].set(expr='f2*0.33189')           #  Lock [O III] 4959, 5007 ratio to Storey & Zeippen2000
        if which_gal == 'S1723' :  
            mypars['f7'].set(expr='f6*0.062')         #  Lock NII/Ha ratio to what was observed for GNIRS
            mypars['f5'].set(expr='f6*0.062/3.071')   #  Lock [N II]  6549, 6585 ratio
        else :  # Any other galaxy
            #mypars['f5'].set(expr='f7/3.071')         #  Lock [N II]  6549, 6585 ratio
            mypars['f5'].set(value=0., vary=False)      # Can't resolve NII from Ha in grism, to set to zero unless we've measured it elsewhere (for S1723)
            mypars['f7'].set(value=0., vary=False)
    if which_grism == 'G102' :
        mypars['f4'].set(expr='f2*0.316555')          # Lock [Ne III] flux ratio to Storey & Zeippen2000
        if  which_gal == 'S2340' :
            mypars['f11'].set(value=0., vary=False)      # These lines fall nar or off the edge of the grism
            mypars['f12'].set(value=0., vary=False) 
            mypars['f13'].set(value=0., vary=False)
        elif which_gal == 'S1723' :
            mypars['f1'].set(expr='f0*1.3699') # Lock [O III] flux ratio to what was measured in ESI
        if which_gal != 'S1723' :
            mypars['f1'].set(value=0., vary=False)  # Cannot resolve OII doublet, so set second trans to zero.
            
    return(mypars)


def finergrid_result(grism_info, fit_result, x, scalefactor, outfile1, outfile2) :
    #The fit result is on the same coarse wavelength array as the data.  Make finer fit for plotting, and save it.
    xfine = np.linspace(grism_info['x1'], grism_info['x2'], 1000)
    fit_df     = pandas.DataFrame(x, columns=('wave',))
    finefit_df = pandas.DataFrame(xfine, columns=('wave',))
    fit_df['bestfit']     = func2fit(x, **fit_result.values) / scalefactor    
    finefit_df['bestfit'] = func2fit(xfine, **fit_result.values) / scalefactor
    fit_df.to_csv(    outfile1, float_format='%.5E', index=False)
    finefit_df.to_csv(outfile2, float_format='%.5E', index=False)
    return(xfine)

def plot_results(df_data, fit_result, func2fit, xfine, title, grism_info, show_initial_fit=False, scalefactor=1.0, units=1.0) :
    print("   Plotting")
    fig, ax  = plt.subplots(figsize=figsize)
    df_data.plot(x='wave', y='flam_contsub_scaled', color='black', linestyle='steps-mid', lw=1.5, legend=False, ax=ax)
    df_data.plot(x='wave', y='flam_u_scaled', color='grey', ax=ax, legend=False)
    if show_initial_fit : plt.plot(df_data['wave'], fit_result.init_fit, color='orange', label='init fit')
    plt.plot(xfine, func2fit(xfine, **fit_result.values), color='blue', label='best fit')
    plt.xlabel(r"observed wavelength ($\AA$)")
    #plt.ylabel(r"continuum-subtracted  $f_\lambda$ ("+str(1./scalefactor) + " " + units + ")")
    plt.ylabel(r"$f_\lambda$ ("+str(1./scalefactor) + " " + units + ")")
    pretty_title = re.sub("wcontMWdr.txt ", "", re.sub("_", " ", title))
    plt.title(pretty_title, position=(0.5, 1.2), fontsize=12)
    plt.xlim(grism_info['x1'], grism_info['x2'])
    (restwaves, linenames, label_df) = get_line_wavelengths(which_grism)
    if 'd0' in list(fit_result.params.keys()) :
        offsets = np.array([fit_result.params[x].value for x in list(fit_result.params.keys()) if re.match('d', x)])
        plt.scatter(restwaves*(1. + fit_result.best_values['zz']) + offsets, np.zeros_like(restwaves)-0.1, marker='+', color='k')
    else :    plt.scatter(restwaves*(1. + fit_result.best_values['zz']), np.zeros_like(restwaves)-0.1, marker='+', color='k')
    ax2 = ax.twiny()
    ax2.set_xlim( grism_info['x1']/(1. + fit_result.best_values['zz']), grism_info['x2']/(1. + fit_result.best_values['zz']))
    ax2.set_xlabel(r"rest-frame wavelength ($\AA$)", labelpad=10)  # vacuum
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    label_df['ydata'] = label_df['yaxcoord'] *  ax2.axes.get_ylim()[1]
    #print "DEBUGGING\n", label_df.head()
    jrr.plot.annotate_from_dataframe(label_df, xcol='restwave', ycol='ydata', text='plotname', xytext=(0,0), ha='center', fontsize=12)
    plt.show()
    return(fig)

def  convert_lmfitresults_2df(LMresult, parnames, sigoff, grism_info, restwaves, linenames):
    # Reformat LMFIT results as a dataframe
    fkeys   = [x for x in parnames if re.match('f', x)]
    dkeys   = [x for x in parnames if re.match('d', x)]
    flux    = [LMresult.params[key].value  for key in fkeys]
    flux_u  = [LMresult.params[key].stderr for key in fkeys]
    if len(dkeys) == 0 :  # In Method 1, no wavelength offsets allowed, so fill w zeros
        doffAng   = np.zeros(shape=len(fkeys))
        doffAng_u = np.zeros(shape=len(fkeys))
    else :
        doffAng  = [LMresult.params[key].value  for key in dkeys]
        doffAng_u= [LMresult.params[key].stderr for key in dkeys]
    df = pandas.DataFrame({ 'restwave':restwaves, 'linename':linenames, 'flux':flux, 'flux_u':flux_u, 'doffAng':doffAng, 'doffAng_u':doffAng_u})
    df = df[['restwave', 'linename', 'flux', 'flux_u', 'doffAng', 'doffAng_u']] # force reorder
    df['doffsig'] = df['doffAng'].abs() / grism_info['wave_unc']
    df['dbad']    = df['doffsig'] > sigoff 
    return(df)

def check4zero_errorbars(LMresult, parnames) :
    fkeys   = [x for x in parnames if re.match('f', x)]
    flux_u  = np.array([LMresult.params[key].stderr for key in fkeys])
    if np.count_nonzero(flux_u) == 0 :
        return(1)     # FAILS,  errorbars are zero
    else : return(0)  # Passes, errorbars are not zero
            
def supplemental_header(LMresult) :
    suphead = "# REDSHIFT " + str(LMresult.params['zz'].value) + "\n# Z_UNCERT " + str(LMresult.params['zz'].stderr)
    suphead += "\n# MORPH_BROAD " +  str(LMresult.params['morph_broad'].value) + "\n# MB_UNCERT "  + str(LMresult.params['morph_broad'].stderr)
    suphead += "\n# How to read this file in python:  pandas.read_csv(INFILE, comment='#')\n#\n";
    return(suphead)

######################################################################

plt.ion()
##  Setup
#infile = 'S1723_G141_grism2process.txt'
infile = 'S1723_G102_grism2process.txt'
#infile = 'S2340_G102_grism2process.txt'  # CHANGE THIS.  Keep format
#infile = 'S2340_G141_grism2process.txt'  

figsize = (12,4)
scalefactor = 1E17 # Scale everything by scalefactor, to avoid numerical weirdness in LMFIT
units = "erg/s/cm^2"
guess_morphbroad = 1.5
sigoff=3 # Warn if the delta wavelengths exceed this
show_initial_fit = False
print_fitreports = False

''' Status report:
This works pretty well.  There are 2 methods, and the user can pick which to use.  
Method 1 moves all the line centroids in unison, by a global redshift. It's run in 
all cases.  Method 2 fixes the redshift as the result from method 1, and then re-fits, 
allowing the wavelengths of groups of lines to  move together.  The input file "infile"
tells the script whether to run Method 2, and scales the initial guess of linestrength.
Method 2 works great for all S1723 spectra, both grisms.  
Method 2 works great for almost all S2340 spectra. For a few at low SNR it stalls, so I used method 1 instead.  
'''

# Input the list of spectra to fit
m = re.search('(\S+)_(\S+)_(\S+)', infile)   # Formart of input file tells us which galaxy, which grism
which_gal   = m.group(1)
which_grism = m.group(2)
(restwaves, linenames, label_df) = get_line_wavelengths(which_grism)
df = pandas.read_table(infile, comment="#", delim_whitespace=True)  # Read the list of files to process

# Gather info needed for fitting
zz = get_redshift(which_gal)
grism_info = jrr.grism.get_grism_info(which_grism) 
pp = PdfPages("grism_fitspectra_"+which_gal+"_"+which_grism+".pdf")

for row in df.itertuples() :
    specfile = basename(row.filename)
    print(row.filename)
    subdir = re.split("/", row.filename)[-3] + "/" # This should be the dir, like 1Dsum
    if not exists(subdir):  makedirs(subdir)
    tweak_wav = row.tweak_wav
    specfile_dict = jrr.grism.parse_filename(specfile)
    outfile0  = subdir + re.sub('.txt', '.fitreport',   specfile)
    outfile1  = subdir + re.sub('.txt', '_meth1.fitdf', specfile)
    outfile2a = subdir + re.sub('.txt', '_meth1.model', specfile)
    outfile2b = subdir + re.sub('.txt', '_meth1.modelfine', specfile)
    outfile3  = subdir + re.sub('.txt', '_meth2.fitdf', specfile)
    outfile4a = subdir + re.sub('.txt', '_meth2.model', specfile)
    outfile4b = subdir + re.sub('.txt', '_meth2.modelfine', specfile)
    f = open(outfile0, 'w')  
    header = "# Fitting HST grism spectra with grism_fit_spectra_v4.py\n# FILENAME "+ specfile
    header += "\n# FILENAME_WPATH " + row.filename + "\n# WHICH_GAL " + which_gal
    header += "\n# WHICH_GRISM " + which_grism + "\n# ROLL " + specfile_dict['roll'] + "\n# DESCRIP " +  specfile_dict['descrip']
    header += "\n# TIME_FIT " + time.strftime("%c") + "\n# SCALEFACTOR " + str(1./scalefactor) + "\n# UNITS " + units
    header += "\n# Note: For S1723 NII/Ha ratios are set from GNIRS, and [O II] 3727/3729 ratio set from ESI.\n"
    header += "# Note: For S2340, NII fluxes were fixed to zero; what is reported as Ha is actually the blend of [N II]+Ha+[N II].\n"
    header += "# Note on fluxing for S1723/1Dsum: In the bothrolls files, Michael had summed the flux over both rolls.  Here those fluxes\n"
    header += "#                   have been divided by 2, and therefore should match the fluxes of the spectra extracted from a single roll.\n"
    header += "# Note on fluxing for S1723/1Dbyclumps, and all S2340: The fluxing is weird, b/c Michael added multiple clumps from whichever\n"
    header += "# rolls were clean.  So the fluxing is super-weird, and should only be used in a relative sense.  Divide by Hbeta and move on.\n#\n"
    f.write(header)
    sp  = pandas.read_csv(row.filename,  comment="#")      ## Read the grism spectrum file
    if which_gal == 'S1723' and 'both' in specfile :  jrr.grism.half_the_flux(sp)   # Michael has summed the flux over both rolls. Compensating
    sp['flam_contsub_scaled'] = (sp['flam'] - sp['cont']) * scalefactor 
    sp['flam_u_scaled'] = sp['flam_u'] * scalefactor
    sp['weight'] = 1.0 / (sp['flam_u_scaled'])**2   # inverse variance weights
    subset = sp.loc[sp['wave'].between(grism_info['x1'], grism_info['x2'])] # subset of spectr, clean wavelength range.
            
    print("   METHOD 1 fit, to determine the best-fit redshift. Wavelengths fixed.")
    (guesses, parnames) = prep_params(which_grism=which_grism, scale=row.scale_guess)   # Make container to hold the parameters
    func2fit = pick_fitting_function(which_grism, waveoff=False)
    mymodel = lmfit.Model(func2fit, independent_vars=('wave',), param_names=parnames,  which_grism=which_grism)  # Set up a model
    mypars = set_params(mymodel, parnames, guesses, zz, morph_broad=guess_morphbroad)  # Set initial parameters for that model
    locked_params = lock_params(mypars, which_grism=which_grism, which_gal=which_gal, limit_wave=False)
    result1  = mymodel.fit(subset['flam_contsub_scaled'], locked_params, wave=subset['wave'])  # fitting is done here
    zz_fit = result1.best_values['zz']
    if print_fitreports : print(result1.fit_report())
    f.write(result1.fit_report())
    df1 = convert_lmfitresults_2df(result1, parnames, sigoff, grism_info, restwaves, linenames)
    df1.to_csv('/tmp/foo', float_format='%.3f', index=False)    
    jrr.util.put_header_on_file('/tmp/foo', header + supplemental_header(result1), outfile1)
    if check4zero_errorbars(result1, parnames) : print("####### WARNING: errorbars were zero! ####### ")
    plot_label = subdir + " " + specfile + " fit 1"
    xfine = finergrid_result(grism_info, result1, subset['wave'], scalefactor, outfile2a, outfile2b)
    fig = plot_results(subset, result1, func2fit, xfine, plot_label, grism_info, show_initial_fit=show_initial_fit, scalefactor=scalefactor, units=units)
    pp.savefig(bbox_inches='tight', pad_inches=0.1, figsize=figsize)
       
    if row.tweak_wav :  # If input file requests a second fit, w line centroids allowed to very:
        #print "   ******************************************************"
        print("   METHOD 2 fit: Redshift fixed, allowing individual line centroids to move, to compensate for relative wavelength uncertanty in grism.")
        (guesses, parnames) = prep_params(which_grism=which_grism, waveoff=True, scale=row.scale_guess)   # Make container to hold the parameters
        func2fit = pick_fitting_function(which_grism, waveoff=True)
        mymodel = lmfit.Model(func2fit, independent_vars=('wave',), param_names=parnames,  which_grism=which_grism)  # Set up a model
        mypars = set_params(mymodel, parnames, guesses, zz, morph_broad=guess_morphbroad)  # Set initial parameters for that model
        locked_params = lock_params(mypars, which_grism=which_grism, which_gal=which_gal, waveoff=True, sigoff=sigoff)
        locked_params['zz'].vary=False  # FIX THE REDSHIFT
        result2  = mymodel.fit(subset['flam_contsub_scaled'], locked_params, wave=subset['wave'], verbose=True)  # fitting is done here
        if print_fitreports : print(result2.fit_report())
        f.write(result2.fit_report())
        df2 = convert_lmfitresults_2df(result2, parnames, sigoff, grism_info, restwaves, linenames)
        df2.to_csv('/tmp/foo', float_format='%.3f', index=False)
        jrr.util.put_header_on_file('/tmp/foo', header + supplemental_header(result2), outfile3)
        if check4zero_errorbars(result2, parnames) : print("####### WARNING: errorbars were zero! ####### ")
        xfine = finergrid_result(grism_info, result2, subset['wave'], scalefactor, outfile4a, outfile4b)
        plot_label = subdir + " " + specfile + " fit 2"
        fig = plot_results(subset, result2, func2fit, xfine, plot_label, grism_info, show_initial_fit=show_initial_fit, scalefactor=scalefactor, units=units)
        pp.savefig(bbox_inches='tight', pad_inches=0.1)
    f.close()



pp.close()
