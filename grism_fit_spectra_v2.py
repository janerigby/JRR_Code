# This is a script to fit the 1D G102 and G141 WFC/HST grism spectra for S1723 and S2340.
# It used to be in IDL (fit_g141.pro and fit_G102.pro), but I am porting it to Python,
# and at present, using LMFIT as the fitting package (though it's slow).
# I have also changed the fitting strategy -- instead of independently varying the
# observed wavelengths of all the lines, a common redshift is instead varied.
# ** This code assumes the continuum has already been removed.**
# Here's the old IDL readme:
#; Fitting the G141 spectrum of S1723.  Methodology adopted from fits
#; to RCS0327 knots E and U, that appeared in Wuyts et al. 2014, and in
#; turn are based on fitting in Rigby et al. 2011.  
# In the python version, I'm trying a custom fitting function, because there aren't as many
# free parameters as the old IDL code assumed.  For example, the old IDL code let the
# central wavelength of each line vary, but in fact, there's only 1 variable, the global redshift.
#; jrigby, feb 2012, march 2012, oct 2016, may 2018

import numpy as np
import pandas
import jrr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import lmfit 
from os.path import expanduser, basename
import glob
import re
import time


def get_filenames(dir, which_grism) :
    myfiles = [ basename(x) for x in glob.glob(dir + "*" + which_grism + "*.txt") ]
    return(myfiles)


def error_unknown_grism(which_grism) :
    raise Exception('ERROR: Cannot understand which_grism', which_grism)
    return(0)

def get_grism_info(which_grism) :
    # R is FWHM, for a pt source.  For extended source will be morphologically broadened
    # wave_unc is relative wavelength uncertainty, in Angstroms, from wfc3 data handbook, S9.4.8
    if   which_grism == 'G141' :   grism_info = {'R':130., 'x1':11000., 'x2':16600., 'wave_unc':20.}
    elif which_grism == 'G102' :   grism_info = {'R':210., 'x1': 8000., 'x2':11500., 'wave_unc':10. }
    else : error_unknown_grism(which_grism)
    return(grism_info)

def fitfunc_G141(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8) : # A custom G141 fitting function
    # Parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, its clearer to have params be FLUXES.
    # This version fits a global redshift zz, and allows no wavelength miscalibration between lines.
    fluxes          = np.array((f0,        f1,        f2,        f3,         f4,          f5,          f6,       f7,   f8))
    restwaves   = np.array((4862.683,  4960.295,  5008.240,  5877.59,       6549.85,  6564.61,  6585.28, 6725,  7137.8))     
    linenames   =         ['Hbeta',   '[O~III]', '[O~III]',  'He~I',     '[N~II]', 'Halpha', '[N~II]', '[S~II]', '[Ar~III]']
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G141'))

def fitfunc_G141_waveoff(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, d0, d1, d2, d3, d4, d5, d6, d7, d8) : # A custom G141 fitting function
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, its clearer to have params be FLUXES.
    # This version lets each line have a dwave offset, to deal w relative wavelength calibration problems in grism spectra
    fluxes          = np.array((f0,        f1,        f2,        f3,      f4,      f5,      f6,       f7,   f8))
    dwaves       = np.array((d0, d1, d2, d3, d4, d5, d6, d7, d8))
    restwaves   = np.array((4862.683,  4960.295,  5008.240,  5877.59,  6549.85,  6564.61,  6585.28, 6725.0,  7137.8))     
    linenames   =         ['Hbeta',   '[O~III]', '[O~III]',  'He~I',  '[N~II]', 'Halpha', '[N~II]', '[S~II]','[Ar~III]']
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G141', dwaves=dwaves))

def fitfunc_G102(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13) : # A custom G102 fitting function
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, it is much easier to have params be FLUXES.  Easier to report, too
    fluxes    = np.array((f0,         f1,        f2,        f3,       f4,       f5,        f6,       f7,       f8,       f9,       f10,      f11,      f12,       f13))
    restwaves = np.array((3727.092,   3729.875,   3869.86,  3890.151, 3971.195, 3968.593,  4102.892, 4341.684, 4364.436, 4687.02,  4712.58,  4714.490, 4741.4489, 4862.683))
    linenames =         ['[O~II]',    '[O~II]',  '[Ne~III]', 'H8+HeI', 'Heps',  '[Ne~III]', 'Hdelt',  'Hgam',  '[O~III]', 'He~II', '[Ar~IV]', 'HeI',   '[Ar~IV]', 'Hbeta']
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G102'))

def fitfunc_G102_waveoff(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13) : 
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, it is much easier to have params be FLUXES.  Easier to report, too
    fluxes    = np.array((f0,         f1,        f2,        f3,       f4,       f5,        f6,       f7,       f8,       f9,       f10,      f11,      f12,       f13))
    dwaves       = np.array((d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13))
    restwaves = np.array((3727.092,   3729.875,   3869.86,  3890.151, 3971.195, 3968.593,  4102.892, 4341.684, 4364.436, 4687.02,  4712.58,  4714.490, 4741.4489, 4862.683))
    linenames =         ['[O~II]',    '[O~II]',  '[Ne~III]', 'H8+HeI', 'Heps',  '[Ne~III]', 'Hdelt',  'Hgam',  '[O~III]', 'He~II', '[Ar~IV]', 'HeI',   '[Ar~IV]', 'Hbeta']
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G102', dwaves=dwaves))

def fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism, dwaves=()) :  # This part is the same for either grism
    if len(dwaves) == 0 :  dwaves = np.zeros_like(restwaves)
    grism_info = get_grism_info(which_grism) 
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

def prep_params(which_grism, scale_guesses=1.0, waveoff=False, scale=1.) : # Create param container for fitfunc_G141 or fitfunc_G102
    if   which_grism == 'G141': 
        guesses =  np.array((             250.,  500., 1500., 25.,  16.,  800.,  50.,  80.,   3.))
        if waveoff :
            parnames = ('zz', 'morph_broad', 'f0',  'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'd0', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8') #fitfunc_G141_waveoff
            guesses  = np.concatenate((guesses, np.zeros(9))) * scale  # initialize the dwavelengths
        else:
            parnames = ('zz', 'morph_broad', 'f0',  'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8')  # Par names from fitfunc_G141

    elif which_grism == 'G102' :
        guesses = np.array((              144,  200, 100,   40,   20,   30,   50,  120,   10,   10,    5,    5,     20,   250))       
        if waveoff : 
            parnames = ('zz', 'morph_broad', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11', 'f12', 'f13', 'd0', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'd10', 'd11', 'd12', 'd13') #fitfunc_G102_waveoff
            guesses  = np.concatenate((guesses, np.zeros(14))) * scale  # initialize the dwavelengths
        else : 
            parnames = ('zz', 'morph_broad', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11', 'f12', 'f13')  # Par names from fitfunc_G102
    else : error_unknown_grism(which_grism)
    return(guesses * scale_guesses, parnames)

def set_params(mymodel, parnames, guesses, zz, morph_broad=1.) :  # Set initial values of params for fitfunc_G141 or fitfunc_G102
    mypars = mymodel.make_params()
    mypars.add('zz', value=zz, vary=True)
    mypars.add('morph_broad', value=morph_broad, vary=True)
    for ii, parname in enumerate(parnames[2:]) :
        mypars.add(parname, value=guesses[ii], vary=True)
    return(mypars)

def lock_params(mypars, which_grism, waveoff=False, sigoff=0.0, S1723=False) :
    grism_info = get_grism_info(which_grism)
    if waveoff :   # need to set bounds on d0...dX, to sigoff*wave uncert
        dwave_keys = [x for x in mypars.keys() if re.match('d', x)]  # find all the d0... dX wavelength offset pars
        for dwave in dwave_keys :
            mypars[dwave].set(min= -1. * sigoff * grism_info['wave_unc'])
            mypars[dwave].set(max=       sigoff * grism_info['wave_unc'])
    if which_grism == 'G141' :
        mypars['f1'].set(expr='f2 * 0.33189')           #  Lock [O III] 4959, 5007 ratio to Storey & Zeippen2000
        if S1723 :  
            mypars['f6'].set(expr='f5 * 0.062')         #  Lock NII/Ha ratio to what was observed for GNIRS
            mypars['f4'].set(expr='f5 * 0.062/3.071')   #  Lock [N II]  6549, 6585 ratio
        else :  # Any other galaxy
            mypars['f4'].set(expr='f6 / 3.071')         #  Lock [N II]  6549, 6585 ratio
    elif which_grism == 'G102' :
        mypars['f5'].set(expr='f2 * 0.316555')          # Lock [Ne III] flux ratio to Storey & Zeippen2000
        if S1723 :  mypars['f1'].set(expr='f0 * 1.3699') # Lock [O III] flux ratio to what was measured in ESI
    else : error_unknown_grism(which_grism)
    return(mypars)

######################################################################



##  Setup 
label = 'SDSSJ1723+3411'
datadir = expanduser("~") + '/Dropbox/Grism_S1723/Grizli_redux/1Dsum/Wcont/'
scalefactor = 1E17 # Scale everything by scalefactor, to avoid numerical weirdness in LMFIT 
guess_morphbroad = 1.5
zz = 1.331366
show_initial_fit = False

# Lets try G141
which_grism = 'G141'
#which_grism = 'G102'
filenames = get_filenames(datadir, which_grism)

#filenames = ('sgas1723_bothrolls_G102_wcontMWdr.txt',)

''' Status report:
After labor day: I was able to implement the new method, where the line centroids are only allowed to move with a global
redshift.  It's nifty, and may work well for other kinds of spectra.  However, for the grism, the relative wavelength 
calibration problems are leading to obvious shifts in the lines, particularly in Hbeta at the red end of G102.  
It's clear that I need to go back to the logic of the old IDL mpfit routines, where the wavelength of each line
is allowed to move by N sigma.   Now, last time I tried this, I couldn't get the fitter to work.  Since then, 
I've gotten use to lmfit, so perhaps I was just doing something dumb.  Let's save changes, and do better tomorrow.
This is really slow in LMFIT, esp if bounds are narrow, but it does converge.  Tomorrow, print out the output, 
and analyze.  Is there any pattern in the wavelength tweaks?  What about if they're bounded or unbounded?  Why is d6 huge? 
Why do the lines still look mis-aligned?    Should I switch to astropy for fitting, since it may be faster, or keep using this?
God, I want to be done w this project...

Other issues and things to do:
#     - what's the extra line at 4481A?
#     - Split HeI/H8? or are they really at the same Wave?
# - Check units of output.  Should be fluxes in erg/s/cm^2
# - put crosses at the wavelngths of each individual line
# - make plots for publication
# - in G102, enforce flux positivity in weak likes bleuward of HBeta?
# - Deal with scaling?  Michael has added multiple positions, may need to scale. **
# - Dump results to files.
# - Document what I've done.

'''

header = "# Fitting HST grism spectra with grism_fit_spectra_v2.py\n"

grism_info = get_grism_info(which_grism) 
for specfile in filenames[0:1] :   # KLUDGE! REMOVE THIS, just fitting 1 spectrum now...
    outfile = re.sub('.txt', '.fit', specfile)
    f = open(outfile, 'w')  # Over-ride the IDL script to fit the continuum
    f.write(header)
    f.write("# FILENAME "+ specfile + "\n")
    f.write("# GRISM " + which_grism + "\n")
    f.write("# TIME_FIT " + time.strftime("%c") + "\n")
    sp  = pandas.read_csv(datadir + specfile, comment="#")      ## Read the grism spectrum file
    sp['flam_contsub_scaled'] = (sp['flam'] - sp['cont']) * scalefactor
    sp['flam_u_scaled'] = sp['flam_u'] * scalefactor
    sp['weight'] = 1.0 / (sp['flam_u_scaled'])**2   # inverse variance weights
    subset = sp.loc[sp['wave'].between(grism_info['x1'], grism_info['x2'])] # subset of spectr, clean wavelength range.
    if 'both' in specfile : scale=1.    # Guesses were done for the bothrolls sum.
    else:                   scale=0.5   # Just adjusting the guesses to be close to right.

    print "INITIAL FIT, to determine the best-fit redshift. Wavelengths fixed."
    (guesses, parnames) = prep_params(which_grism=which_grism, scale=scale)   # Make container to hold the parameters
    func2fit = pick_fitting_function(which_grism, waveoff=False)
    mymodel = lmfit.Model(func2fit, independent_vars=('wave',), param_names=parnames,  which_grism=which_grism)  # Set up a model
    mypars = set_params(mymodel, parnames, guesses, zz, morph_broad=guess_morphbroad)  # Set initial parameters for that model
    locked_params = lock_params(mypars, which_grism=which_grism, S1723=True)
    result1  = mymodel.fit(subset['flam_contsub_scaled'], locked_params, wave=subset['wave'])  # fitting is done here
    print result1.fit_report()
    zz_fit = result1.best_values['zz']

    print " ******************************************************"
    print "SECOND FIT, for good.  Redshift fixed, allow individual line centroids to move, b/c grism wavelength calib is terrible."
    (guesses, parnames) = prep_params(which_grism=which_grism, waveoff=True, scale=scale)   # Make container to hold the parameters
    func2fit = pick_fitting_function(which_grism, waveoff=True)
    mymodel = lmfit.Model(func2fit, independent_vars=('wave',), param_names=parnames,  which_grism=which_grism)  # Set up a model
    mypars = set_params(mymodel, parnames, guesses, zz, morph_broad=guess_morphbroad)  # Set initial parameters for that model
    locked_params = lock_params(mypars, which_grism=which_grism, S1723=True, waveoff=True, sigoff=3.)
    locked_params['zz'].vary=False  # FIX THE REDSHIFT
    #print locked_params
    result2  = mymodel.fit(subset['flam_contsub_scaled'], locked_params, wave=subset['wave'], verbose=True)  # fitting is done here
    print result2.fit_report()
    f.write(result2.params.pretty_print())
    #f.write(result2.fit_report())
    print " ******************************************************"

    print "Done fitting", specfile,  ". Now plotting..."
    ax = subset.plot(x='wave', y='flam_contsub_scaled', color='black', linestyle='steps-mid')
    subset.plot(x='wave', y='flam_u_scaled', color='grey', ax=ax)
    if show_initial_fit : plt.plot(subset['wave'], result2.init_fit, color='orange', label='init fit')
#    plt.plot(subset['wave'], result2.best_fit, color='blue', label='best fit')
    xfine = np.linspace(grism_info['x1'], grism_info['x2'], 1000)
    plt.plot(xfine, func2fit(xfine, **result2.values), color='blue', label='best fit')
    plt.xlabel("observed wavelength (A)")
    plt.ylabel(r"Scaled, continuum-subtracted  $f_\lambda$")
    plt.title(label + "  " + which_grism, position=(0.4, 0.9))
    plt.xlim(grism_info['x1'], grism_info['x2'])
    plt.legend()
    ax2 = ax.twiny()
    ax2.set_xlim( grism_info['x1']/(1. + result2.best_values['zz']), grism_info['x2']/(1. + result2.best_values['zz']))
    ax2.set_xlabel(r"rest-frame vacuum wavelength ($\AA$)")
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))    
    plt.show()
    f.close()



# The uncertainties in the fits are hard to find.  They're here:  result2.params['f1'].stderr
