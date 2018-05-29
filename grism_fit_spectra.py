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
import lmfit 
from os.path import expanduser, basename

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

def fitfunc_G141(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11) : # A custom G141 fitting function
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, it is much easier to have params be FLUXES.  Easier to report, too
    fluxes          = np.array((f0,        f1,        f2,        f3,         f4,          f5,          f6,       f7,       f8,       f9,       f10,       f11))
    restwaves   = np.array((4862.683,  4960.295,  5008.240,  5877.59,    6302.04064,  6313.8086,   6549.85,  6564.61,  6585.28,  6718.29,  6732.674,  7137.8))
    linenames   =         ['Hbeta',   '[O~III]', '[O~III]',  'He~I',    '[O~I]',     '[S~III]',   '[N~II]', 'Halpha', '[N~II]', '[S~II]', '[S~II]',  '[Ar~III]']
    return ( fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G141') )

def fitfunc_G102(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13) : # A custom G102 fitting function
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, it is much easier to have params be FLUXES.  Easier to report, too
    fluxes    = np.array((f0,         f1,        f2,        f3,       f4,       f5,        f6,       f7,       f8,       f9,       f10,      f11,      f12,       f13))
    restwaves = np.array((3727.092,   3729.875,   3869.86,  3890.151, 3971.195, 3968.593,  4102.892, 4341.684, 4364.436, 4687.02,  4712.58,  4714.490, 4741.4489, 4862.683))
    linenames =         ['[O~II]',    '[O~II]',  '[Ne~III]', 'H8+HeI', 'Heps',  '[Ne~III]', 'Hdelt',  'Hgam',  '[O~III]', 'He~II', '[Ar~IV]', 'HeI',   '[Ar~IV]', 'Hbeta']
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G102'))

def fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism) :  # This part is the same for either grism
    grism_info = get_grism_info(which_grism) 
    sigmas      = restwaves * (1+zz) / ( grism_info['R'] *2.35482/morph_broad)
    y = np.zeros_like(wave)
    for ii, restwave in enumerate(restwaves) :
        aa = fluxes[ii] / sigmas[ii] / np.sqrt(2.*np.pi)  #  Convert from gaussian flux to gaussian height
        y += jrr.spec.onegaus(wave, aa, restwaves[ii]*(1+zz), sigmas[ii], 0.0)
    return(y)

def pick_fitting_function(which_grism) :
    if   which_grism == 'G141': return(fitfunc_G141)
    elif which_grism == 'G102': return(fitfunc_G102)
    else : error_unknown_grism(which_grism)

def prep_params(which_grism, scale_guesses=1.0) : # Create param container for fitfunc_G141 or fitfunc_G102
    if   which_grism == 'G141': 
        guesses =  np.array((250., 500.,  1000.,  25.,   3.,    3.,    30.,  500.,    100.,      1.,     1.,     1.))
        parnames = ('zz', 'morph_broad', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11')  # Par names from fitfunc_G141
    elif which_grism == 'G102' :
        guesses = np.array((144,  200,    100,    40,     20,   30,   50,   120,      10,    10,    5,     5,  20,   250))       
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

def lock_params(mypars, which_grism, S1723=False) :
    if which_grism == 'G141' :
        mypars['f1'].set(expr='f2 * 0.33189')           #  Lock [O III] 4959, 5007 ratio to Storey & Zeippen2000
        if S1723 :  
            mypars['f8'].set(expr='f7 * 0.062')         #  Lock NII/Ha ratio to what was observed for GNIRS
            mypars['f6'].set(expr='f7 * 0.062/3.071')   #  Lock [N II]  6549, 6585 ratio
        else :  # Any other galaxy
            mypars['f6'].set(expr='f8 / 3.071')         #  Lock [N II]  6549, 6585 ratio
    elif which_grism == 'G102' :
        mypars['f5'].set(expr='f2 * 0.316555')          # Lock [Ne III] flux ratio to Storey & Zeippen2000
        if S1723 :  mypars['f1'].set(expr='f0 * 1.3699') # Lock [O III] flux ratio to what was measured in ESI
    else : error_unknown_grism(which_grism)
    return(mypars)

######################################################################



##  Setup 
zz = 1.33136
label = 'SDSSJ1723+3411'
datadir = expanduser("~") + '/Dropbox/Grism_S1723/Grizli_redux/1Dsum/Wcont/'
scalefactor = 1E17 # Scale everything by scalefactor, to avoid numerical weirdness in LMFIT 

# Lets try G141
which_grism = 'G141'
specfile = 'sgas1723_bothrolls_G141_wcontMWdr.txt'

which_grism = 'G102'
specfile = 'sgas1723_bothrolls_G102_wcontMWdr.txt'

guess_morphbroad = 1.5

## Read the grism spectrum file
grism_info = get_grism_info(which_grism) 
sp  = pandas.read_csv(datadir + specfile, comment="#")
sp['flam_contsub_scaled'] = (sp['flam'] - sp['cont']) * scalefactor
sp['flam_u_scaled'] = sp['flam_u'] * scalefactor
sp['weight'] = 1.0 / (sp['flam_u_scaled'])**2   # inverse variance weights
subset = sp.loc[sp['wave'].between(grism_info['x1'], grism_info['x2'])] # subset of spectr, clean wavelength range.

##  Fit w a custom function fitfunc_G141 
(guesses, parnames) = prep_params(which_grism=which_grism)   # Make container to hold the parameters
func2fit = pick_fitting_function(which_grism)
mymodel = lmfit.Model(func2fit, independent_vars=('wave',), param_names=parnames,  which_grism=which_grism)  # Set up a model
mypars = set_params(mymodel, parnames, guesses, zz, morph_broad=guess_morphbroad)  # Set initial parameters for that model
locked_params = lock_params(mypars, which_grism=which_grism, S1723=True)
result  = mymodel.fit(subset['flam_contsub_scaled'], locked_params, wave=subset['wave'])  # fitting is done here
print result.fit_report()

# OK, these test cases are working for both G102 and G141.  Swell
# Cool!  Lots of progress here.  Next tasks:
# - Fine-tune the fitting, esp for G102.  Extra lines?
# - What's wrong w Hbeta fit in G102?
#     - what's the extra line at 4481A?
#     - Split HeI/H8? or are they really at the same Wave?
# - Check units of output.  Should be fluxes in erg/s/cm^2
# - Dump results to files.
# - Document what I've done.


print "Done fitting.  Now plotting..."
ax = subset.plot(x='wave', y='flam_contsub_scaled', color='black')
subset.plot(x='wave', y='flam_u', color='grey', ax=ax)
plt.plot(subset['wave'], result.init_fit, color='orange', label='init fit')
plt.plot(subset['wave'], result.best_fit, color='blue', label='best fit')
plt.legend()
plt.show()

