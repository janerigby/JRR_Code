# This is a script to fit the 1D G102 and G141 WFC/HST grism spectra for S1723 and S2340.
# It used to be in IDL (fit_g141.pro and fit_G102.pro), but I am porting it to Python.
# At the moment, the fitting tool I am using is LMFIT.  It seems slooow.
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
    if   which_grism == 'G141' :   grism_info = {'R':130., 'x1':10840., 'x2':16600., 'wave_unc':20.}
    elif which_grism == 'G102' :   grism_info = {'R':210., 'x1': 8000., 'x2':11500., 'wave_unc':10. }
    else : error_unknown_grism(which_grism)
    return(grism_info)

def get_morphological_broadening_factor() :
    return(1.3)  # DUMMY, fill this in later! ********

def fitfunc_G141_small(wave, zz, morph_broad, A0, A1, A2, A3) : #  This is probably obsolete.  Delete when no longer needed.
    theAs       = np.array((A0, A1, A2, A3))
    restwaves   = np.array((4862.683,  4960.295,  5008.240,   6564.61))
    linenames   =  ['Hbeta',   '[O~III]', '[O~III]',  'Halpha']
    grism_info = get_grism_info(which_grism) 
    sigmas      = restwaves * (1+zz) / ( grism_info['R'] *2.35482/morph_broad)
    y = np.zeros_like(wave)
    for ii, restwave in enumerate(restwaves) :
        y += jrr.spec.onegaus(wave, theAs[ii], restwaves[ii]*(1+zz), sigmas[ii], 0.0)
    return(y)  

def fitfunc_G141(wave, zz, morph_broad, A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11) : # A custom G141 fitting function
    theAs       = np.array((A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11))
    restwaves   = np.array((4862.683,  4960.295,  5008.240,  5877.59,    6302.04064,  6313.8086,   6549.85,  6564.61,  6585.28,  6718.29,  6732.674,  7137.8))
    linenames   =         ['Hbeta',   '[O~III]', '[O~III]',  'He~I',    '[O~I]',     '[S~III]',   '[N~II]', 'Halpha', '[N~II]', '[S~II]', '[S~II]',  '[Ar~III]']
    grism_info = get_grism_info(which_grism) 
    sigmas      = restwaves * (1+zz) / ( grism_info['R'] *2.35482/morph_broad)
    y = np.zeros_like(wave)
    for ii, restwave in enumerate(restwaves) :
        y += jrr.spec.onegaus(wave, theAs[ii], restwaves[ii]*(1+zz), sigmas[ii], 0.0)
    return(y)  

def prep_params_G141(scale_guesses=1.0) : # Create param container for fitfunc_G141()
    guesses =  np.array((2.,       5.0,       11.0,       0.3,         0.05,         0.05,         0.03,     5.,       0.1,      0.1,     0.5,        0.2))
    parnames = ('zz', 'morph_broad', 'A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11')  # Par names from fitfunc_G141
    return(guesses * scale_guesses, parnames)

def set_params_G141(mymodel, parnames, guesses) :  # Set initial values of params for fitfunc_G141
    mypars = mymodel.make_params()
    mypars.add('zz', value=zz, vary=True)
    mypars.add('morph_broad', value=1., vary=True)
    for ii, parname in enumerate(parnames[2:]) :
        mypars.add(parname, value=guesses[ii], vary=True)
    # Constraints could go here.
    return(mypars)

def get_emlines(which_grism, zz) :
    if   which_grism == 'G141' :
        #             g0          g1        g2         g3          g4           g5           g6        g7        g8        g9        g10        g11        
        restwave   = [4862.683,  4960.295,  5008.240,  5877.59,    6302.04064,  6313.8086,   6549.85,  6564.61,  6585.28,  6718.29,  6732.674,  7137.8]
        linename =  ['Hbeta',   '[O~III]', '[O~III]',  'He~I',    '[O~I]',     '[S~III]',   '[N~II]', 'Halpha', '[N~II]', '[S~II]', '[S~II]',  '[Ar~III]']
        height_guess =[2.,       5.0,       11.0,       0,         0.0,         0.0,         0.00,     5.,       0.1,      0.0,     0.5,        0.0]

#    if   which_grism == 'G141' :  # TEMP DEBUGGING TEMP DEBUGGING  ****
#        #             g0          g1        g2           g6        g7        g8      
#        restwave   = [4862.683,  4960.295,  5008.240,    6549.85,  6564.61,  6585.28]
#        linename =  ['Hbeta',   '[O~III]', '[O~III]',   '[N~II]', 'Halpha', '[N~II]']
#        height_guess =[2.,       5.0,       11.,         0.03,     5.,       0.1]

        
    elif which_grism == 'G102':
        #            g0         g1          g2         g3         g4         g5          g6         g7         g8          g9        g10         g8          g10          g11       
        restwave  = [3727.092,  3729.875,   3869.86,   3890.151,  3971.195,  3968.593,   4102.892,  4341.684,  4364.436,   4687.02,  4712.58,   4714.490,   4741.4489,   4862.683]
        linenames = ['[O~II]', '[O~II]',   '[Ne~III]', 'H8+HeI', 'Heps',    '[Ne~III]', 'Hdelt',   'Hgam',    '[O~III]',  'He~II',   '[Ar~IV]', 'HeI',      '[Ar~IV]',  'Hbeta'] 
        height_guess = [0.6,    0.6,        0.4,        0.3,      0.1,         0.1,      0.15,      0.35,      0.1,        0.02,      0.02,      0.02,       0.02,        0.77]   
    else : error_unknown_grism(which_grism)
    emlines = pandas.DataFrame({'restwave' : restwave, 'linename' : linename, 'height_guess' :height_guess})  ## Conveniently package emission lines into dataframe
    emlines['obswave'] = emlines['restwave'] * (1.0 + zz)
    grism_info = get_grism_info(which_grism)               # Get resolution, wavelength limits
    morph_broad =  get_morphological_broadening_factor()   # Decrease R by morphological broadening factor. Should be run iteratively, minimize chisq.
    emlines['sigma']   = emlines['obswave'] / ( grism_info['R'] *2.35482/morph_broad) # Expected Gaussian sigma of each line.  Hold constant.
    emlines['amp_guess'] = emlines['height_guess'] * emlines['sigma'] * np.sqrt(2. * np.pi)
    return(emlines)


def lock_lines_together(composite_par, S1723=None)  :  # Tie some amplitudes together by theoretical ratios or other observations
    if which_grism == 'G141' :
        composite_par['g1_amplitude'].set(expr='g2_amplitude * 0.33189')           #  Lock [O III] 4959, 5007 ratio to Storey & Zeippen2000
        if S1723 :
            composite_par['g8_amplitude'].set(expr='g7_amplitude * 0.062')         #  Lock NII/Ha ratio to what was observed for GNIRS
            composite_par['g6_amplitude'].set(expr='g7_amplitude * 0.062/3.071')   #  Lock [N II]  6549, 6585 ratio
        else :  # Any other galaxy
            composite_par['g8_amplitude'].set(expr='g6_amplitude / 3.071')         #  Lock [N II]  6549, 6585 ratio

    if which_grism == 'G102' :
        composite_par['g5_amplitude'].set(expr='g2_amplitude * 0.316555')          # Lock [Ne III] amplitudes by ratio from Storey & Zeippen2000
        if S1723 :  composite_par['g1_amplitude'].set(expr='g0_amplitude * 1.3699') # Lock [O III] ratio to what was measured in ESI
    return(composite_par)  # Do I need to return this, or is it globally readable?  

######################################################################

# This code assumes the continuum has already been removed

## LOTS OF SETUP.
which_grism = 'G141'
datadir = expanduser("~") + '/Dropbox/Grism_S1723/Grizli_redux/1Dsum/Wcont/'
specfile = 'sgas1723_bothrolls_G141_wcontMWdr.txt'
zz = 1.32952   # from GNIRS
zz = 1.3325076 
label = 'SDSSJ1723+3411'
sigoff = 1.0 # *** PLACEHOLDER, REPLACE LATER w something more motivated
scalefactor = 1E17 # Scale everything by scalefactor, to avoid numerical weirdness


emlines = get_emlines(which_grism, zz)  # Read emission lines to fit into data frame
grism_info = get_grism_info(which_grism)               # Get resolution, wavelength limits

## Read the grism spectrum file
sp  = pandas.read_csv(datadir + specfile, comment="#")
sp['flam_contsub_scaled'] = (sp['flam'] - sp['cont']) * scalefactor
sp['flam_u_scaled'] = sp['flam_u'] * scalefactor
sp['weight'] = 1.0 / (sp['flam_u_scaled'])**2   # inverse variance weights
subset = sp.loc[sp['wave'].between(grism_info['x1'], grism_info['x2'])] # subset of spectr, clean wavelength range.

## Set up the fitting
mod = {}  # A dictionary of models
par = {}  # A dictionary of parameters for those models
for row in emlines.itertuples():
    print "DEBUGGING, working on line", row.Index, row.linename, row.restwave, " AA"
    prefix = 'g' + str(row.Index) + '_'
    tempmod = lmfit.models.GaussianModel(prefix=prefix)  # prefix is handy!
    temppar = tempmod.make_params()
    temppar.add(prefix+'sigma',  value=row.sigma, vary=False)  # Fix the linewidths
    temppar.add(prefix+'center', value=row.obswave, vary=True, min=(row.obswave - sigoff*grism_info['wave_unc']), max=(row.obswave + sigoff*grism_info['wave_unc']))
    temppar.add(prefix+'amplitude',  value=row.amp_guess, vary=True)  
#    temppar.add(prefix+'amplitude',  value=10., vary=True)    # TEMP KLUDGE!**** THis shouldn't work but it does w value=1... why???***
    mod[prefix] = tempmod
    par[prefix] = temppar
print "Loaded list of lines to fit"

# Combine all the mods and all the pars into one Composite.   This is super-clunky, may need to refine
for ii, thiskey in enumerate(mod.keys()) :
    if ii == 0 :
        compmod = mod[thiskey]
        comppar = par[thiskey]
    else :
        compmod += mod[thiskey]
        comppar += par[thiskey]


'''
# NEXT LINE IS TEMPORARILY COMMENTED OUT< DEBUGGING!
#locked_comppar = lock_lines_together(comppar, S1723=True)  # Lock together amplitudes of lines with fixed ratios, or w measurements from ESI, GNIRS (for S1723)
#BELOW IS TEMP! Should be handled in lock_lines_together
comppar['g1_amplitude'].set(expr='g2_amplitude * 0.33189')
print "Fitting now..."
## This runs extremely slow!!! Trying alt aproach below.
result = compmod.fit(subset['flam_contsub_scaled'], comppar, x=subset['wave'], weights=subset['weight'])  # Fitting is done here.
print(result.fit_report(min_correl=0.25))
'''

## Trying to fit w custom function fitfunc_G141  *** STOPPED HERE, in progress
#guesses = (2., 5., 11., 5.)  # for the smaller fitfunc_G141.  Delete when obsolete
#parnames =('zz', 'morph_broad', 'A0', 'A1', 'A2', 'A3')
(guesses, parnames) = prep_params_G141()   # Make container to hold the parameters
mymodel = lmfit.Model(fitfunc_G141, independent_vars=('wave',), param_names=parnames)  # Set up a model
mypars = set_params_G141(mymodel, parnames, guesses) # Set initial parameters for that model
#locked_pars = lock_params_G141(mypars, S1723=True)
result  = mymodel.fit(subset['flam_contsub_scaled'], mypars, wave=subset['wave'])  # FITTING DONE HERE
print result.fit_report()
# Cool.  This has all the lines, but it runs fast (1s).  Now, add the constraints. *** Paused here.



print "Done fitting.  Now plotting..."
ax = subset.plot(x='wave', y='flam_contsub_scaled', color='black')
subset.plot(x='wave', y='flam_u', color='grey', ax=ax)
plt.plot(subset['wave'], result.init_fit, color='orange', label='init fit')
plt.plot(subset['wave'], result.best_fit, color='blue', label='best fit')
plt.legend()
plt.show()

