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
from matplotlib.backends.backend_pdf import PdfPages
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

 
def get_line_wavelengths(which_grism) :
    if which_grism == "G141" :
        restwaves   = np.array((4862.683,  4960.295,  5008.240,  5877.59,  6310.,            6549.85,  6564.61,  6585.28, 6725.))
        linenames   =         ['Hbeta',   '[O~III]', '[O~III]',  'He~I',  '[O~I]+[S~III]',  '[N~II]', 'Halpha', '[N~II]', '[S~II]']
        #                       0             1            2         3          4             5           6           7         8  
    elif which_grism == "G102" :
        restwaves = np.array((3727.092,   3729.875,   3869.86,  3890.151,  3968.593,   3971.195,   4025.,        4102.892, 4341.684, 4364.436, 4472.7,  4687.02,   4741.4489, 4862.683))
        linenames =         ['[O~II]',    '[O~II]',  '[Ne~III]', 'H8+HeI', '[Ne~III]', 'Heps',     'HeI+HeII',  'Hdelt',  'Hgam',  '[O~III]', 'He~I', 'He~II',    '[Ar~IV]', 'Hbeta']
        #                       0             1            2         3          4         5           6           7          8        9        10        11         12           13
    else : error_unknown_grism(which_grism)
    return(restwaves, linenames)

# SHould have done the below functions with just a single "params" var that gets unpcked, as in v=params.valuesdict().  See https://lmfit.github.io/lmfit-py/parameters.html
# Also, probably should have done all this in astropy, which has similar functionality. Moving on.

def fitfunc_G141(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8) : # A custom G141 fitting function
    # Parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, its clearer to have params be FLUXES.
    # This version fits a global redshift zz, and allows no wavelength miscalibration between lines.
    fluxes          = np.array((f0,        f1,        f2,        f3,         f4,          f5,          f6,       f7,  f8))
    (restwaves, linenames) = get_line_wavelengths("G141")
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G141'))

def fitfunc_G141_waveoff(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, d0, d1, d2, d3, d4, d5, d6, d7, d8) : # A custom G141 fitting function
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, its clearer to have params be FLUXES.
    # This version lets each line have a dwave offset, to deal w relative wavelength calibration problems in grism spectra
    fluxes          = np.array((f0,        f1,        f2,        f3,      f4,      f5,      f6,       f7, f8))
    dwaves       = np.array((d0, d1, d2, d3, d4, d5, d6, d7, d8))
    (restwaves, linenames) = get_line_wavelengths("G141")
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G141', dwaves=dwaves))

def fitfunc_G102(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13) : # A custom G102 fitting function
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, it is much easier to have params be FLUXES.  Easier to report, too
    fluxes    = np.array((f0,         f1,        f2,        f3,       f4,       f5,        f6,       f7,       f8,       f9,       f10,      f11,      f12, f13))
    (restwaves, linenames) = get_line_wavelengths("G102")
    return (fitfunc_eithergrism(wave, zz, morph_broad, restwaves, fluxes, which_grism='G102'))

def fitfunc_G102_waveoff(wave, zz, morph_broad, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13) : 
    # parameters are the FLUX of each Gaussian, = a*c*sqrt(2pi).  Since sigmas are fixed, it is much easier to have params be FLUXES.  Easier to report, too
    fluxes    = np.array((f0,         f1,        f2,        f3,       f4,       f5,        f6,       f7,       f8,       f9,       f10,      f11,      f12, f13))
    dwaves       = np.array((d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13))
    (restwaves, linenames) = get_line_wavelengths("G102")
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
    grism_info = get_grism_info(which_grism)
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
        dwave_keys = [x for x in mypars.keys() if re.match('d', x)]  # find all the d0... dX wavelength offset pars
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
            mypars['f12'].set(value=0., vary=False)      # These lines fall off the edge of the grism
            mypars['f13'].set(value=0., vary=False)
        elif which_gal == 'S1723' :
            mypars['f1'].set(expr='f0*1.3699') # Lock [O III] flux ratio to what was measured in ESI
        if which_gal != 'S1723' :
            mypars['f1'].set(value=0., vary=False)  # Cannot resolve OII doublet, so set second trans to zero.
            
    return(mypars)

def parse_filename(grism_filename) :
    m = re.search('(\S+)_(\S+)_(\S+)_(\S+)_(\S+)', grism_filename)
    mydict = { 'gname': m.group(1), 'descrip': m.group(2), 'roll': m.group(3),  'grating': m.group(4), 'suffix': m.group(5)}
    return(mydict)
    
def plot_results(df_data, fit_result, func2fit, title, grism_info, show_initial_fit=False) :
    print "Done fitting", title,  ". Now plotting..."
    ax = df_data.plot(x='wave', y='flam_contsub_scaled', color='black', linestyle='steps-mid', lw=1.5, legend=False)
    df_data.plot(x='wave', y='flam_u_scaled', color='grey', ax=ax, legend=False)
    if show_initial_fit : plt.plot(df_data['wave'], fit_result.init_fit, color='orange', label='init fit')
#    plt.plot(subset['wave'], fit_result.best_fit, color='blue', label='best fit')
    xfine = np.linspace(grism_info['x1'], grism_info['x2'], 1000)
    plt.plot(xfine, func2fit(xfine, **fit_result.values), color='blue', label='best fit')
    plt.xlabel("observed wavelength (A)")
    plt.ylabel(r"Scaled, continuum-subtracted  $f_\lambda$")
    plt.title(title + "  ", position=(0.5, 0.9), fontsize=10)
    plt.xlim(grism_info['x1'], grism_info['x2'])
    (restwaves, linenames) = get_line_wavelengths(which_grism)
    if 'd0' in fit_result.params.keys() :
        offsets = np.array([fit_result.params[x].value for x in fit_result.params.keys() if re.match('d', x)])
        plt.scatter(restwaves*(1. + fit_result.best_values['zz']) + offsets, np.zeros_like(restwaves)-0.1, marker='+', color='k')
    else :    plt.scatter(restwaves*(1. + fit_result.best_values['zz']), np.zeros_like(restwaves)-0.1, marker='+', color='k')
    ax2 = ax.twiny()
    ax2.set_xlim( grism_info['x1']/(1. + fit_result.best_values['zz']), grism_info['x2']/(1. + fit_result.best_values['zz']))
    ax2.set_xlabel(r"rest-frame vacuum wavelength ($\AA$)")
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))    
    plt.show()

######################################################################


##  Setup 
scalefactor = 1E17 # Scale everything by scalefactor, to avoid numerical weirdness in LMFIT 
guess_morphbroad = 1.5
sigoff=3 # Warn if the delta wavelengths exceed this
show_initial_fit = True
which_grism = 'G102'
which_gal = 'S2340' #'S1723' # 

if which_gal == 'S1723' :
    zz = 1.331366
    datadir = expanduser("~") + '/Dropbox/Grism_S1723/Grizli_redux/1Dsum/Wcont/'
elif which_gal == 'S2340' :
    zz = 1.4275 
    datadir = expanduser("~") + '/Dropbox/Grism_S2340/Grizli_redux/1Dsum/Wcont/'
else : raise Exception("Unrecognized galaxy:", which_gal)
    
filenames = get_filenames(datadir, which_grism)
figsize = (8,4)
pp = PdfPages("grism_fitspectra_"+which_gal+"_"+which_grism+".pdf")

''' Status report:
This works pretty well, It's a two-step approach.  Step 1 moves all the line centroids in unison, 
by a global redshift.    Step 2 fixes the redshift as the result from step1, then re-fits, allowing 
the wavelengths of groups of lines to  move together.  
Other issues and things to do:
# - make plots for publication 
# - Deal with scaling. **  Michael has added multiple positions, may need to scale. **
# - Dump results to files. My custom output file needs to be output
# - Document what I've done.
'''


grism_info = get_grism_info(which_grism) 
for specfile in filenames :
    specfile_dict = parse_filename(specfile) 
    outfile1 = re.sub('.txt', '.fit', specfile)
    f = open(outfile1, 'w')  
    header = "# Fitting HST grism spectra with grism_fit_spectra_v3.py\n# FILENAME "+ specfile + "\n# GRISM " + which_grism
    header += "\n# TIME_FIT " + time.strftime("%c") + "\n# UNITS " + str(1./scalefactor) + "  erg/s/cm^2"
    header += "# Note: For S1723 NII/Ha ratios are set from GNIRS, and [O II] 3727/3729 ratio set from ESI.\n"
    header += "# Note: For S2340, NII fluxes were fixed to zero; what is reported as Ha is actually the blend of [N II]+Ha+[N II].\n"
    header += "# Note on fluxing for S1723/1Dsum: In the bothrolls files, Michael had summed the flux over both rolls.  Here those fluxes\n"
    header += "#                   have been divided by 2, and therefore should match the fluxes of the spectra extracted from a single roll.\n"
    header += "# Note on fluxing for S1723/1Dbyclumps, and all S2340: The fluxing is weird, b/c Michael added multiple clumps from whichever\n"
    header += "# rolls were clean.  So the fluxing is super-weird, and should only be used in a relative sense.  Divide by Hbeta and move on.
    f.write(header)
    sp  = pandas.read_csv(datadir + specfile, comment="#")      ## Read the grism spectrum file
    if which_gal == 'S1723' and 'both' in specfile : bothrollsfactor = 0.5  # Michael has summed the flux in both rolls.  Compensating here
    else : bothrollsfactor = 1.0
    sp['flam_contsub_scaled'] = (sp['flam'] - sp['cont']) * scalefactor * bothrollsfactor
    sp['flam_u_scaled'] = sp['flam_u'] * scalefactor * bothrollsfactor
    sp['weight'] = 1.0 / (sp['flam_u_scaled'])**2   # inverse variance weights
    subset = sp.loc[sp['wave'].between(grism_info['x1'], grism_info['x2'])] # subset of spectr, clean wavelength range.
    if which_gal == 'S1723': 
        if 'both' in specfile : scale=0.5    # Scale just adjuststhe guesses to be close to right.  Not used for anything else.
        else:                   scale=0.4   # 
    else :
        if 'both' in specfile : scale=0.5  
        else:                   scale=0.05

    print "INITIAL FIT, to determine the best-fit redshift. Wavelengths fixed."
    (guesses, parnames) = prep_params(which_grism=which_grism, scale=scale)   # Make container to hold the parameters
    func2fit = pick_fitting_function(which_grism, waveoff=False)
    mymodel = lmfit.Model(func2fit, independent_vars=('wave',), param_names=parnames,  which_grism=which_grism)  # Set up a model
    mypars = set_params(mymodel, parnames, guesses, zz, morph_broad=guess_morphbroad)  # Set initial parameters for that model
    locked_params = lock_params(mypars, which_grism=which_grism, which_gal=which_gal, limit_wave=False)
    result1  = mymodel.fit(subset['flam_contsub_scaled'], locked_params, wave=subset['wave'])  # fitting is done here
    #print result1.fit_report()
    zz_fit = result1.best_values['zz']
    #plot_results(subset, result1, func2fit, specfile + " FIRST FIT", grism_info, show_initial_fit=show_initial_fit)
    #pp.savefig(bbox_inches='tight', pad_inches=0.1)

    print " ******************************************************"
    print "SECOND FIT, for good.  Redshift fixed, allow individual line centroids to move, b/c grism wavelength calib is terrible."
    (guesses, parnames) = prep_params(which_grism=which_grism, waveoff=True, scale=scale)   # Make container to hold the parameters
    func2fit = pick_fitting_function(which_grism, waveoff=True)
    mymodel = lmfit.Model(func2fit, independent_vars=('wave',), param_names=parnames,  which_grism=which_grism)  # Set up a model
    mypars = set_params(mymodel, parnames, guesses, zz, morph_broad=guess_morphbroad)  # Set initial parameters for that model
    locked_params = lock_params(mypars, which_grism=which_grism, which_gal=which_gal, waveoff=True, sigoff=sigoff)
    locked_params['zz'].vary=False  # FIX THE REDSHIFT
    #print locked_params
    result2  = mymodel.fit(subset['flam_contsub_scaled'], locked_params, wave=subset['wave'], verbose=True)  # fitting is done here
    print result2.fit_report()
    f.write(result2.fit_report())
    print " ******************************************************"

    f.close()
    plot_results(subset, result2, func2fit, specfile + " SECOND FIT", grism_info, show_initial_fit=show_initial_fit)
    #pp.savefig(bbox_inches='tight', pad_inches=0.1)

    (restwaves, linenames) = get_line_wavelengths(which_grism)
    print header
    print "# index restwave  linename   flux     dflux   constraint waveoffset   dwaveoffset"
    print "#  --      (A)       --      ("+str(1/scalefactor)+" erg/s/cm^2)   same   expr    (obsA)     (obsA)"

    # Next up: Reformat result2 as a dataframe, and then print it out
    fkeys   = [x for x in parnames if re.match('f', x)]
    dkeys   = [x for x in parnames if re.match('d', x)]
    flux    = [result2.params[key].value  for key in fkeys]
    flux_u  = [result2.params[key].stderr for key in fkeys]
    Angoff  = [result2.params[key].value  for key in dkeys]
    Angoff_u= [result2.params[key].stderr for key in dkeys]
    df = pandas.DataFrame({ 'restwave':restwaves, 'linenames':linenames, 'flux':flux, 'flux_u':flux_u, 'Angoff':Angoff, 'Angoff_u':Angoff_u})
    df = df[['restwave', 'linenames', 'flux', 'flux_u', 'Angoff', 'Angoff_u']] # force reorder
    df.loc[df['Angoff'] > sigoff*grism_info['wave_unc']]
    
    
    for ii, restwave in enumerate(restwaves):
        
        print ii, restwave, linenames[ii], result2.params['f'+str(ii)].value, result2.params['f'+str(ii)].stderr, result2.params['f'+str(ii)].expr,
        print  result2.params['d'+str(ii)].value, result2.params['d'+str(ii)].stderr
        # Do some error checking: ##
        if np.abs(result2.params['d'+str(ii)].value) > sigoff*grism_info['wave_unc'] : print "# WARNING! Wavelength drift was large!  d"+ str(ii)
        if result2.params['f'+str(ii)].stderr    > result2.params['f'+str(ii)].value : print "# WARNING! FLUX HIGHLY UNCERTAIN!  f" + str(ii)

pp.close()
