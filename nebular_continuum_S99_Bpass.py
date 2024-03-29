from __future__ import print_function
from builtins import str
from os.path import basename, expanduser, exists
from os import chdir, makedirs
from shutil import copyfile
import subprocess
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas
from astropy.table import Table
from numpy import log10, round
import jrr

# Note: run from this dir:  /Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Cloudy_models/Nebular_continuum/

###### I'm gonna make me some functions ################################

def translate_Z_abs2solar(absoluteZ) :
    Zinsolar = {'001': 0.05, '004': 0.2 , '008': 0.4, '02': 1.0, '04': 2.0, '020': 1.0, '040': 2.0}
    # Above is repetitive, bc S99 metallicity format first, then just those that have a different BPASS metallicity syntax
    return( Zinsolar[absoluteZ] )

def normbymed(series) :
    return(series / series.median())

def norminrange(df, wave1, wave2, wavecol='wave', normcol='flam') :  # Normalize by median within a wavelength range
    norm1 = df.loc[df[wavecol].between(wave1,wave2)][normcol].median()
    return(norm1)

def read_JC_S99file(modeldir, S99file) :  # this is packaged as explained in sb99_readme.txt
    tab = Table.read(modeldir + S99file)
    #print "Cols are", tab.colnames
    #for thiscol in tab.colnames :
    #     print type(tab[thiscol].data) ,  tab[thiscol].data.shape
    return(tab)

def make_tab_header(tab, IMF,  style, baseZ, thisage, logU)  :
    tab.meta['IMF']    = IMF
    tab.meta['style']  = style
    tab.meta['absZ']   = baseZ
    tab.meta['ageMyr'] = thisage
    tab.meta['logU']   = logU
    return(0) # acts on global tab
    
def write_JC_fitsfile(tab, outdir, outfile) :  # Output file should be packaged just like John expects.
    tab.write(outdir + outfile, format='fits', overwrite=True)
    return(0)

def read_loresS99file(modeldir, baseZ) : # Read the low-res S99 .spectrum1 file
    infile = modeldir + "low_res_" + str(baseZ) + ".spectrum1"
    colnames = ('time', 'wave_Ang', 'logflam_tot', 'logflam_stellar', 'logflam_neb')
    df = pandas.read_table(infile, comment="#", skiprows=6, names=colnames, delim_whitespace=True)
    df['flam_tot']     = 10**df['logflam_tot']
    df['flam_stellar'] = 10**df['logflam_stellar']
    df['flam_neb']     = 10**df['logflam_neb']
    return(df)

def grab_age_from_lores(lores_df, age) :
    return lores_df.loc[lores_df['time'] == float(age)*1E6]   # return a subset w that age

def process_S99_spectrumfile(tab, Z, age, age_index) :  # These are fits files in which JC has packaged the highres S99 spectra
    table_outfile = name_that_cloudy_file(Z, age) + '.sed' 
    wave = tab['WAVE'].data[0,]
    flam = tab['FLUX'].data[0, age_index, :]
    fnu  = jrr.spec.flam2fnu(wave, flam)  # convert from flam to fnu, which is what Cloudy needs
    df = pandas.DataFrame({'wave' : wave, 'fnu' : fnu, 'flam' : flam})  # For convenience, make a pandas dataframe
    return(wave, flam, fnu, df)

def name_that_cloudy_file(Z, age) :
    return("Z" + Z + "_" + str(age) + "Myr")

def write_cloudy_infile(Z, age, logU, cloudy_template, style, ver='2.1') :  
    solarZ =  translate_Z_abs2solar(Z)
    if   style == 'JC_S99'  :   model = name_that_cloudy_file(Z, age)
    elif 'BPASS' in style   :   model = jrr.bpass.default_filenames(ver=BPASS_ver, style=style)[1]
    else : raise Exception("Unrecognized style")
    prefix = name_that_cloudy_file(Z, age)
    infile = prefix + '.in'  ;  outfile = prefix + '.out'
    logZ_absolute = round(log10(0.02 * solarZ), 2) # log(absolute Z), where 0.02 is solar.
    jrr.util.replace_text_in_file("MODELGOESHERE", model, cloudy_template, infile) 
    jrr.util.replace_text_in_file("AGEGOESHERE", str(age), infile)
    jrr.util.replace_text_in_file("PREFIXGOESHERE", prefix, infile)
    jrr.util.replace_text_in_file("METALLICITYWRTSOLARGOESHERE", str(solarZ), infile) # Z in solar units
    jrr.util.replace_text_in_file("ABSOLUTEMETALLICITYGOESHERE", str(Z), infile) # Z in absolute units (solar is 0.02)
    jrr.util.replace_text_in_file("ABSOLUTELOGMETALLICITYGOESHERE", str(logZ_absolute), infile) # Z in absolute units (solar is 0.02)
    jrr.util.replace_text_in_file("LOGUGOESHERE", str(logU), infile) # ionization parameter
    return(infile, outfile)

def retrieve_cloudy_nebcont(Z, age, directory) :
    nebcontfile = directory + name_that_cloudy_file(Z, age) + '.con'
    colnames = ("wave_um", "incident", "trans", "outward", "transnet", "reflected", "total", "reflectedline", "outwardline", "dum")
    df = pandas.read_table(nebcontfile, comment="#", names=colnames, delim_whitespace=True)
    df['nu'] = jrr.spec.wave2freq(df['wave_um'])
    df['wave_Ang'] =  df['wave_um'] * 1E4 # "nu" col is wave in micron, not nu, as requested in cloudy infile. convert from um to Ang
    df['nebcont']  =  (df['outward'] - df['outwardline'])  # col4 - col9, see Hazy 16.44.8 "Save continuum" also S99/neb_cont_allmodels.gnu
    # Units of file should be nu*Jnu, according to Hazy.  Must divide columns by nu to get units in fnu.  Then convert to flam b/c that's what JC wants.
    df['flam_incident'] = jrr.spec.fnu2flam(df['wave_Ang'], (df['incident'] / df['nu']) )
    df['flam_nebcont']  = jrr.spec.fnu2flam(df['wave_Ang'], (df['nebcont']  / df['nu']) )
    df['flam_nebline']  = jrr.spec.fnu2flam(df['wave_Ang'], (df['outwardline']  / df['nu']) )  # Adding line emission
    return(df)

def add_nebular_to_stellar(stellar_df, cloudy_df) :
    # Cloudy calculates the nebular continuum, and outputs a lowres neb cont spectrum. Need to add it to the hires spectrum.
    # Cloudy applies a weird scaling factor, so I can't just add them.  So, measure the multiplicative offset between lores & highres
    # spectrum, use that to scale the nebular, and then add it directly.  Update: adding line emission, in case we want it later.
    nwave1 = 1760. ;  nwave2 = 1790.
    norm1 =   cloudy_df.loc[cloudy_df['wave_Ang'].between(nwave1, nwave2)]['flam_incident'].median()
    norm2 = stellar_df.loc[stellar_df['wave'].between(  nwave1, nwave2)]['flam'].median()
    stellar_df['flam_nebcont'] =  norm2 / norm1 * jrr.spec.rebin_spec_new(cloudy_df['wave_Ang'], cloudy_df['flam_nebcont'], stellar_df['wave']) #nebular continuum
    stellar_df['flam_nebline'] =  norm2 / norm1 * jrr.spec.rebin_spec_new(cloudy_df['wave_Ang'], cloudy_df['flam_nebline'], stellar_df['wave']) #nebular line emissin
    stellar_df['flam_stellar_nebcont'] = stellar_df['flam'] + stellar_df['flam_nebcont']
    return(0) # acts on stellar_df 

def plot_stellar_nebular(style, stellar_df, cloudy_df, Z, age, logU, lores_df=None, nwave1=1760., nwave2=1790.) :
    # Plot the input spectrum (from S99 or BPASS) against the Cloudy outputs)
    norm1 =  cloudy_df.loc[cloudy_df['wave_Ang'].between(nwave1, nwave2)]['flam_incident'].median()
    norm2 = stellar_df.loc[stellar_df['wave'].between(  nwave1, nwave2)]['flam'].median()
    plt.plot(stellar_df['wave'], stellar_df['flam'] /norm2, label='High-res stellar spectrum')
    plt.plot(cloudy_df['wave_Ang'], cloudy_df['flam_incident']/ norm1, label="Cloudy incident", color='k', linewidth=2)
    plt.plot(stellar_df['wave'], stellar_df['flam_nebcont'] /norm2, label='scaled neb for hires', color='pink', lw=3)
    plt.plot(cloudy_df['wave_Ang'], cloudy_df['flam_nebcont'] / norm1, label="Cloudy nebular continuum", color='red', linewidth=0.7)
    if lores_df is not None :
        norm3 =   lores_df.loc[lores_df['wave_Ang'].between(nwave1, nwave2)]['flam_stellar'].median()
        plt.plot(lores_df['wave_Ang'],   lores_df['flam_stellar']   / norm3, label='S99 lores', linewidth=0.7)
    plt.xlim(800,3000)
    plt.ylim(0,6)
    plt.title(style + " Z" + str(baseZ) + ", age " + str(thisage) + "Myr, logU="+str(logU))
    plt.legend()
    pp.savefig()
    plt.clf()
    return(0)

def get_cloudydir(logU, style='JC_S99') :   
    return(style + "_logU" + str(logU) + "/")

def parse_filename(filename) :
    # format will be Z001_10Myr.con or Z001_10Myr.fits
    myre = re.compile('Z(\S+)_(\S+)Myr.(\S+)')  
    (baseZ, baseage, extension) = myre.match(filename).groups()
    return(baseZ, baseage, extension)
    
def make_label(baseZ, thisage, logU, style='JC_S99'):
    label = style + " Z" + str(baseZ) + ", " + str(thisage) + " Myr, logU=" + str(logU)
    return(label)


########################################################################
#BPASS_ver = '2.1'  # Ran everything for 2.1. 
#BPASS_ver = '2.2'  # But, 2.2 now available. 
BPASS_ver = '2.2.1' # Patch on 2.2, to solve negative fluxes
homedir   = expanduser("~")
modeldir  = homedir + '/Dropbox/MagE_atlas/Contrib/S99/models/'
outnebdir = homedir + '/Dropbox/S99_Bpass2_wnebcont/'
nebcontdir = '/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Cloudy_models/Nebular_continuum/'  # Where cloudy model results live
# Here's what I expect in J. Chisholm's packaged S99 files. Getting age, Z from S99 files themselves. Using same ages for BPASS.
ages = (1, 2, 3, 4, 5, 8, 10, 15, 20, 40) # Myr
absZs = ('001', '004', '008', '020', '040')  # Metallicities as BPASS expects 
logUs = (-2.0, -2.3, -2.5, -2.7, -3.0)  # John C. wants to see impact of neb cont on different logUs. This is a big range
styles = ('JC_S99', 'BPASS_binary', 'BPASS_single')

# Formatting for cloudy
cloudy_exe = "/Volumes/Apps_and_Docs/JRR_Utils/c17.00/source/cloudy.exe"
#cloudy_exe = "/Volumes/Apps_and_Docs/JRR_Utils/c17.01/source/cloudy.exe"
cloudy_script = 'run_cloudy.sh'  # Script to run all the Cloudy models
cloudy_template = { 'JC_S99': '../S99_cloudy_template.in' , 'BPASS_binary': '../bpass_cloudy_template.in', 'BPASS_single': '../bpass_cloudy_template.in'}

#  *** SWITCHES ************
prep_cloudy =    False   # Prepare all the Cloudy input files?
analyze_cloudy = True      # Analyze the cloudy outputs, to grab nebular continuua?
sanity_plots   = True
#  *************************


# Realized that high-res S99 output is inappropriate for Cloudy, bc it has no flux <900A.  Therefore, have compiled the
# lowres S99 spectra from JC into Cloudy format, and am now running that in Cloudy.
if prep_cloudy :  # Part 1: THIS SECTION PREPARES THE CLOUDY BATCH MODE SCRIPT
    for style in styles :
        for logU in logUs:
            cloudydir = get_cloudydir(logU, style=style)
            if not exists(nebcontdir + cloudydir):   makedirs(nebcontdir + cloudydir)
            print("Working in dir", cloudydir)
            f = open(nebcontdir + cloudydir + cloudy_script, 'w')
            chdir(nebcontdir + cloudydir)
            if style == 'JC_S99' :  # S99 specific file wrangling
                filenames = [ basename(x) for x in glob.glob(modeldir + '*fits')]      # For each of J. Chisholm's packaged S99 files:
                for thisfile in filenames:
                    baseZ = re.sub('.fits', '', re.sub('sb99', '', thisfile))
                    hires_tab = read_JC_S99file(modeldir, thisfile)
                    for age_index, thisage in enumerate(hires_tab['AGE'][0,].data) :
                        thisage = "{0:g}".format(np.float(thisage) / 1E6)  # convert to Myr, and drop the trailing .0
                        (infile, outfile) = write_cloudy_infile(baseZ, thisage, logU, cloudy_template[style], style)  # Write Cloudy infile
                        f.write(cloudy_exe + " < " + infile + " > " + outfile + "\n")
                f.close()  # Close the cloudy_script
            elif('BPASS' in style) :  # Make the Cloudy model using BPASS as ionizing spectrum
                for bpassZ in absZs :
                    for thisage in ages :
                        (infile, outfile) = write_cloudy_infile(bpassZ, thisage, logU, cloudy_template[style], style)
                        f.write(cloudy_exe + " < " + infile + " > " + outfile + "\n")
                f.close()  # Close the cloudy_script
            else : raise Exception("Unrecognized style")
        chdir(nebcontdir)


####  Part 2: NOW, GO RUN the Cloudy models in each dir, as:
####      xjobs -j 7 -s run_cloudy.sh


#### Part 3: Cloudy models run, and ready to be processed?  Good.  Pair each Cloudy .con continuum outfile with
##### its corresponding high-res S99 output spectra, interpolate as needed, add the two spectra, and repackage for JC.


if analyze_cloudy :  
    plt.close("all")
    for style in styles :  # Adding BPASS compatitibility now...
        for logU in logUs:
            cloudydir = get_cloudydir(logU, style=style)
            the_pdf = "Sanity_checks_" + style + "_stellar_neb_continuua_logU" + str(logU) + ".pdf"
            pp = PdfPages(the_pdf)  # output  
            chdir(nebcontdir + cloudydir)
            if style == 'JC_S99' :  # S99 specific file wrangling. Read JC's infiles, check that lowres, hires stellar spectra sortof match; check nebcont scaled right
                filenames = [ basename(x) for x in glob.glob(modeldir + '*fits')]  
                for thisfile in filenames:  
                    baseZ = re.sub('.fits', '', re.sub('sb99', '', thisfile))
                    lores_dfall = read_loresS99file(modeldir, baseZ)  # Grab the lores S99 file that was fed into Cloudy        
                    hires_tab = read_JC_S99file(modeldir, thisfile)
                    hires_tab['flam_stellar_nebcont']    = hires_tab['FLUX'] *0.0  # Create this new column w the right size.
                    hires_tab['flam_nebcont'] = hires_tab['FLUX'] *0.0  
                    hires_tab['flam_nebline'] = hires_tab['FLUX'] *0.0  
                    for age_index, thisage in enumerate(hires_tab['AGE'][0,].data) : 
                        thisage = "{0:g}".format(np.float(thisage) / 1E6)  # convert to Myr, and drop the trailing .0
                        (wave, flam, fnu, stellar_df) = process_S99_spectrumfile(hires_tab, baseZ, thisage,age_index)
                        print("Working on: (Z age logU style cloudydir)", baseZ, thisage, logU, style, cloudydir)
                        cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)   # Grab the Cloudy output file w nebular continuum
                        lores_df = grab_age_from_lores(lores_dfall, thisage)
                        add_nebular_to_stellar(stellar_df, cloudy_df)
                        hires_tab['flam_stellar_nebcont'][0, age_index, :] = np.array(stellar_df['flam_stellar_nebcont'])
                        hires_tab['flam_nebline'][0, age_index, :]         = np.array(stellar_df['flam_nebline'])
                        hires_tab['flam_nebcont'][0, age_index, :]         = np.array(stellar_df['flam_nebcont'])
                        # above step repackages stellar_df into fits files for J. Chisholm
                        plot_stellar_nebular(style, stellar_df, cloudy_df, baseZ, thisage, logU, lores_df=lores_df)
                        suffix = "_wnebcont_logU" + str(logU) + ".fits" 
                        outfile = re.sub('.fits', suffix, thisfile)
                    make_tab_header(hires_tab, style,  style, baseZ, str(ages), logU)
                    write_JC_fitsfile(hires_tab, outnebdir, outfile) 
                    newtab = read_JC_S99file(outnebdir, outfile)      # Check that I can read the newly-generated file
                    (nwave, nflam, nfnu, nstellar_df) = process_S99_spectrumfile(newtab, baseZ, thisage, age_index)
            elif('BPASS' in style) :   # Input spectrum is BPASS. 
                for baseZ in absZs :
                    for thisage in ages :
                        print("Working on: (Z age logU style cloudydir)", baseZ, thisage, logU, style, cloudydir)
                        cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)   # Grab the Cloudy output file w nebular continuum
                        cloudy_df.head()  # Need more after this; just want to show that I can read the cloudy file.
                        stellar_df = jrr.bpass.wrap_load_1spectrum(baseZ, thisage*1E6, style, ver=BPASS_ver)   # Get high-res BPASS stellar continuum. cloudy .con too low res.
                        add_nebular_to_stellar(stellar_df, cloudy_df)
                        plot_stellar_nebular(style+BPASS_ver, stellar_df, cloudy_df, baseZ, thisage, logU)  #added BPASS version for sanity-checking
                        outfile = style + "_Z" + baseZ + "_" + str(thisage) + "Myr_logU" + str(logU) + ".fits"
                        tab = Table.from_pandas(stellar_df)
                        make_tab_header(tab, jrr.bpass.default_filenames(ver=BPASS_ver, style=style)[1],  style, baseZ, thisage, logU)  # I think this adds metadata
                        write_JC_fitsfile(tab,  outnebdir, outfile)  # Output file as John expects.  Not bothering to concatenate ages.
            pp.close()
            chdir(nebcontdir)


style = "JC_S99" # Have just run sanity checks for S99
if sanity_plots :  # Part 4: make a bunch of sanity plots
    plt.close("all")
    the_pdf = "behavior_nebular_continua.pdf"
    pp = PdfPages(the_pdf)  # output   
    logU = -2.5
    cloudydir = get_cloudydir(logU,  style=style)
    filenames = [ basename(x) for x in glob.glob(modeldir + '*fits')]  
    for thisfile in filenames: 
        baseZ = re.sub('.fits', '', re.sub('sb99', '', thisfile))
        thisage = '2' # Myr  kludge
        cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)
        subset = cloudy_df.loc[cloudy_df['wave_Ang'].between(1200,3000)]
        plt.plot(cloudy_df['wave_Ang'], cloudy_df['flam_nebcont'], label=baseZ)
    plt.xlim(900,1E4)
    plt.ylim(0, 1E-4)
    plt.legend()
    plt.xlabel("wavelength (Angstrom)")  ; plt.ylabel("flam nebular continuum")
    plt.title("Neb cont intensity should scale w Z.  For age=" + str(thisage) + " Myr, logU=" + str(logU))
    pp.savefig()
    plt.clf()
    for thisfile in filenames: 
        baseZ = re.sub('.fits', '', re.sub('sb99', '', thisfile))
        thisage = '2' # Myr  kludge
        cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)
        subset = cloudy_df.loc[cloudy_df['wave_Ang'].between(1200,3000)]
        norm1 = norminrange(subset, 2000., 2020., wavecol='wave_Ang', normcol='flam_nebcont')
        plt.plot(cloudy_df['wave_Ang'], cloudy_df['flam_nebcont']/norm1, label=baseZ)
    plt.xlim(900,1E4)
    plt.ylim(0, 5)
    plt.legend()
    plt.xlabel("wavelength (Angstrom)")  ; plt.ylabel("flam nebular continuum, normalized at 2000A")
    plt.title("Nebcont shape varies w Z, for given age=" + str(thisage) + " Myr, logU=" + str(logU))
    pp.savefig()
    plt.clf()

    for logU in logUs:
        cloudydir = get_cloudydir(logU,  style=style)
        filenames = [ basename(x) for x in glob.glob(nebcontdir + cloudydir + '*2Myr.con')]
        for thisfile in filenames:   # Make a sanity check that nebcont scales w Z
            (baseZ, thisage, extension) = parse_filename(thisfile)
            cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)
            label = make_label(baseZ, thisage, logU, style=style)
            plt.plot(cloudy_df['wave_Ang'], cloudy_df['flam_nebcont']/cloudy_df['flam_incident'], label=label)
        plt.ylabel("nebular / stellar ratio")
        plt.xlim(900,1E4)
        plt.ylim(0,2)
        plt.title("Nebular cont to stellar cont ratio, vs Z, for a given logU")
        plt.legend()
        pp.savefig()
        plt.clf()


# Want these to be just for a given age.
    colors=("blue", "orange", "green", "red", "purple")
    for logU in logUs:
        cloudydir = get_cloudydir(logU,  style=style)
        filenames = [ basename(x) for x in glob.glob(nebcontdir + cloudydir + '*2Myr.con')]
        for ii, thisfile in enumerate(filenames) :
            (baseZ, thisage, extension) = parse_filename(thisfile)
            cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)
            label = make_label(baseZ, thisage, logU, style=style)
            plt.plot(cloudy_df['wave_Ang'], cloudy_df['flam_nebcont']/cloudy_df['flam_incident'], label=label, color=colors[ii])
    plt.ylabel("nebular / stellar ratio")
    plt.xlim(900,1E4)
    plt.ylim(0,2)
    plt.title("Nebular cont to stellar cont, logU=" + str(logUs) + ", age=2Myr")
    pp.savefig()
    plt.clf()

    for logU in logUs[1:-1] :      # repeat for narrower range of logUs
        cloudydir = get_cloudydir(logU, style=style)
        filenames = [ basename(x) for x in glob.glob(nebcontdir + cloudydir + '*2Myr.con')]
        for ii, thisfile in enumerate(filenames) :
            (baseZ, dummyage, extension) = parse_filename(thisfile)
            cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)
            label = make_label(baseZ, thisage, logU, style=style)
            plt.plot(cloudy_df['wave_Ang'], cloudy_df['flam_nebcont']/cloudy_df['flam_incident'], label=label, color=colors[ii])
    plt.ylabel("nebular / stellar ratio")
    plt.xlim(900,1E4)
    plt.ylim(0,2)
    plt.title("Nebular cont to stellar cont, logU=" + str(logUs[1:-1]) + ", age=2Myr")
    pp.savefig()
    plt.clf()

    
    for ii, thisfile in enumerate(filenames) :
        (baseZ, thisage, extension) = parse_filename(thisfile)
        for logU in logUs :
            cloudydir = get_cloudydir(logU, style=style)
            cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)
            norm1 = norminrange(cloudy_df, 2000., 2020., wavecol='wave_Ang', normcol='flam_nebcont')
            plt.plot(cloudy_df['wave_Ang'], cloudy_df['flam_nebcont']/norm1, label="logU="+str(logU), color=colors[ii])
    plt.xlabel("wavelength (Angstrom)")  ; plt.ylabel("flam nebular continuum, normalized at 2000A")
    plt.xlim(900,1E4)
    plt.ylim(0,5)
    plt.title("How nebcont shape vary w logU, for a given age (2 Myr), all Z")
    pp.savefig()
    pp.close()

