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
import jrr

###### I'm gonna make me some functions ################################

def translate_metallicity(baseZ) :
    Zinsolar = {'001': 0.05, '004': 0.2 , '008': 0.4, '02': 1.0, '04': 2.0}
    return( Zinsolar[baseZ] )  # test this
    
def normbymed(series) :
    return(series / series.median())

def norminrange(df, wave1, wave2, wavecol='wave', normcol='fnu') :  # Normalize by median within a wavelength range
    norm1 = df.loc[df[wavecol].between(wave1,wave2)][normcol].median()
    return(norm1)

def read_JC_S99file(modeldir, S99file) :  # this is packaged as explained in sb99_readme.txt
    tab = Table.read(modeldir + S99file)
    print "Cols are", tab.colnames
    #for thiscol in tab.colnames :
    #     print type(tab[thiscol].data) ,  tab[thiscol].data.shape
    return(tab)

def write_JC_S99file(tab, outdir, outfile) :  # Output file should be packaged just like John expects.
    tab.write(outdir + outfile, format='fits', overwrite=True)
    return(0)

def read_loresS99file(modeldir, baseZ) : # Read the low-res S99 .spectrum1 file
    infile = modeldir + "low_res_" + str(baseZ) + ".spectrum1"
    print "DEBUGGING, lores file should be", infile
    colnames = ('time', 'wave_Ang', 'logflam_tot', 'logflam_stellar', 'logflam_neb')
    df = pandas.read_table(infile, comment="#", skiprows=6, names=colnames, delim_whitespace=True)
    df['fnu_tot']     = jrr.spec.flam2fnu(df['wave_Ang'], 10**df['logflam_tot'])
    df['fnu_stellar'] = jrr.spec.flam2fnu(df['wave_Ang'], 10**df['logflam_stellar'])
    df['fnu_neb']     = jrr.spec.flam2fnu(df['wave_Ang'], 10**df['logflam_neb'])
    return(df)

def grab_age_from_lores(lores_df, age) :
    return lores_df.loc[lores_df['time'] == float(age)*1E6]   # return a subset w that age

def name_that_cloudy_file(Z, age) :
    return("Z" + Z + "_" + str(age) + "Myr")

def process_S99_spectrumfile(tab, Z, age, age_index) :
    table_outfile = name_that_cloudy_file(Z, age) + '.sed' 
    wave = tab['WAVE'].data[0,]
    flam = tab['FLUX'].data[0, age_index, :]
    fnu  = jrr.spec.flam2fnu(wave, flam)  # convert from flam to fnu, which is what Cloudy needs
    df = pandas.DataFrame({'wave' : wave, 'fnu' : fnu, 'flam' : flam})  # For convenience, make a pandas dataframe
    return(wave, flam, fnu, df)

def write_cloudy_continuum(tab, Z, age, age_index) :  # NOT USED AT PRESENT.  Using low-res but larger wavelength range spectrum
    # Write a continuum file (wave, fnu), formatted for Cloudy "Table SED" command.
    tempfile = "temp.txt"
    (temp_wave, temp_flam, temp_fnu, tempdf) = process_S99_spectrumfile(tab, Z, age, age_index)
    np.savetxt(tempfile, np.transpose([temp_wave,  temp_fnu]), "%.7f  %.4E", comments="")
    header_text = "#wave_angs  fnu\n199.000     1.0E-99  units Angstrom\n"  # See Hazy1, Table SED command.  Kludge to state units
    jrr.util.put_header_on_file(tempfile, header_text, table_outfile)
    return(0)

def write_cloudy_infile(Z, age, logU, cloudy_template) :
    solarZ = translate_metallicity(Z)
    contfile = name_that_cloudy_file(Z, age)
    infile   = name_that_cloudy_file(Z, age) + '.in'
    outfile  = name_that_cloudy_file(Z, age) + '.out'
    jrr.util.replace_text_in_file("SEDGOESHERE", contfile, cloudy_template, infile)
    jrr.util.replace_text_in_file("AGEGOESHERE", str(age), infile)
    jrr.util.replace_text_in_file("METALLICITYWRTSOLARGOESHERE", str(solarZ), infile) # Z in solar units
    jrr.util.replace_text_in_file("ABSOLUTEMETALLICITYGOESHERE", str(Z), infile) # Z in absolute units (solar is 0.02)
    jrr.util.replace_text_in_file("LOGUGOESHERE", str(logU), infile) # ionization parameter
    return(infile, outfile)

def retrieve_cloudy_nebcont(Z, age, directory) :
    nebcontfile = directory + name_that_cloudy_file(Z, age) + '.con'
    print "DEBUGGING, nebcontfile is", nebcontfile
    #tempneb = jrr.util.strip_pound_before_colnames(nebcontfile)
    colnames = ("wave_um", "incident", "trans", "outward", "transnet", "reflected", "total", "reflectedline", "outwardline", "dum")
    df = pandas.read_table(nebcontfile, comment="#", names=colnames, delim_whitespace=True)
    df['nu'] = jrr.spec.wave2freq(df['wave_um'])
    df['wave_Ang'] =  df['wave_um'] * 1E4 # "nu" col is wave in micron, not nu, as requested in cloudy infile. convert from um to Ang
    df['nebcont']  =  (df['outward'] - df['outwardline'])  # col4 - col9, see Hazy 16.44.8 "Save continuum" also S99/neb_cont_allmodels.gnu
    # Units of file should be nu*Jnu, according to Hazy.  Must divide columns by nu to get units in fnu
    df['fnu_incident'] = df['incident'] / df['nu'] 
    df['fnu_nebcont']  = df['nebcont']  / df['nu'] 
    return(df)

def add_nebular_to_stellar_deprecated(stellar_df, cloudy_df) :  # deprecated bc this adds artifacts to the nebular cont.
    # fnu(stellar + neb)_hires  =  [ fnu(stellar)_lores + fnu(neb)_lores ] / fnu(stellar)_lores     * fnu(stellar)_hires
    cloudy_df['multby4neb'] =  (cloudy_df['fnu_incident'] + cloudy_df['fnu_nebcont'])  / cloudy_df['fnu_incident']
    stellar_df['multby4neb'] = jrr.spec.rebin_spec_new(cloudy_df['wave_Ang'], cloudy_df['multby4neb'], stellar_df['wave'])
    stellar_df['fnu_wneb'] = stellar_df['fnu'] * stellar_df['multby4neb']
    stellar_df['fnu_nebonly'] = stellar_df['fnu'] * stellar_df['multby4neb'] -  stellar_df['fnu'] 
    return(0)  # acts on stellar_df

def add_nebular_to_stellar(stellar_df, cloudy_df) :
    # Cloudy calculates the nebular continuum for the lores S99 spectrum.  I want to add it to the hires S99 spectrum.
    # Cloudy applies a weird scaling factor, so I can't just add them.  I had done a scaling short-cut, but that
    # introduces artifical abs lines in the nebular spectrum, bc of mismatch between highres, lores stellar spectra.
    # New method: measure the multiplicative offset between lores, highres spectrum, use that to scale the nebular,
    # and then add it directly.
    nwave1 = 1760. ;  nwave2 = 1790.
    norm1 =   cloudy_df.loc[cloudy_df['wave_Ang'].between(nwave1, nwave2)]['fnu_incident'].median()
    norm2 = stellar_df.loc[stellar_df['wave'].between(  nwave1, nwave2)]['fnu'].median()
    print "DEBUGGING, norms are", norm1, norm2
    stellar_df['fnu_nebonly'] =  norm2 / norm1 * jrr.spec.rebin_spec_new(cloudy_df['wave_Ang'], cloudy_df['fnu_nebcont'], stellar_df['wave'])
    stellar_df['fnu_wneb'] = stellar_df['fnu'] + stellar_df['fnu_nebonly']
    return(0) # acts on stellar_df

    
def get_cloudydir(logU, style='JC_S99') :   # later, should add style='BPASS'
    return(style + "_logU" + str(logU) + "/")

def parse_filename(filename) :
    # format will be Z001_10Myr.con or Z001_10Myr.fits
    myre = re.compile('Z(\S+)_(\S+)Myr.(\S+)')  
    (baseZ, baseage, extension) = myre.match(filename).groups()
    return(baseZ, baseage, extension)
    
def make_label(baseZ, thisage, logU, style='JC_S99'):
    label = "Z" + str(baseZ) + ", " + str(thisage) + " Myr, logU=" + str(logU)
    return(label)

########################################################################

homedir   = expanduser("~")
modeldir  = homedir + '/Dropbox/MagE_atlas/Contrib/S99/models/'
nebcontdir = '/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Cloudy_models/Nebular_continuum/'  # Where cloudy model results live
# Here's what I expect in J. Chisholm's packaged S99 files. However, getting age, Z from files themselves.
# ages = (1, 2, 3, 4, 5, 8, 10, 15, 20, 40) # Myr
# metallicities = (0.01, 0.2, 0.4, 1.0, 2.0) # fraction of solar
logUs = (-2.0, -2.3, -2.5, -2.7, -3.0)  # John C. wants to see impact of neb cont on different logUs
style = ('JC_S99')  #, 'BPASS')  # need to add BPASS

# Formatting for cloudy
cloudy_exe = "/Volumes/Apps_and_Docs/JRR_Utils/c17.00/source/cloudy.exe"
cloudy_template = "../S99_cloudy_template.in"  # Cloudy input template
cloudy_script = 'run_cloudy.sh'  # Script to run all the Cloudy models

# Realized that high-res S99 output is inappropriate for Cloudy, bc it has no flux <900A.  Therefore, have compiled the
# lowres S99 spectra from JC into Cloudy format, and am now running that in Cloudy.
# Once Cloudy finishes, I need to read each continuum spectrum, add it to its corresponding S99 highres spectrum,
# and pack it up like JC expects.  Then repeat for BPASS.

#  *** SWITCHES ************
prep_cloudy =    False     # Prepare all the Cloudy input files?
analyze_cloudy = True      # Analyze the cloudy outputs, to grab nebular continuua?
sanity_plots   = True
#  *************************



if prep_cloudy :  # Part 1: THIS SECTION PREPARES THE CLOUDY BATCH MODE SCRIPT
    for logU in logUs:
        cloudydir = get_cloudydir(logU, style=style)
        if not exists(nebcontdir + cloudydir):   makedirs(nebcontdir + cloudydir)
        print "DEBUGGING, working on logU=", logU, "in dir", cloudydir
        f = open(nebcontdir + cloudydir + cloudy_script, 'w')
        chdir(nebcontdir + cloudydir)
        filenames = [ basename(x) for x in glob.glob(modeldir + '*fits')]      # For each of J. Chisholm's packaged S99 files:
        for thisfile in filenames:
            baseZ = re.sub('.fits', '', re.sub('sb99', '', thisfile))
            hires_tab = read_JC_S99file(modeldir, thisfile)
            for age_index, thisage in enumerate(hires_tab['AGE'][0,].data) :
                thisage = "{0:g}".format(np.float(thisage) / 1E6)  # convert to Myr, and drop the trailing .0
                (infile, outfile) = write_cloudy_infile(baseZ, thisage, logU, cloudy_template)  # Write Cloudy infile (ported from wrapper.py)
                f.write(cloudy_exe + " < " + infile + " > " + outfile + "\n")
        f.close()  # Close the cloudy_script
        print "Done this dir", cloudydir
    chdir(nebcontdir)

####  Part 2: NOW, GO RUN the Cloudy models in each dir, as:   xjobs -j 7 -s run_cloudy.sh


#### Part 3: Cloudy models run, and ready to be processed?  Good.  Pair each Cloudy .con continuum outfile with
##### its corresponding high-res S99 output spectra, interpolate as needed, add the two spectra, and repackage for JC.
if analyze_cloudy :  
    plt.close("all")
    for logU in logUs:
        cloudydir = get_cloudydir(logU, style=style)
        the_pdf = "Sanity_checks_S99_stellar_continuua_for_Cloudy_logU" + str(logU) + ".pdf"
        pp = PdfPages(the_pdf)  # output   
        chdir(nebcontdir + cloudydir)
        filenames = [ basename(x) for x in glob.glob(modeldir + '*fits')]  
        for thisfile in filenames:  # REAL USE THIS
            baseZ = re.sub('.fits', '', re.sub('sb99', '', thisfile))
            lores_dfall = read_loresS99file(modeldir, baseZ)  # Grab the lores S99 file that was fed into Cloudy        
            hires_tab = read_JC_S99file(modeldir, thisfile)
            hires_tab['flam_wneb']    = hires_tab['FLUX'] *0.0  # Create this new column w the right size.
            hires_tab['flam_nebonly'] = hires_tab['FLUX'] *0.0  
            for age_index, thisage in enumerate(hires_tab['AGE'][0,].data) :  
                thisage = "{0:g}".format(np.float(thisage) / 1E6)  # convert to Myr, and drop the trailing .0
                (wave, flam, fnu, stellar_df) = process_S99_spectrumfile(hires_tab, baseZ, thisage, age_index)
                cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)   # Grab the Cloudy output file w nebular continuum
                lores_df = grab_age_from_lores(lores_dfall, thisage)
                add_nebular_to_stellar(stellar_df, cloudy_df)
                hires_tab['flam_wneb'][0, age_index, :] = np.array(jrr.spec.fnu2flam(stellar_df['wave'], stellar_df['fnu_wneb']))
                hires_tab['flam_nebonly'][0, age_index, :] = np.array(jrr.spec.fnu2flam(stellar_df['wave'], stellar_df['fnu_nebonly']))
                # above step repackages stellar_df into fits files for john. 
                
                nwave1 = 2000. ;  nwave2 = 2020.  # Sanity check plots
                norm1 =   cloudy_df.loc[cloudy_df['wave_Ang'].between(nwave1, nwave2)]['fnu_incident'].median()
                norm2 = stellar_df.loc[stellar_df['wave'].between(  nwave1, nwave2)]['fnu'].median()
                norm3 =   lores_df.loc[lores_df['wave_Ang'].between(nwave1, nwave2)]['fnu_stellar'].median()
                #print "Norms are", norm1, norm2, norm3
                plt.plot(stellar_df['wave'], stellar_df['fnu'] /norm2, label='S99 hires spectrum')
                plt.plot(cloudy_df['wave_Ang'], cloudy_df['fnu_incident']/ norm1, label="Cloudy incident", color='k', linewidth=2)
                plt.plot(stellar_df['wave'], stellar_df['fnu_nebonly'] /norm2, label='scaled neb for hires', color='pink', lw=3)
                plt.plot(cloudy_df['wave_Ang'], cloudy_df['fnu_nebcont'] / norm1, label="Cloudy nebular continuum", color='red', linewidth=0.7)
                plt.plot(lores_df['wave_Ang'],   lores_df['fnu_stellar']   / norm3, label='S99 lores', linewidth=0.7)
                #plt.plot(lores_df['wave_Ang'],   lores_df['fnu_neb']       / norm3, label='S99-predicted nebcont') # this is zero. not sure why
                plt.xlim(800,3000)
                plt.ylim(0,2)
                plt.title("Z" + baseZ + ", age " + thisage + "Myr, logU="+str(logU))
                plt.legend()
                pp.savefig()
                plt.clf()
                # Check that I can read the newly-generated file
                suffix = "_wnebcont_logU" + str(logU) + ".fits" 
                outfile = re.sub('.fits', suffix, thisfile)
            write_JC_S99file(hires_tab, modeldir + "Wnebcont/", outfile) 
            newtab = read_JC_S99file(modeldir  + "Wnebcont/", outfile)
            (nwave, nflam, nfnu, nstellar_df) = process_S99_spectrumfile(newtab, baseZ, thisage, age_index)
        pp.close()
        chdir(nebcontdir)


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
        plt.plot(cloudy_df['wave_Ang'], cloudy_df['fnu_nebcont'], label=baseZ)
    plt.xlim(900,1E4)
    plt.ylim(0, 8E-17)
    plt.legend()
    plt.xlabel("wavelength (Angstrom)")  ; plt.ylabel("fnu nebular continuum")
    plt.title("Neb cont intensity should scale w Z.  For age=" + str(thisage) + " Myr, logU=" + str(logU))
    pp.savefig()
    plt.clf()
    for thisfile in filenames: 
        baseZ = re.sub('.fits', '', re.sub('sb99', '', thisfile))
        thisage = '2' # Myr  kludge
        cloudy_df = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir + cloudydir)
        subset = cloudy_df.loc[cloudy_df['wave_Ang'].between(1200,3000)]
        norm1 = norminrange(subset, 2000., 2020., wavecol='wave_Ang', normcol='fnu_nebcont')
        plt.plot(cloudy_df['wave_Ang'], cloudy_df['fnu_nebcont']/norm1, label=baseZ)
    plt.xlim(900,1E4)
    plt.ylim(0, 5)
    plt.legend()
    plt.xlabel("wavelength (Angstrom)")  ; plt.ylabel("fnu nebular continuum, normalized at 2000A")
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
            plt.plot(cloudy_df['wave_Ang'], cloudy_df['fnu_nebcont']/cloudy_df['fnu_incident'], label=label)
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
            plt.plot(cloudy_df['wave_Ang'], cloudy_df['fnu_nebcont']/cloudy_df['fnu_incident'], label=label, color=colors[ii])
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
            plt.plot(cloudy_df['wave_Ang'], cloudy_df['fnu_nebcont']/cloudy_df['fnu_incident'], label=label, color=colors[ii])
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
            norm1 = norminrange(cloudy_df, 2000., 2020., wavecol='wave_Ang', normcol='fnu_nebcont')
            plt.plot(cloudy_df['wave_Ang'], cloudy_df['fnu_nebcont']/norm1, label="logU="+str(logU), color=colors[ii])
    plt.xlabel("wavelength (Angstrom)")  ; plt.ylabel("fnu nebular continuum, normalized at 2000A")
    plt.xlim(900,1E4)
    plt.ylim(0,5)
    plt.title("How nebcont shape vary w logU, for a given age (2 Myr), all Z")
    pp.savefig()
    pp.close()

