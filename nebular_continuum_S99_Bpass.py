from os.path import expanduser
from os.path import basename
from os import chdir
from shutil import copyfile
import subprocess
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas
from astropy.table import Table
import jrr

###### I'm gonna make me some functions ################################

def normbymed(series) :
    return(series / series.median())

def read_JC_S99file(modeldir, S99file) :  # this is packaged as explained in sb99_readme.txt
    tab = Table.read(modeldir + S99file)
    print "Cols are", tab.colnames
    #for thiscol in tab.colnames :
    #     print type(tab[thiscol].data) ,  tab[thiscol].data.shape
    return(tab)

def write_JC_S99file(tab, outdir, outfile) :  # Output file should be packaged just like John expects.
    tab.write(outdir + outfile, format='fits', overwrite=True)
    return(0)

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

def write_cloudy_infile(Z, age, cloudy_template) :
    contfile = name_that_cloudy_file(Z, age)
    infile   = name_that_cloudy_file(Z, age) + '.in'
    outfile  = name_that_cloudy_file(Z, age) + '.out'
    jrr.util.replace_text_in_file("SEDGOESHERE", contfile, cloudy_template, infile)
    jrr.util.replace_text_in_file("AGEGOESHERE", str(age), infile)
    jrr.util.replace_text_in_file("METALLICITYGOESHERE", str(Z), infile)
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

def add_stellar_and_nebular(stellar_tab, cloudy_df) :
    foo = jrr.spec.rebin_spec_new(cloudy_df['wave_Ang'], cloudy_df['nebcont'], stellar_tab['WAVE'].data[0,])
    # This is not done yet.  First, need to check units.
    return(foo)
    

########################################################################
homedir   = expanduser("~")
modeldir  = homedir + '/Dropbox/MagE_atlas/Contrib/S99/models/'
nebcontdir = '/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Cloudy_models/Nebular_continuum/JC_S99/'

# THIS SECTION PREPARES THE CLOUDY BATCH MODE SCRIPT

# Here's what I expect in J. Chisholm's packaged S99 files. However, getting age, Z from files themselves.
# ages = (1, 2, 3, 4, 5, 8, 10, 15, 20, 40) # Myr
# metallicities = (0.01, 0.2, 0.4, 1.0, 2.0) # fraction of solar

# Formatting for cloudy
cloudy_exe = "/Volumes/Apps_and_Docs/JRR_Utils/c17.00/source/cloudy.exe"
cloudy_template = "template.in"  # Cloudy input template
cloudy_script = 'run_cloudy.sh'  # Script to run all the Cloudy models
f = open(cloudy_script, 'w')

# Realized that high-res S99 output is inappropriate for Cloudy, bc it has no flux <900A.  Therefore, have compiled the
# lowres S99 spectra from JC into Cloudy format, and am now running that in Cloudy.
# Once Cloudy finishes, I need to read each continuum spectrum, add it to its corresponding S99 highres spectrum,
# and pack it up like JC expects.  Then repeat for BPASS.

prep_cloudy = False      # Prepare all the Cloudy input files?
analyze_cloudy = True    # Analyze the cloudy outputs, to grab nebular continuua?

if prep_cloudy :
    # For each of J. Chisholm's packaged S99 files:
    chdir(nebcontdir)
    filenames = [ basename(x) for x in glob.glob(modeldir + '*fits')]
    for thisfile in filenames:  
        baseZ = re.sub('.fits', '', re.sub('sb99', '', thisfile))
        hires_tab = read_JC_S99file(modeldir, thisfile)
        for age_index, thisage in enumerate(tab['AGE'][0,].data) :
            thisage = "{0:g}".format(np.float(thisage) / 1E6)  # convert to Myr, and drop the trailing .0
            #write_cloudy_continuum(tab, baseZ, thisage, age_index)  # The continuum file that Cloudy expects
            # Above is commented out, bc highres spectrum doesn't have hard ionizing spectrum that Cloudy needs.
            (infile, outfile) = write_cloudy_infile(baseZ, thisage, cloudy_template)  # Write Cloudy infile (ported from wrapper.py)
            f.write(cloudy_exe + " < " + infile + " > " + outfile + "\n")
    f.close()  # Close the cloudy_script

####  NOW, GO RUN the Cloudy models as:   xjobs -j 7 -s run_cloudy.sh

    
#### Cloudy models run, and ready to be processed?  Good.  We need to pair each Cloudy .con continuum outfile with
##### its corresponding high-res S99 output spectra, interpolate as needed, add the two spectra, and repackage for JC.
if analyze_cloudy :
    chdir(nebcontdir)
    filenames = [ basename(x) for x in glob.glob(modeldir + '*fits')]  
    #for thisfile in filenames:  # REAL USE THIS
    for thisfile in filenames[0:1] : # KLUDGE TO ONLY RUN ONCE 
        baseZ = re.sub('.fits', '', re.sub('sb99', '', thisfile))
        hires_tab = read_JC_S99file(modeldir, thisfile)
        #for age_index, thisage in enumerate(hires_tab['AGE'][0,].data) :  # REAL USE THIS
        for age_index, thisage in enumerate(hires_tab['AGE'][0,].data[0:1]) :
            thisage = "{0:g}".format(np.float(thisage) / 1E6)  # convert to Myr, and drop the trailing .0
            clouddf = retrieve_cloudy_nebcont(baseZ, thisage, nebcontdir)   # Grab the Cloudy output file w nebular continuum
            subset_clouddf = clouddf.loc[clouddf['wave_Ang'].between(900,3000)]
            plt.plot(subset_clouddf['wave_Ang'], normbymed(subset_clouddf['incident']), label="Cloudy incident")
            plt.plot(subset_clouddf['wave_Ang'], normbymed(subset_clouddf['nebcont']),  label="Cloudy nebular continuum")
            plt.xlim(800,3000)
            (wave, flam, fnu, stellar_df) = process_S99_spectrumfile(hires_tab, baseZ, thisage, age_index)
            plt.plot(stellar_df['wave'], normbymed(stellar_df['fnu']), label='S99 hires spectrum')
            # These are scaled oddly wrt each other.  Also, are units (fnu) correct?
            plt.legend()
            # not ready for this yet# nebcont = add_stellar_and_nebular(hires_tab, clouddf)


#    # mess with this file, just to prove I know how.
#    outfile = re.sub('.fits', '_wnebcont.fits', thisfile)
#    tab['FLUX'][0, 0, :] += 1E5   # arbitrarily modify the file  **KLUDGE*
#    write_JC_S99file(tab, './', outfile) 
#    newtab = read_JC_S99file("./", outfile)   # Good, I can read the new version of the file.

