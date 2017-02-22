import shutil
import jrr
import glob
from os.path import basename
import pandas
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy import constants, units
# Needed for barycentric correction
from astropy import coordinates
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation

# Combine (straight sum), with clipping, the ESI/Keck spectra.
# jrigby, Feb 2017


def get_the_spectra(filenames, obslog, colwav='obswave') :
    df = {}
    for thisfile in filenames :
        print "loading file ", thisfile
        df[thisfile] = pandas.read_table(thisfile, delim_whitespace=True, comment="#")
        # Apply the Barycentric correction
        thisobs = obslog.loc[thisfile]
        keck = EarthLocation.of_site('keck')
        my_target = SkyCoord(thisobs['RA'], thisobs['DEC'], unit=(units.hourangle, units.deg), frame='icrs')
        my_start_time = Time( thisobs['DATE-OBS'] + "T" + thisobs['UT_START'] , format='isot', scale='utc')
        midpt = TimeDelta(thisobs['EXPTIME'] / 2.0, format='sec')
        my_time = my_start_time + midpt  # time at middle of observation
        barycor_vel = jrr.barycen.compute_barycentric_correction(my_time, my_target, location=keck)
        #print "DEBUGGING", my_target, thisobs, my_target, my_start_time
        print "FYI, the barycentric correction factor for", thisfile,  "was", barycor_vel
        jrr.barycen.apply_barycentric_correction(df[thisfile], barycor_vel, colwav='obswave', colwavnew='wave') # testing   
        df[thisfile]['Nfiles'] = 1 # N of exposures that went into this spectrum
    return(df)  # return a dictionary of dataframes of spectra

def jrr_filter(out, thresh=2.5, mintoreplace=1E-16, low_wave_cut=4600.) :
    ''' Filter the spectrum.  First, remove the low-wavelength stuff b/c it's junk.
    Second, filter out errant sky lines, without altering the line ratios, as follows:
    Make new column, ['flam_sum_jrr'], which is mostly straight sum ([flam][sum]),
    but replaces extreme values (sum > thresh*median*N) w the median*N    '''
    out2 = out.copy(deep=True)
    out2 = out2[out2.index > low_wave_cut]
    thefilter =  (out2['flam'] > thresh * out2['flam_median_xN']) & (out2['flam'] > mintoreplace)
    out2['flam_sum_jrr']  = np.where(thefilter, out2['flam_median_xN'], out2['flam'])
    out2['replaced']      = np.where(thefilter, True, False)
    return(out2) 

def jrr_filter_old(out, thresh=2.5, mintoreplace=1E-16, low_wave_cut=4600.) :
    ''' Filter the spectrum.  First, remove the low-wavelength stuff b/c it's junk.
    Second, filter out errant sky lines, without altering the line ratios, as follows:
    Make new column, ['flam_sum_jrr'], which is mostly straight sum ([flam][sum]),
    but replaces extreme values (sum > thresh*median*N) w the median*N    '''
    out2 = out.copy(deep=True)
    out2 = out2[out2.index > low_wave_cut]
    thefilter =  (out2['flam']['sum'] > thresh * out2['flam_median_xN']) & (out2['flam']['sum'] > mintoreplace)
    out2['flam_sum_jrr']  = np.where(thefilter, out2['flam_median_xN'], out2['flam']['sum'])
    out2['replaced']      = np.where(thefilter, True, False)
    return(out2) 
    
def plot_the_results(filenames, df, out2, groupby=False) :
    plt.clf()
    print "Plotting...."
    xlim = (8000, 9000)
    fig1 = plt.figure(1, figsize=(24,6))
    for thisfile in filenames :
        plt.plot(df[thisfile]["wave"], df[thisfile]["flam"], label=thisfile)
    plt.legend()
    plt.ylim(0,1E-16)
    plt.xlim(*xlim)
    plt.show()
    fig2 = plt.figure(2, figsize=(20,5))
    plt.ylim(0,7E-16)
    plt.xlim(*xlim)
    if groupby :
        plt.plot(out2.wave, out2['flam']['sum'], color='purple', label="flam")
        plt.plot(out2.wave, out2['flam_median_xN'], color='green', label="median*N")
        plt.plot(out2.wave, out2['flam_u']['add_in_quad'], color='orange', label="flam_u")
        plt.scatter(out2.wave, out2['replaced']*out2['flam_sum_jrr'], label="replaced")
        plt.plot(out2.wave, out2['flam_sum_jrr'], color='black')
    else :
        plt.plot(out2.wave, out2['flam'], color='purple', label="flam")
        plt.plot(out2.wave, out2['flam_median_xN'], color='green', label="median*N")
        plt.plot(out2.wave, out2['flam_u'], color='orange', label="flam_u")
        plt.scatter(out2.wave, out2['replaced']*out2['flam_sum_jrr'], label="replaced")
        plt.plot(out2.wave, out2['flam_sum_jrr'], color='black')
    plt.legend()
    plt.show()
    return(0)
            
def write_spectrum_to_file(out2, outfile, prefix, filenames, header, thresh) :
    '''Make a simplifed outfile '''
    print "Writing to file"
    cols = ['wave', 'flam_sum_jrr', 'flam_u', 'flam_median_xN', 'Nfiles', 'replaced']
    simple_df = out2[cols] 
    shutil.copy(header, outfile)
    with open(outfile, 'a') as f:  simple_df.to_csv(f, sep='\t')
    jrr.util.replace_text_in_file("OBJECTNAME", prefix, outfile)
    jrr.util.replace_text_in_file("THRESHTHRESH", str(thresh), outfile)
    jrr.util.replace_text_in_file("INDYSPECTRANAMES", ''.join(filenames), outfile)
    return(0)

def make_a_stack(filenames, obslog, prefix, thresh, mintoreplace, low_wave_cutoff) :     # Put it all together
    df = get_the_spectra(filenames, obslog)
    out = jrr.spec.stack_observed_spectra(df, do_sum=True, colwav='wave', colf='flam', colfu='flam_u')
    out2 = jrr_filter(out, thresh, mintoreplace, low_wave_cutoff)    
    plot_the_results(filenames, df, out2, groupby=False)
    outfile = prefix + "_ESI_JRR_sum.txt"
    write_spectrum_to_file(out2, outfile, prefix, filenames, "../JRR_header.txt", thresh)
    return(out2)


def get_observation_datetimes() :
    infile = "../ESI_target_datetime_observed.txt"
    obslog = pandas.read_table(infile, delim_whitespace=True, comment="#")
    obslog['galaxy_name'] = obslog['galaxy_name'].str.replace(".txt", "_esi.txt")  # Ayan filename inconsistencies
    obslog.set_index("galaxy_name", inplace=True, drop=False)  # change index
    return(obslog)
        
#############################
#### Actually run things ####

def run_2016Aug() :
    obslog = get_observation_datetimes()
    prefix = "s1723" 
    filenames = [ basename(x) for x in glob.glob(prefix+"*_esi.txt") ]
    out2 = make_a_stack(filenames, obslog, prefix, thresh=2.5, mintoreplace=1E-16, low_wave_cutoff=4600)
    prefix = "s2340"
    filenames = [ basename(x) for x in glob.glob(prefix+"*_esi.txt") ]
    out2 = make_a_stack(filenames, obslog, prefix, thresh=2.5, mintoreplace=0.5E-16, low_wave_cutoff=4600)
    return(0)

def run_2016Apr() :
    obslog = get_observation_datetimes()
    prefix = "J1050"
    filenames = [ basename(x) for x in glob.glob(prefix+"*_esi.txt") ]
    out2 = make_a_stack(filenames, obslog, prefix, thresh=2.5, mintoreplace=1E-16, low_wave_cutoff=4600)
    prefix = "J1458"
    filenames = [ basename(x) for x in glob.glob(prefix+"*_esi.txt") ]
    out2 = make_a_stack(filenames, obslog, prefix, thresh=2.5, mintoreplace=1E-16, low_wave_cutoff=4600)
    return(0)
    

#############################

