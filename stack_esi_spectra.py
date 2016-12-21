import shutil
import jrr
import glob
from os.path import basename
import pandas
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip, sigma_clipped_stats

# Given a list of spectra, want a function that does a weighted average, and
# a straight average, and a median.  W error spectra.    Basically, generalize stacker


def get_the_spectra(filenames) :
    df = {}
    for thisfile in filenames :
        print "loading file ", thisfile
        df[thisfile] = pandas.read_table(thisfile, delim_whitespace=True, comment="#")
        df[thisfile]['Nfiles'] = 1 # N of exposures that went into this spectrum
    return(df)  # return a dictionary of dataframes of spectra

def sum_the_spectra(df) :
    print "Summing the spectra with custom functions.  This is slow."
    df_all = pandas.concat(df)   # Life is easy b/c same rest_wave array for all spectra
    sigclip =3
    f = {'obswave':np.mean, 'flam': [np.sum, np.median], 'flam_u':jrr.util.add_in_quad, 'restwave':np.mean, 'badmask':np.sum, 'Nfiles':np.sum}
    grouped = df_all.groupby(by="obswave")    # GROUPBY!!! I love you!
    out = grouped.agg(f)  # Aggregate the big dataframe by wavelength, and then aggregate by a bunch of functions
    out['flam_median_xN'] = out['flam']['median'] * out['Nfiles'].median().values
    return(out)

def jrr_filter(out, thresh=2.5, mintoreplace=1E-16, low_wave_cut=4600.) :
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

def plot_the_results(filenames, df, out2) :
    plt.clf()
    print "Plotting...."
    xlim = (8000, 9900)
    fig1 = plt.figure(1, figsize=(24,6))
    for thisfile in filenames :
        plt.plot(df[thisfile]["obswave"], df[thisfile]["flam"], label=thisfile)
    plt.legend()
    plt.ylim(0,1E-16)
    plt.xlim(*xlim)
    plt.show()
    fig2 = plt.figure(2, figsize=(20,5))
    plt.ylim(0,7E-16)
    plt.xlim(*xlim)
    plt.plot(out2.obswave, out2['flam']['sum'], color='purple', label="sum")
    #plt.plot(out2.obswave, out2['flam_median_xN'], color='green', label="median*N")
    plt.plot(out2.obswave, out2['flam_u'], color='orange', label="sigma(flam)")
    plt.scatter(out2.obswave, out2['replaced']*out2['flam_sum_jrr'], label="replaced")
    plt.plot(out2.obswave, out2['flam_sum_jrr'], color='black')
    plt.legend()
    plt.show()
    return(0)
            
def write_spectrum_to_file(out2, outfile, prefix, filenames, header, thresh) :
    '''There are way too many columns in out2, and it's MultiIndexed.  I don't want to give
    that mess to collaborators.  Instead, make a simpler new df, and write it to file.
    Kludgy, but I can't figure out syntax for MultiIndexed in .to_csv()'''
    print "Writing to file"
    simple_df = pandas.DataFrame(data=out2['flam_sum_jrr'])
    simple_df['flam_u']     =  out2['flam_u']['add_in_quad']
    simple_df['flam_med*N'] = out2['flam_median_xN']
    simple_df['Nfiles']     = out2['Nfiles']
    simple_df['replaced']   = out2['replaced']
    shutil.copy(header, outfile)
    with open(outfile, 'a') as f:
        simple_df.to_csv(f, sep='\t')
    jrr.util.replace_text_in_file("OBJECTNAME", prefix, outfile)
    jrr.util.replace_text_in_file("THRESHTHRESH", str(thresh), outfile)
    jrr.util.replace_text_in_file("INDYSPECTRANAMES", ''.join(filenames), outfile)
    return(0)

def make_a_stack(filenames, prefix, thresh, mintoreplace, low_wave_cutoff) :     # Put it all together
    df = get_the_spectra(filenames)
    out = sum_the_spectra(df)  # This step can be slow.
    out2 = jrr_filter(out, thresh, mintoreplace, low_wave_cutoff)    
    plot_the_results(filenames, df, out2)
    outfile = prefix + "_ESI_JRR_sum.txt"
    write_spectrum_to_file(out2, outfile, prefix, filenames, "../JRR_header.txt", thresh)
    return(out2)

#############################
#### Actually run things ####

def run_2016Aug() :
    prefix = "s1723" 
    filenames = [ basename(x) for x in glob.glob(prefix+"*_esi.txt") ]
    df = make_a_stack(filenames, prefix, thresh=2.5, mintoreplace=1E-16, low_wave_cutoff=4600)
    prefix = "s2340"
    filenames = [ basename(x) for x in glob.glob(prefix+"*_esi.txt") ]
    df = make_a_stack(filenames, prefix, thresh=2.5, mintoreplace=0.5E-16, low_wave_cutoff=4600)
    return(0)

def run_2016Apr() :
    prefix = "J1050"
    filenames = [ basename(x) for x in glob.glob(prefix+"*_esi.txt") ]
    df = make_a_stack(filenames, prefix, thresh=2.5, mintoreplace=1E-16, low_wave_cutoff=4600)
    prefix = "J1458"
    filenames = [ basename(x) for x in glob.glob(prefix+"*_esi.txt") ]
    df = make_a_stack(filenames, prefix, thresh=2.5, mintoreplace=1E-16, low_wave_cutoff=4600)
    return(0)
    


#############################

