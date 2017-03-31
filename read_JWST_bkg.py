''' Writing a function to read STScIs precompiled backgrounds.
Working prototype in  precompiled_JWST_bkgs2.py Now, rewriting as a function.
Goal is to make bathtub curves of zody vs time.
To read binary files, following tutorial at http://vislab-ccom.unh.edu/~schwehr/rt/python-binary-files.html

This is the right schema.
  Verified against source code generate_stray_light_with_threads.c
#double RA
#double DEC
#double pos[3]
#double nonzodi_bg[SL_NWAVE]
#int[366] date_map  # This maps dates to indices
#for each day in FOR:
#  double zodi_bg[SL_NWAVE]
#  double stray_light_bg[SL_NWAVE]
'''

import jrr
import glob
import re
from os.path import basename
import struct
import numpy as np
import pandas
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
sns.set(font_scale=2)
sns.set_style("white")

def rebin_spec_new(wave, specin, new_wave, fill=np.nan):
    f = interp1d(wave, specin, bounds_error=False, fill_value=fill)  # With these settings, writes NaN to extrapolated regions
    new_spec = f(new_wave)
    return(new_spec)

def read_JWST_precompiled_bkg(infile, base_dir, bkg_dir, showplot=False, verbose=False) :
    wave_file = base_dir + "updated_std_spectrum_wavelengths.txt"  # Standard wavelength array.  Should be SL_NWave=108 long
    wave_array = np.loadtxt(wave_file)    
    SL_NWAVE = len(wave_array)  # should be 108.  Size of wavelength array

    thermal_file = "thermal_curve_jwst_jrigby_1.1.csv"  # The constant (not time variable) thermal self-emission curve
    temp_thermal = np.genfromtxt(base_dir + thermal_file, delimiter=',')
    thermal = rebin_spec_new(temp_thermal[:, 0], temp_thermal[:,1],  wave_array, fill=0.0)  # rebin to same wavelength_array as others.
    
    sbet_file = open(bkg_dir + myfile)
    sbet_data = sbet_file.read()
    # Unpack the constant first part
    if verbose: print "File has", len(sbet_data), "bytes, which is", len(sbet_data)/8., "doubles"
    size_calendar = struct.calcsize("366i") # bytes, not doubles
    partA = struct.unpack(str(5 + SL_NWAVE)+'d', sbet_data[0: (5 + SL_NWAVE)*8])
    RA = partA[0]
    DEC = partA[1]
    pos = partA[2:5]
    nonzodi_bg = np.array(partA[5:5+SL_NWAVE])

    # Unpack the calendar dates      # code goes from 0 to 365 days.
    date_map = np.array(struct.unpack('366i', sbet_data[(5 + SL_NWAVE)*8  : (5 + SL_NWAVE)*8 + size_calendar]))
    if verbose: print "Out of", len(date_map), "days, these many are legal:", np.sum(date_map >=0)
    #print "indices of days:", date_map[date_map>=0]
    calendar = np.where(date_map >=0)[0]
    #print "calendar date:", calendar
    # So, the index dd in zodi_bg[dd, : ]  corresponds to the calendar day lookup[dd]
    Ndays = len(calendar) 
    if verbose: print len(date_map), Ndays

    # Unpack part B, the time-variable part
    zodi_bg        = np.zeros((Ndays,SL_NWAVE))
    stray_light_bg = np.zeros((Ndays,SL_NWAVE))
    perday = SL_NWAVE*2
    partB= struct.unpack(str((len(calendar))*SL_NWAVE*2)+'d', sbet_data[perday*Ndays*-8 : ])

    for dd in range(0, int(Ndays)):
        br1 = dd*perday
        br2 = br1 + SL_NWAVE
        br3 = br2 + SL_NWAVE
        #print "Breaking at:", br1, br2, br3
        zodi_bg[dd, ]        = partB[br1 : br2]
        stray_light_bg[dd, ] = partB[br2 : br3]
    expand = np.ones((Ndays,SL_NWAVE))  # same shape as zodi_bg
    total = nonzodi_bg * expand + thermal * expand + stray_light_bg + zodi_bg

    if showplot :
        plt.clf()
        thisday = 100
        plt.plot(wave_array, nonzodi_bg, label="ISM")
        plt.plot(wave_array, zodi_bg[thisday, :], label="Zodi")
        plt.plot(wave_array, stray_light_bg[thisday, :], label="Stray light")
        plt.plot(wave_array, thermal, label = "Thermal")
        plt.plot(wave_array, total[thisday, :], label = "Total", color='black', lw=3)
        plt.xlim(0.6,31)
        plt.legend()
        plt.yscale('log')
        plt.show()
    return((calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total))  #send it as a tuple

def index_of_wavelength(wave_array, desired_wavelength) :  # look up index of wavelength array corresponding to desired wavelength
    the_index = np.where(wave_array == desired_wavelength)
    return(the_index[0][0])

def make_bathtub(results, wavelength_desired, thresh, showplot=False) :
    # thresh is threshold above minimum background, to calculate number of good days
    (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results # show how to break it up
    the_index = index_of_wavelength(wave_array, wavelength_desired)
    total_thiswave = total[ :, the_index]
    themin = np.min(total_thiswave)
    allgood =  np.sum(total_thiswave < themin * thresh)*1.0
    if showplot: 
        plt.scatter(calendar, total_thiswave)
        percentiles = (themin, themin*thresh)
        plt.hlines(percentiles, 0, 365, color='green')
        plt.xlabel("Day of the year")
        plt.xlim(0,366)
        plt.ylabel("bkg at " + str(wave_array[the_index]) + " um (MJy/SR)")
        plt.show()
    return(allgood)  # Returns the number of days in the FOR with a background below 

    
###################################################
# Setup

base_dir  = "/Volumes/Apps_and_Docs/MISSIONS/JWST/Zody_bathtubs/"  # Satchmo
bkg_dir   = base_dir + "sl_cache.v1.0/"
whichwaves = [1.0, 2.0, 5.0, 10.1, 15.1, 20.5]
whichthresh = [1.05, 1.1, 1.3, 1.5, 2.0]

Bathtub_tutorial = False    # Tutorial, make one example bathtub
if Bathtub_tutorial :
    myfile = "1464/sl_pix_146447.bin"
    thresh=1.1;  wavelength_desired = whichwaves[-1]
    results = read_JWST_precompiled_bkg(myfile, base_dir, bkg_dir, showplot=True)
    (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results 
    allgood = make_bathtub(results, wavelength_desired, thresh, showplot=True)
    print "Ran", myfile, wavelength_desired, "micron", thresh, "threshold"
    #print allgood, "good days out of", len(calendar)

Calc_Gooddays = False   # Loop through every position on the sky, and calculate how many days in the FOR
if Calc_Gooddays :     # have a background below a threshold, for several thresholds, at several wavelengths 
    healpix_dirs = glob.glob(bkg_dir + "*/")
    dirs_to_run = healpix_dirs   
    len(glob.glob(bkg_dir + "*/*bin"))
    allGood   = np.zeros(shape=(len(whichwaves), len(whichthresh), 100*len(dirs_to_run)))  # Big array
    allNday   =  []
    whichfile  = []
    allRA = []
    allDEC = []
    ii = 0  # file_index
    for thisdir in dirs_to_run :
        myfiles = [ basename(x) for x in glob.glob(thisdir + "*.bin") ]
        print len(myfiles), thisdir
        for myfile in myfiles:
            results = read_JWST_precompiled_bkg(thisdir + myfile, base_dir, thisdir, showplot=False)
            (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results 
            allNday.append(len(calendar))
            for jj, wavelength_desired in enumerate(whichwaves) :
                for kk, thresh in enumerate(whichthresh):
                    this_allgood = (make_bathtub(results, wavelength_desired, thresh, showplot=False))
                    allGood[jj,kk,ii] = this_allgood
            whichfile.append(myfile)
            allRA.append(RA)
            allDEC.append(DEC)
            ii += 1  # increment monotonic file counter
    # Done running everything.  Now, make pretty df and output
    df = pandas.DataFrame( { 'file':whichfile, 'RA':allRA, 'DEC':allDEC, 'Nday':allNday})
    df = df[['file', 'RA', 'DEC', 'Nday']]  # reorder
    
    for jj, wavelength_desired in enumerate(whichwaves) :
        for kk, thresh in enumerate(whichthresh):
            goodcol = "Good"+str(wavelength_desired)+ "_" + str(thresh)
            df[goodcol] = allGood[jj,kk, :ii]   # does this fix the problem?
    df.to_csv('tmp')
    header = "#Number of Days in FOR, and number of days with good background, for waves" + str(whichwaves) + " micron, and thresholds" + str(whichthresh) + "\n"
    jrr.util.put_header_on_file('tmp', header, "gooddays.txt") 

    
Analyze_Bathtubs = True   # Once rereun finished, remove whichwaves and whichtresh below, and plot for all
if Analyze_Bathtubs :
    infile = "gooddays_worked_31mar2017.txt"
    pp = PdfPages(re.sub(".txt",  ".pdf", infile))  # the outfile
    df2 = pandas.read_csv(base_dir + infile, comment="#")
    x1 = 90; x2=366; y1=0; y2=1.05

    whichwaves = [1.0, 2.0, 10.1]  
    whichthresh = [1.1, 1.3]
    for wavelength_desired in whichwaves :
        for thresh in whichthresh :
            lab1 = "Ndays in FOR"
            lab2 = "fraction of days w " + str(wavelength_desired) + " um bkg <" + str(thresh) + " of minimum"
            goodcol =  "Good"+str(wavelength_desired)+ "_" + str(thresh)
            print wavelength_desired, thresh, goodcol
            df2['good_frac'] = df2[goodcol] / df2['Nday']
            if False :
                plt.scatter(df2['Nday'], df2['good_frac'], s=1)
                plt.xlabel(lab1)
                plt.ylabel(lab2)
                plt.xlim(x1,x2)
                plt.ylim(y1,y2)
                pp.savefig()
            # Density plot, with histograms on margins.  From NB example,
            cmap = sns.cubehelix_palette(n_colors=10,  start=0, rot=0.0, gamma=2.0, hue=1, light=1, dark=0.4, reverse=False, as_cmap=True)
            g = sns.JointGrid(df2['Nday'], df2['good_frac'], size=8)
            g.set_axis_labels(lab1, lab2)
            g.ax_marg_x.hist(df2['Nday'], bins=np.arange(x1, x2, 8))
            g.ax_marg_y.hist(df2['good_frac'], bins=np.arange(y1,y2,0.03), orientation="horizontal")
            g.plot_joint(plt.hexbin, gridsize=100, extent=[x1, x2, y1, y2], cmap=cmap, mincnt=1, bins='log')
            pp.savefig()
    pp.close()
    
# Since the index of the df is the healpix value, should *should* be easy to read in w healpy,
# and then plot good_frac in healpy as a projection.  Bet I'll rediscover the galaxy, the Ecliptic
