''' These functions and scripts read the JWST backgrounds that have been 
precompiled by STScI.  From prototype in precompiled_JWST_bkgs2.py.
Goals are to make bathtub curves of background versus visibility window,
and then analyze how much different thresholds on background level affect
target visibility.

This is the right schema of the precompiled backgrounds. They are binary files, 
read  following tutorial at http://vislab-ccom.unh.edu/~schwehr/rt/python-binary-files.html
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
import struct
import re
from os.path import basename
import healpy
import numpy as np
import pandas
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from cycler import cycler
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

#### Setup #############################################
sns.set(font_scale=2)
sns.set_style("white")
nside = 128 # from generate_backgroundmodel_cache.c .  
base_dir  = "/Volumes/Apps_and_Docs/MISSIONS/JWST/Zody_bathtubs/"  # Satchmo
bkg_dir   = base_dir + "sl_cache.v1.0/"
whichwaves = [1.0, 2.0, 5.0, 10.1, 15.1, 20.5]  # wavelengths (micron) to calc bkg
whichthresh = [1.05, 1.1, 1.3, 1.5, 2.0]        # thresholds over minimum to consider a good background
########################################################


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
    return(the_idex[0][0])

def myfile_from_healpix(healpix) :
    return ( str(healpix)[0:4] + "/sl_pix_" + str(healpix) + ".bin")

def calc_the_healpix(df):  # For each row of a dtaframe, take RA, DEC in degrees and calculate the healpix number
    df['healpix'] = df.apply(lambda row : format(healpy.pixelfunc.ang2pix(nside, row.RA_deg, row.DEC_deg, nest=False, lonlat=True), '06d'), axis=1)
    return(0)

def make_bathtub(results, wavelength_desired, thresh, showthresh=True, showplot=False, showall=False, title=False, label=False) :
    # thresh is threshold above minimum background, to calculate number of good days
    (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results # show how to break it up
    the_index = index_of_wavelength(wave_array, wavelength_desired)
    total_thiswave = total[ :, the_index]
    stray_thiswave = stray_light_bg[ :, the_index]
    zodi_thiswave =  zodi_bg[ :, the_index]
    
    themin = np.min(total_thiswave)
    allgood =  np.sum(total_thiswave < themin * thresh)*1.0
    if showplot:
        annotation = str(allgood) + " good days out of " + str(len(calendar)) + " days observable, for thresh " + str(thresh)
        plt.annotate(annotation, (0.05,0.05), xycoords="axes fraction", fontsize=12)
        plt.scatter(calendar, total_thiswave, s=20, label=label)
        if showall :
            plt.scatter(calendar, stray_thiswave, s=20, label=label)
            plt.scatter(calendar, zodi_thiswave, s=20, label=label)
        percentiles = (themin, themin*thresh)
        if showthresh : plt.hlines(percentiles, 0, 365, color='green')
        plt.xlabel("Day of the year")
        plt.xlim(0,366)
        plt.ylabel("bkg at " + str(wave_array[the_index]) + " um (MJy/SR)")
        if title : plt.title(title)
        #plt.show()
    return(allgood)  # Returns the number of days in the FOR with a background below 

    
###################################################

##########  This is the heart.  This reads the precompiled backgrounds, and calculates N good days.  Hrs to run.
Calc_Gooddays = False   # Loop through *every* position on the sky, and calculate how many days in the FOR
if Calc_Gooddays :      # have a background below a threshold, for several thresholds, at several wavelengths. 
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

def read_gooddays(base_dir=base_dir, infile="gooddays.txt", addGalEcl=False) :
    #infile = "gooddays.txt"    # "gooddays_worked_31mar2017.txt"
    df2 = pandas.read_csv(base_dir + infile, comment="#")
    df2.rename(columns={"file":"healfile", "RA":"RA_deg", "DEC":"DEC_deg"}, inplace=True)
    if addGalEcl :
        print "Converting from Equatorial Coords to Galactic and Ecliptic.  Takes a minute"
        jrr.util.convert_RADEC_GalEclip_df(df2, colra='RA_deg', coldec='DEC_deg')    # Compute Galactic, Eclptic coords
    return(df2)

def read_deepfields(addGalEcl=False) :
    #coords = "/Volumes/Apps_and_Docs/MISSIONS/JWST/Stray_light_commissioning/coordinates.txt" #satch
    coords = "coordinates.txt" # Milk
    deep_fields = pandas.read_table(coords, delim_whitespace=True, comment="#")
    jrr.util.convert_RADEC_segidecimal_df(deep_fields) # Convert RADEC to format healpy can read
    calc_the_healpix(deep_fields)
    if addGalEcl : jrr.util.convert_RADEC_GalEclip_df(deep_fields, colra='RA_deg', coldec='DEC_deg')
    return(deep_fields)

    
##############################################################################
# Below assumes that Calc_Gooddays above has already been run, that it wrote
# a "gooddays.txt" or similar output file, which we are now analyzing.
Analyze_Bathtubs = True 
if Analyze_Bathtubs :
    df2= read_gooddays(base_dir="", addGalEcl=False)
    pp = PdfPages("gooddays_new.pdf")
    x1 = 90; x2=366; y1=0; y2=1.05
    for wavelength_desired in whichwaves :
        for thresh in whichthresh :
            lab1 = "Ndays in FOR"
            lab2 = "fraction of days w " + str(wavelength_desired) + " um bkg <" + str(thresh) + " of minimum"
            lab3 = "Ndays w " + str(wavelength_desired) + " um bkg <" + str(thresh) + " of minimum"
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
            healpy.mollview(df2['good_frac'], min=0, max=1, cmap=cmap, flip='astro', title=lab2, coord='CE')
            pp.savefig()
            healpy.mollview(df2[goodcol], min=0, max=365, cmap=cmap, flip='astro', title=lab3, coord='CE')
            pp.savefig()
            
            plt.clf()
    pp.close()
    

healpix = 195257
Healpy_tutorial = False
if Healpy_tutorial :
    print "Convert from healpix number to RA, DEC:"
    print healpy.pixelfunc.pix2ang(nside, healpix, nest=False, lonlat=True)
    benchmark = (261.68333333, -73.33222222) # RA, DEC in decimal degrees of 1.2 min zody 
    print "Convert from RA, DEC of 1.2 min zody to healpix numbers:"
    healpix = healpy.pixelfunc.ang2pix(nside, benchmark[0], benchmark[1], nest=False, lonlat=True)
    print "Healpix for", benchmark, "was", healpix

Bathtub_tutorial = False
if Bathtub_tutorial :
    # got healpix from Healpy tutorial above. Should be 1.2 min zody benchmark
    myfile = myfile_from_healpix(healpix)
    thresh=1.1;  wavelength_desired = whichwaves[1]
    results = read_JWST_precompiled_bkg(myfile, base_dir, bkg_dir, showplot=True)
    (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results 
    allgood = make_bathtub(results, wavelength_desired, thresh, showplot=True)
    print "Ran", myfile, wavelength_desired, "micron", thresh, "threshold"
    #print allgood, "good days out of", len(calendar)
    
# Now, calculate bathtubs for commissioning stray light positions
Bathtubs_deepfields = False
if Bathtubs_deepfields :
    deep_fields = read_deepfields()
    pp = PdfPages("deepfields_bathtubs.pdf")
    for row in deep_fields.itertuples() :
        myfile = myfile_from_healpix(row.healpix)
        results = read_JWST_precompiled_bkg(myfile, base_dir, bkg_dir, showplot=False)
        (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results
        print "Running for", row.FIELDS, row.RA_deg, row.DEC_deg, row.healpix, len(calendar)
        for wavelength_desired in whichwaves[0:3] :
            plt.clf()
            title = re.sub("_", " ", row.FIELDS) + " at " + str(wavelength_desired) + " micron"
            allgood = make_bathtub(results, wavelength_desired, 1.1, showplot=True, title=title)
            pp.savefig()
    pp.close()
    plt.close("all")
# Repeat, but all fields on one plot, for each wavelength
    sns.set_palette("hls", 10)
    pp = PdfPages("deepfields_bathtubs_2.pdf")
    for wavelength_desired in whichwaves :
        title = "Stray light calib fields at " + str(wavelength_desired) + " micron"
        plt.clf()
        for row in deep_fields.itertuples() :
            print "Running for", row.FIELDS, row.RA_deg, row.DEC_deg, row.healpix
            myfile = myfile_from_healpix(row.healpix)
            results = read_JWST_precompiled_bkg(myfile, base_dir, bkg_dir, showplot=False)
            (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results
            allgood = make_bathtub(results, wavelength_desired, 1.1, showthresh=False, showplot=True, title=title, label=re.sub("_", " ", row.FIELDS))
        plt.legend(fontsize=10, frameon=False, labelspacing=0)
        pp.savefig()
    pp.close()


# Now, make 2D maps on the sky, and do math by Ecliptic latitude
Plot_on_sky = False
if Plot_on_sky :
    plt.clf()
    df2= read_gooddays(base_dir="", addGalEcl=False)
    cmap = sns.cubehelix_palette(n_colors=10,  start=0, rot=0.0, gamma=2.0, hue=1, light=1, dark=0.4, reverse=False, as_cmap=True)
    healpy.mollview(df2['Nday'], cmap=cmap, flip='astro', title="Ndays observable", coord='CE')
    plt.show()   
    df2= read_gooddays(base_dir="", addGalEcl=True)
    binlatby = 1
    latbins = np.arange(-90,90,binlatby)
    groups = df2.groupby(pandas.cut(df2['Ecl_lat'], latbins))
    plt.plot(latbins[:-1], groups['Good1.0_1.1'].median().values, label="1 micron", color='red')
    plt.plot(latbins[:-1], groups['Good2.0_1.1'].median().values, label="2 micron", color='orange')
    plt.plot(latbins[:-1], groups['Good5.0_1.1'].median().values, label="5 micron", color='green')
    plt.plot(latbins[:-1], groups['Good10.1_1.1'].median().values, label="10 micron", color='blue')
    plt.plot(latbins[:-1], groups['Nday'].median().values, label="Nday", color='black')
    plt.xlabel("Ecliptic Latitude (deg)")
    plt.ylabel("N days with background <1.1 of min")
    plt.xlim(-90,90)
    plt.locator_params(axis='x', nbins=6)
    plt.ylim(0,365)
    plt.legend(labelspacing=0.3)
    plt.grid()
    plt.show()
    # pseudo code:  groupby Ecliptic Latitude, and compute statistics
