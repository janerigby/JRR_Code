''' These functions and scripts read the JWST background cache that was 
precompiled by STScI. Goals are to make bathtub curves of background versus 
visibility window, and to analyze how much background thresholds affect 
target visibility.   Jane.Rigby@nasa.gov, Apr 2017.
#
# Here is the schema for the precompiled background cache.  (I verified the schema
# against the source code, generate_stray_light_with_threads.c)
# The cache uses a Healpix RING tesselation, with NSIDE=128.  Every point on the sky 
# (tesselated tile) corresponds to one binary file, whose name includes its healpix 
# pixel number, in a directory corresponding to the first 4 digits of the healpix number.  
# I used this tutorial to read the binary files: 
#      http://vislab-ccom.unh.edu/~schwehr/rt/python-binary-files.html
# Anyway, here's the schema of each binary file:
#double RA
#double DEC
#double pos[3]
#double nonzodi_bg[SL_NWAVE]
#int[366] date_map  # This maps dates to indices.  **There are 366 days, not 365!**
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
nside = 128  # Healpy parameter, from generate_backgroundmodel_cache.c .  
base_dir  = "/Volumes/Apps_and_Docs/MISSIONS/JWST/Zody_bathtubs/"  # Satchmo
coords_dir = "/Volumes/Apps_and_Docs/MISSIONS/JWST/Stray_light_commissioning/"
bkg_dir   = base_dir + "sl_cache.v1.0/"

trywaves  =  [1.0, 2.0, 5.0, 10.1, 15.1, 20.5]  # wavelengths (micron) to calc bkg
trythresh = [1.05, 1.1, 1.3, 1.5, 2.0]        # thresholds over minimum to consider a good background

# Above was OK for initial exploration.  Now, switching to pivot wavelengths of broad-band filters.
# Or rather, their closest wavelengths in the ETC wavelength array.
nircam_broad = [0.7, 0.9, 1.1, 1.5, 2.0, 2.8, 3.5, 4.5]
miri_broad   = [5.5, 7.7, 10.1, 12.7, 15.1, 17.5, 21.5, 25.5]
good_nircammiri = "gooddays_nircam_miri.txt"


########################################################

whichwaves = nircam_broad + miri_broad
whichthresh = [1.1]

def rebin_spec_new(wave, specin, new_wave, fill=np.nan):
    f = interp1d(wave, specin, bounds_error=False, fill_value=fill)  # With these settings, writes NaN to extrapolated regions
    new_spec = f(new_wave)
    return(new_spec)

def read_JWST_precompiled_bkg(infile, base_dir, bkg_dir, showplot=False, verbose=False) :
    # Reads one JWST background file, in binary format.
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
        plt.xlabel("wavelength (micron)")
        plt.ylabel("Equivalent in-field radiance (MJy/SR)")
        plt.legend()
        plt.yscale('log')
        plt.show()
    return((calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total))  #pack it up as a tuple

def index_of_wavelength(wave_array, desired_wavelength) :  # look up index of wavelength array corresponding to desired wavelength
    the_index = np.where(wave_array == desired_wavelength)
    return(the_index[0][0])

def myfile_from_healpix(healpix) :
    return ( str(healpix)[0:4] + "/sl_pix_" + str(healpix) + ".bin")

def calc_the_healpix(df):  # For each row of a dataframe, take RA, DEC in degrees and calculate the healpix number
    df['healpix'] = df.apply(lambda row : format(healpy.pixelfunc.ang2pix(nside, row.RA_deg, row.DEC_deg, nest=False, lonlat=True), '06d'), axis=1)
    df['healpix'] = df['healpix'].astype('str')
    df['healpix_asindex'] = df['healpix'].str.lstrip('0').astype('int')
    return(0)

def make_bathtub(results, wavelength_desired, thresh, showthresh=True, showplot=False, showsubbkgs=False, showannotate=True, title=False, label=False) :
    # Once binary file was read w  read_JWST_precompiled_bkg(), compute the bathtub, and optionally, make a plot.
    # thresh is threshold above minimum background, to calculate number of good days
    (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results # show how to break it up
    the_index = index_of_wavelength(wave_array, wavelength_desired)
    total_thiswave = total[ :, the_index]
    stray_thiswave = stray_light_bg[ :, the_index]
    zodi_thiswave =  zodi_bg[ :, the_index]
    
    themin = np.min(total_thiswave)
    allgood =  np.sum(total_thiswave < themin * thresh)*1.0
    if showplot:
        if showannotate:
            annotation = str(allgood) + " good days out of " + str(len(calendar)) + " days observable, for thresh " + str(thresh)
            print str(wavelength_desired) + " " + annotation
            plt.annotate(annotation, (0.05,0.05), xycoords="axes fraction", fontsize=12)
        if not label : label="Total"
        plt.scatter(calendar, total_thiswave, s=20, label=label)
        if showsubbkgs :
            plt.scatter(calendar, zodi_thiswave, s=20, label="Zodiacal")
            plt.scatter(calendar, stray_thiswave, s=20, label="Stray light")
            plt.legend()
        percentiles = (themin, themin*thresh)
        if showthresh : plt.hlines(percentiles, 0, 365, color='green')
        plt.xlabel("Day of the year")
        plt.xlim(0,366)
        plt.ylabel("bkg at " + str(wave_array[the_index]) + " um (MJy/SR)")
        if title : plt.title(title)
        #plt.show()
    return(allgood)  # Returns the number of days in the FOR with a background below 

def read_gooddays(base_dir=base_dir, infile="gooddays.txt", addGalEcl=False) :
    df2 = pandas.read_csv(base_dir + infile, comment="#")
    df2.rename(columns={"file":"healfile", "RA":"RA_deg", "DEC":"DEC_deg"}, inplace=True)
    if addGalEcl :
        print "Converting from Equatorial Coords to Galactic and Ecliptic.  Takes a minute"
        jrr.util.convert_RADEC_GalEclip_df(df2, colra='RA_deg', coldec='DEC_deg')    # Compute Galactic, Eclptic coords
    return(df2)

def read_deepfields(filename, addGalEcl=False) :
    deep_fields = pandas.read_table(filename, delim_whitespace=True, comment="#")
    jrr.util.convert_RADEC_segidecimal_df(deep_fields) # Convert RADEC to format healpy can read
    calc_the_healpix(deep_fields)    
    if addGalEcl : jrr.util.convert_RADEC_GalEclip_df(deep_fields, colra='RA_deg', coldec='DEC_deg')
    return(deep_fields)

###################################################


## This is the heart of my analysis.  Reads the precompiled background files for every
# point on the sky, and calculates the number of low-background ("good") days as a function
# of wavelength and threshhold.  Takes a few hours to run b/c of input/output of 2E5 files.
# Write output to file named gooddays_XX.txt, to be analyzed by other parts below.
Calc_Gooddays = False  
if Calc_Gooddays :     
    outfile =  good_nircammiri
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
    jrr.util.put_header_on_file('tmp', header, outfile) 


##############################################################################
# Below assumes that Calc_Gooddays above has already been run, that it wrote
# a "gooddays.txt" or similar output file, which we are now analyzing.
Analyze_Bathtubs = False 
if Analyze_Bathtubs :
    df2= read_gooddays(infile=good_nircammiri, base_dir="", addGalEcl=False)
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
Healpy_tutorial = True
if Healpy_tutorial :
    print "Convert from healpix number to RA, DEC:"
    print healpy.pixelfunc.pix2ang(nside, healpix, nest=False, lonlat=True)
    benchmark = (261.68333333, -73.33222222) # RA, DEC in decimal degrees of 1.2 min zody 
    print "Convert from RA, DEC of 1.2 min zody to healpix numbers:"
    healpix = healpy.pixelfunc.ang2pix(nside, benchmark[0], benchmark[1], nest=False, lonlat=True)
    print "Healpix for", benchmark, "was", healpix

Bathtub_tutorial = True
if Bathtub_tutorial :
    # got healpix from Healpy tutorial above. Should be 1.2 min zody benchmark
    myfile = myfile_from_healpix(healpix)
    thresh=1.1;  wavelength_desired = whichwaves[4]
    results = read_JWST_precompiled_bkg(myfile, base_dir, bkg_dir, showplot=True)
    (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results 
    allgood = make_bathtub(results, wavelength_desired, thresh, showplot=True, showsubbkgs=True)
    print "Ran", myfile, wavelength_desired, "micron", thresh, "threshold"
    print allgood, "good days out of", len(calendar)
    
# Now, make plots of bathtubs for famous deep fields, and fields used in commissioning stray light.
Bathtubs_deepfields = False
if Bathtubs_deepfields :
    thefiles = (coords_dir + "deep_fields.txt", coords_dir + "commissioning_targets.txt")
    outpdf = ("deepfields_bathtubs.pdf", "commissioning_bathtubs.pdf")
    for ii, infile in enumerate(thefiles) :
        deep_fields = read_deepfields(infile)
        pp = PdfPages(outpdf[ii])
        for row in deep_fields.itertuples() :
            myfile = myfile_from_healpix(row.healpix)
            results = read_JWST_precompiled_bkg(myfile, base_dir, bkg_dir, showplot=False)
            (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results
            print "Running for", row.FIELD, row.RA_deg, row.DEC_deg, row.healpix, len(calendar)
            for wavelength_desired in whichwaves :
                plt.clf()
                title = re.sub("_", " ", row.FIELD) + " at " + str(wavelength_desired) + " micron"
                allgood = make_bathtub(results, wavelength_desired, 1.1, showplot=True, showannotate=True, title=title)
                plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
                pp.savefig()
        pp.close()
        plt.close("all")
        # Repeat, but all fields on one plot, for each wavelength
        sns.set_palette("hls", 12)
        pp = PdfPages(re.sub(".pdf", "_2.pdf", outpdf[ii]))
        for wavelength_desired in whichwaves :
            title = "Fields at " + str(wavelength_desired) + " micron"
            plt.clf()
            for row in deep_fields.itertuples() :
                print "Running for", row.FIELD, row.RA_deg, row.DEC_deg, row.healpix
                myfile = myfile_from_healpix(row.healpix)
                results = read_JWST_precompiled_bkg(myfile, base_dir, bkg_dir, showplot=False)
                (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results
                allgood = make_bathtub(results, wavelength_desired, 1.1, showthresh=False, showplot=True, title=title, showannotate=False, label=re.sub("_", " ", row.FIELD))
            plt.legend(fontsize=10, frameon=False, labelspacing=0)
            plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
            pp.savefig()
        pp.close()

# Make a table of good days versus wavelength for the deep fields.
Gooddays_deepfields = False
if Gooddays_deepfields :
    df2= read_gooddays(infile=good_nircammiri, base_dir="", addGalEcl=True)
    thefiles = ("commissioning_targets.txt", "deep_fields.txt")
    for infile in thefiles :
        deep_fields = read_deepfields(coords_dir + infile, addGalEcl=True)
        deep_fields.set_index('healpix_asindex', inplace=True)
        subset = df2.ix[deep_fields.index]  # Grab subset of big dataframe for the healpix vals of the deep_fields
        subset.insert(0, 'FIELD', deep_fields['FIELD'])
        subset.reset_index(inplace=True)
        subset2 = subset.reindex(subset.Ecl_lat.abs().order().index)  # sorts by abs val of Ecliptic latitude
        subset2.drop('Unnamed: 0', axis=1, inplace=True)
        subset2.drop('healpix_asindex', axis=1, inplace=True)
        subset2.drop('healfile', axis=1, inplace=True)
        subset2.drop(['RA_deg', 'DEC_deg'], axis=1, inplace=True)
        subset2.transpose().to_csv(re.sub(".txt", '_gooddays.csv', infile))  # dump to csv
        subset2.drop('Ecl_lon', axis=1, inplace=True)
        print subset2.transpose().to_clipboard()  # write to clipboard, so I can copy to the memo

# Now, make 2D maps on the sky, and bin by Ecliptic latitude
Plot_on_sky = False
if Plot_on_sky :
    plt.clf()
    df2= read_gooddays(infile=good_nircammiri, base_dir="", addGalEcl=False)
    cmap = sns.cubehelix_palette(n_colors=10,  start=0, rot=0.0, gamma=2.0, hue=1, light=1, dark=0.4, reverse=False, as_cmap=True)
    healpy.mollview(df2['Nday'], cmap=cmap, flip='astro', title="Ndays observable", coord='CE')
    plt.show()   
    df2= read_gooddays(infile=good_nircammiri, base_dir="", addGalEcl=True)
    binlatby = 1
    latbins = np.arange(-90,90,binlatby)
    groups = df2.groupby(pandas.cut(df2['Ecl_lat'], latbins))
    groups.median().to_csv('gooddays_binned_by_ecliptic_latitude.csv')
    sns.set_palette("hls", 11)
    plt.plot(latbins[:-1], groups['Good21.5_1.1'].median().values, label="21.5 micron")
    plt.plot(latbins[:-1], groups['Good17.5_1.1'].median().values, label="17.5 micron")
    plt.plot(latbins[:-1], groups['Good15.1_1.1'].median().values, label="15.1 micron")
    plt.plot(latbins[:-1], groups['Good12.7_1.1'].median().values, label="12.7 micron")
    plt.plot(latbins[:-1], groups['Good10.1_1.1'].median().values, label="10.1 micron")
    plt.plot(latbins[:-1], groups['Good5.5_1.1'].median().values, label=" 5.1 micron") 
    plt.plot(latbins[:-1], groups['Good4.5_1.1'].median().values, label="4.5 micron") 
    plt.plot(latbins[:-1], groups['Good2.0_1.1'].median().values, label="2.0 micron") 
    plt.plot(latbins[:-1], groups['Good1.1_1.1'].median().values, label="1.1 micron") 
    plt.plot(latbins[:-1], groups['Good0.7_1.1'].median().values, label="0.7 micron")
    plt.plot(latbins[:-1], groups['Nday'].median().values, label="Nday", color='black')
    plt.xlabel("Ecliptic Latitude (deg)")
    plt.ylabel("N days with background <1.1 of min")
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plt.xlim(-90,90)
    plt.locator_params(axis='x', nbins=6)
    plt.ylim(0,365)
    plt.legend(labelspacing=0.2, fontsize=12)
    plt.grid()
    plt.savefig("gooddays_vs_eclipticlatitude.png")
    plt.show()

