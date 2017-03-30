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

import glob
from os.path import basename
import struct
import numpy as np
import pandas
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

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
        plt.plot(wave_array, total[thisday, :], label = "Total")
        plt.xlim(0.6,31)
        plt.legend()
        plt.show()
    return((calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total))  #send it as a tuple

def index_of_wavelength(wave_array, desired_wavelength) :  # look up index of wavelength array corresponding to desired wavelength
    the_index = np.where(wave_array == desired_wavelength)
    return(the_index[0][0])

def make_bathtub(results, wavelength_desired=2.0, thresh=1.1, showplot=False) :
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
thresh = 1.1  # 10 percent  # Threshhold above minimum background, to calculate Ndays
wavelength_desired = 2.0
base_dir  = "/Volumes/Apps_and_Docs/MISSIONS/JWST/Zody_bathtubs/"  # Satchmo
bkg_dir   = base_dir + "sl_cache.v1.0/"


# Tutorial, make one example bathtub 
myfile = "1464/sl_pix_146447.bin"
results = read_JWST_precompiled_bkg(myfile, base_dir, bkg_dir, showplot=False)
(calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results 
allgood = make_bathtub(results, wavelength_desired, thresh, showplot=False)
print "Ran", myfile, wavelength_desired, "micron", thresh, "threshold"
print allgood, "good days out of", len(calendar)


# Loop through all the directories
healpix_dirs = glob.glob(bkg_dir + "*/")
allNday   =  []
allGood   = []
whichfile  = []

for thisdir in healpix_dirs[0:3] : 
    myfiles = [ basename(x) for x in glob.glob(thisdir + "*.bin") ]
    print len(myfiles)
    for ii, myfile in enumerate(myfiles) :
        results = read_JWST_precompiled_bkg(thisdir + myfile, base_dir, thisdir, showplot=False)
        (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results  # break up results tuple
        allNday.append(len(calendar))
        allGood.append(make_bathtub(results, wavelength_desired, thresh, showplot=False))
        whichfile.append(myfile)
allNday   = np.array(allNday)
allGood   = np.array(allGood)
whichfile = np.array(whichfile, dtype='str')
plt.scatter(allNday, allGood)
plt.show()

df = pandas.DataFrame(

header = "Number of Days in FOR, and number of days with good background, for thresh" + str(thresh) + "at wavelength" + str(wavelength_desired)
temp =  np.transpose([whichfile, allNday, allGood], dtype=(np.int, np.int, np.str))
np.savetxt("gooddays.txt", np.transpose([whichfile, allNday, allGood]), "%25s  %d  %d", header=header, comments="")
