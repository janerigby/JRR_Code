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

import struct
import numpy as np
from matplotlib import pyplot as plt

myfiles = ('sl_pix_000593.bin','sl_pix_041303.bin','sl_pix_099422.bin','sl_pix_157139.bin','sl_pix_196512.bin')
myfile =  myfiles[0]

def read_JWST_precompiled_bkg(infile, showplots=False) :
    thedir = "/Volumes/Apps_and_Docs/MISSIONS/JWST/Zody_bathtubs/bg_samples_for_jane/"  # Satchmo
    thedir = "/Users/jrrigby1/MISSIONS/JWST/Zody_bathtubs/bg_samples_for_jane/" # Milk
    wave_file = "updated_std_spectrum_wavelengths.txt"  # Standard wavelength array.  Should be SL_NWave=108 long
    wave_array = np.loadtxt(thedir + wave_file)    
    SL_NWAVE = len(wave_array)  # should be 108.  Size of wavelength array
    sbet_file = open(thedir + myfile)
    sbet_data = sbet_file.read()

    # Unpack the constant first part
    print "File has", len(sbet_data), "bytes, which is", len(sbet_data)/8., "doubles"
    size_calendar = struct.calcsize("366i") # bytes, not doubles
    partA = struct.unpack(str(5 + SL_NWAVE)+'d', sbet_data[0: (5 + SL_NWAVE)*8])
    RA = partA[0]
    DEC = partA[1]
    pos = partA[2:5]
    nonzodi_bg = partA[5:5+SL_NWAVE]

    # Unpack the calendar dates      # code goes from 0 to 365 days.
    date_map = np.array(struct.unpack('366i', sbet_data[(5 + SL_NWAVE)*8  : (5 + SL_NWAVE)*8 + size_calendar]))
    print "Out of", len(date_map), "days, these many are legal:", np.sum(date_map >=0)
    #print "indices of days:", date_map[date_map>=0]
    calendar = np.where(date_map >=0)[0]
    #print "calendar date:", calendar
    # So, the index dd in zodi_bg[dd, : ]  corresponds to the calendar day lookup[dd]
    Ndays = len(calendar) 
    print len(date_map), Ndays

    # Unpack part B, the time-variable part
    zodi_bg        = np.zeros((Ndays,SL_NWAVE))
    stray_light_bg = np.zeros((Ndays,SL_NWAVE))
    perday = SL_NWAVE*2
    partB= struct.unpack(str((len(calendar))*SL_NWAVE*2)+'d', sbet_data[perday*Ndays*-8 : ])

    print "Ndays:", Ndays
    for dd in range(0, int(Ndays)):
        br1 = dd*perday
        br2 = br1 + SL_NWAVE
        br3 = br2 + SL_NWAVE
        #print "Breaking at:", br1, br2, br3
        zodi_bg[dd, ]        = partB[br1 : br2]
        stray_light_bg[dd, ] = partB[br2 : br3]
    #print br3, "was last double" 
    #print "DEBUGGING", zodi_bg
    #print "DEBUGGING", stray_light_bg

    if showplots :
        plt.clf()
        thisday = 100
        plt.plot(wave_array, nonzodi_bg, label="ISM")
        plt.plot(wave_array, zodi_bg[thisday, :], label="Zodi")
        plt.plot(wave_array, stray_light_bg[thisday, :], label="Stray light")
        plt.xlim(0.6,31)
        plt.legend()
        plt.show()
        
    twomicron = 16 # 2.0 micron is index 16 in wave_array
    total_2mic = nonzodi_bg[twomicron]*np.ones(Ndays) + stray_light_bg[ : , twomicron] + zodi_bg[ : , twomicron]
    themin = np.min(total_2mic)
    print "Fraction of FOR  <110% of min bkg", np.sum(total_2mic < themin*1.1)*1.0/ Ndays  # What fraction of days is below 10% threshold
    if showplots:
        plt.clf()      # Plot bathtub curve for 2um
        plt.scatter(calendar, total_2mic)
        percentiles = (themin, themin*1.1)
        plt.hlines(percentiles, 0, 365, color='green')
        plt.xlabel("Day of the year")
        plt.ylabel("bkg at " + str(wave_array[twomicron]) + " um (MJy/SR)")
        plt.show()
    
    return(calendar, RA, DEC, pos, nonzodi_bg, zodi_bg, stray_light_bg)
    

# actually run and test this.
for myfile in myfiles :
    print myfile
    (calendar, RA, DEC, pos, nonzodi_bg, zodi_bg, stray_light_bg) = read_JWST_precompiled_bkg(("bg_samples_for_jane/" + myfile),  showplots=False)
