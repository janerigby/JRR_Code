''' quick test to get a bpass model loaded and fit
jrigby oct 2015 '''

import numpy as np
from  astropy.io import ascii
import matplotlib.pyplot as plt
import jrr

CIV_range  = (1551.5, 1560.0)  # this is new; remake all plots
# What should the He II range be?
HeII_range = (1630.0, 1650.0)
scale = 1E17  # for plotting, multiply by scale to get y axis ~ unity



def measure_the_EWs(rest_wave, rest_flam, rest_flam_u, rest_cont, rest_cont_u, short_label, zz, counter, out) :
    civ  = np.where(  (rest_wave > CIV_range[0])  & (rest_wave < CIV_range[1])  )   # indices 
    heii = np.where(  (rest_wave > HeII_range[0]) & (rest_wave < HeII_range[1])  )           
    # Calculate dispersion.  Kludgy but works
    dciv  = np.average(np.array([  (rest_wave[civ][1] - rest_wave[civ][0]),   (rest_wave[civ][-1] - rest_wave[civ][-2])  ]))
    dheii = np.average(np.array([  (rest_wave[heii][1] - rest_wave[heii][0]),   (rest_wave[heii][-1] - rest_wave[heii][-2])  ]))
    # Calculate EW
    (EW_CIV, uncert_EWCIV)   = jrr.calc_EW( rest_flam[civ],  rest_flam_u[civ],  rest_cont[civ],  rest_cont_u[civ],  dciv, 0.0)
    (EW_HeII, uncert_EWHeII) = jrr.calc_EW( rest_flam[heii], rest_flam_u[heii], rest_cont[heii], rest_cont_u[heii], dheii, 0.0)
    out.write("%30s    %5.3f      %5.3f      %5.3f      %5.3f\n" % (short_label, EW_CIV, uncert_EWCIV, EW_HeII, uncert_EWHeII))
    return(0)


def measure_EWs_for_a_Bpass_model(infile, outfile) :
    out = open(outfile, 'w')
    out.write("#t(yr)    EW(CIV,A)   dEW     EW(HeII,A),  dEW\n")
    print "Reading bigass bpass model", infile
    s = ascii.read(infile, comment="#")   # input has units of **flam**
    wave = s['col1']  # wave for all times

    for nn in range(2,40) :
        colname = "col" + str(nn)
        flam = s[colname]  
        flam_u = np.zeros(shape=flam.shape) 
        time = 10**(6+0.1*(nn -2))
        # Make a simple linear continuum fit from clean regions for CIV and HeII
        x1 = 1475.
        x2 = 1690.
        dx = 10.
        ind1 = np.where((wave > (x1-dx)) & (wave < (x1+dx)) )  # clean continuum on blue end
        ind2 = np.where((wave > (x2-dx)) & (wave < (x2+dx)) )  # clean continuum on red end
        ind3 = np.where((wave > (x1-dx)) & (wave < (x2+dx)) )  # span over CIV + HeII + continuum
        y1 = np.average(flam[ind1])
        y2 = np.average(flam[ind2])
        slope = (y1-y2)/(x1-x2)
        b = y1 - slope*x1
        plt.plot(wave[ind3], flam[ind3]/y1, color="black")
        plt.plot(wave[ind3], (slope*wave[ind3] + b)/y1, color="green")
        plt.show()
        pretty_time = "%.3E" % (time)
        measure_the_EWs(wave[ind3], flam[ind3], flam_u[ind3], slope*wave[ind3]+b, 0.0*wave[ind3], pretty_time, 0.0, 1, out)
    out.close()
    return(0)
    
bpassfile = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Mage/Models/BPASS/BPASSv2_imf135_100/OUTPUT_CONT/spectra.z020.dat.gz"
outfile = "BPASSv2_imf135_100_CONT_z020_EWs.dat"
measure_EWs_for_a_Bpass_model(bpassfile, outfile)
