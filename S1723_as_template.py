import numpy as np
import pandas
#from astropy.io.ascii import read
from astropy import constants
import re
import jrr  # Jane's routines
import matplotlib.pyplot as plt

def pretty_plot() :
    plt.xlabel("wavelength (micron)", fontsize=22)
    #plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()

def plot_synthetic_spectrum_flam(zz, sim_df, objname) :
    fig, ax1 = plt.subplots()
    plt.title("Synthetic spectrum for " + objname)
    ax1.set_ylabel("observed flambda (erg/s/cm^2)")
    ax1.set_xlabel("observed wavelength (microns)")
    ax1.plot( sim_df['wave_um'] , sim_df['flam'])
    xlims = (0.6, 5) #observed
    ax1.set_xlim(xlims[0], xlims[1])
    ax2 = ax1.twiny()
    ax2.set_xlim( xlims[0] / (1+zz), xlims[1]/(1+zz))
    ax2.set_xlabel('rest wavelength (micron)')
    ax1.set_ylim(0,3E-17)
    return(0)

def plot_linefluxes_and_sens(zz, detected_lines, objname) :
    df_sens = pandas.read_json('pandeia_lineflux_sensitivities_v1.5.0.json')
    # Column PSSL is the limiting lineflux sensitivity for a spectrally unresovled line in  a point source at SNR=10 in 1E4  From Pandeia 1.5
    fig3, ax3 = plt.subplots()
    ax3.stem(detected_lines['obs_wave']/1E4, detected_lines['scaled_flux'], use_line_collection=True)
    # add the curves from pandeia
    for index, row in df_sens.iterrows() :     # Go row by row, for each grating, order, etc
        plt.plot(row.wavelengths, row.PSSL)
    plt.yscale('log')
    plt.ylabel('line flux, erg/s/cm^2')
    plt.xlabel("observed wavelength (micron)")
    xlims = (0.6, 5) #observed
    plt.xlim( xlims[0], xlims[1])
    plt.title(objname + " line fluxes. Stems above sensitivity curves are detected at SNR>10 in 1E4s")
    ax4 = ax3.twiny()
    ax4.set_xlim( xlims[0] / (1+zz), xlims[1]/(1+zz))
    ax4.set_xlabel('rest wavelength (micron)')
    return(0)

def test_flux_conservation(df_sim) :
    print("Total flux in lines in rest-frame: ", '{:0.3e}'.format((df_sim.rest_flam * df_sim.d_restwave).sum()), "erg/s/cm^2")
    print("Total flux in lines in obs-frame: ", '{:0.3e}'.format((df_sim.flam * df_sim.d_wave).sum()), "erg/s/cm^2")
    return(0)

def write_simulated_spectrum(outfilename, df_sim, zz, scaleby) :
    df_sim.to_csv('temp', columns=('wave_um', 'fnu_mJy'), sep=' ', index=False)
    header_text = "# Simulated spectrum at redshift zz=" + str(np.round(zz, 5)) + ", scaling from S1723 fluxes by a factor of " + str(np.round(scaleby, 5)) + "\n"
    jrr.util.put_header_on_file('temp', header_text, outfilename)
    return(0)

##################################################

plt.close("all")

# Simulate the z=5 COOL LAMPS arc
zz_z5 = 5.0
# For reference,    make_synthetic_spectrum_based_on_S1723 returns ( df_s1723, detected_lines, df_sim)
(dum1, dum2, dum3) = make_synthetic_spectrum_based_on_S1723(zz_z5, fwhm_kms=10, scaleby=1.0)  # Run this once, just to get the df_s1723 spectrum to compute the real scaleby
O2flux_s1723 = dum2[dum2['wave'].between(3725,3732)]['flux'].sum()   # Scale the line fluxes for the z=5 arc, linearly from flux of [O II] sum
O2flux_z5  = 4.5E-16  # From Gourav, for z=5 source
# Jane: the 16--84 percentile range is 2.48 to 6.57 X 10^-16 erg /s/cm^2/ w the median 4.5E-16
scaleby_z5 = O2flux_z5 / O2flux_s1723 / 30.  # Scale to flux of z=5 source, and split into 30 pieces
(df_s1723, z5_detected_lines, z5_sim) = make_synthetic_spectrum_based_on_S1723(zz_z5, fwhm_kms=10, scaleby=scaleby_z5)  # run again, but this time with scaleby
outfilename = "sim_spectrum_z" + str(np.round(zz_z5, 2)) + "_v1.1.txt"
write_simulated_spectrum(outfilename, z5_sim, zz_z5, scaleby_z5)

# Simulate S1110 at z=2.481.  Simply scale S1110 line fluxes from apparent magnitudes
traci_tab4 =  jrr.sgas.load_s1110_Johnson_Tab4()  # load table 4 of Traci's paper.  Somethings weird about the fluxes. Ignore for now.
zz_s1110 = 2.481                 # redshift from Traci's paper
mAB_606 = 27.5 # 27. 28.  # AB magnitude apparent 606, from Clmp_brightness/notes: mAB(606) stats:  median 27.5, quartiles are 27.0 and 28.0, mode is ~27.25


Lnu = jrr.util.mAB_to_Lnu(mAB_606, zz_s1110)
SFR = jrr.util.Kennicutt_LUV_to_SFR(Lnu)
fHa  = jrr.util.Kennicutt_SFR_to_LHa(SFR) / (4. * np.pi * jrr.util.luminosity_distance(zz_s1110)**2)
fHB = fHa / 3.0
scaleby_s1110 =  fHB / df_s1723[df_s1723['wave'].between(4860, 4870)]['flux'].sum()
(dummy3, s1110_detected_lines, s1110_sim) =  make_synthetic_spectrum_based_on_S1723(zz_s1110, fwhm_kms=10, scaleby=scaleby_s1110)  
outfilename = "sim_spectrum_z" + str(np.round(zz_s1110, 2)) + "_v1.1_median.txt"
write_simulated_spectrum(outfilename, s1110_sim, zz_s1110, scaleby_s1110)


#plot_synthetic_spectrum_flam(zz_z5, z5_sim, "z5")
plot_linefluxes_and_sens(zz_z5, z5_detected_lines, 'z5')

#plot_synthetic_spectrum_flam(zz_s1110, s1110_sim, 's1110')
plot_linefluxes_and_sens(zz_s1110, s1110_detected_lines, 's1110')

print("Testing flux conservation.  Total lineflux in obs, rest frames should match:")
test_flux_conservation(z5_sim)
test_flux_conservation(s1110_sim)



plt.show()
