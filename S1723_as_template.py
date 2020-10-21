import numpy as np
import pandas
from astropy.io.ascii import read
from astropy.stats import gaussian_fwhm_to_sigma
from astropy import units
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

def load_s1723_lineflux_template(infile='tab_s1723_measured_emission_v6.csv') :
    # The infile is in CSV format.  It was created by converting Table 1 of Rigby et al. 2020 Table 1 (in latex format)
    # to csv, using JRR_Code/s1723_convert_latex_table_csv.py.
    # JRR copied it from the S1723 paper to here (in Dropbox/SGAS/) so team has access.
    #
    # Need to change the table a bit before it can be used -- it has S_II blended rather than the indy lines,
    # and it double-counts C~III] [O III]
    # The SII doublet was blended at grism resoln. Split it, with the flux evenly distributed.
    names = ['lineID', 'wave', 'spectrograph', 'Wr', 'Wr_u', 'significance', 'flux', 'flux_u', 'dr_flux', 'dr_flux_u']  #column names
    df3 = pandas.read_csv(infile, comment='#', index_col=0)
    df3.drop(index=68, inplace=True)  # Drop the SII row
    df_SII = pandas.read_csv('fake_SII.txt', names= names + ['detected',])
    df_SII['detected'] =  df_SII['detected'].astype(bool)
    df4 = pandas.concat([df3, df_SII])
    df4.reset_index(inplace=True)
    df4.drop('index', axis=1, inplace=True)
    df4 = df4.reindex()

    # Now that df is well-behaved, convert all eligible columns to float64
    floatcols = ['wave', 'Wr', 'Wr_u', 'significance', 'flux', 'flux_u', 'dr_flux', 'dr_flux_u']
    for thiscol in floatcols :
        df4[thiscol] = df4[thiscol].astype(float)

    flux_scaling = 1.0E-17  # The S1723 table's fluxes are in units of 10^-17 erg/s/cm^2
    fluxcols = ['flux', 'flux_u', 'dr_flux', 'dr_flux_u']
    for thiscol in fluxcols :
        df4[thiscol] *= flux_scaling
    return(df4)

def load_s1110_Johnson_Tab4() :
    zz_s1110 = 2.481  
    johnson_tab4 = pandas.read_csv("Johnson_2017a_Table4.txt", delim_whitespace=True, comment="#")
    # The magnitudes in Tab 4 of Johnson et al 2017a are too faint by 7.5 mag (factor of 1000, due to a scaling factor.)
    fix_mag_factor =  -2.5*np.log10(1000.)      # Not sure this is right.  Resulting mags are crazy
    fix_mag_keys = ('m_F606W', 'm_F390W')
    for thiskey in fix_mag_keys :   johnson_tab4[ thiskey] += fix_mag_factor
    # convert (no longer screwewd up) m_AB(F606W) to flux density, correct for 1+z, and convert UV fnu to SFR following Kennicutt
    johnson_tab4['fnu606'] = johnson_tab4['m_F606W'].map(lambda x: jrr.util.mAB_to_fnu(x))    # convert AB mags to observed-frame fnu
    johnson_tab4['Lnu']    =  johnson_tab4['fnu606'].map(lambda x: jrr.util.fnu_to_Lnu(x, zz_s1110))  # convert obs fnu to rest flambda
    #johnson_tab4['fnu606'] = 10**((johnson_tab4['m_F606W'] + 48.60)/-2.5)  # cgs units
    #johnson_tab4['Lnu'] = johnson_tab4['fnu606'] / (1 + zz_s1110)  * 4. * np.pi  *  dL_s1110**2
    johnson_tab4['SFR'] = johnson_tab4['Lnu'].map(lambda x: jrr.util.Kennicutt_LUV_to_SFR(x))  # Eqn 1 of Kennicutt 1998, in Msol/yr
    return(johnson_tab4)

def make_synthetic_spectrum_based_on_S1723(zz, fwhm_kms=100, scaleby=1.0) :  # fwhm_kms is assumed line width, fwhm, in km/s.  scaleby is factor to scale fluxes
    # Convenient wrapper function to load the S1723 template, and use it to generate a simulated spectrum
    # Load the S1723 lineflux measurements into a well-behaved pandas dataframe.  Uses the local copy which has duplicate line measurements commented out
    s1723_df = load_s1723_lineflux_template() 
    c_kms = constants.c.to('km/s').value   # Someday I'll trust astropy.units
    s1723_df['sigma'] = (fwhm_kms  / c_kms ) * gaussian_fwhm_to_sigma * s1723_df['wave']  # linewidths to use
    detected_lines =  s1723_df[s1723_df['detected']]             
    print("Will scale the fluxes of S1723 by a factor of", np.round(scaleby, 4))

    # Make a synthetic spectrum based on the fluxes of S1723, redshifted and scaled
    startwave = 1300.  # Angstroms
    endwave   = 7000.  # Angstroms
    Npix      = 4.0E4   #kludge
    df_sim = pandas.DataFrame(data=pandas.Series(np.geomspace(startwave, endwave, Npix)), columns=('rest_wave',))
    df_sim['rest_flam'] = 0.0     # no continuum to start
    df_sim['rest_flam_u'] = 0.0   # no uncertainty either
    df_sim['d_restwave'] = df_sim['rest_wave'].diff()  # dlambda array, the width of 1 pixel
    df_sim['d_restwave'].iloc[0] = df_sim['d_restwave'].iloc[1] # first value is bogus
    for row in detected_lines.itertuples() :  # Make rest frame synthetic spectrum.  Add a gaussian for each detected line, at the right flux
        amplitude = scaleby * row.flux / row.sigma / np.sqrt(2.*np.pi)  
        y_thisgauss = jrr.spec.onegaus(df_sim['rest_wave'],  amplitude, row.wave, row.sigma, 0.0)  # array to fill, AMPLITDE, central wavelength, sigma, continuum
        df_sim['rest_flam'] +=  (y_thisgauss / df_sim['d_restwave'])  # per Angstrom

    detected_lines['obs_wave']    = detected_lines['wave'] * (1.0 + zz)
    detected_lines['scaled_flux'] = detected_lines['flux'] * scaleby   # Scale lineflux from above.  No (1+z), b/c flux is conserved
    jrr.spec.convert2obsframe_df(df_sim, zz, units='flam', colrw='rest_wave', colrf='rest_flam', colrf_u='rest_flam_u')
    df_sim['d_wave'] = df_sim['d_restwave'] * (1+zz)
    df_sim['wave_um'] = df_sim['wave']/1E4
    df_sim['fnu_mJy'] = jrr.spec.flam2fnu(df_sim['wave'], df_sim['flam']) /1E-23 *1E3 # convert cgs fnu to mJy
    return( s1723_df, detected_lines, df_sim)

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

##################################################
plt.close("all")

# Simulate the z=5 COOL LAMPS arc
zz_z5 = 5.0
(s1723_df, dummy, dummy) = make_synthetic_spectrum_based_on_S1723(zz_z5, fwhm_kms=100, scaleby=1.0)  # Run this once, just to get the s1723_df spectrum to compute the real scaleby
O2flux_s1723 = s1723_df[s1723_df['wave'].between(3725,3732)]['flux'].sum()   # Scale the line fluxes for the z=5 arc, linearly from flux of [O II] sum
O2flux_z5  = 4.5E-16  # From Gourav, for z=5 source
# Jane: the 16--84 percentile range is 2.48 to 6.57 X 10^-16 erg /s/cm^2/ w the median 4.5E-16
scaleby_z5 = O2flux_z5 / O2flux_s1723 / 30.  # Scale to flux of z=5 source, and split into 30 pieces
(s1723_df, z5_detected_lines, z5_sim) = make_synthetic_spectrum_based_on_S1723(5.0, fwhm_kms=100, scaleby=scaleby_z5)  # run again, but this time with scaleby
outfilename = "sim_spectrum_z" + str(np.round(5.0, 2)) + "_v1.1.txt"
z5_sim.to_csv(outfilename, sep=' ', columns=('wave_um', 'fnu_mJy'), index=False)


# Simulate S1110 at z=2.481.  Simply scale S1110 line fluxes from apparent magnitudes
s1110_df =  load_s1110_Johnson_Tab4()  # load table 4 of Traci's paper.  Somethings weird about the fluxes.
zz_s1110 = 2.481                 # redshift from Traci's paper
mAB_606 = 26.6  # AB magnitude apparent 606, from Matt on Slack
Lnu = jrr.util.mAB_to_Lnu(mAB_606, zz_s1110)
SFR = jrr.util.Kennicutt_LUV_to_SFR(Lnu)
fHa  = jrr.util.Kennicutt_SFR_to_LHa(SFR) / (4. * np.pi * jrr.util.luminosity_distance(zz_s1110)**2)
fHB = fHa / 3.0
scaleby_s1110 =  fHB / s1723_df[s1723_df['wave'].between(4850, 4870)]['flux'].sum()
(dummy2, s1110_detected_lines, s1110_sim) =  make_synthetic_spectrum_based_on_S1723(zz_s1110, fwhm_kms=100, scaleby=scaleby_s1110)  # run again, but this time with scaleby
outfilename = "sim_spectrum_z" + str(np.round(zz_s1110, 2)) + "_v1.1.txt"
s1110_sim.to_csv(outfilename, sep=' ', columns=('wave_um', 'fnu_mJy'), index=False)



plot_synthetic_spectrum_flam(zz_z5, z5_sim, "z5")
plot_linefluxes_and_sens(zz_z5, z5_detected_lines, 'z5')

plot_synthetic_spectrum_flam(zz_s1110, s1110_sim, 's1110')
plot_linefluxes_and_sens(zz_s1110, s1110_detected_lines, 's1110')




test_flux_conservation = True
if test_flux_conservation :
    print("Total flux in lines in rest-frame was", '{:0.3e}'.format((z5_sim.rest_flam * z5_sim.d_restwave).sum()), "erg/s/cm^2")
    print("Total flux in lines in obs-frame  was", '{:0.3e}'.format((z5_sim.flam * z5_sim.d_wave).sum()), "erg/s/cm^2")
    print("These should match because flux conservation.")

plt.show()
