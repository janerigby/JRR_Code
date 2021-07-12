from astropy.cosmology import WMAP9 as wmap_cosmo  # example in the astropy documentation
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import pandas
import jrr

'''Given deflection maps from Keren Sharon's gravitational lensing models, compute the separations of the 
multiple sight-lines to the intervening absorber, at the redshift of the intervening absorber.  Using
code developed and tested in lensing_deflections_A2744_RCS0327.py.
Using jrr/lens.py (lensing calculations) and jrr/util.py (sky coords).
jrigby, March 2021
'''

def load_RCS0327_pos():
    infile = 'RCS0327_test_table.txt'
    names = ('position', 'RA_seg', 'DEC_seg')
    df = pandas.read_csv(infile, delim_whitespace=True, names=names, comment='#')
    jrr.util.convert_RADEC_segidecimal_df(df, colra='RA_seg', coldec='DEC_seg', newRA='RA', newDEC='DEC')
    return(df)

def get_Megasaura_coordinates(infile=None) :
    if infile==None:
        infile = '/Users/jrrigby1/SCIENCE/Lensed-LBGs/Observing/Where-pointed-MagE/MagE_master_coordinates/megasaura_andfriends_coordinates.txt'
    df = pandas.read_csv(infile, delim_whitespace=True, comment='#')
    return(df)

def get_intervening_absorber_multsightlines():
    infile='Absorber_list_multiple_sightline_v2.dat'
    df = pandas.read_csv(infile, delim_whitespace=True, comment='#')
    return(df)

def get_lensmodel_info(label) :
    if label == 'RCSGA032727-132609' :
        datadir='/Users/jrrigby1/SCIENCE/Lensed-LBGs/RCS0327/Lens_Model/'
        deflection_map= { 'x' : 'RCS0327_dplx_1.701.fits', 'y' : 'RCS0327_dply_1.701.fits'}
        z_lens = 0.564 # Wuyts et al. 2010
        # For RCS0327, the deflection maps are scaled for the giant arc(z_arc = 1.70), per KS email 01/08/2021
        deflection_norm = jrr.lens.scale_lensing_deflection(cosmo, z_lens, 1.701)

    elif label == "SGASJ122651.3+215220_for0.77" :
        datadir = '/Volumes/Lens_models/SDSS1226/'
        deflection_map = { 'x' : 'dplx_0.77_abs.fits',   'y' : 'dply_0.77_abs.fits'}
        z_lens = 0.43  # Koester et al. 2010
        deflection_norm = jrr.lens.scale_lensing_deflection(cosmo, z_lens, 0.77)  #  Should be 0.77, checked w Keren on Slack 

    elif label == "SGASJ122651.3+215220" :
        datadir = '/Volumes/Lens_models/SDSS1226/'
        deflection_map = { 'x' : 'dplx_2.9233.fits',   'y' : 'dply_2.9233.fits'}
        z_lens = 0.43  # Koester et al. 2010
        deflection_norm = jrr.lens.scale_lensing_deflection(cosmo, z_lens, 2.9233)
    
    elif label == "SGASJ152745.1+065219" :
        datadir = '/Volumes/Lens_models/Sharon_2020/MAST-SGAS-UPLOAD/sdssj1527p0652/'
        deflection_map = { 'x' : 's1527_dplx_2.761.fits',   'y' : 's1527_dply_2.761.fits'}       
        z_lens = 0.39  # Koester et al. 2010
        deflection_norm = 1.0  # most models are scaled to ds/dL=1

    elif label == "SunburstArc" :
        datadir = '/Volumes/Lens_models/Sunburst_Arc/'
        deflection_map = {'x' : 'dplx1.fits', 'y' : 'dply1.fits'}
        z_lens = 0.443  #Dahle et al. 2016
        deflection_norm = 1.0  # most models are scaled to ds/dL=1
        
    else : raise Exception("ERROR: Unrecognized lens model label", label)
    return(datadir, deflection_map, z_lens, deflection_norm)


# GENERAL SETUP
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)  # Keren's favorite cosmology.
coords_df = get_Megasaura_coordinates() 
abs_df =  get_intervening_absorber_multsightlines()
unit=(u.hourangle, u.deg)

for index, row in abs_df.iterrows() : # Loop over each intervening absorber
    print(row.cluster, row.redshift)
    (datadir, deflection_map, z_lens, deflection_norm) = get_lensmodel_info(row.cluster)  # Get the lens model info
    coords_df['zabs'] = row.redshift
    subset_coords = coords_df.loc[coords_df['system'] == row.cluster].copy(deep=True)              # get the positions within that arc
    subset_coords['zabs'] = row.redshift # bookkeeping    #change below
    prefix = row.cluster + '_z' + str(row.redshift) + '_intervabs'
    jrr.lens.compute_sourceplane_positions(subset_coords, prefix, z_lens, datadir, deflection_map, cosmo, deflection_norm=deflection_norm, debug=False, unit=unit)
    #results get dumped to subset_coords
    print(subset_coords.head(100))

    # Need to do this relative to the first item in the table, grab it  Import it
    #jrr.lens.compute_sourceplane_offsets(abs_df, fiducialRA, fidicualDEC)
    
'''    
# RCS0327 as a test case
arcs_df       = load_RCS0327_pos()
abs_zz        = 0.982920
arcs_df['zz'] = abs_zz 
(datadir, deflection_map, z_lens, deflection_norm) = get_lensmodel_info('RCS0327')
RCS0327_df = jrr.lens.compute_sourceplane_positions(arcs_df, 'RCS0327', z_lens, datadir, deflection_map, cosmo, deflection_norm=deflection_norm, debug=False)
# Now, need to interpret these distances.  RA_srcplane DEC_srcplane are the src plane of the *intervening absorber*.
RCS0327_df['angscale'] = cosmo.kpc_proper_per_arcmin(abs_zz).to(u.kpc / u.arcsec).value
# Now, calculate the distance from fiducial point G1 from Lopez et al. paper, to each position:
jrr.util.distancefrom_coords_df(arcs_df, 51.865374, -13.438874, 'RA_srcplane', 'DEC_srcplane', newcol='offset_arcsec')
# Checking against Lopez paper.  Agrees pretty well with their figures.  I'm happy with RCS0327
RCS0327_df['offset_kpc'] = RCS0327_df['angscale'] * RCS0327_df['offset_arcsec']
print(RCS0327_df.head())
'''



''' OK, I think RCS0327 is working above.  Need to repeat for other multiple sight-lines.  Breaking up into small tasks:
DONE - make file with coordinates for every multiple sightline in each system   megasura_andfriends_coordinates.txt
DONE - make file of redshifts to check for each system, Absorber_list_multiple_sightline_v2.dat
- loop through each system  (There are only 5.  Have no lens model for PSZ0441; RCS0327 check, S1226 need, S1527 need, Sunburst check)
    - find the right lens model
    - Get the missing lens models from Keren.
    - loop through each redshift
    - calculate the offsets relative to a fiducial (the first position listed?)
    - print the output
    - send the output to Rongmon to make a plot.
'''



# Should also calculate, for all intervening absorbers, the size of the object in the slit, and how big that is in the intervening absorber frame.
