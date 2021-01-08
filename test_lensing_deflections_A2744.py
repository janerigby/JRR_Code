from astropy.cosmology import WMAP9 as wmap_cosmo  # example in the astropy documentation
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import pandas
import jrr

'''Given deflection maps from Keren Sharon's gravitational lensing models, compute the source-plane positions
of the lensed arcs.  The main test that this code works is that each multiply imaged family gets pushed
to the same source-plane position within a few 0.1".  Checked with Keren, it appears to work.  jrigby, Jan 2021
The functions I wrote to solve this problem are in jrr/lens.py (lensing calculations) and jrr/util.py (sky coords)
'''

def open_lensfamily_table(datadir) :
    # Read Table 3, a list of lens family coords and redshifts, from Keren's Readme about Abell 2744 lens model
    infile = datadir + "tab3.txt"
    df = pandas.read_csv(infile, delim_whitespace=True, comment='#')
    return(df)

def make_lensed_regions_file(df_in, racol, deccol, color='green', suffix='reg') :  # Helper function stays in this file
    # Make a regions file for the image plane positions.  Input is a dataframe containing RA, DEC, ID of arcs
    df_reg = df_in.copy(deep=True)
    df_reg['text'] = 'text={' + df_reg['ID'] + ' ' + suffix  + '}'
    df_reg['color'] = 'color=' + color
    df_reg['radius'] = '1.5'
    jrr.util.make_ds9_regions_file('A2744_' + suffix + '.reg', df_reg, racol=racol, deccol=deccol)
    return(df_reg)


#### Let's compute some stuff!
debug = False
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)  # Keren's favorite cosmology.
# test out with some defaults
z_lens = 0.5
z_src = 1.7
jrr.lens.lensing_distances(cosmo, z_lens, z_src, verbose=True)
sigma_critical = jrr.lens.critical_density(cosmo, z_lens, z_src)
print("Critical density was", sigma_critical.decompose())
print("Done testing functions.  Let's try an example image family in A2744.\n")

datadir = '/Users/jrrigby1/SCIENCE/Lensed-LBGs/How_lensing_works/FF_A2744_Sharon_v4cor/' 
deflection_map = { 'x' : "hlsp_frontier_model_abell2744_sharon_v4cor_x-arcsec-deflect.fits", 'y' : "hlsp_frontier_model_abell2744_sharon_v4cor_y-arcsec-deflect.fits"}    
outregfile = datadir + 'img_src_plane.reg'


#Running for all arcs in A2744.  The first distance in each family will be nonsense.
z_lens = 0.3080  # Abell 2744
prev_imgpos = SkyCoord(0.0, 0.0, unit='deg') # initalize
prev_srcpos = SkyCoord(0.0, 0.0, unit='deg') # initalize
print("Coordinates                                 distance from previous")
print("ID   image plane    src plane               d_img    d_src")
arcs_df = open_lensfamily_table(datadir)
df_imgplane_reg =  make_lensed_regions_file(arcs_df, 'RA', 'DEC', color='green', suffix='img') 

src_pos_list = []
for row in arcs_df.itertuples() :
    img_pos = SkyCoord(row.RA, row.DEC, unit='deg')
    z_src = row.zz
    thisdeflect = {}
    for key in deflection_map : #Loop through the X and Y deflection maps
        # Retrieve the deflection, and scale it
        defscale = jrr.lens.scale_lensing_deflection(cosmo, z_lens, z_src)
        imval, xy = jrr.util.imval_at_coords(datadir + deflection_map[key], img_pos)
        #print("Using", key, "deflection map", deflection_map[key])
        thisdeflect[key] = imval*defscale
        if debug:
            print("  ", row.ID, "  Deflection in", key, "at", xy, "(stupid numpy yx index):", np.round(imval, 2), end='  ')
            print("scale by", np.round(defscale,3), "to get ", np.round(thisdeflect[key],2))
    src_pos =  jrr.util.offset_coords(img_pos, delta_RA=thisdeflect['x'] *1, delta_DEC=thisdeflect['y']*-1) # minus sign here, per lensing eqn?
    src_pos_list.append(jrr.util.print_skycoords_decimal(src_pos)) # save this to list, then put it in the data frame after the loop.
    print(row.ID, jrr.util.print_skycoords_decimal(img_pos), jrr.util.print_skycoords_decimal(src_pos), end=' ')
    print(prev_imgpos.separation(img_pos).to_string(unit=u.arcsec, precision=2),  prev_srcpos.separation(src_pos).to_string(unit=u.arcsec, precision=2))
    prev_imgpos = img_pos
    prev_srcpos = src_pos

# Store the source-plane positions in the dataframe.
arcs_df['temp'] = src_pos_list
arcs_df[['RA_srcplane', 'DEC_srcplane']] = arcs_df['temp'].str.split(expand=True).astype('float')
arcs_df.drop(['temp',],  inplace=True, axis=1)
df_srcplane_reg =  make_lensed_regions_file(arcs_df, 'RA_srcplane', 'DEC_srcplane', color='red', suffix='src') 
    
        
# According to Keren's FF readme, the deflection maps are scaled to d_LS/d_S = 1
# Pseudo code for work to be done:
# - DONE    Convert from RA, DEC to pixel coordinates in the fits images.
# - DONE    Grab the values from the x, y deflection maps. 
# - DONE    Compute d_LS/d_S for the redshifts in question, and scale the deflections.
# - DONE   Go from image plane coords and deflection to source plane coords. 
# - DONE   Test that this works on all the multiply imaged galaxies in A2744.
#   DONE   Debugged, Keren found the problem, X was swapped (RA). It was the format of Keren's deflection maps.
# NEXT  Get Keren's lens models for the intervening absorbers
# NEXT  Compute source-plane coordinates for each of multiple sight-lines, at the redshift of each intervening absorber
