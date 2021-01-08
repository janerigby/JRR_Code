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

def compute_lensing_test() :
    # test out with some defaults
    z_lens = 0.5
    z_src = 1.7
    jrr.lens.lensing_distances(cosmo, z_lens, z_src, verbose=True)
    sigma_critical = jrr.lens.critical_density(cosmo, z_lens, z_src)
    print("Critical density was", sigma_critical.decompose())
    print("Done testing functions.  Let's try an example image family in A2744.\n")
    return(0)

def open_A2744_lensfamily_table(datadir) :
    # Read Table 3, a list of lens family coords and redshifts, from Keren's Readme about Abell 2744 lens model
    infile = datadir + "tab3.txt"
    df = pandas.read_csv(infile, delim_whitespace=True, comment='#')
    return(df)

def open_RCS0327_knotfamilies(datadir) :
    # Read Keren's list of knot families for RCS0327
    infile = datadir + 'arcs1.dat'
    names = ('ID', 'RA', 'DEC', 'a', 'b', 'c', 'zz', 'd') 
    df = pandas.read_csv(infile, delim_whitespace=True, names=names, comment='#')
    return(df)




#### Actually run stuff
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)  # Keren's favorite cosmology.
compute_lensing_test()

# Configure for A2744
datadir = '/Users/jrrigby1/SCIENCE/Lensed-LBGs/How_lensing_works/FF_A2744_Sharon_v4cor/' 
deflection_map = { 'x' : "hlsp_frontier_model_abell2744_sharon_v4cor_x-arcsec-deflect.fits", 'y' : "hlsp_frontier_model_abell2744_sharon_v4cor_y-arcsec-deflect.fits"}    
z_lens = 0.3080  # Abell 2744
arcs_df = open_A2744_lensfamily_table(datadir)
# According to Keren's FF readme, for A2744 the deflection maps are scaled to d_LS/d_S = 1.
A2744_df = jrr.lens.compute_sourceplane_positions(arcs_df, 'A2744', z_lens, datadir, deflection_map, cosmo, deflection_norm=1.0, debug=False)

datadir='/Users/jrrigby1/SCIENCE/Lensed-LBGs/RCS0327/Lens_Model/'
deflection_map= { 'x' : 'RCS0327_dplx_1.701.fits', 'y' : 'RCS0327_dply_1.701.fits'}
z_lens = 0.564 # Wuyts et al. 2010
arcs_df = open_RCS0327_knotfamilies(datadir)
# For RCS0327, the deflection maps are scaled for the giant arc(z_arc = 1.70), per KS email 01/08/2021
deflection_norm = jrr.lens.scale_lensing_deflection(cosmo, z_lens, 1.701)
RCS0327_df = jrr.lens.compute_sourceplane_positions(arcs_df, 'RCS0327', z_lens, datadir, deflection_map, cosmo, deflection_norm=deflection_norm, debug=False)

# Pseudo code for work to be done:
# - DONE    Convert from RA, DEC to pixel coordinates in the fits images.
# - DONE    Grab the values from the x, y deflection maps. 
# - DONE    Compute d_LS/d_S for the redshifts in question, and scale the deflections.
# - DONE   Go from image plane coords and deflection to source plane coords. 
# - DONE   Test that this works on all the multiply imaged galaxies in A2744.
#   DONE   Debugged, Keren found the problem, X was swapped (RA). It was the format of Keren's deflection maps.
# NEXT  Get Keren's lens models for the intervening absorbers
# NEXT  Compute source-plane coordinates for each of multiple sight-lines, at the redshift of each intervening absorber
