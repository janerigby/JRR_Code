''' Compute the relative velocity offsets between the nebular, photospheric, and ISM redshifts
measured for the MagE sample.  jrigby 4/2016'''

import jrr
mage_mode = "reduction"
mage_mode = "released"
(ss) =  jrr.Mage_getlist(mage_mode)

print "#Columns:  1) spectrum label   2) Nebular redshift   3) uncertainty  "
print "#4) velocity offset between measured redshifts for ISM and nebular lines."
print "# Negative sign means blueshift of ISM.  5) Uncertainty"
for ii in range(0, len(ss)) :                  #nfnu_stack[ii] will be ii spectrum
    label     =  ss['short_label'][ii]
    filename  =  ss['filename'][ii]
    z_neb   =   ss['z_neb'][ii]
    z_stars =   ss['z_stars'][ii]
    z_ISM =     ss['z_ISM'][ii]
    sig_neb   =  ss['sig_neb'][ii]
    sig_st    =  ss['sig_st'][ii]
    sig_ISM =    ss['sig_ISM'][ii]

#    if(ss['fl_st'][ii] or ss['fl_neb'][ii]) :
#        (vel_starsneb, dvel_starsneb) = (-99, -99)
#    else :
#        (vel_starsneb, dvel_starsneb) = jrr.velocity_offset(z_stars, sig_st, z_neb, sig_neb)

    if(ss['fl_neb'][ii] or ss['fl_ISM'][ii]) :
        (vel_nebISM, dvel_nebISM) = (-99, -99)
    else :
        (vel_nebISM, dvel_nebISM) = jrr.velocity_offset(z_neb,  sig_neb, z_ISM, sig_ISM)

    print label, z_neb, sig_neb, 
    print("%.2f   %.2f  " % (vel_nebISM, dvel_nebISM))
print "#DONE"

# Note: Now that Mage_getlist uses Pandas, not Astropy.tables, there are far more
# elegant ways to do this math.  However, it worked, so let's move on.
