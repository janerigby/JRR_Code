# from https://gist.github.com/StuartLittlefair/5aaf476c5d7b52d20aa9544cfaa936a1
# This uses an experimental correction that is not yet implemented in Astropy.
# Need to check against jskycalc.

from astropy.time import Time
from astropy.coordinates import SkyCoord, solar_system, EarthLocation, ICRS
from astropy import units
from jrr import barycen

my_target = SkyCoord("04:40:00.00", "-09:45:00", unit=(units.hourangle, units.deg), frame='icrs')
my_target = SkyCoord("18:00:00.00", "-09:00:00", unit=(units.hourangle, units.deg), frame='icrs')

# EarthLocation.get_site_names()  # Print names of all sites astropy knows about
keck = EarthLocation.of_site('keck')
lco = EarthLocation.of_site('Las Campanas Observatory')
mytime = Time('2017-02-25T08:00:00.00000', format='isot', scale='utc')
mytime = Time('2017-07-25T01:00:00.00000', format='isot', scale='utc')
result = barycen.velcorr(mytime, my_target, location=lco)
print result

# Conclusion:  This works great.  Verified for LCO, Keck, for a few coordinates and times.
