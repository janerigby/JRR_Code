import numpy as np

names = ("SGAS1723", "RCS0327", "PSZ1G311.65-18.48", "SGAS1110+6459", "SGAS1226p2152", "SGAS1050p0017")
zz = np.array((1.329279, 1.7034, 2.369,                    2.481, 2.925257, 3.6253))

waves = np.array((0.4861, 0.6563))

#NIRSPEC GRATING & FILTER CHOICES
G140_F100 = np.array((0.97, 1.89))
G235_F170 = np.array((1.66, 3.17))
G396_F290 = np.array((2.87, 5.27))
the_gratings  = (G140_F100, G235_F170, G396_F290)
grating_names = ("G140_F100", "G235_F170", "G396_F290")

for ii, name in enumerate(names) :
    print name, "at z =",zz[ii]
    print "HB, Ha:",  (1.0 +zz[ii]) * waves
    for jj, grating in enumerate(the_gratings) :
        print "   ", grating / (1.0+zz[ii]) , grating_names[jj]

print "Ideal redshifts for each grating:"
for jj, grating in enumerate(the_gratings) :
    print grating_names[jj], np.average(grating) / np.average(waves) - 1.0
    
