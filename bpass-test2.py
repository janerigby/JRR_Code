import pandas
import jrr.bpass
import matplotlib.pyplot as plt


# I wrote some simple routines (jrr.bpass) to read in bpass2.1 spectra.  Testing them here.

######    
print "Loading a BPASS2 file, all ages (slow)"
filename = 'BPASSv2.1_imf135_100/spectra.z020.dat.gz'
df_bpass = jrr.bpass.load_spectra(filename)

# Find just one age of that model
age_to_find = 8E6

# Show that the age setup is working
(closest_age) = jrr.bpass.find_closest_age(age_to_find)
print "Closest age was", closest_age, ", so the age column is", closest_age['colname'].values[0]

print "Load a BPASS spectrum for one age (faster than loading all the ages)"
df_bpass2 = jrr.bpass.load_1spectrum(filename, age_to_find)
print df_bpass2.head()
ax = df_bpass2.plot(x='wave', y='flam')
plt.xlim(100,5000)
plt.show()
