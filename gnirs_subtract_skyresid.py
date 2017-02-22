import jrr
import re
import glob
import numpy as np
import pandas
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
import datetime
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units
import pyds9

''' Steph LaMassa reduced and extracted the GNIRS spectra, using the IRAF pipeline.
Frankly, the results are pretty bad -- there are huge sky residual left over.  I had
her (boxcar) extract the sky residials in the same ways as the object, in a region
right between the A and B images.  This script applies those sky residuals.  This
script loads her extracted sky residual, checks that the wavelength arrays are the
same in both files, and subtracts the sky residuals from the extracted spectra.
It then saves the residual-sky-subtracted spectra to a big new dataframe newdf,
does some simple bad-pixel rejection, and writes the mean, median, and error
in the mean to a file.   jrigby, oct 2016'''

s1723_dir = "/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/SDSSJ1723+3411/Gemini_GNIRS_highres/" # satchmo
#s1723_dir = '/Users/jrrigby1/SCIENCE/Gemini_GNIRS_highres/'  #milk

# PART 1:  These are the 1D spectra I will actually use.
in_dir    = s1723_dir + "Small_ap/"
resid_dir = s1723_dir + "sky_residuals/"
out_dir   = s1723_dir + "Small_ap_minus_skyresiduals/"
Bs = ["310m309_B_small_ap.txt", "311m312_B_small_ap.txt", "314m313_B_small_ap.txt", "315m316_B_small_ap.txt", "318m317_B_small_ap.txt"]
# 309m310_A_small_ap.txt  weird, commenting out until Steph fixes it
As = ["312m311_A_small_ap.txt", "313m314_A_small_ap.txt", "316m315_A_small_ap.txt", "317m318_A_small_ap.txt", "320m319_A_small_ap.txt"]

extracted_spectra_files  = Bs + As
names = ('oldwave', 'cts')

# Open the first spectrum, and set up a dataframe to put successive spectra
newdf =  pandas.read_table(in_dir + extracted_spectra_files[0], delim_whitespace=True, comment="#", skiprows=71, names=names)

for extract_file in extracted_spectra_files :
    df1 = pandas.read_table(in_dir + extract_file, delim_whitespace=True, comment="#", skiprows=71, names=names)
    skyresid_file = re.sub("small_ap.txt", "mid.txt", extract_file)
    df2 = pandas.read_table(resid_dir + skyresid_file, delim_whitespace=True, comment="#", skiprows=71, names=names)
    newdf[extract_file] = jrr.spec.rebin_spec_new(df1['oldwave'], df1.cts - df2.cts, newdf['oldwave'], fill=np.nan)
    
    sum_wave = (df1.oldwave - df2.oldwave).sum()  # this should be zero, if wavelength arrays are identical
    if sum_wave > 0 : print "WARNING, nonzero sum_wave", sum_wave
#    plt.step(df1.wave, df1.cts, color="blue")
#    plt.step(df2.wave, df2.cts, color="red")
    plt.step(df1.oldwave, df1.cts - df2.cts, linewidth=0.3)
    plt.ylim(-100,200)
    plt.xlim(1.5E4,1.6E4)
    # For B nods, subtraction of the two dataframes looks right.
#    plt.show()
Nspectra = len(extracted_spectra_files)*1.0
newdf.drop('cts', axis=1, inplace=True)  # remove cts, it's stale

# Applying barycentric correction to wavelengths.  All observations at about the same time.
thistime =  Time('2016-09-07T06:09:00', format='isot', scale='utc')
thisradec = SkyCoord("17:23:37.23", "+34:11:59.07", unit=(units.hourangle, units.deg), frame='icrs')
keck = EarthLocation.of_site('keck') # Gemini_N not in catalog, so use Keck.  Needs internet connection 
barycor_vel = jrr.barycen.compute_barycentric_correction(thistime, thisradec, location=keck)
print "FYI, the barycentric correction factor for s1723 was", barycor_vel
jrr.barycen.apply_barycentric_correction(newdf, barycor_vel, colwav='oldwave', colwavnew='wave') 

newdf.set_index('wave', inplace=True)  # set wavelength as the axis
badval_hi = 150.
badval_lo = -30.
newdf[newdf.gt(badval_hi)] = np.nan  # remove bad pixels from consideration
newdf[newdf.lt(badval_lo)] = np.nan
newdf['mean'] = newdf.mean(axis=1)
newdf['errinmean']  = newdf.std(axis=1) / np.sqrt(Nspectra)
newdf['std']  = newdf.std(axis=1)
newdf['median']  = newdf.median(axis=1)
plt.step(newdf.index, newdf['mean'], color='black', lw=2)
plt.step(newdf.index, newdf['errinmean'], color='red', lw=1)
plt.step(newdf.index, newdf['median'], color='green', lw=2)
plt.show()
myheader = "# Combined high-resolution GNIRS spectra for S1723.  Reduced by S. LaMassa.\n# Cleaned up and combined by jrigby on " + str(datetime.datetime.today()) + "\n# Wave is vacuum barycentric-corrected wavelength in Angstroms.\n# Other columns are counts (given as median and mean, with error in the mean.)  Spectra have NOT been fluxed.\n"
newdf.to_pickle(s1723_dir + "s1723_gnirs_residual_subtracted_spectrum_JRR.p")
newdf.to_csv("/tmp/tmpspec", columns=('mean', 'median', 'errinmean'))
jrr.util.put_header_on_file("/tmp/tmpspec", myheader,  s1723_dir + "s1723_gnirs_residual_subtracted_spectrum_JRR.csv")

# Part 2: I have promised myself I won't publish 1D spectra without publishing the
# corresponding 2D spectra as well.  Will have to hand-do this, and it's for display
# only (in other words, I'm not separately extracting it), but let's take a few min
# and generate a decent combined 2D spectrum with more-or-less the same residual sky
# subtraction.
indir_2D = s1723_dir + "2D_images/"
TwoD = ["309m310_Asci2.fits", "312m311_Asci2.fits", "313m314_Asci2.fits", "316m315_Asci2.fits", "317m318_Asci2.fits", "320m319_Asci2.fits"]

Axpos = np.array([24,   27.,  29.,  36.,  37., 37.])
Bxpos = np.array([47.,  49.,  53.,  56.,  60., 60.])
midpos = (Bxpos + Axpos)/2.0 # middle
boxcar_width =  7 # pixels
boxcar2 = 31
d = pyds9.DS9('foo1') 
for ii, image in enumerate(TwoD) :
    (data_in, header_in) = fits.getdata(indir_2D + image, header=True)
    lo = np.int(midpos[ii] - np.floor(boxcar_width/2.0))
    hi = np.int(midpos[ii] +  np.ceil(boxcar_width/2.0))
    midcut = data_in[:, lo:hi]
    thing_to_subtract = np.reshape(np.repeat(np.mean(midcut, axis=1), data_in.shape[1]), data_in.shape)
    subtr_im = data_in - thing_to_subtract
    #d.set("frame " + str(ii))
    #d.set_np2arr(subtr_im)
    if ii == 0 :
        output_sum = np.zeros_like(subtr_im)
        sumA    = np.zeros(shape=(subtr_im.shape[0], boxcar2))
        sumB    = np.zeros(shape=(subtr_im.shape[0], boxcar2))

    output_sum += subtr_im

    # extract just A nods, and sum them
    lo = np.int(Axpos[ii] - np.floor(boxcar2/2.0))
    hi = np.int(Axpos[ii] +  np.ceil(boxcar2/2.0))
    thisA = subtr_im[:, lo:hi]
    sumA += thisA

    # extract just B nods, and sum them
    lo = np.int(Bxpos[ii] - np.floor(boxcar2/2.0))
    hi = np.int(Bxpos[ii] +  np.ceil(boxcar2/2.0))
    thisB = subtr_im[:, lo:hi]
    sumB += thisB
    
d.set("frame " + str(ii+1))
d.set_np2arr(output_sum)
d.set("frame " + str(ii+2))
d.set_np2arr(sumA)
d.set("frame " + str(ii+3))
d.set_np2arr(sumB)
fits.writeto(s1723_dir + 'Asum_JRR.fits', sumA, header=header_in, clobber=True)
fits.writeto(s1723_dir + 'Bsum_JRR.fits', sumB, header=header_in, clobber=True)
fits.writeto(s1723_dir + 'ugly_all_sum_JRR.fits', output_sum, header=header_in, clobber=True)


# Part 3:  Read in the individual extractions for S2340, to get average and error in mean
# (Steph already calculated average, but from A and B, so no useful uncertainty.
# Doing it this way, uncertainty falls out.
s2340_dir = "/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/S2340/Gemini_GNIRS_highres/"  # satchmo

in_dir = s2340_dir + "2340+2947_indiv_spec/"
# no resid_dir, no need to remove residual skylines
indy_files = glob.glob(in_dir + "*_small_ap.txt")

newdf =  pandas.read_table(indy_files[0], delim_whitespace=True, comment="#", skiprows=71, names=names)
plt.clf()
for infile in indy_files :
    df1 =  pandas.read_table(infile, delim_whitespace=True, comment="#", skiprows=71, names=names)
    newdf[infile] = jrr.spec.rebin_spec_new(df1['oldwave'], df1.cts, newdf['oldwave'], fill=np.nan)
    plt.step(df1['oldwave'], df1['cts'], linewidth=0.3)

# Apply barycentric correction to wavelengths
thistime =  Time('2016-09-07T07:45:00', format='isot', scale='utc')
thisradec = SkyCoord("23:40:29.27", "+29:48:01.16", unit=(units.hourangle, units.deg), frame='icrs')
keck = EarthLocation.of_site('keck') # Gemini_N not in catalog, so use Keck
barycor_vel = jrr.barycen.compute_barycentric_correction(thistime, thisradec, location=keck)
print "FYI, the barycentric correction factor for s2340 was", barycor_vel
jrr.barycen.apply_barycentric_correction(newdf, barycor_vel, colwav='oldwave', colwavnew='wave') # testing       

newdf.drop('cts', axis=1, inplace=True)  # remove cts, it's stale
newdf.set_index('wave', inplace=True)  # set wavelength as the axis
newdf['mean'] = newdf.mean(axis=1)
newdf['errinmean']  = newdf.std(axis=1) / np.sqrt(len(indy_files))
newdf['median']  = newdf.median(axis=1)
plt.step(newdf.index, newdf['mean'], color='black', lw=2)
plt.ylim(-50,100)
plt.xlim(1.52E4,1.65E4)
plt.show()
myheader = "# Combined high-resolution GNIRS spectra for S2340.  Reduced by S. LaMassa.\n# Cleaned up and combined by jrigby on " + str(datetime.datetime.today()) + "\n# Wave is vacuum barycentric-corrected wavelength in Angstroms.\n# Other columns are counts (given as median and mean, with error in the mean.)  Spectra have NOT been fluxed.\n"
newdf.to_pickle(s2340_dir + "s2340_gnirs_spectrum_JRR.p")
newdf.to_csv("/tmp/tmpspec2", columns=('mean', 'median', 'errinmean'))
jrr.util.put_header_on_file("/tmp/tmpspec2", myheader,  s2340_dir + "s2340_gnirs_spectrum_JRR.csv")
