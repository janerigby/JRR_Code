import numpy as np
import stsynphot as stsyn
from astropy import constants as c

def freq(wave):   # c over lambda, that's what's nu.  wave in micron
    cc = c.c.to('micron/s')
    return(cc/wave)

def octave(w1, w2) : # input is wave1 wave2 in microns
    return(np.abs(np.log2( freq(w1) / freq(w2))))

print("Making raw comparison of number of megapixels")
Mp_nircam  = 2048**2 * 10
Mp_nirspec = 2048**2 * 2
Mp_niriss  = 2048**2 * 1.
Mp_miri    = 1024**2 * 3
results = np.array((Mp_nircam, Mp_nirspec, Mp_niriss, Mp_miri))
print("JWST megapixels:", results/1E6, np.sum(results)/1E6)

D_hst =  np.sqrt(stsyn.config.conf.area / np.pi )* 2. / 100.  # in m
D_jwst = np.sqrt(331830.72404           / np.pi )* 2. / 100.  # in m
Mp_wfc3_uvis = 2051 * 4096 * 2.
Mp_wfc3_ir   = 1024**2 * 1
results_hst = np.array((Mp_wfc3_uvis, Mp_wfc3_ir))
print("WFC3 megapixels:", results_hst/1E6, np.sum(results_hst/1E6))

# N_Megapixels comparison of JWST (57) to HST (18 on WFC3, haven't counted up ACS yet)
# is not that impressive a comparison.  So much easier to make large-format CCDs than
#HgTeCd devices.  It is true that HST has only 1 megapixel of SIMILAR detectors
#to JWST, but that's not a fair comparison, b/c IR is not HST's strength.
# Let's skip the comparison to # Mpc HST, and compare to an iphone instead.

print("\nComparing number of spatial resoln elements per field of view.")
FOV_UVIS = 162*162 #square arcsec
PSF_UVIS = 0.555E-6 / D_hst  * 206265   #in arcsec, at F555W
Apsf_UVIS = (PSF_UVIS/2.)**2 * np.pi

FOV_nircam = (2.2*60)**2 * 2.   # square arcsec, summing over both modules, only SW
PSF_nircam = 1.75E-6 / D_jwst * 206265  # in arcsec, at 2 micron
Apsf_nircam   = (PSF_nircam / 2.)**2 * np.pi
print("      fov(sqarcsec), PSF(arcsec),  Npsf/fov(in millions)")
print("UVIS:",   FOV_UVIS,  PSF_UVIS,    FOV_UVIS / Apsf_UVIS /1E6)
print("Nircam:", FOV_nircam, PSF_nircam, FOV_nircam / Apsf_nircam / 1E6)
print("Result: both UVIS and NIRCam have about 10M resolution elements per FOV.")
print("After noodling w diffraction limit, if JWST is diffraction limited at 1.75 micron,")
print("then it will equal the WFC3 figure of merit.  If 2um, a bit worse (11 vs 15).")

print("\nOctaves:")
  # Human hearing range is 10 octaves
#https://en.wikipedia.org/wiki/Piano_key_frequencies
middle_C = 261.6  #Hz
piano_names = (1, "middleC", 88) 
piano_freq  = np.array((27.5, 261.6, 4186.0))  # Hz
piano_octave = octave(piano_freq[-1],  piano_freq[0])


octave_jwst = octave(28.3, 0.7)
octave_hst = octave(2.5, 0.115) 
octave_hst_that_overlaps_jwst = octave(2.5, 0.7)
octave_eyeball = octave( 0.7200, 0.3800)

print("JWST  HST   HST/JWST_overlap   eyeball  piano")
print(np.round((octave_jwst, octave_hst, octave_hst_that_overlaps_jwst, octave_eyeball, piano_octave), 2))

