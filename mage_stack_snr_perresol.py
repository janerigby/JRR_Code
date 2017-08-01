import jrr
import numpy as np
import  matplotlib.pyplot as plt
(sp, LL) = jrr.mage.open_stacked_spectrum('released')
plt.ion()
plt.plot(sp['rest_wave'], sp['X_avg']/sp['X_jack_std'], color='k', lw=0.5)
plt.ylim(0,100)
RR = 3300
sp['pixperfwhm'] = sp['rest_wave'] / RR / 0.1
plt.grid()
sp['SNRsmooth'] = (sp['X_avg']/sp['X_jack_std']).rolling(window=11, center=True).median()
plt.plot(sp['rest_wave'], sp['SNRsmooth'], color='b')
plt.xlim(1400,1500)

lo = np.array((1440., 1460., 1680., 1760.))
hi = np.array((1450., 1470., 1700., 1800.))
print "wave_lo wave_hi   SNR_perpix_median  SNR_perresolnelement_median  over_resolved"
for ii, low in enumerate(lo) :
    subset =  sp.loc[sp['rest_wave'].between(lo[ii],hi[ii])]
    over_resolved = subset['pixperfwhm'].median()
    SNR_med_pp = (subset['X_avg']/subset['X_jack_std']).median()
    #print SNR_med_pp, SNR_med_pp * np.sqrt( over_resolved ), over_resolved
    sp['temp'] = (sp['X_avg']/sp['X_jack_std']).rolling(window=int(np.floor(over_resolved)), center=True).apply(jrr.util.add_in_quad)
    SNR_method1 = np.round(SNR_med_pp * np.sqrt( over_resolved ), 0)
    SNR_method2 = np.round((sp['temp'].loc[sp['rest_wave'].between(lo[ii],hi[ii])]).median(), 0)
    print "", SNR_method1, SNR_method2, lo[ii], hi[ii]

