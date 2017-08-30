import jrr
import numpy as np
import  matplotlib.pyplot as plt

def calc_smoothSNR(sp, resoln, smooth_window, colwave='wave', colf='fnu', colfu='fnu_u', coldisp='disp') :
    # For a mage spectrum, calculate the SNR *per resoln element*, and smooth it.
    sp['SNRperres'] = sp[colf]/sp[colfu] * np.sqrt(sp[colwave] / resoln / sp[coldisp])
    sp['smoothSNRperres'] = sp['SNRperres'].rolling(window=smooth_window, center=True).median()
    return(0)

mage_mode = 'released'
smooth_window = 151

lab1 = [ 'rcs0327-E', 'rcs0327-U', 'rcs0327-B', 'rcs0327-G', 'rcs0327-counterarc',]
lab2 = [ 'S1527+0652', 'S1527+0652-fnt']
lab3 = [ 'S0004-0103', 'S0900+2234', 'S1226+2152']
lab4 = [ 'S0108+0624', 'Cosmic~Eye', 'S1429+1202']
lab5 = [ 'S0033+0242', 'S2243-0935', 'S2111-0114', 'S0957+0509']
lab6 = ['S1458-0023',  'Horseshoe',  'S1050+0017']

labs = (lab1, lab2, lab3, lab4, lab5, lab6)
#labs = (lab1,)

for labels in labs:
    (df, resoln, dresoln, LL, zz_sys, speclist) = jrr.mage.open_many_spectra(mage_mode, which_list="labels", labels=labels, verbose=True, zchoice='stars', addS99=False, MWdr=False)
    for label in labels:
        sp = df[label]
        calc_smoothSNR(sp, resoln[label], smooth_window)
        #plt.plot(sp['wave'], sp['fnu']/sp['fnu_u'], label=label+' unsmoothed')
        plt.plot(sp['wave'], sp['smoothSNRperres'], label=label)
    plt.legend()
    plt.ylim(0,30)
    plt.show()

    
'''
#Below works on the stack
(sp, LL) = jrr.mage.open_stacked_spectrum(mage_mode)
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
'''
