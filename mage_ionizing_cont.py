import jrr
import numpy as np
import scipy  
import scikits.bootstrap as bootstrap  
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from time import sleep
from matplotlib.backends.backend_pdf import PdfPages
mage_mode = "released"

the_pdf = "Lyman_continuum.pdf"
pp = PdfPages(the_pdf) # output

Lyc =  (850., 910.)  # region for Lyman continuum (Angstroms)
red = (1480.,1520.)  # region for comparion (Angstroms)  Following Vasei+ 2016)

def avg_errorbar(f_this) :
#    return( np.mean(np.abs(f_this[1:3])))
    return( np.mean(np.abs(f_this)))

    
def measure_Lyman_continuum(df, wavcol="rest_wave", fnucol="rest_fnu"):
    # Need to check whether the wavelength regimes are valid
    if df[wavcol][0] > Lyc[0] :
        #print "Warning: Lyman continuum not covered for this object"
        packitup = (-99, -99, -99, -99, -99, -99)
    else :
        f_Lyc = jrr.util.bootstrap_val_confint( df[df[wavcol].between(*Lyc)][fnucol], np.mean) # Median flux in Ly cont
        f_red = jrr.util.bootstrap_val_confint( df[df[wavcol].between(*red)][fnucol], np.mean)  # and at 1500A
        # This should be median, not mean, but it is breaking.  Google when off a plane ******
        ratio = f_Lyc[0] / f_red[0]
        ratio_u =  jrr.util.sigma_adivb(f_Lyc[0], avg_errorbar(f_Lyc[1:3]), f_red[0], avg_errorbar(f_red[1:3]))
        #packitup = (f_Lyc, f_Lyc_u, f_red, f_red_u, ratio, ratio_u)
        packitup = (f_Lyc[0], f_Lyc[1], f_Lyc[2], f_red[0], f_red[1], f_red[2], ratio, ratio_u)
    return(packitup)

def plot_the_measurement(df, result, label, zz=0, wavcol="rest_wave", fnucol="rest_fnu") :
    plt.clf()
    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.step(df[wavcol], df[fnucol], color='black', lw=0.5)
    if zz >0 : # Show where extraction may be dodgy
        toobluewave = 3200.
        tooblue = df[df['wave'].lt(toobluewave)]
        ax1.step(tooblue[wavcol], tooblue[fnucol], color='yellow', lw=0.5)
    ax1.plot( Lyc, np.ones_like(Lyc)*result[0], color='blue', lw=2)
    ax1.plot( red, np.ones_like(red)*result[3], color='red',  lw=2)         
    ax1.errorbar( np.mean(Lyc), result[0], xerr=None, yerr=avg_errorbar(result[1:3]), lw=2, color='blue', capthick=2)
    ax1.errorbar( np.mean(red), result[3], xerr=None, yerr=avg_errorbar(result[4:6]), lw=2, color='red',  capthick=2)
    ax1.set_xlim(Lyc[0], red[1])
    ax2.set_xlim(Lyc[0]*(1.0+zz), red[1]*(1.0+zz))
    tick_spacing = 200.
    loc = ticker.MultipleLocator(base=tick_spacing) # this locator puts ticks at regular intervals
    ax2.xaxis.set_major_locator(loc)
    #ax1.xaxis.set_major_locator(loc)
#    ax1.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
#    ax2.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    plt.ylim((df[fnucol].median()*-0.1), jrr.util.robust_max(df[fnucol]))
    ax1.plot( (Lyc[0], red[1]), (0,0), color="green")
#    plt.title(label)
    plt.annotate(label, (0.5,0.9), xycoords="axes fraction", fontsize=14)
    pp.savefig()
    return(0)

# Measure this for the stack
chuck  = jrr.mage.read_chuck_UVspec(addS99=True, autofitcont=True)
(stack, LL) = jrr.mage.open_stacked_spectrum(mage_mode, which_stack="Stack-A", addS99=True)
wht_stack_result = measure_Lyman_continuum(stack, fnucol="X_avg")
med_stack_result = measure_Lyman_continuum(stack, fnucol="X_median")
print "Stack-A wtd_avg: ", wht_stack_result
print "Stack-A median:  ", med_stack_result
plot_the_measurement(stack, wht_stack_result, "Stack-A wtd avg", zz=0)
plot_the_measurement(stack, med_stack_result, "Stack-A median ", fnucol="X_median", zz=0)

## Measure for individual MagE spectra, where Lycont is covered
## Comment out reloading spectra while debugging
(sp, resoln, dresoln, LL, zz_sys, speclist) = jrr.mage.open_many_spectra(mage_mode, verbose=False)
result = {}
for label in speclist['short_label'] :
    z_syst = speclist.ix[label]['z_neb']
    (result[label]) = measure_Lyman_continuum(sp[label])
    if result[label][5] != -99 :
        print label, result[label][2], "pm", result[label][3]
        print "    ", label, result[label]
        plot_the_measurement(sp[label], result[label], label, zz=z_syst)        
pp.close()
plt.close("all")
# Next, use google, and figure out why np.median breaks the boostrap.ci confidence estimator.
#,Use median, and move on.
