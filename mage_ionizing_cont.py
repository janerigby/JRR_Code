import jrr
import numpy as np
import scipy  
import scikits.bootstrap as bootstrap  
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from time import sleep
from matplotlib.backends.backend_pdf import PdfPages
import pandas
mage_mode = "released"

pandas.set_option('precision', 5)
np.set_printoptions(precision=5)

the_pdf = "Lyman_continuum.pdf"
pp = PdfPages(the_pdf) # output

Lyc =  (850., 910.)  # region for Lyman continuum (Angstroms)
red = (1481.,1520.)  # region for comparion (Angstroms)  Following Vasei+ 2016)

def avg_errorbar(f_this) :
#    return( np.mean(np.abs(f_this[1:3])))
    return( np.mean(np.abs(f_this)))

    
def measure_Lyman_continuum(df, colwave="rest_wave", colfnu="rest_fnu"):
    # Need to check whether the wavelength regimes are valid
    if df[colwave][0] > Lyc[0] :
        #print "Warning: Lyman continuum not covered for this object"
        packitup = (-99, -99, -99, -99, -99, -99, -99, -99)
    else :
        f_Lyc = jrr.util.bootstrap_val_confint( df[df[colwave].between(*Lyc)][colfnu], np.median, alpha=0.05) # Median flux in Ly cont
        f_red = jrr.util.bootstrap_val_confint( df[df[colwave].between(*red)][colfnu], np.median, alpha=0.05)  # and at 1500A
        ratio = f_Lyc[0] / f_red[0]
        ratio_u =  jrr.util.sigma_adivb(f_Lyc[0], avg_errorbar(f_Lyc[1:3]), f_red[0], avg_errorbar(f_red[1:3]))
        packitup = (f_Lyc[0], f_Lyc[1], f_Lyc[2], f_red[0], f_red[1], f_red[2], ratio, ratio_u)
    return(packitup)

def plot_the_measurement(df, result, label, zz=0, colwave="rest_wave", colfnu="rest_fnu") :
    plt.clf()
    fig = plt.figure(figsize=(14,4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.step(df[colwave], df[colfnu], color='black', lw=0.5)
    if zz >0 : # Show where extraction may be dodgy
        toobluewave = 3200.
        tooblue = df[df['wave'].lt(toobluewave)]
        ax1.step(tooblue[colwave], tooblue[colfnu], color='yellow', lw=0.5)
    ax1.plot( Lyc, np.ones_like(Lyc)*result[0], color='blue', lw=2)
    ax1.plot( red, np.ones_like(red)*result[3], color='red',  lw=2)         
    ax1.errorbar( np.mean(Lyc), result[0], xerr=None, yerr=avg_errorbar(result[1:3]), lw=2, color='blue', capthick=2)
    ax1.errorbar( np.mean(red), result[3], xerr=None, yerr=avg_errorbar(result[4:6]), lw=2, color='red',  capthick=2)
    x1=800 ; x2=red[1]+100.
    ax1.set_xlim(x1, x2)
    ax2.set_xlim(x1*(1.0+zz), x2*(1.0+zz))
    tick_spacing = 200.
    loc = ticker.MultipleLocator(base=tick_spacing) # this locator puts ticks at regular intervals
    ax2.xaxis.set_major_locator(loc)
    #ax1.xaxis.set_major_locator(loc)
#    ax1.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
#    ax2.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    plt.ylim((df[colfnu].median()*-0.2), jrr.util.robust_max(df[colfnu]))
    ax1.plot( (Lyc[0], red[1]), (0,0), color="green")
#    plt.title(label)
    plt.annotate(label, (0.5,0.9), xycoords="axes fraction", fontsize=14)
    pp.savefig()
    return(0)

# Measure this for the stack
print "#Currently printing 90% confidence intervals, bc util.bootstrap_val_confint(alpha=0.05)"
print "Object,  fnu_Lyc_median,  fLyc_conflo,  fLyc_confhi,  fnu_1500A,  f1500_conflo, f1500_confhi, fLycF1500rat, fLycF1500uncert"
chuck  = jrr.mage.read_chuck_UVspec(addS99=True, autofitcont=True)
(stack, LL) = jrr.mage.open_stacked_spectrum(mage_mode, which_stack="Stack-A", addS99=True)
wht_stack_result = measure_Lyman_continuum(stack, colfnu="fnu")
med_stack_result = measure_Lyman_continuum(stack, colfnu="fmedian")
print "StackA_wtdavg, ", str(wht_stack_result)[1:-1]
print "StackA_median, ", str(med_stack_result)[1:-1]
plot_the_measurement(stack, wht_stack_result, "Stack-A wtd avg", zz=0)
plot_the_measurement(stack, med_stack_result, "Stack-A median ", colfnu="fmedian", zz=0)

## Measure for individual MagE spectra, where Lycont is covered
## Comment out reloading spectra while debugging
debuglist = ("Cosmic~Eye", "rcs0327-E")
(sp, resoln, dresoln, LL, zz_sys, speclist) = jrr.mage.open_many_spectra(mage_mode, verbose=False, MWdr=False, silent=True)
result = {}
for label in speclist['short_label'] :
    z_syst = speclist.ix[label]['z_neb']
    (result[label]) = measure_Lyman_continuum(sp[label].interpolate())
    if result[label][5] != -99 :
        print label, ",", str(result[label])[1:-1]
        plot_the_measurement(sp[label], result[label], label, zz=z_syst)        

# COPY THE RESULTS INTO DATAFRAMES, TO PLOT AND ANALYZE
temp1 = pandas.DataFrame.from_dict(result, orient='index')
df_indy = temp1.loc[temp1[0] > -99]  # Drop those where Lycont not covered
temp2 = pandas.DataFrame.from_records([wht_stack_result,], index=('StackA_wtdavg',))
temp3 = pandas.DataFrame.from_records([med_stack_result,], index=('StackA_median',))
df_stack = pandas.concat([temp2, temp3])
## The dataframes I will want to plot are df_ind and df_stack.
df_indy.columns =     ['fnu_Lyc', 'fLyc_unclo', 'fLyc_unchi', 'fnu_1500', 'f1500_unclo', 'f1500_unchi', 'fLycF1500rat', 'fLycF1500uncert']
df_stack.columns  = ['fnu_Lyc', 'fLyc_unclo', 'fLyc_unchi', 'fnu_1500', 'f1500_unclo', 'f1500_unchi', 'fLycF1500rat', 'fLycF1500uncert']
df_indy.to_csv("Lycont_measurements_individual_mage_spectra.csv", na_rep='NaN')  # output the results
df_stack.to_csv("Lycont_measurements_stacked_mage_spectra.csv", na_rep='NaN')

#Plot fnu at Lyc versus fnu_1500A.
fig = plt.figure(figsize=(8,8))
plt.scatter( df_indy['fnu_1500'], df_indy['fnu_Lyc'], color='k')
for row in df_indy.itertuples():
    plt.annotate(row.Index, xy=(row.fnu_1500, row.fnu_Lyc), xycoords='data', xytext=(4,3), textcoords="offset points")
yerrors = [df_indy['fLyc_unclo'].values,  df_indy['fLyc_unchi'].values]
xerrors = [df_indy['f1500_unclo'].values, df_indy['f1500_unchi'].values]
plt.errorbar(df_indy['fnu_1500'], df_indy['fnu_Lyc'], yerr=yerrors, xerr=xerrors, fmt='none', ecolor='k', lw=2, elinewidth=2, capthick=2, capsize=3)  # need to solve plotting..
plt.xlabel ("median fnu at rest-frame 1500A")
plt.ylabel ("median fnu at Lyman continuum")
plt.xscale('log')
plt.yscale('log')
plt.xlim(1E-30, 1E-28)
plt.ylim(4E-32, 6E-28)
pp.savefig()

# Plot flux ratio (Lyc to 1500A vers 1500A)
fig = plt.figure(figsize=(8,8))
plt.clf()
df_indy['fLyc_unc'] = (df_indy['fLyc_unclo']   + df_indy['fLyc_unchi']) / 2.0
df_indy['f1500_unc'] = (df_indy['f1500_unclo'] + df_indy['f1500_unchi'])/ 2.0
yerrors = jrr.util.sigma_adivb_df(df_indy, 'fnu_Lyc', 'fLyc_unc', 'fnu_1500', 'f1500_unc')
plt.scatter(df_indy['fnu_1500'],   df_indy['fnu_Lyc']/df_indy['fnu_1500'], color='k')
for row in df_indy.itertuples():
    plt.annotate(row.Index, xy=(row.fnu_1500, row.fnu_Lyc/row.fnu_1500), xycoords='data', xytext=(4,3), textcoords="offset points")
plt.errorbar(df_indy['fnu_1500'],   df_indy['fnu_Lyc']/df_indy['fnu_1500'], xerr=xerrors,  yerr=yerrors, fmt='none', ecolor='k', lw=2, elinewidth=2, capthick=2, capsize=3)
plt.scatter( (0.9E-28,1E-28),   df_stack['fnu_Lyc']/df_stack['fnu_1500'], color='r', label='stacks')   # add the stacks

plt.xscale('log')
plt.xlim(1E-30, 1E-28)
plt.ylim(0.,0.5)
plt.xlabel("median fnu at rest-frame 1500A")
plt.ylabel("fnu(Lyman cont.) / fnu(1500A)")
pp.savefig()

pp.close()
plt.close("all")


#plt.show()
