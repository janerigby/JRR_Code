''' This is a check of interfaces. JC's S99 fits use spectra whose continuum has been divided
by a constant (the median in a set interval), converted to flambda, converted to rest-frame. 
JRR tends to plot spectra in observed units, and in fnu.  Need to be able to convert back 
and forth.  Therefore, wrote jrr.mage.redo_open_S99_spectrum(), which opens the spectrum,
figures out what normalization was applied, and writes a pandas data frame that has different
columns, but otherwise can be treated like a megasaura spectrum.  
jrigby, 9Dec2016
'''

import jrr
import  matplotlib.pyplot  as plt
from matplotlib.backends.backend_pdf import PdfPages
xlo=1000.
xhi=2000.

def check_reading_in_S99_rest(rootname) :
    (sp, LL) = jrr.mage.redo_open_S99_spectrum(rootname)
    (orig_sp, orig_resoln, orig_dresoln, orig_LL, orig_zz_syst)  =  jrr.mage.wrap_open_spectrum(rootname, "released")
    sp.rest_fnu_data.median()
    orig_sp.rest_fnu.median()
    plt.plot(orig_sp.rest_wave, orig_sp.rest_fnu, linewidth=2, color='black')
    plt.plot(sp.rest_wave, sp.rest_fnu_data, linewidth=0.5, color="red")
    plt.plot(sp.rest_wave, sp.rest_fnu_s99, linewidth=1, color="green")
    #plt.xlim(1400,1460)
    plt.xlim(xlo, xhi)
    plt.ylim(0.0, 1E-28)
    plt.title(rootname)
    #plt.show()
    return(0)

def check_reading_in_S99_obs(rootname) :
    (sp, LL) = jrr.mage.redo_open_S99_spectrum(rootname)
    (orig_sp, orig_resoln, orig_dresoln, orig_LL, orig_zz_syst)  =  jrr.mage.wrap_open_spectrum(rootname, "released")
    plt.plot(orig_sp.wave, orig_sp.rest_fnu, linewidth=2, color='black')
    plt.plot(sp.wave, sp.rest_fnu_data, linewidth=0.5, color="red")
    plt.plot(sp.wave, sp.rest_fnu_s99, linewidth=1, color="green")
    plt.xlim(xlo*(1+orig_zz_syst), xhi*(1+orig_zz_syst))
    plt.ylim(0.0, 1E-28)
    plt.title(rootname)
    #plt.show()
    return(0)

the_pdf = "test_readinS99_rest.pdf"
pp = PdfPages(the_pdf)  # output   

rootname = "chuck"       #doublecheck that John's S99 spectrum matches Chuck's 
chuck =  jrr.mage.read_chuck_UVspec()
(sp, LL) = jrr.mage.redo_open_S99_spectrum(rootname, denorm=True)
sp.rest_fnu_data.median()
chuck.rest_fnu.median()
fig = plt.figure(0, figsize=(8,5))
plt.plot(chuck.rest_wave, chuck.rest_fnu, linewidth=2, color='black')
plt.plot(sp.rest_wave, sp.rest_fnu_data, linewidth=0.5, color="red")
plt.plot(sp.rest_wave, sp.rest_fnu_s99, linewidth=1, color="green")
plt.xlim(xlo, xhi)
plt.ylim(0.0, 2)
plt.title(rootname)
pp.savefig()  #yes, it does match
plt.clf()

rootname = "ChisholmstackA"   # doublecheck that John's S99 spectrum matches my stack-A
alt_file = "magestack_byneb_ChisholmstackA_spectrum.txt"  # this is the spectrum that JChisholm fit
stack = jrr.mage.open_stacked_spectrum("reduction", alt_file)
(sp, LL) = jrr.mage.redo_open_S99_spectrum("Stack-A", denorm=True)
sp.rest_fnu_data.median()
stack.rest_fnu.median()
fig = plt.figure(1, figsize=(8,5))
plt.plot(stack.rest_wave, stack.rest_fnu, linewidth=2, color='black')
plt.plot(sp.rest_wave, sp.rest_fnu_data, linewidth=0.5, color="red")
plt.plot(sp.rest_wave, sp.rest_fnu_s99, linewidth=1, color="green")
plt.xlim(xlo, xhi)
plt.ylim(0.0, 2)
plt.title(rootname)
pp.savefig()  


S99files = ("S0004-0103", "S0033+0242", "S0108+0624", "S0900+2234", "S0957+0509", "S1226+2152", "S1429+1202", "S1458-0023", "S1527+0652", "S2111-0114", "Cosmic~Eye", "Horseshoe", "rcs0327-E", "rcs0327-G", "rcs0327-U") 
print "My spectrum should be in black, thinner overlapping red spectrum to verify normalization, and a green S99 fit"
for ii, rootname in enumerate(S99files) :
    print "Testing", rootname
    plt.clf()
    fig = plt.figure(ii+2, figsize=(8,5))
    check_reading_in_S99_rest(rootname)
    pp.savefig()
pp.close()

the_pdf = "test_readinS99_obs.pdf"
pp = PdfPages(the_pdf)  # output   
for ii, rootname in enumerate(S99files) :
    print "Testing", rootname
    plt.clf()
    fig = plt.figure(ii, figsize=(8,5))
    check_reading_in_S99_obs(rootname)
    pp.savefig()
pp.close()

