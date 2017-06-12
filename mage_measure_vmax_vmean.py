import operator
import jrr
import numpy as np
import pandas
import re
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

#pandas.set_option('display.float_format', lambda x: '%.1f' % x)  # suppress scientific notation in output

# FRIDAY TOMORROW, PICK UP HERE WHEN FRESH. ***
# Good, I am now getting the uncertainty in vmean, vmax by varying the continuum via scalecont.  Hokey but realistic.
# Measure Si 1260 by hand, continuum is screwed up for both Mage, COS lowres.
# Wait, why is Si 1260 screwed up?  Appears to be a continuum fitting problem...
# Use COS lowres rather than COS R2E4.  It has better noise properties, fairer comparison to MagE.  
def wrap_measure_vmaxvmean(sp, colwave, colf, colcont, Nover_bluemax, Nover_red, thecenters, thelabels, IParray, pdfout) :
    scalecont = np.array((1.0, 1.05, 0.995, 1.01, 0.99, 1.02, 0.98))
    # Measure the velocities and put them in a dataframe
    v_df = pandas.DataFrame(data=np.array(thelabels), columns=('linelab',))
    v_df['linecen'] = thecenters
    v_df.set_index('linelab', drop=True, inplace=True)
    vmean_ar = np.zeros_like(thecenters)    ;   vmean_std = np.zeros_like(thecenters)
    vmax_ar  = np.zeros_like(thecenters)    ;   vmax_std  = np.zeros_like(thecenters)
    pp = PdfPages(pdfout)
    for ii, linelab in enumerate(thelabels) :      # Measure vmean, vmax for the real stacked MagE spectrum
        if jrr.spec.test_wave_in_spectrum(sp, thecenters[ii], colwave) :
            plt.clf()
            vmean     = np.zeros_like(scalecont)
            vmax_blue = np.zeros_like(scalecont)
            vmax_red  = np.zeros_like(scalecont)
            for jj, sc in  enumerate(scalecont) :
                (vmean[jj], vmax_blue[jj], vmax_red[jj]) = jrr.spec.calc_vmean_vmax(sp, colwave, colf, colcont, 1, 1, thecenters[ii], scalecont=sc, isabs=True, plotit=True, label= linelab)
                # Now that I'm using the R3500 COS spectrum, don't need Norder_blue>1.  Should vary continuum by N% and take avg.
            vmean_ar[ii]  =  np.average(vmean)
            vmax_ar[ii]   =  np.average(vmax_blue)
            vmean_std[ii] =  np.std(vmean)
            vmax_std[ii]  =  np.std(vmax_blue)
            #print "DEBUG", linelab, vmean, vmean_ar[ii], vmean_std[ii], vmax_blue, vmax_ar[ii], vmax_std[ii] #, vmax_red
            plt.grid()
            pp.savefig()
    pp.close()
    v_df['IP'] = IParray
    v_df['vmean'] = np.round(vmean_ar)
    v_df['vmean_std'] = np.round(vmean_std)
    v_df['vmax']  = np.round(vmax_ar)
    v_df['vmax_std']  = np.round(vmax_std)
    v_df['vlowlim'] = np.nan
    v_df['comment'] = ""
    return(v_df)

# This would be nice, butoptional.  Deprioritizing.
def open_mage_jackknife() :
    jackfile = "/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Mage/Analysis/Stacked_spectra/magestack_byneb_standard_jackknife.txt"
    sp =  pandas.read_table(jackfile, delim_whitespace=True, comment="#", header=0)
    jrr.spec.calc_dispersion(sp, 'restwave', 'rest_disp')
    sp['badmask'] = False
    sp['linemask'] = False
    boxcar = jrr.spec.get_boxcar4autocont(sp)
    print "DEBUGGING, boxcar is ", boxcar
    # ** NEED AN LL
    jrr.spec.fit_autocont(sp, LL, zz=0.0, vmask=1000, boxcar=3001)
    # Want to step through each jackknife spectrum, fit fautocont, then loop through vmax,vmean procedure, get mean, std
    return(0)

mage_mode = 'reduction'

# Not to self.  Originally I measured vmax manually in interactive plots,
# using hack_to_get_vmax.py.  It was clunky.
# Now, actually measuring it the right way.

##### Define a bunch of lines. Copied from mage_winds.py
line_label_a         = ('O I 1302', 'Si II 1260', 'Si II 1526',  'Al II 1670', 'C II 1334')
line_center_a = np.array((1302.1685, 1260.4221,  1526.7066,  1670.7874,   1334.5323))      
line_label_b  =  ('Mg II 2796', 'Fe II 2344', 'Fe II 2383')
line_center_b = np.array((2796.352,  2344.214,  2382.765))
line_label_c     = ("Al III 1854", "Si IV 1393", "C IV 1548", "N V 1238")
line_center_c = np.array((1854.72, 1393.76, 1548.19,  1238.82))
line_label_all  = line_label_a  + line_label_b  + line_label_c# + line_label_extra
line_center_all = np.concatenate((line_center_a, line_center_b, line_center_c))#, line_center_extra))
# Add an array of ionization potentials, so I can plot velocities versus IP
IP = np.array((13.62, 16.35, 16.35, 18.83, 24.38,15.035, 16.1878, 16.1878, 28.44765, 45.14181, 64.4939, 97.8902))
##############################

# Pick which lines to measure
(thelabels, thecenters) = (line_label_all, line_center_all)

# Read spectra
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(sp, dumLL) = jrr.mage.open_stacked_spectrum(mage_mode, which_stack='standard', addS99=True)
sp['unity' ] = 1.0  # continuum
cos_df = jrr.mage.read_our_COS_stack(resoln="full")
cos_df2 = jrr.mage.read_our_COS_stack(resoln="matched_mage")

plt.ioff()
print "COMPUTING vmax vmean for the MagE stacked spectrum"
print "#  (velocities are in systemic rest frame, in km/s)"
Nover_blue = 1 # Nth pixel over the continuum counts as the edge of the line for calculation of vmax, vmean
Nover_red = 1
vmage_df = wrap_measure_vmaxvmean(sp, 'wave', 'X_avg', 'unity', Nover_blue, Nover_red, thecenters, thelabels, IP, 'mage_vmaxvmean.pdf')

# Manually enter lower limits for blends
vmage_df.loc['O I 1302', 'vlowlim']   = -510.
vmage_df.loc['O I 1302', 'comment']   =  "blend w photospheric abs"
vmage_df.loc['Al II 1670', 'vlowlim'] = -760.
vmage_df.loc['Al II 1670', 'comment'] = "blend with O III] emission"
vmage_df.loc['Fe II 2383', 'vlowlim'] = -840.
vmage_df.loc['Fe II 2383', 'comment'] = "blend with O III] emission"
vmage_df.loc['N V 1238', 'vlowlim']   = -10.
vmage_df.loc['N V 1238', 'comment']   = "Not clearly detected"

print vmage_df.head(100)
vmage_df.drop('linecen', axis=1).to_latex('mage_vmaxvmean.tex')               
plt.close("all")

print "COMPUTING vmax vmean for the COS stacked spectrum"
vcos_df = wrap_measure_vmaxvmean(cos_df, 'rest_wave', 'fweightavg', 'unity', Nover_blue, Nover_red, thecenters, thelabels, IP, 'cos_vmaxvmean_R2E4.pdf')
# Manually enter lower limits for blends
vcos_df.loc['O I 1302', 'vlowlim'] = -500.
vcos_df.loc['O I 1302', 'comment']  =  "blend w photospheric abs"
print vcos_df.head(100)
vcos_df.drop('linecen', axis=1).to_latex('cos_vmaxvmean_R2E4.tex')
plt.close("all")

print "COMPUTING vmax vmean for the COS stacked spectrum"
vcos_df2 = wrap_measure_vmaxvmean(cos_df2, 'rest_wave', 'fweightavg', 'unity', Nover_blue, Nover_red, thecenters, thelabels, IP, 'cos_vmaxvmean_R3500.pdf')
# Manually enter lower limits for blends
vcos_df2.loc['O I 1302', 'vlowlim'] = -500.
vcos_df2.loc['O I 1302', 'comment']  =  "blend w photospheric abs"
print vcos_df2.head(100)
vcos_df2.drop('linecen', axis=1).to_latex('cos_vmaxvmean_R3500.tex')
plt.close("all")

plt.ion()
## Plot the results versus ionization potential
matplotlib.rcParams.update({'font.size': 16})
s=60
vmage_notlim = vmage_df[vmage_df['vlowlim'].isnull()]
ax1 = vmage_notlim.plot(x='IP', y='vmax', kind='scatter', label=r'MagE $v_{max}$', color='red', s=s)
ax1.errorbar(vmage_notlim['IP'], vmage_notlim['vmax'], yerr=vmage_notlim['vmax_std'], ls='none', color='k', lw=1.5, label=None)
vmage_notlim.plot(x='IP', y='vmean', kind='scatter', label=r'MagE $v_{mean}$', color='blue', s=s, ax=ax1)
ax1.errorbar(vmage_notlim['IP'], vmage_notlim['vmean'], yerr=vmage_notlim['vmean_std'], ls='none', color='k', lw=1.5, label=None)
plt.quiver(vmage_df['IP'], vmage_df['vlowlim'], np.zeros(shape=vmage_df.shape[0]), np.ones(shape=vmage_df.shape[0])*100, color='r')
plt.xlabel("IP (eV)") ; plt.ylabel("v (km/s)")
plt.xlim(10,70) ; plt.ylim(0,-2900)
ax1.legend(loc='upper left', labelspacing=0.2, borderpad=0.1)
plt.tight_layout()
plt.savefig('mage_vmaxvmean_vsIP.pdf')

vcos_notlim = vcos_df[vcos_df['vlowlim'].isnull()]
ax2 = vcos_notlim.plot(x='IP', y='vmax', kind='scatter', label=r'COS $v_{max}$', color='red', s=s)
ax2.errorbar(vcos_notlim['IP'], vcos_notlim['vmax'], yerr=vcos_notlim['vmax_std'], ls='none', color='k', lw=1.5, label=None)
vcos_notlim.plot(x='IP', y='vmean', kind='scatter', label=r'COS $v_{mean}$', color='blue', s=s, ax=ax2)
ax2.errorbar(vcos_notlim['IP'], vcos_notlim['vmean'], yerr=vcos_notlim['vmean_std'], ls='none', color='k', lw=1.5, label=None)
plt.quiver(vcos_df['IP'], vcos_df['vlowlim'], np.zeros(shape=vcos_df.shape[0]), np.ones(shape=vcos_df.shape[0])*100, color='r')
plt.xlim(10,70) ; plt.ylim(0,-2900)
plt.xlabel("IP (eV)") ; plt.ylabel("v (km/s)")
ax2.legend(loc='upper left', labelspacing=0.2, borderpad=0.1)
plt.tight_layout()
plt.savefig('cos_vmaxvmean_vsIP.pdf')



##### set up EXAMPLE data frame, w wavelength, flux, continuum.  Was useful for debugging
if False :  
    colrwave = 'rest_wave' ; colf = 'f' ; colcont = 'cont'
    linecen = 1260.4221
    wavelo = 1100. ; wavehi = 1400. ; dwave = 0.2
    df = pandas.DataFrame(data=np.arange(wavelo, wavehi, dwave), columns=(colrwave,))  # New data frame
    df[colcont] = 1.0  # the continuum
    df[colf] = jrr.spec.onegaus(df[colrwave], -0.7,  linecen, 4., 1.0) + jrr.spec.onegaus(df[colrwave], -0.3,  linecen - 15., 9., 0.0) + np.random.normal(0.0, 0.03, size=len(df[colrwave]))  # the flux
    # Measure vmean, vmax 
    (vmean, vmax_blue, vmax_red) = jrr.spec.calc_vmean_vmax(df, colrwave, colf, colcont, Nover, linecen, isabs=True, plotit=True)
    print "Vmean, vmax_blue, vmax_red:", vmean, vmax_blue, vmax_red



