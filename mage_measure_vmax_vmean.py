import operator
import jrr
import numpy as np
import pandas
import re
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

# Originally I measured vmax manually in interactive plots, using hack_to_get_vmax.py.  It was clunky.
# Now, actually measuring it the right way, using v_max and absorption-weighted mean velocity v_mean

#pandas.set_option('display.float_format', lambda x: '%.1f' % x)  # suppress scientific notation in output

def mage_mark_blends(df) :  # Manually enter lower limits for blends
    df.loc['O I 1302', 'vlowlim']   = -510.
    df.loc['O I 1302', 'comment']   =  "blend w photospheric abs."
    df.loc['Fe II 2383', 'vlowlim'] = -840.
    df.loc['Fe II 2383', 'comment'] = "blend with Fe II 2374 abs"
    df.loc['N V 1238', 'vlowlim']   = -10.
    df.loc['N V 1238', 'comment']   = "Not clearly detected"
    return(0)

def cos_mark_blends(df) :
    df.loc['O I 1302', 'vlowlim'] = -500.
    df.loc['O I 1302', 'comment']  =  "blend w photospheric abs"
    return(0)

# Good, I am now getting the uncertainty in vmean, vmax by varying the continuum via scalecont.  Hokey but realistic.
# Use COS lowres rather than COS R1.4E4.  It has better noise properties, fairer comparison to MagE.  
def wrap_measure_vmaxvmean(sp, colwave, colf, colcont, Nover_bluemax, Nover_red, thecenters, thelabels, IParray, pdfout) :
    scalecont = np.array((1.0, 0.98, 0.985, 0.99, 0.995, 1.01, 1.01, 1.015, 1.02)) 
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
                # Now that I'm using the R3300 COS spectrum, don't need Norder_blue>1.  Should vary continuum by N% and take avg.
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

def plot_vs_IP(df, color1, color2, ax, label, marker='o', s=60) :  # Plotting procedure gets run on several data frames, so make a function
    df_notlim = df[df['vlowlim'].isnull()]
    df_notlim.plot(x='IP', y='vmax',  kind='scatter', label=label[0], color=color2, s=s, ax=ax, marker=marker)
    df_notlim.plot(x='IP', y='vmean', kind='scatter', label=label[1], color=color1, s=s, ax=ax, marker=marker)
    ax.errorbar(df_notlim['IP'], df_notlim['vmean'], yerr=df_notlim['vmean_std'], ls='none', color='k', lw=1.5, label=None)
    ax.errorbar(df_notlim['IP'], df_notlim['vmax'],  yerr=df_notlim['vmax_std'],  ls='none', color='k', lw=1.5, label=None)
    plt.quiver(df['IP'], df['vlowlim'], np.zeros(shape=df.shape[0]), np.ones(shape=df.shape[0])*100, color=color1)
    return(0)

def finish_IPplot(ax, pdfout) :
    plt.xlabel("IP (eV)") ; plt.ylabel(r'v ($km$ $s^{-1}$)') 
    plt.xlim(10,70) ; plt.ylim(50,-2900)
    ax.legend(loc='upper left', labelspacing=0.2, borderpad=0.1)
    plt.tight_layout()
    plt.savefig(pdfout)
    return(0)



#############################################
mage_mode = 'reduction'
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
Nover_blue = 1 # Nth pixel over the continuum counts as the edge of the line for calculation of vmax, vmean
Nover_red = 1
plt.ioff()

# UGH what follows below is ugly, lots of copy and paste rather than loops.  Sorry, future self.

#New names:
#vmage_whtavg  # use the weighted avg stack of MagE
#vmage_median  # Use the median stack of MagE
#vcos_loR_whtavg
#vcos_hiR_whtavg
#vcos_loR_whtavg
#vcos_hiR_median

print "COMPUTING vmax vmean for the MagE stacked spectrum"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(sp, dumLL) = jrr.mage.open_stacked_spectrum(mage_mode, which_stack='standard', addS99=True)
sp['unity' ] = 1.0  # continuum
print "#  (velocities are in systemic rest frame, in km/s)"
#vmage_df is the shape-normalized MagE stack, weighted avg
vmage_whtavg = wrap_measure_vmaxvmean(sp, 'wave', 'fweightavg',  'unity', Nover_blue, Nover_red, thecenters, thelabels, IP, 'mage_vmaxvmean_whtavg.pdf')
vmage_median = wrap_measure_vmaxvmean(sp, 'wave', 'fmedian',     'unity', Nover_blue, Nover_red, thecenters, thelabels, IP, 'mage_vmaxvmean_median.pdf')
mage_mark_blends(vmage_whtavg)  ;   mage_mark_blends(vmage_median)   # manually enter lower limits for blends
vmage_whtavg.drop('linecen', axis=1).to_latex('mage_vmaxvmean_whtavg.tex')
vmage_median.drop('linecen', axis=1).to_latex('mage_vmaxvmean_median.tex')
plt.close("all")

print "COMPUTING vmax vmean for the R=1.4E4 COS stacked spectrum"
cos_hiR = jrr.mage.read_our_COS_stack(resoln="full")
cos_loR = jrr.mage.read_our_COS_stack(resoln="matched_mage")
vcos_hiR_whtavg = wrap_measure_vmaxvmean(cos_hiR, 'rest_wave', 'fweightavg', 'unity', Nover_blue, Nover_red, thecenters, thelabels, IP, 'cos_vmaxvmean_R1.4E4_whtavg.pdf')
vcos_hiR_median = wrap_measure_vmaxvmean(cos_hiR, 'rest_wave', 'fmedian',    'unity', Nover_blue, Nover_red, thecenters, thelabels, IP, 'cos_vmaxvmean_R1.4E4_median.pdf')
cos_mark_blends(vcos_hiR_whtavg)  ;  cos_mark_blends(vcos_hiR_median)  # Manually enter lower limits for blends
vcos_hiR_whtavg.drop('linecen', axis=1).to_latex('cos_vmaxvmean_R1.4E4_whtavg.tex')
vcos_hiR_median.drop('linecen', axis=1).to_latex('cos_vmaxvmean_R1.4E4_median.tex')
plt.close("all")

print "COMPUTING vmax vmean for the R=3300 COS stacked spectrum"
vcos_loR_whtavg = wrap_measure_vmaxvmean(cos_loR, 'rest_wave', 'fweightavg', 'unity', Nover_blue, Nover_red, thecenters, thelabels, IP, 'cos_vmaxvmean_R3300_whtavg.pdf')
vcos_loR_median = wrap_measure_vmaxvmean(cos_loR, 'rest_wave', 'fmedian',    'unity', Nover_blue, Nover_red, thecenters, thelabels, IP, 'cos_vmaxvmean_R3300_median.pdf')
cos_mark_blends(vcos_loR_whtavg)  ;  cos_mark_blends(vcos_loR_median) # Manually enter lower limits for blends
vcos_loR_whtavg.drop('linecen', axis=1).to_latex('cos_vmaxvmean_R3300_whtavg.tex')
vcos_loR_median.drop('linecen', axis=1).to_latex('cos_vmaxvmean_R3300_median.tex')
plt.close("all")


#########   Make plots of vmean, vmax vs IP    ###########################
label1 = (r'MagE $v_{max}$', 'MagE $v_{mean}$')
label2 = (r'COS $v_{max}$',  'COS $v_{mean}$')
label3 = ('_nolegend_',      '_nolegend_')

matplotlib.rcParams.update({'font.size': 16})
#plt.ion()   #  Plot velocities versus IP for the MagE stack
fig = plt.figure(1)
ax1 = fig.add_subplot(111)
plot_vs_IP(vmage_whtavg, 'red', 'blue', ax1, label1, marker='o')
plot_vs_IP(vmage_median,  'red', 'blue', ax1, label3, marker='s')  # temp commented out, debuggung
finish_IPplot(ax1, 'mage_vmaxvmean_vsIP.pdf')

plt.close("all")
fig = plt.figure(2)
ax2 = fig.add_subplot(111)  # Plot velocities for the COS stack
plot_vs_IP(vcos_hiR_whtavg, 'red', 'blue', ax2, label2, marker='o')
plot_vs_IP(vcos_hiR_median, 'red', 'blue', ax2, label3, marker='s')
finish_IPplot(ax2, 'cos_R1.4E4_vmaxvmean_vsIP.pdf')

plt.close("all")
fig = plt.figure(3)
ax3 = fig.add_subplot(111)  # Plot velocities for the COS stack
plot_vs_IP(vcos_loR_whtavg, 'red', 'blue', ax3, label2, marker='o')
plot_vs_IP(vcos_loR_median, 'red', 'blue', ax3, label3, marker='s')
finish_IPplot(ax3, 'cos_R3300_vmaxvmean_vsIP.pdf')



################################################################
##### set up EXAMPLE data frame, w wavelength, flux, continuum.  Was useful for debugging
if False :  
    colrwave = 'rest_wave' ; colf = 'f' ; colcont = 'cont'
    linecen = 1260.4221
    wavelo = 1100. ; wavehi = 1400. ; dwave = 0.2
    df = pandas.DataFrame(data=np.arange(wavelo, wavehi, dwave), columns=(colrwave,))  # New data frame
    df[colcont] = 1.0  # the continuum
    df[colf] = jrr.spec.onegaus(df[colrwave], -0.7,  linecen, 4., 1.0) + jrr.spec.onegaus(df[colrwave], -0.3,  linecen - 15., 9., 0.0) + np.random.normal(0.0, 0.03, size=len(df[colrwave]))  # the flux
    # Measure vmean, vmax 
    (vmean, vmax_blue, vmax_red) = jrr.spec.calc_vmean_vmax(df, colrwave, colf, colcont, Nover_blue, Nover_red, linecen, isabs=True, plotit=True)
    print "Vmean, vmax_blue, vmax_red:", vmean, vmax_blue, vmax_red



