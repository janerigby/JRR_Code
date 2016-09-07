''' Same as multipanel-spectrum-ids.py, but plot multiple spectra.  Using for sub-stacks (highZ, loZ),
young, middle-aged, old.
jrigby, Sept 2016. '''

import jrr
import matplotlib.pyplot as plt
from numpy import ones_like, zeros_like, median
from astropy.io import ascii
from commands import getoutput
import sys
from re import split, sub, search
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages

mage_mode = "reduction"  
#mage_mode = "released"

plotsize = (11,16)
Npanels = 24  # N panels total
Npages   = 4   # divided into this many pages of plots
lincolor = 'k'
errcolor = '0.55'
contcolor = '0.75'

stacks = ['magestack_bystars_highZ_spectrum.txt', 'magestack_bystars_lowZ_spectrum.txt']
sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
the_dfs = [sp1, sp2]
the_zzs = [0.0, 0.0]  # Stacked spectra are already in rest frame.
the_pdf =  "multipanel_stacked_byZ.pdf"
(spec_path, line_path) = getpath(mage_mode)
(LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.
plot_multipanel_spectrum(the_dfs, thezzs, LL, Npages, Npanels, the_pdf)  # Move this down below once the function is written

def plot_multipanel_spectrum(the_dfs, the_zzs, LL, Npages, Npanels, outfile) :
    ''' Plot an echelle(ette) spectrum with several pages and many panels for page.  Same as multipanel-spectrum-ids.py,
    but modular.  Inputs:
    the_dfs:  an array of dataframes containing spectra, to plot. Usually this is just one dataframe, one spectrum,
              but want to allow overplotting of multiple spectra.  First spectrum will be used to set the wavelength range.
    the_zzs:  an array of redshifts
    LL:       a linelist to plot features.  At present, takes redshift from the_zzs[0]
    Npages:   pages of plots in the output PDF.  Will split evenly.
    Npanels:  Total number of panels to plot, across all pages
    outfile:  Name of output PDF file.
    '''
    pp = PdfPages(the_pdf)  # output
    for jj in range(0,Npages):
        print "Drawing page ", jj, " of ", Npages-1
        fig    = plt.figure(num=jj+1, figsize=plotsize)

        # Set up the wavelength ranges
        begin_wave = the_dfs[0].wave.iloc[0:1]
        end_wave   = the_dfs[0].wave.iloc[-2:-1]
        waves_per_panel = len(the_dfs[0].wave) / Npanels
        start = waves_per_panel * np.linspace(0, Npanels) + begin_wave
        end   = waves_per_panel * np.linspace(0, Npanels) + begin_wave
        print "DEBUGGING", first_wave, last_wave
        print "This should be an array of start wavelengths", start

        
        
        for kk in range(0, Npanels/Npages):
            format = (Npanels/Npages)*100 + 11  # makes label: 311 for 3 rows, 1 col, start at 1
            subit = host_subplot(format + kk)
            for ss, df in enumerate(the_dfs) :
                if ss == 0 :  # If the first df, set the plot ranges
                    start =  jj*(len(df.wave)/Npages) +  kk*(len(df.wave)/Npanels)
                    end   =  start + len(df.wave)/Npanels
                    top = median(df.fnu[start:end])*1.1 + jrr.util.mad(df.fnu[start:end])*5
                    plt.ylim(0, top)  # trying to fix weird autoscaling from bad pixels
            plt.plot(wave[start:end], fnu[start:end], lincolor, linestyle='steps')   # plot spectrum
            plt.plot(wave[start:end], sig[start:end], errcolor, linestyle='steps')   # plot 1 sigma uncert spectrum
            plt.autoscale(axis='x', tight=True)
            if not skipcont :
                plt.plot(wave[start:end], cont[start:end], contcolor, linestyle='steps', zorder=1, linewidth=2) # plot the continuum

            jrr.mage.plot_linelist(LL, zz)   # Plot the line IDs.

            upper = subit.twiny()  # make upper x-axis in rest wave (systemic, or of abs lines)
            upper.set_xlim(wave[start]/(1.0+zz), wave[end]/(1.0+zz))

            majorLocator   = MultipleLocator(50)  # put subticks on lower x-axis
            minorLocator   = MultipleLocator(10)
            subit.xaxis.set_major_locator(majorLocator)
            subit.xaxis.set_minor_locator(minorLocator)
            subit.xaxis.tick_bottom()  # don't let lower ticks be mirrored  on upper axis

            subit.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=False, Symmetric=False))
            if(kk == Npanels/Npages-1):
                subit.set_xlabel(ur"observed-frame vacuum wavelength (\u00c5)")
                plt.ylabel('fnu') # fnu in cgs units: erg/s/cm^2/Hz
            if(kk == 0):
                upper.set_xlabel(ur"rest-frame vacuum wavelength (\u00c5)")
                plt.annotate("f_nu * "+str(factor) + " erg/s/cm^2/Hz",  xy=(0.05,0.055), color="black", xycoords="figure fraction")
                plt.annotate("Line color coding:", xy=(0.05,0.043), color="black", xycoords="figure fraction")
                plt.annotate("Photosphere", xy=(0.05,0.03), color="blue", xycoords="figure fraction")            
                plt.annotate("Emission",    xy=(0.17,0.03), color="red", xycoords="figure fraction")
                plt.annotate("ISM",         xy=(0.25,0.03), color="green", xycoords="figure fraction")
                plt.annotate("Wind",        xy=(0.32,0.03), color="cyan", xycoords="figure fraction")
                plt.annotate("Fine structure", xy=(0.4,0.03), color="purple", xycoords="figure fraction")
                plt.annotate("Intervening", xy=(0.55,0.03), color="orange", xycoords="figure fraction")
        plt.suptitle(rootname + " z=" +  str(zz))  # global title
        pp.savefig()    
        fig.canvas.draw()

    pp.close()
    print "   Generated PDF:  ", the_pdf
    return(0)
