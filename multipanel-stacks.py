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
Nfigs   = 4   # divided into this many pages of plots
lincolor = 'k'
errcolor = '0.55'
contcolor = '0.75'

stacks = ['magestack_bystars_highZ_spectrum.txt', 'magestack_bystars_lowZ_spectrum.txt']
sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[0])
sp1 = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=stacks[1])
the_pdf =  "multipanel_stacked_byZ.pdf"
                               

def plot_one_multipanel(sp1, sp2) :
    pp = PdfPages(the_pdf)  # output
    (spec_path, line_path) = getpath(mage_mode)
    (LL, zz) = jrr.mage.get_linelist(line_path + "stacked.linelist")  #z_syst should be zero here.

for j in range(0,Nfigs):
    print "Drawing page ", j, " of ", Nfigs-1
    fig    = plt.figure(num=j+1, figsize=plotsize)  #from matplotlib
    for k in range(0, Npanels/Nfigs):
        format = (Npanels/Nfigs)*100 + 11  # makes label: 311 for 3 rows, 1 col, start at 1
        subit = host_subplot(format + k)
        start =  j*(len(wave)/Nfigs) +  k*(len(wave)/Npanels)
        end   =  start + len(wave)/Npanels
        top = median(fnu[start:end])*1.1 + jrr.util.mad(fnu[start:end])*5
        plt.ylim(0, top)  # trying to fix weird autoscaling from bad pixels
#        plt.autoscale(enable=False, axis='y', tight=True)

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
        if(k == Npanels/Nfigs-1):
            subit.set_xlabel(ur"observed-frame vacuum wavelength (\u00c5)")
            plt.ylabel('fnu') # fnu in cgs units: erg/s/cm^2/Hz
        if(k == 0):
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
print "Ran ", this_script
print "   Plotted spectrum: ", inspec
print "   Used linelist: ", linelist
print "   Generated PDF:  ", the_pdf
print "   DONE."
