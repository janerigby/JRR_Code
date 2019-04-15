from __future__ import print_function
# plot a figure with multiple panels
from builtins import range
import matplotlib
import pylab
import numpy
import asciitable
import sys
import os
import re # regular expressions
import lineid_plot
from subprocess import call
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('multipage.pdf')

milk = False  # Running on milk?
plotU = False  # Which spectrum to plot?  Knot E or Knot U
plotE = True

if(milk):
    inE = 'Data/rcs0327-knotE-allre.txt'  # local on milk
    inU = 'Data/rcs0327-knotU.txt'  # local on milk
    plotsize = (10,13)
else:
    inE=  "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Mage/Combined-spectra/RCS0327/rcs0327-knotE-allre.txt"
    inU = "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Mage/Combined-spectra/RCS0327/rcs0327-knotU.txt"
    plotsize = (11,16)
        
Npanels = 24  # N panels total
Nfigs   = 4   # divided into this many pages of plots
zabs = 1.702 # bulk redshift of ISM, for "rest wave" axis

if(plotE):
    contfile  = re.sub(".txt", ".cont", inE)  # Read continuum files.  Same dir, rootname; extension is .cont.  
    infile = inE
    lincolor = 'k'
    errcolor = '0.5'

if(plotU):
    contfile = re.sub(".txt", ".cont", inU)
    infile = inU
    lincolor = 'b'
    errcolor = 'c'

def onclick(event):  # Setting up interactive clicking.  Right now, just prints location.  Need to add gaussian
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        event.button, event.x, event.y, event.xdata, event.ydata))

# Need a linelist.  concat-linelist.pl will grab several different linelists and merge, w appropriate redshifts for different components
call("/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Mage/RCS0327/Plot-all/concat-linelist.pl")  # Sort the linelists, as plotIDs requires

lines = asciitable.read("plotted.linelist", Reader=asciitable.NoHeader, guess=False)
Lwav = lines.col1 * (1.0 + lines.col7)  # redshift the lines
Lcolor = lines.col6
LID = []  # Grab nicely named IDs.  Can drop the waves if they're in the way
for i in range(0,len(lines)):
    LID.append( lines.col2[i]) # + " " + str(lines.col3[i]))
#sys.exit()

# loadtext actually comes from numpy, which pylab has.
wave,fnu,sig = pylab.loadtxt(infile, skiprows=3, usecols=(0,1,2), unpack=True)

# for plotting purposes, multiply by scale factor to bring to order unity
factor = 1E29
fnu *= factor
sig *= factor


# If continuum file exists, grab it, for later plotting, fitting of EW.
if(os.path.exists(contfile)):
    pixcont, wavcont, concont = numpy.loadtxt(contfile,comments="#", usecols=(0,1,2), unpack=True)
    concont *=factor

# Hard-code the pairs of resonant absorption and fine structure emission I want to plot
# Just need to give an approx central wavelength that will capture the group
#            SiII   SiII   SiII  #FeII   
locations = [1285., 1530., 1812.] #, 1616.] 
win=70

#locations = [2365.,  2605.]  # FeII
#win = 80.  # size of wavelength window to plot, in Angstroms, rest-frame

for p in range(0,len(locations)):
    wstart = (locations[p] - win/2.)*(1.0+zabs)  # plot first in observed frame
    wend   = (locations[p] + win/2.)*(1.0+zabs) 
    
    fig    = pylab.figure(num=1, figsize=plotsize)  #from matplotlib    
    format = len(locations)*100 + 11   # makes label: 311 for 3 rows, 1 col, start at 1
    subit = host_subplot(format + p)
    inwin = numpy.where ( (wave > wstart) & (wave < wend))
    pylab.plot(wave[inwin], fnu[inwin], lincolor, linestyle='steps')   # plot spectrum
    pylab.plot(wave[inwin], sig[inwin], errcolor, linestyle='steps')   # plot 1 sigma uncert spectrum
    pylab.autoscale(axis='x', tight=True)

    # If continuum file exists, plot the continuum
    if(os.path.exists(contfile)):
        this = numpy.where((wavcont > wstart) & (wavcont < wend))
        #print "Debug ", wave[start], wave[end]
        #print "Debug", this[0]
#jrr        pylab.plot(wavcont[this[0]], concont[this[0]], 'y', linestyle='steps') # plot this portion of the continuum

    # lineid changes indexes, so can't reliably change color.  Instead, let's try my own simpleline IDer
    for i in range(0, len(LID)):
        temp = pylab.ylim()[1] * 0.93
        if ((Lwav[i] > pylab.xlim()[0]) and (Lwav[i] < pylab.xlim()[1])):  #don't plot if in the margins
            #pylab.text(Lwav[i], temp, LID[i], rotation=90, ha='center', va='center', color=Lcolor[i])
            pylab.annotate(LID[i], xy=(Lwav[i], temp*0.7), xytext=(Lwav[i], temp), color=Lcolor[i], rotation=90, ha='center', va='center', arrowprops=dict(color=Lcolor[i], shrink=0.00,width=0.5, headwidth=0))
            

    upper = subit.twiny()  # make upper x-axis in rest wave of abs lines
    upper.set_xlim(wstart/(1.0+zabs), wend/(1.0+zabs))

    majorLocator   = MultipleLocator(50)  # put subticks on lower x-axis
    minorLocator   = MultipleLocator(10)
    subit.xaxis.set_major_locator(majorLocator)
    subit.xaxis.set_minor_locator(minorLocator)
    subit.xaxis.tick_bottom()  # don't let lower ticks be mirrored  on upper axis

    subit.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=False, Symmetric=False))
    if(p == len(locations)-1):
        subit.set_xlabel(ur"observed-frame vacuum wavelength (\u00c5)", fontsize=24)
        pylab.ylabel('fnu') # fnu in cgs units: erg/s/cm^2/Hz
    if(p == 0):
        upper.set_xlabel(ur"rest-frame vacuum wavelength (\u00c5)", fontsize=24)
            
pp.savefig()    
fig.canvas.draw()
fig.canvas.mpl_connect('button_press_event', onclick)  # This is called an event connection.
#   pylab.suptitle('global title')

pp.close()
