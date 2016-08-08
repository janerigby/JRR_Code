# plot a figure with multiple panels
import matplotlib
#matplotlib.use('MacOSX')
import pylab
import numpy
import asciitable
import sys
import os
import re # regular expressions
import lineid_plot

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

 # JRR simplify, just make 1 plot
start =  2600
end = 3700
pylab.plot(wave[start:end], fnu[start:end], lincolor, linestyle='steps')   # plot spectrum
pylab.plot(wave[start:end], sig[start:end], errcolor, linestyle='steps')   # plot 1 sigma uncert spectrum
pylab.autoscale(axis='x', tight=True)

# If continuum file exists, plot the continuum
if(os.path.exists(contfile)):
    this = numpy.where((wavcont > wave[start]) & (wavcont < wave[end]))
    pylab.plot(wavcont[this[0]], concont[this[0]], 'y', linestyle='steps') # plot this portion of the continuum

# Here is where simple linear interpolation of continuum should go, to have continuum values at every wavelength.
# Next, look up the EW formula.  
# Next, define boundaries, where EW goes positive again on either side.
# Next, sum up the EW for every pixel within the boundaries
# Finally, print it up.
