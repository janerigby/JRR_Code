import asciitable
import pylab

inE = '/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Combined-spectra/RCS0327/rcs0327-knotE-allre.txt'

data = asciitable.read(inE)
pylab.plot(data.col1, data.col2)  #plot comes from matplotlib
show()

