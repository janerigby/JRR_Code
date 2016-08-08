import asciitable
import pylab
import lineid_plot
import matplotlib
from matplotlib import pyplot as plt

inE= "/Volumes/Apps_and_Docs/WORK/Lensed-LBGs/Combined-spectra/RCS0327/rcs0327-knotE-allre.txt"
#inE = 'Data/rcs0327-knotE-allre.txt' MILK
data = asciitable.read(inE)
pylab.plot(data.col1, data.col2)  #plot comes from matplotlib

# plot an example linelist
linelist = 'Linelists/Christy/wind.xidl'
lines    = asciitable.read(linelist)
ztemp = 1.7016
obslines = lines.col1 * (1.0+ztemp)
lineid_plot.plot_line_ids(data.col1,data.col2, obslines, lines.col2, ax=plt.gca())

ax = pylab.gca()
a = ax.findobj(matplotlib.text.Annotation) 
for i in a:
    i.set_color("red")
ax.figure.canvas.draw() # update

for i in range(0, len(obslines)):
    temp = plt.ylim()[1] * 0.93
    plt.text(obslines[i], temp, lines.col2[i], rotation=90, ha='center', va='center') 
