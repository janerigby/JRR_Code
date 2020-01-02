from  matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter, MultipleLocator
import pandas
import numpy as np

''' This is a quick python script to make a graphical visualization of the JWST commissioning timeline.
It reads a short list of events from an excel file, and plots them.  jrigby, Jan 2020'''

def stem_and_label(xs, ys, labels, labelsize, bottom=0.7):
    plt.stem(xs, ys, use_line_collection=True, bottom=bottom, basefmt=' ')
    for ii, label in enumerate(labels):
        xoff = 0 ; yoff=0.1
        plt.annotate(label, xy=(xs[ii]+xoff, ys[ii]+yoff), fontsize=labelsize[ii])
    return(0)

plt.ion()
fig, ax = plt.subplots(figsize=(18,6))
fs0 = 14 ; fs1 = 16 ; fs2 = 20 ; fs3=24; fs4=30  # shortcuts for Font sizes
ys = [1, 2, 3, 4, 5]
y_label = 0.2

day  = np.arange(0,181,1)
plt.plot(day, np.zeros_like(day), lw=4, color='k')

# Major stages.  Dates are hardcoded.
pheight = 0.7  ; lheight=0.1
patch1 = patches.Rectangle(xy=(0, 0), width=25.3, height=pheight, fill=True, color='C0')  # Deploy
patch0 = patches.Rectangle(xy=(25.3, 0), width=35-25.3, height=pheight, fill=True, color='grey')  # Cooldown
patch2 = patches.Rectangle(xy=(35, 0), width=89, height=pheight, fill=True, color='C1')  # OTE commis
patch3 = patches.Rectangle(xy=(35+89, 0), width=180-35-89, height=pheight, fill=True, color='C2')  # OTE commis

ax.get_xaxis().set_major_locator(MultipleLocator(20))
plt.xlabel("Days after launch", fontsize=fs3)
plt.xticks(fontsize=fs2)
ax.xaxis.set_tick_params(width=2, length=10)
for patch in (patch1, patch2, patch3, patch0) :
    ax.add_patch(patch)
plt.annotate("Deploy", xy=(5, y_label), fontsize=fs4)
plt.annotate("cooling", xy=(25.5, y_label), fontsize=fs1)
plt.annotate("Telescope commissioning", xy=(50, y_label), fontsize=fs4)
plt.annotate("SI commissioning", xy=(130, y_label), fontsize=fs4)

df = pandas.read_excel("events.xlsx")  # Grab an excel file with a short list of events, time, and label sizes
df['fontsize'] = df['labsize']
stem_and_label(df['day'], df['height'], df['event'], df['labsize'], bottom=pheight)

# hot/cold thermal slews from 132.5d to 153d.  Hatch this region
plt.rcParams['hatch.linewidth'] = 3
bar = ax.bar(x=132, height=pheight, width=21, bottom=0, align='edge', fill=False, hatch='//', edgecolor='darkgrey')

plt.ylim(-0.2, 3.5)
plt.tight_layout()
# Hide the frame of the plot and the y axis
ax.set_frame_on(False)
plt.yticks([])
fig.savefig("jwst_visual_commis_timeline_jrigby.pdf")
