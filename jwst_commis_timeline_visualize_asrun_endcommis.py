from  matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter, MultipleLocator
import pandas
import numpy as np

''' This python script makes a graphical visualization of the JWST commissioning timeline.
It reads a short list of events from an excel file, and plots them.  jrigby, Jan 2020.  
jrigby updating May 2022 to as-run in commissioning'''


def stem_and_label(xs, ys, labels, labelsize, bottom=0.7):
    plt.stem(xs, ys, use_line_collection=True, bottom=bottom, basefmt=' ')
    for ii, label in enumerate(labels):
        xoff = 0 ; yoff=0.05
        plt.annotate(label, xy=(xs[ii]+xoff, ys[ii]+yoff), fontsize=labelsize[ii])
    return(0)

plt.ion()
fig, ax = plt.subplots(figsize=(18,6))
fs0 = 14 ; fs1 = 14 ; fs2 = 20 ; fs3=24; fs4=30  # shortcuts for Font sizes
ys = [1, 2, 3, 4, 5]
y_label = 0.2


t_deploy  = 25.2
t_startOTE = 33.7
t_endOTE   = 116.0
t_end_SI   = 197. # commissioning ended on July 10 2022, when NIRCam coronagraphy was approved. # 180.

day  = np.arange(0,t_end_SI +0.01 ,1)
plt.plot(day, np.zeros_like(day), lw=4, color='k')

# Major stages.  Dates are hardcoded.
pheight = 0.7  ; lheight=0.1
patch1 = patches.Rectangle(xy=(0, 0), width=t_deploy, height=pheight, fill=True, color='C0')  # Deploy  #blue
patch0 = patches.Rectangle(xy=(t_deploy, 0), width= t_startOTE - t_deploy, height=pheight, fill=True, color='grey')  # Cooldown, grey
patch2 = patches.Rectangle(xy=(t_startOTE, 0), width=t_endOTE - t_startOTE, height=pheight, fill=True, color='C1')  # OTE commis
patch3 = patches.Rectangle(xy=(t_endOTE, 0), width=t_end_SI - t_endOTE, height=pheight, fill=True, color='C2')  # SI commis

ax.get_xaxis().set_major_locator(MultipleLocator(20))
plt.xlabel("Days after launch", fontsize=fs3)
plt.xticks(fontsize=fs2)
ax.xaxis.set_tick_params(width=2, length=10)
for patch in (patch1, patch2, patch3, patch0) :
    ax.add_patch(patch)
plt.annotate("Deploy", xy=(2, y_label), fontsize=fs4)
plt.annotate("cooling", xy=(25.1, y_label), fontsize=fs1)
plt.annotate("Telescope commissioning", xy=(37, y_label), fontsize=fs4)
plt.annotate("SI commissioning", xy=(130, y_label), fontsize=fs4)

plt.plot((132.2, 146.2), (2,2), ls='-', color='k') #draw line for thermal characterization

df = pandas.read_excel("events_endcommis.xlsx")  # Grab an excel file with a short list of events, time, and label sizes
df['fontsize'] = df['labsize']
#stem_and_label(df['day'], df['height'], df['event'], df['labsize'], bottom=pheight)   # PLANNED
stem_and_label(df['elapsed_day_actual'], df['height'], df['event'], df['labsize'], bottom=pheight)  # ACTUAL
#### Will need to update these to day_actual as events unfold***

# hot/cold thermal slews from 132.5d to 153d.  Hatch this region
plt.rcParams['hatch.linewidth'] = 3
thermalslew_begin = 132.2 # hardcoded  # started 5/6/2022
thermalslew_end   = 146.2 # hardcoded  # ended May 20 2022
thermalslew_width = thermalslew_end - thermalslew_begin
bar = ax.bar(x=thermalslew_begin, height=pheight, width=thermalslew_width, bottom=0, align='edge', fill=False, hatch='//', edgecolor='darkgrey')

plt.ylim(-0.2, 3.5)
plt.tight_layout()
# Hide the frame of the plot and the y axis
ax.set_frame_on(False)
plt.yticks([])

plt.annotate("Jane.Rigby@nasa.gov, at end of commissioning" , xycoords='figure fraction', xy=(0.7,0.03), fontsize=fs0)
fig.savefig("jwst_visual_commis_timeline_jrigby_asrun_endcommis.pdf")

#plt.annotate("Jane.Rigby@nasa.gov" , xycoords='figure fraction', xy=(0.7,0.03), fontsize=fs0)
#fig.savefig("jwst_visual_commis_timeline_jrigby_nodate.pdf")
