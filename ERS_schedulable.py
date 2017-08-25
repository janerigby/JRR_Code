import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

def midpoint_date(dateA, dateB) :
    return (dateA + (dateB - dateA)/2)

def draw_ERS_box(ax, ERSwin, ymin, ymax, color='grey') :
    ax.fill( [ERSwin[0], ERSwin[1], ERSwin[1], ERSwin[0]], [ymin, ymin, ymax, ymax], color, alpha=0.2, edgecolor='k')
    return(0)

def adjust_date_axis(ax) :
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y')
    ax.xaxis.set_major_locator(months)
    ax.xaxis.set_major_formatter(yearsFmt)
#    ax.xaxis.set_minor_locator(months)
    return(0)

def adjust_date_axis2(ax) :
    years = mdates.YearLocator()
    months = mdates.MonthLocator()
    monthsFmt = mdates.DateFormatter('%b') 
    yearsFmt = mdates.DateFormatter('\n\n%Y')  # add some space for the year label
    ax.xaxis.set_minor_locator(months)
    ax.xaxis.set_minor_formatter(monthsFmt)
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=90)
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    return(0)

noL = '_nolegend_'
# The 4 high-priority targets.  
win1 = (datetime.date(2019, 2,25), datetime.date(2019, 9,19))
win2 = (datetime.date(2020, 2,25), datetime.date(2020, 9,19))
win3 = (datetime.date(2019, 4,27), datetime.date(2019, 6,24))
win4 = (datetime.date(2019,12,13), datetime.date(2020, 2, 6))
win5 = (datetime.date(2019, 7, 3), datetime.date(2020, 1,29))
win7 = (datetime.date(2019, 7,24), datetime.date(2020, 2,16))
#

# Alternate SPT targets in case of later start to Cycle 1
awin1 = (datetime.date(2019, 4,22), datetime.date(2019, 6,28))
awin2 = (datetime.date(2019, 9, 3), datetime.date(2019,11, 7))
awin3 = (datetime.date(2019, 4,25), datetime.date(2019, 7, 1))
awin4 = (datetime.date(2019, 9, 5), datetime.date(2019,11,10))
awin5 = (datetime.date(2019, 1, 1), datetime.date(2019, 3,23))
awin6 = (datetime.date(2019, 8,21), datetime.date(2020, 3,23))
#

# RCSO327
awin7 = (datetime.date(2019, 1, 1), datetime.date(2019, 2,11))
awin8 = (datetime.date(2019, 8, 2), datetime.date(2019,10, 4))
awin9 = (datetime.date(2019,12,13), datetime.date(2020, 2,11))
# Planck arc
awin10 =(datetime.date(2019, 2,28), datetime.date(2019, 9,21))
awin11 =(datetime.date(2020, 2,28), datetime.date(2020, 4, 1))
# S1110
awin12 =(datetime.date(2019, 1, 1), datetime.date(2019, 5,12))
awin13 =(datetime.date(2019,10,28), datetime.date(2020, 4, 1))
# S1050
awin14 =(datetime.date(2019, 1, 1), datetime.date(2019, 1,19))
awin15 =(datetime.date(2019, 4,19), datetime.date(2019, 6,10))
awin16 =(datetime.date(2019,12, 1), datetime.date(2020, 1,20))


# New primary targets, after SPT folks shifted to lower redshift
windows = (win1,                win2, win3,                win4,  awin1,                 awin2, awin3,                 awin4)
namesz  = ("S1723+34 (z=1.32)", noL,  "S1226+21 (z=2.92)", noL,   "SPT2134-50 (z=2.78)", noL,   "SPT2147-50 (z=3.76)", noL)   
offsets = (0,                   0,    2,                   2,     1,                     1,     3,                     3)   
color   = ('b',                'b',   'b',                 'b',   'r',                  'r',    'r',                   'r')    

# New alternate targets, moving two highest-z sources from primary to alternate
awindows = ( win5,                win7,               awin5,              awin6,               awin7,              awin8,             awin9,  awin10,               awin11, awin12,           awin13, awin14,           awin15,           awin16)  
anamesz  = ( "SPT0346 (z=5.66)",  "SPT0418 (z=4.22)", 'SPT0532 (z=3.40)', 'SPT0532 (z=3.40)', 'RCS0327 (z=1.70)', 'RCS0327 (z=1.70)', noL,    'PSZ G311.65-18 (z=2.37)', noL,    "S1110 (z=2.48)", noL,    "S1050 (z=3.63)", "S1050 (z=3.63)", noL)     
aoffsets = ( 3,                   2.5,                1.5,                1.5,                 0.0,                0.0,               0.0,    0.5,                  0.5,      1.0,            1.0,     2.0,              2.0,             2.0)       
acolor   = ( 'r',                 'r',                'r',                'r',                 'b',                'b',               'b',    'b',                  'b',    'b',              'b',    'b',              'b',              'b')     


ERSwin = (datetime.date(2019, 4, 1), datetime.date(2019, 8, 31))
annotate_xpos1 = datetime.date(2019,5,1)
annotate_xpos2 = datetime.date(2019,9,10)
ymin = -0.5
ymax = 3.3
xmin = datetime.date(2019, 3,1)
xmax = datetime.date(2020, 3,1)

fig = plt.figure(figsize=(12,3))
ax1 = fig.add_subplot(1, 2, 1)
for ii, thiswin in enumerate(windows) :
    ax1.plot(thiswin, (offsets[ii], offsets[ii]), label=namesz[ii], color=color[ii], lw=3)
    print "Primary target", namesz[ii], "window duration:", thiswin[1] - thiswin[0]
    if namesz[ii] != noL :
        #ax1.annotate(namesz[ii], xy=(midpoint_date(*thiswin), 0.05+offsets[ii]), xycoords='data', ha='center')
        ax1.annotate(namesz[ii], xy=(annotate_xpos1, 0.05+offsets[ii]), xycoords='data', ha='left')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, ymax)
draw_ERS_box(ax1, ERSwin, ymin, ymax)


# Draw some alternative boxes, temporary
#draw_ERS_box(ax1, (datetime.date(2019, 7, 1), datetime.date(2019, 11, 30)), ymin, ymax, color='yellow')
draw_ERS_box(ax1, (datetime.date(2019, 10, 1), datetime.date(2020, 2, 28)), ymin, ymax, color='pink')

ax1.set_title("Primary targets")
adjust_date_axis2(ax1) 

ax2 = fig.add_subplot(1, 2, 2)
for ii, thiswin in enumerate(awindows) :
    print "Alt target", anamesz[ii], "window duration:", thiswin[1] - thiswin[0]
    ax2.plot(thiswin, (aoffsets[ii], aoffsets[ii]), label=anamesz[ii], color=acolor[ii], lw=3)
    if anamesz[ii] != noL :
        #ax2.annotate(anamesz[ii], xy=(midpoint_date(*thiswin), 0.05+aoffsets[ii]), xycoords='data', ha='center')
        ax2.annotate(anamesz[ii], xy=(annotate_xpos2, 0.05+aoffsets[ii]), xycoords='data', ha='left')
ax2.set_title("Alternate targets")
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(ymin, ymax)
draw_ERS_box(ax2, ERSwin, ymin, ymax)
adjust_date_axis2(ax2)

ax1.get_yaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)

fig.subplots_adjust(left=0.03, bottom=0.21, right=0.98, top=0.88, wspace=0.08, hspace=0.2)

plt.show()

