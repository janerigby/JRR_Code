import pandas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#import seaborn as sns

# Before calling this, ran idl> .run .run print_ages_Z.pro, print_ages_Z, and then dumped output to print_ages_Z.out
# Had to use IDL because I can't figure out how to filter binary fits tables in astropy 


figsize = (8,10)
the_pdf = "print_ages_Z.pdf"
pp = PdfPages(the_pdf)  # output

infile = "print_ages_Z.out"
df = pandas.read_table(infile, delim_whitespace=True, comment="#", header=0, names=('filename', 'frac_light', 'uncert', 'age', 'metallicity'))

#df['filename'] = df['filename'].str.replace("/Volumes/Apps_and_Docs/jrrigby1/Dropbox/MagE_atlas/Contrib/S99/", "")
df['filename'] = df['filename'].str.replace("-continuum-properties.fits", "")
files = pandas.Series.unique(df.filename)

Npages = 2 # Number of pages
Ncol = 2 #  metallicity and age
Nrow = int(np.ceil(len(files) / (Npages*1.0)))
plotnum = 1 # initialize

fig = plt.figure(figsize=figsize)
for file in files :
    subset = df[df['filename'].eq(file)]
    #print subset.shape
    print subset.head(2)
    # now need to bin by age, metallicity
    Zbins = [0.01, 0.2, 0.4, 1.0, 2.0] # fraction of solar
    groupbyZ = subset.groupby(['metallicity'])['frac_light'].sum()
    ax = fig.add_subplot(Nrow, Ncol, plotnum)
    #plt.scatter(subset.metallicity, subset.frac_light)
    #plt.scatter(groupbyZ.index, groupbyZ.values, color='r')
    width = 0.1
    plt.bar(groupbyZ.index - width/2., groupbyZ.values, width=0.1) 
    plt.xlabel("metallicity (Fraction of solar)")
    plt.ylabel("light fraction")
    plt.ylim(0,1)
    plt.xlim(-0.05,2.1)
    plt.annotate(file, (0.4,0.8), xycoords="axes fraction", fontsize=12)
    plotnum += 1

    ax = fig.add_subplot(Nrow, Ncol, plotnum)
    groupbyt = subset.groupby(['age'])['frac_light'].sum()
    #plt.scatter(subset.age, subset.frac_light)
    #plt.scatter(groupbyt.index, groupbyt.values, color='r')
    width = 1.0
    plt.bar(groupbyt.index / 1.0E6, groupbyt.values, width=1) 
    plt.xlabel("age (Myr)")
    plt.ylabel("light fraction")
    plt.ylim(0,1)
    plt.xlim(-0.05,42)
    plt.annotate(file, (0.4,0.8), xycoords="axes fraction", fontsize=12)
    plotnum +=1
    print plotnum

    if plotnum == 2*Nrow+1 :
        print "Reached the end of a page", plotnum
        pp.savefig()  # save the plot
        plotnum = 1   # reset
        fig = plt.figure(figsize=figsize)  # start a new figure

if plotnum != 1 :  pp.savefig() # save the last plot
pp.close()

