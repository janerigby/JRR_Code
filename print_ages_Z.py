import pandas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
#import seaborn as sns
import re

# Before calling this, ran idl> .run .run print_ages_Z.pro, print_ages_Z, and then dumped output to print_ages_Z.out
# Had to use IDL because I can't figure out how to filter binary fits tables in astropy 

figsize = (8,10)
the_pdf = "print_ages_Z.pdf"
pp = PdfPages(the_pdf)  # output

infile = "print_ages_Z.out"
df = pandas.read_table(infile, delim_whitespace=True, comment="#", names=('filename', 'frac_light', 'uncert', 'age', 'metallicity'))


#df['filename'] = df['filename'].str.replace("/Volumes/Apps_and_Docs/jrrigby1/Dropbox/MagE_atlas/Contrib/S99/", "")
df['filename'] = df['filename'].str.replace("-continuum-properties.fits", "")
df['age'] /= 1E6  # convert from yr to Myr
files = pandas.Series.unique(df.filename)

# Make a table out of the fits
aslatex = df.to_latex()
f = open('print_ages_Z.tex', 'w')
f.write(aslatex)  # Not sure I will use this in the paper, it's long and i don't know that we trust fits that far
f.close()


ticked =[0.0,0.3, 0.6,0.9]

Npages = 3 # Number of pages
Ncol = 2 #  metallicity and age
#Nrow = int(np.ceil(len(files) / (Npages*1.0)))
Nrow = 8 # JRR KLUDGE
plotnum = 1 # initialize
print "DEBUGGING Nrow Ncol Npage", Nrow, Ncol, Npages

fig = plt.figure(figsize=figsize)
for file in files :
    subset = df[df['filename'].eq(file)]
    #print subset.shape
    pretty_label = re.sub("_", " ", file)

    #print to a file
    print subset.head(2)
    # now need to bin by age, metallicity
    Zbins = [0.01, 0.2, 0.4, 1.0, 2.0] # fraction of solar
    groupbyZ = subset.groupby(['metallicity'])['frac_light'].sum()
    ax1 = fig.add_subplot(Nrow, Ncol, plotnum)
    #plt.scatter(subset.metallicity, subset.frac_light)
    #plt.scatter(groupbyZ.index, groupbyZ.values, color='r')
    width = 0.1
    plt.xlabel("metallicity (fraction of solar)")
    plt.bar(groupbyZ.index - width/2., groupbyZ.values, width=0.1) 
    plt.ylabel("light frac")
        
    plt.ylim(0,1)
    plt.xlim(-0.05,2.1)
    plt.annotate(pretty_label, (0.4,0.8), xycoords="axes fraction", fontsize=12)
    plt.yticks(ticked, ticked)
    plotnum += 1

    ax2 = fig.add_subplot(Nrow, Ncol, plotnum)
    groupbyt = subset.groupby(['age'])['frac_light'].sum()
    #plt.scatter(subset.age, subset.frac_light)
    #plt.scatter(groupbyt.index, groupbyt.values, color='r')
    width = 1.0
    plt.bar(groupbyt.index, groupbyt.values, width=1) 
    plt.xlabel("age (Myr)")
    #plt.ylabel("light frac")
    plt.ylim(0,1)
    plt.xlim(-0.05,42)
    plt.annotate(pretty_label, (0.4,0.8), xycoords="axes fraction", fontsize=12)
#    if plotnum < 2*Nrow-3 :  
#        ax1.xaxis.set_major_formatter(plt.NullFormatter())
#        ax2.xaxis.set_major_formatter(plt.NullFormatter())
    plt.yticks(ticked, ticked)
    plotnum +=1
    print plotnum

    if plotnum == 2*Nrow+1 or file == "S1226+2152":
        fig.subplots_adjust(hspace=0)
        print "Reached the end of a page", plotnum
        pp.savefig()  # save the plot
        plotnum = 1   # reset
        fig = plt.figure(figsize=figsize)  # start a new figure

if plotnum != 1 :
    fig.subplots_adjust(hspace=0)
    pp.savefig() # save the last plot

    
# Manually add the rollup of all individual galaxies
fig = plt.figure(figsize=figsize)  # start a new figure
not_stack = df[~df['filename'].str.contains('tack') & ~df['filename'].str.contains('teidel')]  # First, exclude stack
groupbyZ = not_stack.groupby(['metallicity'])['frac_light'].sum()
ax1 = fig.add_subplot(Nrow, Ncol, plotnum)
width = 0.1
plt.xlabel("metallicity (fraction of solar)")
plt.bar(groupbyZ.index - width/2., groupbyZ.values/not_stack.frac_light.sum(), width=0.1) 
plt.ylabel("light frac")
plt.ylim(0,1)
plt.xlim(-0.05,2.1)
plt.annotate("Full sample", (0.4,0.8), xycoords="axes fraction", fontsize=12)
plt.yticks(ticked, ticked)
plotnum += 1

ax2 = fig.add_subplot(Nrow, Ncol, plotnum)
groupbyt = not_stack.groupby(['age'])['frac_light'].sum()
width = 1.0
plt.bar(groupbyt.index, groupbyt.values/ not_stack.frac_light.sum(), width=1) 
plt.xlabel("age (Myr)")
#plt.ylabel("light frac")
plt.ylim(0,1)
plt.xlim(-0.05,42)
plt.annotate("Full sample", (0.4,0.8), xycoords="axes fraction", fontsize=12)
plt.yticks(ticked, ticked)
fig.subplots_adjust(hspace=0)
pp.savefig() # save the last plot
    
    
# Let's measure some properties of the different components, using pandas
print "These are the components that are old (40 Myr)"
oldcomps = df[df['age'].eq(40)].sort_values(by='metallicity')
print oldcomps.sort_values("frac_light")

age_breakdown = (8, 16)
young_comps = df[df['age'].lt(age_breakdown[0]) & ~df['filename'].str.contains('tack')].sort_values(by='metallicity')
print "Light-weighted metallicity for young, lt ", age_breakdown[0], "yr: ",
print (young_comps['frac_light'] * young_comps['metallicity']).sum() / young_comps['frac_light'].sum() 

midage_comps = df[df['age'].between(age_breakdown[0], age_breakdown[1])  & ~df['filename'].str.contains('tack')].sort_values(by='metallicity')
print "Light-weighted metallicity for middle-age components, between", age_breakdown, "yr: ",
print (midage_comps['frac_light'] * midage_comps['metallicity']).sum() / midage_comps['frac_light'].sum()

old_comps = df[df['age'].gt(age_breakdown[1]) & ~df['filename'].str.contains('tack')].sort_values(by='metallicity')
print "Light-weighted metallicity for old, gt ", age_breakdown[1], "yr: ",
print (old_comps['frac_light'] * old_comps['metallicity']).sum() / old_comps['frac_light'].sum() 


# Want a plot that shows relative contribution of both Z and age models.
groupby_tZ = not_stack.groupby(by=('metallicity', 'age'))
bytZ = groupby_tZ.sum()  # Sum over same age, metallicity
sum_to_norm = bytZ.frac_light.sum()
print "sum of frac_light:", sum_to_norm, "should be ~14"
pandas.set_option('display.multi_sparse', False)
print bytZ
#                   age                         metallicity
scatscale = 8000
fig = plt.figure(figsize=(8,8))
plt.scatter(bytZ.index.get_level_values(1), bytZ.index.get_level_values(0), s=((bytZ.frac_light.values/sum_to_norm)**2)*scatscale)
plt.xlabel("age (Myr)", fontsize=20)
plt.ylabel("Metallicity (fraction of solar)", fontsize=20)
plt.xlim(-0.5,41)
plt.ylim(-0.05,2.1)
plt.xticks(size=18)
plt.yticks(size=18)
# Show what size of scatter points means:
ptsize =np.array((0.14, 0.7, 1.4, 2.8))   # this is 1%, 5%, 10%, 20% of the sample
plt.scatter( np.ones_like(ptsize)*30, np.linspace(1.7,2,num=len(ptsize)), s=((ptsize/sum_to_norm)**2*scatscale) )
plt.annotate(" 1%", (31,1.68), xycoords="data", fontsize=14)
plt.annotate(" 5%", (31,1.78), xycoords="data", fontsize=14)
plt.annotate("10%", (31,1.88), xycoords="data", fontsize=14)
plt.annotate("20%", (31,1.98), xycoords="data", fontsize=14)
plt.plot( (27,35,35,27,27), (1.6,1.6,2.1,2.1,1.6), color='k', linewidth=2)    
pp.savefig()
pp.close()

# Do simple math that Rongmon suggested, to get numbers for the paper.
Zcuts = (0.01, 0.1, 0.5)
for Zcut in Zcuts :
    print "Fraction of sample's stellar light fit by >", Zcut, "solar metallicity: ",
    print not_stack[not_stack['metallicity'] > Zcut].frac_light.sum()  / not_stack.frac_light.sum()
agecuts = (10., 20., 30.)
for agecut in agecuts:
    print "Fraction of sample's stellar light fit by >", agecut, "Myr age:",
    print not_stack[not_stack['age'] > agecut].frac_light.sum()  / not_stack.frac_light.sum()
