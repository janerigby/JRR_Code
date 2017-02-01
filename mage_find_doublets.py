import jrr
import ayan
import numpy as np
import pandas
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
mage_mode = "released"

''' Code to find intervening doublets within the Megasaura MagE spectra.
1st, conducts a blind search of the spectra for significant absorption lines,
following methodology of Schneider et al. 1993 (HST QAL key project).
2nd, for each peak, tests whether there is a partner line at the expected position
were it a doublet.  Writes to a data frame.
jrigby Jan 2017
'''

def get_doublet_waves() :  # vacuum barycentric NIST
    MgII = np.array((2796.352, 2803.531))
    CIV  = np.array((1548.195, 1550.770))
    SiIV = np.array((1393.755, 1402.770))
    return(MgII, CIV, SiIV)

def make_empty_doublet_dataframe() :
    df = pandas.DataFrame(columns=('flag', 'gal', 'zz', 'doubname', 'wave1', 'wave2', 'EW1', 'EW2', 'snr1', 'snr2', 'in_forest'))
    return(df)
    
def find_lines_Schneider(sp, resoln, siglim=3., abs=True, delta=0.15) :
    # Blind search for absorption lines, following Schneider et al. 1993
    # Delta seems pretty damned arbitary, may bite me later.
    # Significant peaks identified as sp['peak']=True
    ayan.mage.calc_schneider_EW(sp, resoln, plotit=False)  # Calculate EW and EW limits
    sp['temp'] = False  # 1st pass, found a peak
    sp['peak'] = False  # 2nd pass, peak is significant
    maxtab, mintab = jrr.peakdet.peakdet(sp.W_interp,delta)  # Find peaks.
    if abs:   peak_ind =  [np.int(p[0]) for p in mintab] # The minima
    else:     peak_ind =  [np.int(p[0]) for p in maxtab] # The maxima
    sp['temp'].iloc[peak_ind] = True  # Convert back into pandas style
    # Choose only peaks that are greater than siglim significant
    if abs :  # If looking for absorption lines 
        subset = sp['temp']  &  (sp['W_interp'] < sp['W_u_interp'] * siglim * -1)
    else :    # If looking for emission lines
        subset = sp['temp']  &  (sp['W_interp'] > sp['W_u_interp'] * siglim)
    sp['peak'].ix[subset] = True  # Peaks now marked 
    print "FINDING PEAKS (this is slow), N peaks: ", sp['temp'].sum(),  "  significant peaks: ", sp['peak'].sum()
    return(0)

def plot_peakfinding(sp) :
    # Turned off, used to write/debug find_lines_Schneider
    plt.plot(sp.wave, sp.fnu/sp.fnu_autocont, color='black', label="fnu contnorm")
    plt.plot(sp.wave, sp.W_interp,   color='blue',  label="W_interp")
    plt.plot(sp.wave, sp.W_u_interp*siglim, color='green')
    plt.plot(sp.wave, sp.W_u_interp*siglim*-1, color='green')
    subset = sp['peak']
    plt.plot(sp['wave'].ix[subset], sp['W_interp'].ix[subset], 'o', label="try", color="pink")
    plt.xlim(4000,7000)
    plt.title(thisgal)
    plt.ylim(-3,3)
    plt.legend()
    plt.show()
    return(0)


def test_candidate_doublets(sp, zz_syst, resoln, doublet, doubname, ylims=(-2,2)) :
    '''For a given doublet, step through each peak, check whether there is a second peak
    at the expected doublet separation.  Outputs a dataframe of candidates that user marked yes or possible.'''
    plotwin = 20.  #plotting width, A
    df = make_empty_doublet_dataframe()
    counter = 0
    subset = sp[sp['peak']]   # subset of sp where find_lines_Schneider() found a peak
    for candidate in subset.index :   # For each peak, test whether there is a 2nd transition
        testz = sp.ix[candidate]['wave'] / doublet[0] - 1.0
        if testz <= zz_syst :
            dwave = doublet[1]*(1.+testz) / 2. / resoln  # Search +- 1 HWHM
            searchreg = sp['wave'].between((doublet[1]*(1.0+testz)-dwave), (doublet[1]*(1.0+testz)+dwave)) & sp['peak']
            if sp[searchreg]['peak'].sum() :
                plt.clf()
                in_forest = sp.ix[candidate]['wave'] < zz_syst * 1216.  #  In Lya Forest
                print "Possible ", doubname, " doublet at", testz, doublet[0]*(1+testz), doublet[1]*(1+testz)
                plt.step(sp.wave, sp.W_interp, color='blue')
                plt.title(thisgal + " " + doubname + "?")
                plt.step(sp.wave, sp.fnu/sp.fnu_autocont, color='black', label="fnu contnorm")
                plt.step(sp.wave, sp.fnu_u/sp.fnu_autocont, color='grey', label="fnu_u contnorm")
                plt.hlines( 1.0, 3000, 8500, colors='grey')
                plt.vlines( (doublet[0]*(1+testz), doublet[1]*(1+testz)), *ylims, colors='red')
                plt.step(sp.wave, sp.W_u_interp*siglim*-1, color='green')
                plt.xlim(doublet[0]*(1+testz) - plotwin, doublet[1]*(1+testz) + plotwin)
                plt.ylim(*ylims)
                jrr.mage.plot_linelist(LL, zz_syst)
                plt.draw()
                plt.pause(1)
                EW1  = sp.ix[candidate].W_interp  # Using Schneider method to estimate EW
                EW2 = float(np.min(sp[searchreg]['W_interp']))
                snr1 = sp.ix[candidate].W_interp / sp.ix[candidate].W_u_interp *-1.
                snr2 = float(np.max(sp[searchreg]['W_interp'] / sp[searchreg]['W_u_interp'] * -1.))
                user_flag = (raw_input("Is this real? y for yes, n for no, p for possibly, b for yes but blended:")) 
                if user_flag == "y" or user_flag == "p" or user_flag == "b":
                    df.loc[counter] = (user_flag, thisgal, testz, doubname, doublet[0]*(1+testz), doublet[1]*(1+testz), EW1, EW2, snr1, snr2, in_forest)
                    pp.savefig()
                    counter += 1
    return(df)

# Actually run things
(MgII, CIV, SiIV) = get_doublet_waves()
siglim=5
ylims = (-2,2)
the_pdf = "found_doublets_SNR" + str(siglim) + "cleanup.pdf"
outfile = "found_doublets_SNR" + str(siglim) + "cleanup.txt"
pp = PdfPages(the_pdf)
plt.clf()
df = make_empty_doublet_dataframe()
#speclist = jrr.mage.wrap_getlist(mage_mode, which_list="wcont")
speclist = jrr.mage.wrap_getlist(mage_mode, which_list="labels", labels=("Cosmic~Eye",))
for thisgal in speclist.short_label :
    (sp, resoln, dresoln, LL, zz_syst) = jrr.mage.wrap_open_spectrum(thisgal, mage_mode, addS99=True)  # load spectrum
    find_lines_Schneider(sp, resoln, siglim=siglim, abs=True)  #identify all the >N sigma peaks in spectrum
    #plot_peakfinding(sp)  # for debugging
    df1 = test_candidate_doublets(sp, zz_syst, resoln, CIV,  "CIV",  ylims=ylims)
    df2 = test_candidate_doublets(sp, zz_syst, resoln, MgII, "MgII", ylims=ylims)
    df3 = test_candidate_doublets(sp, zz_syst, resoln, SiIV, "SiIV", ylims=ylims)
    df  = pandas.concat([df, df1, df2, df3]).reset_index(drop=True)
print df.head(20)
with open(outfile, 'w') as f:
    df.to_csv(f, sep='\t')
pp.close()


# Wrote a wrapper for Rongom to go analyze the doublet ratio
plt.close("all")
def doublet_ratio_partialcoverage_analyzer(doublet_file) :
    doublet_list = pandas.read_table(doublet_file, delim_whitespace=True, comment="#")
    for row in doublet_list.itertuples():
        plt.clf()
        (sp, resoln, dresoln, LL, zz_syst) = jrr.mage.wrap_open_spectrum(row.gal, mage_mode, addS99=True)  # load spectrum
        print "debugging, working on this row of doublet list: ", row
        plt.step(sp.wave, sp.fnu/sp.fnu_autocont, color='black')
        plt.step(sp.wave, sp.fnu_u/sp.fnu_autocont, color='grey')
        plt.step(sp.wave, sp.fnu/sp.fnu_s99model, color='blue')
        plt.xlim(row.wave1 - 30., row.wave2+30.)
        plt.ylim(-1,2)
        plt.vlines( (row.wave1, row.wave2), -3, 3, colors='red')
        # And Rongmon does magical fitting here.
        plt.pause(1)
    return(0)
    