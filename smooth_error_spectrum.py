import jrr
import numpy as np
import pandas
from matplotlib import pyplot as plt
mage_mode = "reduction"
#mage_mode = "released"
import astropy.convolution
from astropy.stats import gaussian_sigma_to_fwhm

labels = ['rcs0327-E'] #['S1226+2152'] #'rcs0327-B', 'S0957+0509', 'S2111-0114']
(specs) = jrr.mage.getlist_labels(mage_mode, labels)
#(specs) = jrr.mage.getlist_wcont(mage_mode)
print "Grabbed this many spectra:", len(specs)
(spec_path, line_path) = jrr.mage.getpath(mage_mode)

for ii in range(0, len(specs)) :                  #nfnu_stack[ii] will be ii spectrum
    label     = specs['short_label'][ii]
    filename  = specs['filename'][ii]
    zz =  specs['z_neb'][ii]
    (sp, resoln, dresoln)  = jrr.mage.open_spectrum(filename, zz, mage_mode)
    linelist = jrr.mage.get_linelist_name(filename, line_path)   # convenience function
    (LL, zz1) = jrr.mage.get_linelist(linelist)  # Load the linelist into Pandas data frame

    # Automatically fit continuum.  results written to sp.fnu_autocont, sp.flam_autocont.
    jrr.mage.auto_fit_cont(sp, LL, zz)  # Fit the continuum.
    # Will want to run this on S99, and on the MagE data, for the photospehric lines, in the same way.
    good = sp[~sp.badmask & ~sp.linemask]
    wavelo = good.wave.iloc[0]
    wavehi = good.wave.iloc[-1]
    
    Nfwhm = 1. * gaussian_sigma_to_fwhm    # Number of resolution FWHMs in each bin.
    #Using this b/c Ayan is using +-2 sigma to measure direct sum FWHMs
    nmax = int(np.ceil(np.log10(wavehi/wavelo) / np.log10(1. + Nfwhm / resoln)))
    nn = np.arange(0,nmax,1)
    binleft = wavelo*(1.0+Nfwhm/resoln)**nn  # wavelength array to bin to.  Binned by Nfwhm resoln elements.
    bincen = pandas.Series(binleft[0:-1] * (1.0 + Nfwhm/resoln/2.0))
    good['cut'] = pandas.cut(good.wave, binleft, right=False)
    good['temp'] = good['fnu'] - good['fnu_autocont']
    group1 = good.groupby(['cut'], )['temp'].sum()  # sum over bins.
    plt.plot(bincen, group1.values, color='r') # sums
    plt.plot(good.wave, good.fnu)

    group2 = group1.reset_index()    
    Nformad = 11 # number of measurements to consider std
    
    mad     = group2.groupby(group2.index / Nformad).mad()
    outbins = bincen.groupby(bincen.index / Nformad).median()
    #print outbins, mad
    plt.scatter(outbins, mad)
    rolling = mad.rolling(window=10).median()
    plt.plot(outbins, rolling, color='k')
    plt.ylim(0,10E-29)
    

    # NOTE:  stoppped working on this b/c Ayan has imported Schneider et al. 1993's EW significance algorithm
    # In principal, this should work, but I'm not going to keep working on it.  Better to cite Schneider et al.
    # and move on w our lives.  jrigby, july 2016
