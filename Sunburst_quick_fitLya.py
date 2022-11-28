# Quick estimate of Lya flux in the Sunburst arc, knot 1 (MagE slit M-7)
import jrr
import matplotlib.pyplot as plt
import warnings
from astropy.stats import gaussian_fwhm_to_sigma
warnings.filterwarnings('ignore')

mage_mode = 'released'
scalefactor = 1E-17 # optimizer does not like small numbers
outfilename = "/Users/jrrigby1/SCIENCE/Lensed-LBGs/Mage/Analysis/Fit-lines/Lya_Mage_Sunburst.txt"
outfile = open(outfilename, "w")
outfile.write("#label  flux_Lya_directsum    flux_Lya_gaussfit  gauss0 gauss1 gauss2 gauss3\n")

for ii in (0, 2, 3, 4, 5, 6, 7, 8 ,9) : # there's no M-1
    label = 'sunburst_M-'+str(ii)
    (sp, resoln, dresoln, LL, zz_syst) = jrr.mage.wrap_open_spectrum(label, mage_mode)
    sp['flam_big'] = sp['flam']  / scalefactor
    sp['flam_u_big'] = sp['flam_u'] / scalefactor
    sp['flam_autocont_big'] = sp['flam_autocont']  / scalefactor
    plt.ion()
    ax = sp.plot(x='wave', y='flam_big', color='k')
    Lya_region = sp['wave'].between(4091., 4114.)
    plt.xlim(3900,4200)
    plt.ylim(-1, 80)
    plt.ylabel("flam in units of" + str(scalefactor) + "erg/s/cm^2/A")
    guesscont = (sp['flam_autocont'].loc[sp['wave'].between(3900.,4075.)]).median()/scalefactor
    guess_pars = (1E-16/scalefactor, 4100., 3., guesscont)
    jrr.spec.onegaus(sp['wave'], *guess_pars)
    guess_y = jrr.spec.onegaus(sp['wave'], *guess_pars)
    plt.plot(sp['wave'], guess_y, color='green')
    plt.plot(sp['wave'],  sp['flam_autocont']/scalefactor, color='blue')
    popt, fit_y = jrr.spec.fit_gaussian_fixedcont(sp, guess_pars[0:-1], contlevel=guesscont, colf='flam_big', method='lm', sigma=sp['flam_u_big'])
    plt.plot(sp['wave'], fit_y, color='red')
    plt.title(label)
    net_flux_Lya = (sp['flam'].loc[Lya_region] - sp['flam_autocont'].loc[Lya_region]).sum() * (sp['disp'].loc[Lya_region]).median()
    print("M-", ii, net_flux_Lya)
    outfile.write(   " ".join((label, str(net_flux_Lya), str(jrr.spec.sum_of_gaussian(popt)*scalefactor), str(popt[0]*scalefactor), str(popt[1]), str(popt[2]), str(popt[3]*scalefactor), "\n")))
outfile.close()

