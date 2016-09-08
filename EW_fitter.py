''' Equivalent width and flux fitter.  Given a spectrum, a continuum, a linelist, and a redshift,
fit a bunch of emission or absorption lines.
Started july 2016, Jane Rigby and Ayan Acharyya

README
Usage: python EW_fitter.py --<options>
<options> are:
--fcen FCEN ; FCEN = 1 if initial guess for gaussian centers are to be fixed, default 0(False)
--fcon FCON ; FCON = 1 if initial guess for continuum is to be fixed, default 1(True)
--dx DX ; DX = how many angstroms of wavelength should be plotted in one panel in the final plot, default = 300A
--only ONLY ; ONLY = 1 if final plot should have only those panels(patches of wavelength axis) where a line was fitted,
                else 0 to display all frames, default = 1
--vmax VMAX ; VMAX = in km/s, to set the maximum FWHM that can be fit as a line, default = 300km/s
--frame FRAME ; FRAME = index of the panel if you want to see specific panels, the panel indices can be found in top 
                left of each panel in the final plot, default plots all frames
--nbin NBIN ; NBIN = # of wavelength points to be binned together to calculate median and MAD binned fluxes and errors
--lines LINES ; LINES = 'emission' OR 'photospheric' depending on which lines you want to be fit, linelists can be 
                found in files labframe.shortlinelist_emission and labframe.shortlinelist_photospheric, respectively.
                default = emission
--fout FOUT ; FOUT = filename you want the output ASCII file to have, default = fitted_line_list.txt
--keepprev ; boolean option, if present then doesn't kill the previous matplotlib plots
--silent ; boolean option, if present then doesn't print whole bunch of info on console while running
--mymask ; boolean option, if present then MODIFIES the masking around skylines (as opposed to jrr.mage.flag_skylines), 
            using thefunction flag_skylines in ayan.mage
--check ; boolean option, if present then prints at the end the root mean square deviation of 1 sigma error between EWs
        computed via fitting and via summing
--allspec ; boolean option, if present then run the code over all the spectra files
--savepdf ; boolean option, if present then save the plots as pdfs
--noplot ; boolean option, if present then does not create any plot
--hide ; boolean option, if present then does not show plots at the end
--showbin ; boolean option, if present then show the binning (on either side of line/s used to calculate median and
            MAD errors) on the resultant plot
--fullmad ; boolean option, if present then calculate the MAD at every point on spectrum and add a column to dataframe
--showerr ; boolean option, if present then plots the Schneider precription EWs with errors
'''
import sys
sys.path.append('../')
import ayan.mage as m
import jrr
import ayan.splot_util as s
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
from  astropy.io import ascii
from matplotlib import pyplot as plt
mage_mode = "released"
import argparse as ap
from matplotlib.backends.backend_pdf import PdfPages

#-----------Main function starts------------------
parser = ap.ArgumentParser(description="Mage spectra fitting tool")
parser.add_argument("--shortlabel")
parser.add_argument("--fcen")
parser.add_argument("--fcon")
parser.add_argument("--dx")
parser.add_argument("--only")
parser.add_argument("--vmax")
parser.add_argument("--frame")
parser.add_argument("--nbin")
parser.add_argument("--lines")
parser.add_argument("--fout")
parser.add_argument('--keepprev', dest='keepprev', action='store_true')
parser.set_defaults(keepprev=False)
parser.add_argument('--silent', dest='silent', action='store_true')
parser.set_defaults(silent=False)
parser.add_argument('--mymask', dest='mymask', action='store_true')
parser.set_defaults(mymask=False)
parser.add_argument('--check', dest='check', action='store_true')
parser.set_defaults(check=False)
parser.add_argument('--allspec', dest='allspec', action='store_true')
parser.set_defaults(allspec=False)
parser.add_argument('--savepdf', dest='savepdf', action='store_true')
parser.set_defaults(savepdf=False)
parser.add_argument('--hide', dest='hide', action='store_true')
parser.set_defaults(hide=False)
parser.add_argument('--noplot', dest='noplot', action='store_true')
parser.set_defaults(noplot=False)
parser.add_argument('--showbin', dest='showbin', action='store_true')
parser.set_defaults(showbin=False)
parser.add_argument('--fullmad', dest='fullmad', action='store_true')
parser.set_defaults(fullmad=False)
parser.add_argument('--showerr', dest='showerr', action='store_true')
parser.set_defaults(showerr=False)
args, leftovers = parser.parse_known_args()
if args.dx is not None:
    dx = float(args.dx)
else:
    dx = 310.
if args.only is not None:
   display_only_success = args.only
else:
    display_only_success = 1
if args.frame is not None:
    frame = args.frame
else:
    frame = None
if args.shortlabel is not None:
    labels = [args.shortlabel]
else:
    labels = ['rcs0327-E']
if args.lines is not None:
    listname = str(args.lines)
else:
    listname = 'emission'
if args.fout is not None:
    fout = str(args.fout)+'.txt'
else:
    fout = 'fitted_line_list.txt'
if args.allspec:
    labels = [
    'rcs0327-B',\
    'rcs0327-E',\
    #'rcs0327-Ehires',\
    #'rcs0327-Elores',\
    'rcs0327-G',\
    'rcs0327-U',\
    #'rcs0327-BDEFim1',\
    #'rcs0327-counterarc',\
    'S0004-0103',\
    #'S0004-0103alongslit',\
    #'S0004-0103otherPA',\
    'S0033+0242',\
    'S0108+0624',\
    'S0900+2234',\
    'S0957+0509',\
    'S1050+0017',\
    'Horseshoe',\
    'S1226+2152',\
    #'S1226+2152hires',\
    #'S1226+2152lores',\
    'S1429+1202',\
    'S1458-0023',\
    'S1527+0652',\
    #'S1527+0652-fnt',\
    'S2111-0114',\
    'Cosmic~Eye',\
    #'S2243-0935',\
    ]
if not args.keepprev:
    plt.close('all')

#-------------------------------------------------------------------------
(specs) = jrr.mage.getlist_labels(mage_mode, labels)
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
line_table = pd.DataFrame(columns=['label', 'line_lab', 'obs_wav', 'rest_wave', 'type','EWr_fit','EWr_fit_u', 'EWr_sum', \
'EWr_sum_u', 'f_line','f_line_u', \
#'wt_mn_flux', 'onesig_err_wt_mn_flux', \
'med_bin_flux', 'mad_bin_flux', 'MAD_significance', 'EWr_3sig_lim_Schneider', 'fl_3sig_lim_Sndr', 'fit_cont','fit_fl','fit_cen', 'fit_cen_u', \
'fit_sig','zz','zz_u'])

for ii in range(0, len(specs)) :                  #nfnu_stack[ii] will be ii spectrum
    shortlabel     = specs['short_label'][ii]
    print 'Spectrum', (ii+1), 'of', len(specs),':', shortlabel #Debugging
    filename  = specs['filename'][ii]
    zz_sys = specs['z_syst'][ii] # from the new z_syst column in spectra_filename file
    zz_dic = {'EMISSION':specs['z_neb'][ii], 'FINESTR':specs['z_neb'][ii], 'PHOTOSPHERE': specs['z_stars'][ii] if specs['fl_st'][ii]==0 else specs['z_neb'][ii], 'ISM':specs['z_ISM'][ii], 'WIND':specs['z_ISM'][ii]}
    zz_err_dic = {'EMISSION':specs['sig_neb'][ii] if specs['fl_neb'][ii]==0 else specs['sig_ISM'][ii], 'FINESTR':specs['sig_neb'][ii] if specs['fl_neb'][ii]==0 else specs['sig_ISM'][ii], 'PHOTOSPHERE': specs['sig_st'][ii] if specs['fl_st'][ii]==0 else specs['sig_neb'][ii], 'ISM':specs['sig_ISM'][ii], 'WIND':specs['sig_ISM'][ii]}    
    if shortlabel != 'stack':
        (sp_orig, resoln, dresoln)  = jrr.mage.open_spectrum(filename, zz_sys, mage_mode)
    else:
        sp_orig = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=filename)
        resoln = 2884. # from that of rcs0327
        dresoln = 81. 
    #-----------fitting continuum----------------------------
    m.fit_autocont(sp_orig, zz_sys, line_path,filename)
    #-------masking sky lines-----------------
    if args.mymask:
        m.flag_skylines(sp_orig) #modified masking for skylines, as compared to jrr.mage.flag_skylines
    sp_orig = sp_orig[~sp_orig['badmask']].copy(deep=True)
    #-----calculating MAD error over entire spectrum--------
    if args.fullmad:
        m.calc_mad(sp_orig, resoln, 5)
        continue
    #------calculating the EW limits at every point following Schneider et al. 1993---------
    m.calc_schneider_EW(sp_orig, resoln, plotit=args.showerr)
    #makelist(linelist) #required if you need to make a new labframe.shortlinelist file
    line_full = m.getlist('labframe.shortlinelist_'+listname, zz_dic, zz_err_dic)
    #------------Preparing to plot----------------------------------------
    xstart = max(np.min(line_full.wave) - 50.,np.min(sp_orig.wave))
    xlast = min(np.max(line_full.wave) + 50.,np.max(sp_orig.wave))
    if frame is None:
        n_arr = np.arange(int(np.ceil((xlast-xstart)/dx))).tolist()
    else:
        n_arr = [int(ar) for ar in frame.split(',')] #Use this to display selected frame/s
    name = '/Users/acharyya/Dropbox/MagE_atlas/Contrib/EWs/'+listname+'/'+shortlabel+'-'+listname+'_fit'
    if args.savepdf:
        pdf = PdfPages(name+'.pdf')
    #---------pre check in which frames lines are available if display_only_success = 1---------
    #---------------------------just a visualisation thing-----------------------------------
    if display_only_success:
        for jj in n_arr:
            xmin = xstart + jj*dx
            xmax = min(xmin + dx, xlast)
            sp = sp_orig[sp_orig['wave'].between(xmin,xmax)]
            try:
                line = line_full[line_full['wave'].between(xmin*(1.+5./resoln), xmax*(1.-5./resoln))]
            except IndexError:
                continue
            if not len(line) > 0 or not line['wave'].between(np.min(sp.wave),np.max(sp.wave)).all():
                n_arr[jj] = np.ma.masked
        n_arr = np.ma.compressed(n_arr)
    #------------------------------------------------------------
    n = len(n_arr) #number of frames that would be displayed
    if not args.noplot:
        fig = plt.figure(figsize=(18+10/(n+1),(12 if n > 2 else n*3)))
        #fig = plt.figure(figsize=(14+8/(n+1),(9 if n > 2 else n*3)))
        plt.title(shortlabel + "  z=" + str(zz_sys)+'. Vertical lines legend: Blue=initial guess of center,'+\
        ' Red=fitted center, Green=no detection(< 3sigma), Black=unable to fit gaussian', y=1.02)
    for fc, jj in enumerate(n_arr):
        xmin = xstart + jj*dx
        xmax = min(xmin + dx, xlast)
        if not args.noplot:
            ax1 = fig.add_subplot(n,1,fc+1)
        sp = sp_orig[sp_orig['wave'].between(xmin,xmax)]
        try:
            line = line_full[line_full['wave'].between(xmin*(1.+5./resoln), xmax*(1.-5./resoln))]
        except IndexError:
            continue
        #------------Plot the results------------
        if not args.noplot:
            plt.step(sp.wave, sp.fnu, color='b')
            plt.step(sp.wave, sp.fnu_u, color='gray')
            plt.step(sp.wave, sp.fnu_cont, color='y')
            plt.plot(sp.wave, sp.fnu_autocont, color='k')
            plt.ylim(0, 1.2E-28)
            plt.xlim(xmin, xmax)
            plt.text(xmin+dx*0.05, ax1.get_ylim()[1]*0.8, 'Frame '+str(int(jj)))
        if not args.fullmad:
            m.fit_some_EWs(line, sp, resoln, shortlabel, line_table, dresoln, sp_orig, args=args) #calling line fitter
    
        if not args.noplot:
            ax2 = ax1.twiny()
            ax2.set_xlim(ax1.get_xlim())       
            ax2.set_xticklabels(np.round(np.divide(ax1.get_xticks(),(1.+zz_sys)),decimals=0))        
            labels2 = [item.get_text() for item in ax2.get_xticklabels()]
            ax2.set_xticks(np.concatenate([ax2.get_xticks(), line.wave*(1.+zz_sys)/(1.+line.zz)]))
            ax2.set_xticklabels(np.concatenate([labels2,np.array(line.label.values).tolist()]), rotation = 45, ha='left', fontsize='small')
            fig.subplots_adjust(hspace=0.7, top=0.95, bottom=0.05)
    if not args.hide:
        #fig.tight_layout()
        plt.show(block=False)
    if args.savepdf:
        pdf.savefig(fig)
        pdf.close()
    #fig.savefig(name+'.png')
#------------changing data types------------------------------
line_table.obs_wav = line_table.obs_wav.astype(np.float64)
line_table.rest_wave = line_table.rest_wave.astype(np.float64)
line_table.EWr_fit = line_table.EWr_fit.astype(np.float64)
line_table.EWr_fit_u = line_table.EWr_fit_u.astype(np.float64)
line_table.EWr_sum = line_table.EWr_sum.astype(np.float64)
line_table.EWr_sum_u = line_table.EWr_sum_u.astype(np.float64)
line_table.f_line = line_table.f_line.astype(np.float64)
line_table.f_line_u = line_table.f_line_u.astype(np.float64)
#line_table.wt_mn_flux = line_table.wt_mn_flux.astype(np.float64)
#line_table.onesig_err_wt_mn_flux = line_table.onesig_err_wt_mn_flux.astype(np.float64)
line_table.med_bin_flux = line_table.med_bin_flux.astype(np.float64)
line_table.mad_bin_flux = line_table.mad_bin_flux.astype(np.float64)
line_table.MAD_significance = line_table.MAD_significance.astype(np.float64)
line_table.EWr_3sig_lim_Schneider = line_table.EWr_3sig_lim_Schneider.astype(np.float64)
line_table.fl_3sig_lim_Sndr = line_table.fl_3sig_lim_Sndr.astype(np.float64)
line_table.fit_cont = line_table.fit_cont.astype(np.float64)
line_table.fit_fl = line_table.fit_fl.astype(np.float64)
line_table.fit_cen = line_table.fit_cen.astype(np.float64)
line_table.fit_cen_u = line_table.fit_cen_u.astype(np.float64)
line_table.fit_sig = line_table.fit_sig.astype(np.float64)
line_table.zz = line_table.zz.astype(np.float64)
line_table.zz_u = line_table.zz_u.astype(np.float64)
line_table['EW_signi']=3.*line_table['EWr_fit']/line_table['EWr_3sig_lim_Schneider']
line_table['EW_signi']=line_table['EW_signi'].map('{:,.3f}'.format)
line_table['fl_signi']=3.*line_table['f_line']/line_table['fl_3sig_lim_Sndr']
line_table['fl_signi']=line_table['fl_signi'].map('{:,.3f}'.format)
#----------------Saving dataframe to ASCII file--------------------------------------------
head = 'This file contains the measurements of lines in the MagE sample. Generated by EW_fitter.py.\n\
Columns are:\n\
label:        Short, unique label of the object\n\
line_lab:     Label of the line\n\
obs_wav:      Observed-frame wavelength of the line (A)\n\
rest_wave:    Rest-frame wavelength of the line (A)\n\
type:         Type of line, emission, photospheric absorption, ISM, etc.\n\
EWr_fit:      Equivalent width in the rest-frame [RIGHT AYAN?], as calculated from Gaussian fit (A)\n\
EWr__fit_u:   Error in above qty. (A)\n\
EWr_sum:      Equivalent width in the rest-frame [RIGHT AYAN?], as calculated by directly summing the flux (A)\n\
EWr_sum_u:    Error in above qty. (A)\n\
f_line:       Flux i.e. area under Gaussian fit (erg/s/cm^2) \n\
f_line_u:     Error in above qty. (erg/s/cm^2)\n\
med_bin_flux: Median of binned fluxes. Each bin is +/-2 sigma wide. There are 5 bins\n\ 
              on either side of a group of line asked to be fit, beyond the +/-5 sigma\n\ 
              window. This is to take into account the wiggles in the spectrum\n\
mad_bin_flux: Median absolute deviation of the above binned fluxes\n\
MAD_significance: (f_line - med_bin_flux)/mad_bin_flux. Probably should NOT be USED anymore\n\
EWr_3sig_lim_Schneider: 3sigma upper limit for unresolved OR detection criteria for resolved EWs,\n\ 
                     as determined using Schneider et al. 1993 prescription\n\
fl_3sig_lim_Sndr: 3sigma upper limit for unresolved OR detection criteria for resolved FLUXES following above prescription\n\
EW_signi:     Ratio of EWr_fit to (EWr_3sig_lim_Schneider/3). Probably use this for SIGNIFICANCE\n\
fl_signi:     Ratio of f_line to (fl_3sig_lim_Sndr/3).\n\
fit_cont:     Continuum, as from the fit (continuum normalised fit)\n\
fit_fl:       Amplitude, as from the fit (continuum normalised fit)\n\
fit_cen:      Center, as from the fit (continuum normalised fit)\n\
fit_cen_u:    Error in above qty. (A)\n\
fit_sig:      Sigma, as from the fit (A)\n\
zz:           Corrected redshift of this line, from the fitted center\n\
zz_u:         Error in above qty.\n\
EW_significance: How many sigma is the fitted EWr_fit as compared to the one calculated by Schneider et al. 1993 prescription\n\
NaN means the code was unable to fit the line.\n\
'
np.savetxt(fout, [], header=head, comments='#')
line_table.to_csv(fout, sep='\t',mode ='a', index=None)
print 'Full table saved to', fout
#----------Displaying part of dataframe if asked to---------------------------
line_table['f_SNR']=np.abs(line_table['f_line'])/line_table['f_line_u']
short_table = line_table[['line_lab','EWr_fit','EWr_fit_u','EWr_3sig_lim_Schneider','EW_signi','f_line','f_line_u','f_SNR','fl_3sig_lim_Sndr','fl_signi']]
print 'Some columns of the table are:'
print short_table
#----------------Sanity check: comparing 2 differently computed EWs------------------
if args.check:
    err_sum, es, n = 0., 0., 0
    for p in range(0,len(line_table)):
        EWr_fit = float(line_table.iloc[p].EWr_fit)
        EWr_sum = float(line_table.iloc[p].EWr_sum)
        EWr_sum_u = float(line_table.iloc[p].EWr_sum_u)
        #print line_table.iloc[p].line_lab, EWr_fit, EWr_sum, EWr_sum_u #Debugging
        if np.abs(EWr_sum) > EWr_sum_u and EWr_fit > 0.:
            err_sum += EWr_fit - EWr_sum
            es += (EWr_fit - EWr_sum)**2.
            n += 1
    if n > 0:
        print err_sum/n, np.sqrt(es/n)
    else:
        print 'No lines detected.'
#------------------------------------------End of main function------------------------------------------------
