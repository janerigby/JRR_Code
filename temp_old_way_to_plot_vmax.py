### This works, but is kludgy as hell.  Rewriting for pretty.
## Plot the results versus ionization potential
matplotlib.rcParams.update({'font.size': 16})
s=60
vmage_whtavg_notlim = vmage_whtavg[vmage_whtavg['vlowlim'].isnull()]
vmage_median_notlim = vmage_median[vmage_median['vlowlim'].isnull()]
ax1 = vmage_whtavg_notlim.plot(x='IP', y='vmax', kind='scatter', label=r'MagE $v_{max}$', color='red', s=s)
vmage_median_notlim.plot(x='IP', y='vmax', kind='scatter', label=r'MagE $v_{max}$', color='orange', s=s, ax=ax1)
ax1.errorbar(vmage_whtavg_notlim['IP'], vmage_whtavg_notlim['vmax'], yerr=vmage_whtavg_notlim['vmax_std'], ls='none', color='k', lw=1.5, label=None)
ax1.errorbar(vmage_median_notlim['IP'], vmage_median_notlim['vmax'], yerr=vmage_median_notlim['vmax_std'], ls='none', color='k', lw=1.5, label=None)
vmage_whtavg_notlim.plot(x='IP', y='vmean', kind='scatter', label=r'MagE $v_{mean}$', color='blue', s=s, ax=ax1)
vmage_median_notlim.plot(x='IP', y='vmean', kind='scatter', label=r'MagE $v_{mean}$', color='purple', s=s, ax=ax1)
ax1.errorbar(vmage_whtavg_notlim['IP'], vmage_whtavg_notlim['vmean'], yerr=vmage_whtavg_notlim['vmean_std'], ls='none', color='k', lw=1.5, label=None)
ax1.errorbar(vmage_median_notlim['IP'], vmage_median_notlim['vmean'], yerr=vmage_median_notlim['vmean_std'], ls='none', color='k', lw=1.5, label=None)
plt.quiver(vmage_whtavg['IP'], vmage_whtavg['vlowlim'], np.zeros(shape=vmage_whtavg.shape[0]), np.ones(shape=vmage_whtavg.shape[0])*100, color='red')
plt.quiver(vmage_median['IP'], vmage_median['vlowlim'], np.zeros(shape=vmage_median.shape[0]), np.ones(shape=vmage_median.shape[0])*100, color='orange')
plt.xlabel("IP (eV)") ; plt.ylabel("v (km/s)")
plt.xlim(10,70) ; plt.ylim(50,-2900)
ax1.legend(loc='upper left', labelspacing=0.2, borderpad=0.1)
plt.tight_layout()
plt.savefig('mage_vmaxvmean_vsIP.pdf')




indf = (vcos_hiR_whtavg, vcos_loR_whtavg, vcos_hiR_median, vcos_loR_median)   # Make vs ionization potential plot for both version of COS stack
outpdf = ('cos_R2E4_vmaxvmean_whtavg_vsIP.pdf', 'cos_R3500_vmaxvmean_whtavg_vsIP.pdf', 'cos_R2E4_vmaxvmean_median_vsIP.pdf', 'cos_R3500_vmaxvmean_median_vsIP.pdf')
for ii, thedf in enumerate(indf) :
    notlim = thedf[thedf['vlowlim'].isnull()]
    ax2 = notlim.plot(x='IP', y='vmax', kind='scatter', label=r'COS $v_{max}$', color='red', s=s)
    ax2.errorbar(notlim['IP'], notlim['vmax'], yerr=notlim['vmax_std'], ls='none', color='k', lw=1.5, label=None)
    notlim.plot(x='IP', y='vmean', kind='scatter', label=r'COS $v_{mean}$', color='blue', s=s, ax=ax2)
    ax2.errorbar(notlim['IP'], notlim['vmean'], yerr=notlim['vmean_std'], ls='none', color='k', lw=1.5, label=None)
    plt.quiver(thedf['IP'], thedf['vlowlim'], np.zeros(shape=thedf.shape[0]), np.ones(shape=thedf.shape[0])*100, color='r')
    plt.xlim(10,70) ; plt.ylim(50,-2900)
    plt.xlabel("IP (eV)") ; plt.ylabel("v (km/s)")
    ax2.legend(loc='upper left', labelspacing=0.2, borderpad=0.1)
    plt.tight_layout()
    plt.savefig(outpdf[ii])
