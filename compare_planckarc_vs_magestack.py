import jrr
import numpy as np
import pandas
from matplotlib import pyplot as plt

# Load the MagE stack
mage_mode = "reduction"
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
(sp1, LL_magestack)   = jrr.mage.open_stacked_spectrum(mage_mode, colfnu='X_avg', colfnuu='X_jack_std')
sp1['rest_fnu_norm'] = sp1['rest_fnu']
sp1['rest_fnu_unorm'] = sp1['rest_fnu_u']


zz = 2.36988
(pf, LL_planck) = jrr.mage.open_planckarc_sum(zz, vmask1=300., vmask2=3000., smooth_length=50.)
pf['rest_fnu_norm'] = pf['rest_fnu'] / pf['rest_fnu_autocont']
pf['rest_fnu_unorm'] = pf['rest_fnu_u'] / pf['rest_fnu_autocont']

ax = pf.plot(x='rest_wave', y='rest_fnu')
pf.plot(x='rest_wave', y='rest_fnu_autocont', ax=ax)

# May 2017, compare our COS stack to the MagE stack.
the_dfs = [sp1, pf]
the_zzs = [0.0, 0.0]
colortab = ('black', 'blue')
jrr.plot.echelle_spectrum(the_dfs, the_zzs, LL_magestack, Npages=3, outfile="magestack_vs_planckarc.pdf", plot_cont=True, norm_by_cont=False, apply_bad=False, colwave='rest_wave', colfnu='rest_fnu_norm', colfnu_u='rest_fnu_unorm', colcont='unity', title="magestack_vs_planck", waverange=(1400,2760), colortab=colortab, plotx2=False)
plt.close("all")
