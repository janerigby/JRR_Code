import jrr
import matplotlib.pyplot as plt
import pysynphot as S

mage_mode = "released"
label = "Cosmic~Eye"
(sp, resoln, dresoln, LL, zz_syst) = jrr.mage.wrap_open_spectrum(label, mage_mode, addS99=True, zchoice='stars', MWdr=True) # Load the spectrum
sp_orig = sp.copy(deep=True)

jrr.mage.deredden_internal_extinction(sp, 0.2, colf="rest_fnu_cont", deredden_uncert=False)  # Try dereddening it
jrr.spec.calz_unred_df(sp_orig, 0.2, colf='rest_fnu_cont')  # Try with IDL routine ported over
# Good, these 2 approaches agree.  I like mine better b/c it zeros out numbers below 912 A.

ax = sp.plot(x='rest_wave', y='rest_fnu_cont', label='orig')
sp.plot(x='rest_wave', y='rest_fnu_cont_dered', label='dered_extinction', ax=ax)
sp_orig.plot(x='rest_wave', y='rest_fnu_cont_dered', label='dered_IDL', ax=ax)

plt.xlim(900,1400)
plt.ylim(-1E-28, 1E-27)
plt.show()

