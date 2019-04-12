from __future__ import print_function
import numpy as np
import mpfit
import jrr
import matplotlib.pyplot as plt
from lmfit.models import GaussianModel, ConstantModel  # test-driving lmfit

# Following example in mpfit.py
def myfunct(p, fjac=None, x=None, y=None, err=None) :
    model =  p[0] * np.exp((x - p[1])**2 / (-2*p[2]**2))  + p[3]
    status = 0
    return([status, (y-model)/err])

plt.close("all")

xmax = 200   # size of x array
aa = 100.    # Amplitude of gaussian
bb = xmax/2. # Center of gaussian
cc = 3.0     # sigma of gaussian
cont = 3.0   # continuum level
randamp = 3. # amplitude of gaussian noise

x   = np.arange(0,xmax)
y   = randamp * np.random.randn(xmax) + jrr.spec.onegaus(x, aa, bb, cc, cont)
err = randamp * np.random.randn(xmax)
plt.plot(x, y, color='black')

# Let's try fitting that gaussian w the python version of MPFIT
p0 = (50., xmax/2, 5, 2)
fa = {'x':x, 'y':y, 'err':err}
m = mpfit.mpfit(myfunct, p0, functkw=fa)
# There has got to be a less kludgy way to do line below
bestfit = jrr.spec.onegaus(x, m.params[0], m.params[1], m.params[2], m.params[3])
plt.plot(x, bestfit, color='blue')


# The same, but with LMFIT
mod1 = GaussianModel()
mod2 = ConstantModel()  # The continuum level
mod = mod1 + mod2
pars = mod1.guess(y, x=x) + mod2.guess(y, x=x)  # This is cool.  It made rough guesses for us.
pars['c'].min = 0  # Set bounds on continuum
pars['c'].max = 10 # Set bounds on continuum
#pars['amplitude'].vary = False  # Fix a parameter
out = mod.fit(y, pars, x=x, weights=1/err**2)  # Fitting is done here.
plt.plot(x, out.best_fit, color='orange')
print(out.fit_report(min_correl=0.25))
plt.show()
# This is actually more elegant.  I think I should learn LMFIT and use it...
