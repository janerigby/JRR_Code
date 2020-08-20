# Evaluate the polynomial metallicity fits in Byler et al. 2020, given UV line ratios.

import numpy as np
import pandas
import matplotlib.pyplot as plt
import  astropy.modeling
filename_tab4 = 'apjab7ea9t4_ascii.txt'
filename_tab5 = 'apjab7ea9t5_ascii.txt'

def evaluate_Byler_poly(df, diagnostic, degree, x, y) :
    coeffs = df[str(degree)].to_dict()[diagnostic]
    foo = astropy.modeling.models.Polynomial2D(degree, **coeffs)  
    #print("   Evaluating using x=", x, "y=", y, "yields 12+log(O/H)=", foo(x, y))
    #print("   Evaluating yields 12+log(O/H)=", foo(x, y))
    return(foo(x,y))

        
# Footnote 13 in Nell Byler's 2020 paper:
# Quick-start for Python users: Create a dictionary of coefficient names and values from either Tables 4 or 5, “coeffs.”
# Input this dictionary to Polynomial2D(degree, ∗∗coeffs) from astropy.modeling. models.

df = {}
df['3'] = pandas.read_csv(filename_tab4, comment="#", delim_whitespace=True)
df['4'] = pandas.read_csv(filename_tab5, comment="#", delim_whitespace=True)
degrees = (3, 4)
x = -0.1 ;  y = 0.0  # Need to populate this w values from S1723

# dereddened fluxes and uncertainties from table 1 of S1723 integrated paper (Rigby et al. 2020 about to resubmit)
o3_1666   = 6.1      ; pm_1666 = 1.5  
si3_1883  = 8.6      ; pm_1883 = 1.0
he2_1640  = 7.8      ; pm_1640 = 1.8
c3sum = 21.3 + 14.85 ; pm_c3   = 1.4

# Error propogation through dumb repetition. # Wow, I'm so hip
size = 10000
o3_ar    = np.random.normal(loc=o3_1666,  scale=pm_1666, size=size)
si3_ar   = np.random.normal(loc=si3_1883, scale=pm_1883, size=size)
he2_ar   = np.random.normal(loc=he2_1640, scale=pm_1640, size=size)
c3sum_ar = np.random.normal(loc=c3sum,    scale=pm_c3,   size=size)

# For Si3-O3C3, x=log10(OIII 1666  / CIII 1906,1909)   From S5.4 of Nell's paper.
# For He2-O3C3, x=log10(O III 1666 / C III 1906, 1909)
x = {'Si3-O3C3' : np.log10(o3_1666 / c3sum),    'He2-O3C3' : np.log10(o3_1666 / c3sum) }
# For Si3-O3C3  y=log10(SiIII 1883 / CIII 1906, 1909)
# For He2-O3C3, y=log10(HeII 1640 /  C III 1906, 1909)
y = {'Si3-O3C3' : np.log10(si3_1883 / c3sum),  'He2-O3C3' : np.log10(he2_1640 / c3sum) }

x_ar = {'Si3-O3C3' : np.log10(o3_ar  / c3sum_ar),   'He2-O3C3' : np.log10(o3_ar  / c3sum_ar) }
y_ar = {'Si3-O3C3' : np.log10(si3_ar / c3sum_ar),   'He2-O3C3' : np.log10(he2_ar / c3sum_ar) }

Z = {} ; Z_median = {} ; Z_std = {}
diagnostics = ('Si3-O3C3', 'He2-O3C3')
for diagnostic in diagnostics :
    for degree in degrees:
        print("Diagnostic is", diagnostic, ", degree is", degree)
        label = diagnostic + "_deg" + str(degree)
        Z[label] = evaluate_Byler_poly(df, diagnostic, degree, x_ar[diagnostic], y_ar[diagnostic])

# Plot the results
for diagnostic in diagnostics:
    plt.figure()
    for degree in degrees:
        thisone = diagnostic + '_deg' + str(degree)
        plt.hist(Z[thisone], label='degree=' + str(degree), bins=100)
        Z_median[thisone] = np.median(Z[thisone])
        Z_std[thisone]    = np.std(Z[thisone])
        print(thisone, Z_median[thisone], Z_std[thisone])
        plt.axvline(Z_median[thisone], color='black')
    plt.title(diagnostic)
    plt.xlabel("12 + log(O/H)")
    plt.xlim(7.5, 8.5)
    plt.legend()

plt.show()
