maindir = '/Users/jrrigby1/SCIENCE/Lensed-LBGs/Planck_Arc/IR_grisms/WFC3_fit_1Dspec/'
results_dir = maindir + '1dspec/'
proccess_file = maindir + "/Sunburst_G141_grism2process.txt"  # lists whether to use method 1 or 2

# knot4 only has fit1.  For others, fit2 has the right results.  See Sunburst_G141_grism2process.txt

''' Pseudo code
- loop through the process file
- grab all the meth2.fitdf files
- loop through them, grab all the OIII/OII ratios.
- compute uncertainties on the ratio
- make a summary outfile.
- plot it
- repeat for NeIII/OII
'''

# SHould this just be in the main line_ratio NB? **
