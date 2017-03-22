import pandas
import jrr

file_ESI  = "../../../MagE_atlas/Contrib/ESI_spectra/2016Aug/s1723_ESI_JRR_sum.txt"
spec_ESI = pandas.read_table(file_ESI, delim_whitespace=True, comment="#")
spec_ESI.rename(columns={'flam_sum_jrr' : 'flam'}, inplace=True)
