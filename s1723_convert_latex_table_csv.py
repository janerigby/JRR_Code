from astropy.io.ascii import read
import re
import numpy as np
import jrr  # Jane's routines

''' For reasons of back-and-forth collaboration with ANU, the file of record of the line fluxes for S1723 
    is a latex table. Converting it to a .csv text file for publication as an MRT, and for general re-use.
    jrigby, Aug 2020'''


def convert_s1723_lineflux_table_latex2csv(infile="tab_s1723_measured_emission_v6.tex") :
    # Read in the line fluxes of S1723 as a template.  Should move this to jrr.mage, for re-use
    # So much data munging here, to get this Latex table (the file of record) into a machine readable format.
    names = ['lineID', 'wave', 'spectrograph', 'Wr', 'Wr_u', 'significance', 'flux', 'flux_u', 'dr_flux', 'dr_flux_u']  #column names
    df3 = read(infile).to_pandas()  # Latex -->astropy tables --> pandas.  Kinda disreputable
    df3.columns = names                                    # Replace Astropy tables' stupid column names.
    df3['detected'] = ~df3['flux'].str.contains('<').astype(bool)  # Is line detected (True) or a limit (False)
    df3 = df3.applymap(lambda x: re.sub('\$', '', x) if isinstance(x, str) else x)    # kill dollar signs
    df3 = df3.applymap(lambda x: re.sub('>', '', x) if isinstance(x, str)  else x)    # Kill greaterthan symbol
    df3 = df3.applymap(lambda x: re.sub('<', '', x) if isinstance(x, str)  else x)    # Kill lessthan symbol
    df3 = df3.applymap(lambda x: re.sub('\\\\tablenotemark{a}', '', x) if isinstance(x, str) else x)  # remove table notes
    df3.replace('-', np.nan, inplace=True)
    outfile = re.sub('\.tex', '.csv', infile)  
    print("DEBUG!", outfile)
    df3.to_csv(outfile)
    header_text = '#Line fluxes for lensed galaxy SDSS~J1723+3411, from Rigby et al. 2020.  Fluxes are in units of 10^-17 erg/s/cm^2. See latex table for full notes.\n'
    jrr.util.put_header_on_file(outfile, header_text, outfile)
    return(df3)

# Take the original latex table from the paper, and make a CSV as the MRT of record.
not_used = convert_s1723_lineflux_table_latex2csv('/Users/jrrigby1/Dropbox/Grism_S1723/Latex/Resubmit/tab_s1723_measured_emission_v6.tex')

