from shutil import copyfile
import jrr
# Update the old "Planck arc" mage names in a table to the new M-X format.

thedir = "/Users/jrrigby1/Dropbox/MagE_atlas/Contrib/Intervening/Latex_asof_Jun232020/"
oldfile = "tab_doublet_table_JRRedit_v2.tex"
newfile = "tab_doublet_table_JRRedit_v3.tex"
translate_dict = jrr.mage.sunburst_translate_names(longnames_old=False, longnames_new=False)

jrr.util.replace_text_in_file('h5 &',     'M-2 (h5or6) &',  thedir+oldfile, thedir+newfile)  # JRR typo a long time ago
for key, value in translate_dict.items():
    print(key, value)
    jrr.util.replace_text_in_file(key + ' &', value + ' ('+key+') &', thedir+newfile)


