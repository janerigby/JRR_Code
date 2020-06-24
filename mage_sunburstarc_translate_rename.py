import jrr
import shutil

# Copy over the old Planck Arc (old naming convention) files to the new Sunburst_M-X naming convention.
new_dir = "Planckarc/M-Xnames/"

# grab the dictionary of old and new names for the Sunburst Arc (Planck Arc)
(sunburst_dict, reversed_dict) =  jrr.mage.sunburst_translate_names()
example = 'planckarc_h1' 
print(sunburst_dict[example])

speclist = jrr.mage.wrap_getlist('reduction', which_list='all')
for key, value in sunburst_dict.items() :
    print("#        Copying from" , key, "to ", value)
    oldfilename = speclist.loc[key]['filename']
    newfilename = new_dir + value + "-comb1_MWdr.txt"
    print(oldfilename, newfilename)
    new_header_text = "# New filename: " + newfilename + " following Sunburst Arc M-X naming convention of MagE slits for publication.\n"
    jrr.util.put_header_on_file(oldfilename, new_header_text, newfilename)  # Copy file to new naming convention, adding 1 line to header

oldfilename =   "Planckarc/All/planckarc_all-tel-comb_MWdr.txt"
newfilename =  new_dir + "sunburst_all-tel-comb_MWdr.txt"
new_header_text = "# New filename: " + newfilename + " following Sunburst Arc M-X naming convention of MagE slits for publication.\n"
jrr.util.put_header_on_file(oldfilename, new_header_text, newfilename)  # Copy file to new naming convention, adding 1 line to header

# Now, copy over the linelists.  Used     /Users/jrrigby1/SCIENCE/Lensed-LBGs/Mage/Analysis/Plot-all/Lines/copy_sunburst_linelists.sh
