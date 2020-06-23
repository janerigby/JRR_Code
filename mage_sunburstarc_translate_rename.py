import jrr
import shutil

# Use case:  given a new name (for example, sunburst_M-7), go find the right file
# with the old naming convention.

# grab the dictionary of old and new names for the Sunburst Arc (Planck Arc)
sunburst_dict =  jrr.mage.sunburst_translate_names()
example = 'planckarc_h1' 
print(sunburst_dict[example])

speclist = jrr.mage.wrap_getlist('reduction', which_list='all')

for key, value in sunburst_dict.items() :
    print("#        Copying from" , key, "to ", value)
    oldfilename = speclist.loc[key]['filename']
    newfilename = "Planckarc/M-Xnames/" + value + "-comb1_MWdr.txt"
    print(oldfilename, newfilename)
    new_header_text = "# New filename: " + newfilename + " following Sunburst Arc M-X naming convention of MagE slits for publication.\n"
    jrr.util.put_header_on_file(oldfilename, new_header_text, newfilename)  # Copy file to new maning convention, adding 1 line to header
