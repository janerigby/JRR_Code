import jrr
import shutil
import glob
from os.path import basename
import re
from shutil import copyfile

# Copy over the old Planck Arc (old naming convention) files to the new Sunburst_M-X naming convention.
new_dir = "Planckarc/M-Xnames/"

# grab the dictionary of old and new names for the Sunburst Arc (Planck Arc)
(sunburst_dict, reversed_dict) =  jrr.mage.sunburst_translate_names()
# example = 'planckarc_h1' 
# print(sunburst_dict[example])

# Doing this in chunks
move_spectra_files = False
move_S99_bpass_fits = True


if move_spectra_files :
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

# Rename John's starburst99 and Bpass fits to new sunburst arc naming convention
if move_S99_bpass_fits :
    main_dir = "/Users/jrrigby1/Dropbox/MagE_atlas/Contrib/S99/Sunburst-Arc/"
    new_dir = "S99_Bpass_newnames/"
    old_dir = "S99_Bpass_oldnames/"
    myfiles = [ basename(x) for x in glob.glob(main_dir + old_dir + "*.*") ]
    for thisfile in myfiles :
        splitit = re.split("-|_", thisfile)  # Split on hyphen or underscore
        if  splitit[0] == "planck" : splitit[0] = "planckarc"
        if (splitit[0] + '_' + splitit[1]) == 'planckarc_fire' : new_prefix = "sunburst_fire"
        else :            
            new_prefix = sunburst_dict[ splitit[0] + "_" + splitit[1] ]
            print("Splitting ", thisfile, "into ", splitit[0], splitit[1], "    --> ", new_prefix)
        new_name = new_prefix + '-' + '-'.join(splitit[2:])
        print("renaming to", new_name)
        # Next line corrupts the fits files, so get rid of it!
        #new_header_text = "# Renaming on 6/25/2020 by jrigby, from " + thisfile + " to the new naming convention: " + new_name  
        #jrr.util.put_header_on_file(main_dir + old_dir + thisfile,   new_header_text,   main_dir + new_dir + new_name)
        copyfile(main_dir + old_dir + thisfile, main_dir + new_dir + new_name)
