from __future__ import print_function
import jrr
import glob
from os.path import basename, expanduser, exists

# Set this up.  Kludgy
thisdir = "/Volumes/Apps_and_Docs/jrrigby1/Dropbox/Grism_S1723/WFC3_fit_1Dspec/1Dsum/"

myfiles = [ basename(x) for x in glob.glob(thisdir + "*.fitdf") ]
for infile in myfiles :
    head = jrr.util.read_header_from_file(infile, comment="#")
    print(head)




#jrr.util.put_header_on_file(infile, head, outfile)
