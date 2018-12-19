import glob
import subprocess
from os.path import expanduser, basename, normpath, exists
from os import makedirs

def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line


atlasdir = expanduser('~') + '/Dropbox/MagE_2D/'  # output dir
reduxdir = '/Volumes/Apps_and_Docs/SCIENCE/Lensed-LBGs/Mage/Redux/'  # reduction dir
infile  = '2D-rectify.dirs'  # list of all the 2D dirs to copy over

thisdir = 'foo/'  # initalize

with open(reduxdir + infile) as f_in:
    for line in nonblank_lines(f_in):  # strip blank lines
        if '#' in line : next          # strip comments
        elif 'NAME' in line :
            (first, subdir) = line.split()
            print "Working on object   ", subdir
            outdir = atlasdir + subdir + '/'
            if not exists(outdir):  makedirs(outdir)
        else :
            print line
            subsubdir = basename(normpath(line))  # grab the last dir in the path
            goeshere = outdir  + subsubdir + '/'
            if not exists( goeshere ) : makedirs( goeshere)
            dothis = '/bin/cp -pr ' + normpath(line) + '/*sum*fits  ' + goeshere
            print "Planning to ", dothis
            temp = subprocess.check_output(dothis, shell=True).split()


dothis = '/bin/cp -pr  ' + reduxdir + infile + ' ' + atlasdir
temp = subprocess.check_output(dothis, shell=True).split()
