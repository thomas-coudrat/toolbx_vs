#!/vlsci/VR0024/tcoudrat/bin/python275/bin/python

#-------------------------------------------------------------
#
#   Create an index file (.inx) of a given .sdf ligand library
#   file, by submitting that .sdf as argument. The .inx file
#   is created in the same directory as .sdf, with the same
#   name but with the .inx extension. It can then be moved to
#   the VS directory.
#   A blank .icm file is copied temporarly, modified temp
#   and exectuted, then deleted.
#
#   Thomas Coudrat, February 2014
#
#-------------------------------------------------------------

import argparse
import os
import shutil
import sys
from subprocess import check_output, STDOUT, CalledProcessError

# Path to the blank .icm file
indexFile = "/vlsci/VR0024/tcoudrat/Scripts/VS/index.icm"
icmExec = "/vlsci/VR0024/tcoudrat/bin/icm-3.7-3b/icm64"

# Argument parsing stuff
descr = "Create a .inx file of a .sdf file for VS"
descr_sdf = "Provide a .sdf library to create a .inx file for"
parser = argparse.ArgumentParser(description=descr)
parser.add_argument("sdf", help=descr_sdf)
args = parser.parse_args()
sdfFile = args.sdf

print
print "CREATING .INX FOR:", sdfFile
print

# Copy the temp .icm script, modify it, run it, delete it
# Copy
shutil.copy(indexFile, "./temp.icm")
inxFile = sdfFile.replace(".sdf", ".inx")
# Modify
os.system("sed -e 's|SDF_LIB|" + sdfFile + "|g' ./temp.icm -i")
os.system("sed -e 's|INX_FILE|" + inxFile + "|g' ./temp.icm -i")
# Execute
try:
    check_output(icmExec + " -s ./temp.icm", stderr=STDOUT, shell=True)
except CalledProcessError, e:
    print "\n Error executin the ICM script"
    print e.output
    sys.exit()
# Delete
os.remove("./temp.icm")
