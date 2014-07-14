#!/usr/bin/env python

# -------------------------------------------------------------
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
# -------------------------------------------------------------

import argparse
import os
import shutil
import sys
from subprocess import check_output, STDOUT, CalledProcessError
import socket


def main():
    """
    Run the script
    """

    # Get the right paths
    icm, script = setPaths()

    # Extract the arguments
    sdfFile = getArgs()

    # Execute the script
    executeScript(icm, script, sdfFile)


def setPaths():
    """
    Determine which platform this is ran on based on the access to the scripts
    """
    # Path to the .icm script and icm executable
    scriptVlsci = "/vlsci/VR0024/tcoudrat/Scripts/toolbx_vs/index.icm"
    icmVlsci = "/vlsci/VR0024/tcoudrat/bin/icm-3.7-3b/icm64"

    scriptLocal = "/home/thomas/Copy/toolbx_vs/index.icm"
    icmDesktop = "/usr/icm-3.7-3b/icm64"
    icmLaptop = "/home/thomas/bin/icm-3.8-0/icm64"

    scriptMcc = "/nfs/home/hpcpharm/tcoudrat/Scripts/toolbx_vs/index.icm"
    icmMcc = "/nfs/home/hpcpharm/tcoudrat/bin/icm-3.7-3b/icm64"

    # Get the hostname to know which computer this is executed on
    hostname = socket.gethostname()

    # Select the proper scripts
    if hostname == "barcoo":
        script = scriptVlsci
        icm = icmVlsci
    elif hostname == "linux-T1650":
        script = scriptLocal
        icm = icmDesktop
    elif hostname == "Ideapad":
        script = scriptLocal
        icm = icmLaptop
    elif hostname == "msgln6.its.monash.edu.au":
        script = scriptMcc
        icm = icmMcc
    else:
        print "Script or icm executable not found"
        sys.exit()

    return icm, script


def getArgs():
    """
    Argument parsing stuff
    """

    descr = "Create a .inx file of a .sdf file for VS"
    descr_sdf = "Provide a .sdf library to create a .inx file for"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("sdf", help=descr_sdf)
    args = parser.parse_args()
    sdfFile = args.sdf

    return sdfFile


def executeScript(icm, script, sdfFile):
    """
    Execute the index.icm script
    """

    print
    print "CREATING .INX FOR:", sdfFile
    print

    # Get working directory
    workDir = os.getcwd()
    sdfPath = workDir + "/" + sdfFile
    inxPath = sdfPath.replace(".sdf", ".inx")

    # Copy the temp .icm script, modify it, run it, delete it
    # Copy
    shutil.copy(script, "./temp.icm")
    # Modify
    os.system("sed -e 's|SDF_LIB|" + sdfPath + "|g' ./temp.icm -i")
    os.system("sed -e 's|INX_FILE|" + inxPath + "|g' ./temp.icm -i")
    # Execute
    try:
        check_output(icm + " -s ./temp.icm", stderr=STDOUT, shell=True)
    except CalledProcessError, e:
        print "\n Error executin the ICM script"
        print e.output
        sys.exit()
    # Delete
    os.remove("./temp.icm")


if __name__ == "__main__":
    main()
