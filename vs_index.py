#!/usr/bin/env python

# Create an index file (.inx) of a given .sdf ligand library
# file, by submitting that .sdf as argument. The .inx file
# is created in the same directory as .sdf, with the same
# name but with the .inx extension. It can then be moved to
# the VS directory.
# An .icm script is created temporarly, modified and exectuted, then deleted.
#
# https://github.com/thomas-coudrat/toolbx_vs
# Thomas Coudrat <thomas.coudrat@gmail.com>

import argparse
import os
import sys
from subprocess import check_output, STDOUT, CalledProcessError
import socket
import json


def main():
    """
    Run script
    """

    # Get the path from the Json file
    icm = getPath()

    # Extract the arguments
    sdfFile, suffix = parseArgs()

    # Generate the ICM script
    scriptPath = generateScript(icm, sdfFile, suffix)

    # Execute the script
    executeScript(icm, scriptPath, sdfFile, suffix)


def getPath():
    """
    Get the ICM executable from the ICMHOME environment variable
    """

    # Get environment variable. It returns None if not defined on the system.
    icmHome = os.environ.get('ICMHOME')

    # Return path to executable if the environment variable was found
    if icmHome is None:
        "The ICMHOME environment variable must be set for your system. Exiting."
        sys.exit()
    else:
        icm = icmHome + "/icm64"

    return icm


def parseArgs():
    """
    Argument parsing stuff
    """

    descr = "Create a .inx file of a .sdf file for VS"
    descr_sdf = "Provide a .sdf library to create a .inx file for"
    descr_suffix = "Provide a suffix that will identify the cluster for " \
                   "this index is designed"

    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("sdf", help=descr_sdf)
    parser.add_argument("suffix", help=descr_suffix)

    try:
        args = parser.parse_args()
    except SystemExit:
        print("\n***************************************")
        parser.print_help()
        print("***************************************\n")
        sys.exit()

    sdfFile = args.sdf
    suffix = args.suffix

    return sdfFile, suffix


def generateScript(icm, sdfFile, suffix):
    """
    Generate the ICM script in the current working directory
    Modify its content to apply to the .sdf and .inx files
    """

    # Script base
    scr_string = """#!ICM_EXEC
call "_startup"

# Create the .inx file, that indexes this database
makeIndexChemDb "SDF_LIB" "INX_FILE" "mol" { "ID" }

quit
"""

    # Modify the script
    workDir = os.getcwd()
    sdfPath = workDir + "/" + sdfFile
    inxPath = sdfPath.replace(".sdf", "_" + suffix + ".inx")
    scr_string = scr_string.replace("SDF_LIB", sdfPath)
    scr_string = scr_string.replace("INX_FILE", inxPath)
    scr_string = scr_string.replace("ICM_EXEC", icm)

    # print scr_string

    # Write it to file
    scr_path = workDir + "/temp.icm"
    with open(scr_path, "w") as scr_file:
        scr_file.write(scr_string)

    return scr_path


def executeScript(icm, scriptPath, sdfFile, suffix):
    """
    Execute the index.icm script
    """

    print("\nCreating .inx for: \t" + sdfFile)
    print("Will work on: \t" + suffix + "\n")

    # Execute
    try:
        check_output(icm + " -s " + scriptPath, stderr=STDOUT, shell=True)
    except CalledProcessError as e:
        print("\n Error executing the ICM script")
        print(e.output)
        sys.exit()
    # Delete
    os.remove("./temp.icm")


if __name__ == "__main__":
    main()
