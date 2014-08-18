#!/usr/bin/env python

# -------------------------------------------------------------
#
#   Create an index file (.inx) of a given .sdf ligand library
#   file, by submitting that .sdf as argument. The .inx file
#   is created in the same directory as .sdf, with the same
#   name but with the .inx extension. It can then be moved to
#   the VS directory.
#   An .icm script is created temporarly, modified and exectuted, then deleted.
#
#   Thomas Coudrat, August 2014
#
# -------------------------------------------------------------

import argparse
import os
import sys
from subprocess import check_output, STDOUT, CalledProcessError
import socket
import json


def main():
    """
    Run the script
    """

    # Get the path from the Json file
    icm = getPath()

    # Extract the arguments
    sdfFile = parseArgs()

    # Generate the ICM script
    scriptPath = generateScript(icm, sdfFile)

    # Execute the script
    executeScript(icm, scriptPath, sdfFile)


def getPath():
    """
    Read in the json file containing the ICM executable paths for each machine
    Return the ICM executable corresponding to the hostname of the machine being
    used currently
    """

    # This Json file stores the ICM executable locations for each platform
    icmExecJson = os.path.dirname(os.path.realpath(__file__)) + "/icm_exec.json"

    # Read content of .json file
    with open(icmExecJson, "r") as jsonFile:
        icmExec = json.load(jsonFile)

    # Get the hostname to know which computer this is executed on
    hostname = socket.gethostname()

    # Assign the ICM executable path corresponding to the hostname, if it is not
    # defined then stop the execution
    if hostname in icmExec.keys():
        icm = icmExec[hostname]
    else:
        print("The ICM executable is not defined for this machine, please edit\
            the icm_exec.json file")
        sys.exit()

    return icm


def parseArgs():
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


def generateScript(icm, sdfFile):
    """
    Generate the ICM script in the current working directory
    Modify its content to apply to the .sdf and .inx files
    """

    # Write the script
    script_string = """#!ICM_EXEC
call "_startup"

# Create the .inx file, that indexes this database
makeIndexChemDb "SDF_LIB" "INX_FILE" "mol" { "ID" }

quit
"""
    workDir = os.getcwd()
    scriptPath = workDir + "/temp.icm"
    with open(scriptPath, "w") as script_file:
        script_file.write(script_string)

    # Modify the script
    sdfPath = workDir + "/" + sdfFile
    inxPath = sdfPath.replace(".sdf", ".inx")
    os.system("sed -e 's|SDF_LIB|" + sdfPath + "|g' ./temp.icm -i")
    os.system("sed -e 's|INX_FILE|" + inxPath + "|g' ./temp.icm -i")
    os.system("sed -e 's|ICM_EXEC|" + icm + "|g' ./temp.icm -i")

    return scriptPath


def executeScript(icm, scriptPath, sdfFile):
    """
    Execute the index.icm script
    """

    print
    print "CREATING .INX FOR:", sdfFile
    print
    # Execute
    try:
        check_output(icm + " -s " + scriptPath, stderr=STDOUT, shell=True)
    except CalledProcessError, e:
        print "\n Error executing the ICM script"
        print e.output
        sys.exit()
    # Delete
    os.remove("./temp.icm")


if __name__ == "__main__":
    main()
