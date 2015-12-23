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
    sdfFile, suffix = parseArgs()

    # Generate the ICM script
    scriptPath = generateScript(icm, sdfFile, suffix)

    # Execute the script
    executeScript(icm, scriptPath, sdfFile, suffix)


def getPath():
    """
    Read in the json file containing the ICM executable paths for each machine
    Return the ICM executable corresponding to the hostname of the machine being
    used currently
    """

    # This Json file stores the ICM executable locations for each platform
    icmExecJson = os.path.dirname(os.path.realpath(__file__)) + "/../icm_exec.json"

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
