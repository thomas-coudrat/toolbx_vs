#!/usr/bin/env python

#------------------------------------------------------------
#
#   This creates maps for a VS from a .ob object containing
#   the target protein only. The .inx index file pointing
#   to the library to be used must also be present in the
#   directory.
#   It uses an .icm script to create the ICM docking maps
#
#   Thomas Coudrat, February 2014
#
#------------------------------------------------------------

import os
import shutil
import sys
import argparse
from subprocess import check_output, STDOUT, CalledProcessError


def main():

    # Get the paths to the scripts to be executed
    icm, script = setPaths()

    # Get the arguments
    obPath, inxPath, pocket = parseArgs()

    # Running the .icm script
    runScript(icm, script, obPath, pocket)

    # Set of calls to modify the .dtb file
    modifyDtb("i_maxHdonor", "  15", obPath)
    modifyDtb("i_maxLigSize", "  1000", obPath)
    modifyDtb("i_maxNO", "  20", obPath)
    modifyDtb("i_maxTorsion", "  20", obPath)
    modifyDtb("i_ringFlexLevel", "  1", obPath)
    modifyDtb("l_sampleRacemic", "  yes", obPath)
    modifyDtb("r_ScoreThreshold", "  -25", obPath)
    modifyDtb("r_maxPk", "  15", obPath)
    modifyDtb("r_minPk", "  -10", obPath)
    modifyDtb("s_chargeGroups", "  auto", obPath)
    modifyDtb("s_dbIndex", inxPath, obPath)


def setPaths():
    """
    Devide which platform this is executed on, and select the proper script
    paths
    """

    # Path to the .icm script
    scriptVlsci = "/vlsci/VR0024/tcoudrat/Scripts/vs_scripts/maps.icm"
    scriptDesktop = "/home/thomas/Copy/Tools/vs_scripts/maps.icm"
    icmVlsci = "/vlsci/VR0024/tcoudrat/bin/icm-3.7-3b/icm64"
    icmDesktop = "/usr/icm-3.7-3b/icm64"

    # Select the right path, depending on the platform where this is executed
    if os.path.exists(scriptVlsci):
        script = scriptVlsci
        icm = icmVlsci
    elif os.path.exists(scriptDesktop):
        script = scriptDesktop
        icm = icmDesktop
    else:
        print "maps.icm script or icm executable not found"
        sys.exit()

    return icm, script


def parseArgs():
    """
    Get the arguments, define the script's help
    """

    descr = "Creates maps for VS with ICM, requires target in .ob and " \
            "library index in .inx format"
    descr_obPath = "Provide the path to the .ob file containing the target " \
                   "for the VS"
    descr_inxPath = "Provide the path to the .inx index to the ligand " \
                    "to be used for this VS"
    descr_pocket = "Define the pocket number to use (determined by " \
                   "ICMpocket finder). Default is pocket #1."
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("obPath", help=descr_obPath)
    parser.add_argument("inxPath", help=descr_inxPath)
    parser.add_argument("--pocket", help=descr_pocket)

    args = parser.parse_args()

    # Deal with presence of absence of the argument pocket
    if args.pocket:
        pocket = args.pocket
    else:
        pocket = "1"

    return args.obPath, args.inxPath, pocket


def runScript(icm, script, obPath, pocket):
    """
    Copy a temp copy of the .icm script, modify it, and run it
    """
    projName = obPath.replace(".ob", "")

    # Copy
    shutil.copy(script, "./temp.icm")
    # Modify
    os.system("sed -e 's|VS_PROJECT|" + projName + "|g' ./temp.icm -i")
    os.system("sed -e 's|POCKET_NUM|" + pocket + "|g' ./temp.icm -i")
    # Execute
    try:
        check_output(icm + " -s ./temp.icm", stderr=STDOUT, shell=True)
    except CalledProcessError, e:
        print e.output
        sys.exit()
    # Delete temp script
    os.remove("./temp.icm")


def modifyDtb(keyword, value, obPath):
    """
    Make modifications to the .dtb file to fine tune parameters
    """

    dtbPath = obPath.replace(".ob", ".dtb")
    dtbFile = open(dtbPath, "r")
    dtbLines = dtbFile.readlines()
    dtbFile.close()

    totalLines = len(dtbLines)
    newDtbLines = []
    i = 0

    #print dtbLines[i]

    while i < totalLines:
        if keyword in dtbLines[i]:
            newDtbLines.append(dtbLines[i])
            newDtbLines.append(value + "\n")
            i += 2
            continue
        else:
            newDtbLines.append(dtbLines[i])
            i += 1

    newDtbFile = open(dtbPath, "w")
    for line in newDtbLines:
        newDtbFile.write(line)


if __name__ == "__main__":
    main()
