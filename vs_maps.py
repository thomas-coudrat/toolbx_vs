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
import socket
from subprocess import check_output, STDOUT, CalledProcessError


def main():

    # Get the paths to the scripts to be executed
    icm, pocketScript, ligScript = setPaths()

    # Get the arguments
    obPath, inxPath, mapMode, pocket = parseArgs()

    # Running the .icm script creating the maps based on ICMpocket finder
    if mapMode == "pocket":
        runPocketScript(icm, pocketScript, obPath, pocket)

    # Running the .icm script creating the maps based on the bound ligand
    if mapMode == "ligand":
        runLigScript(icm, ligScript, obPath)

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

    # Paths to the .icm scripts
    pocketVlsci = "/vlsci/VR0024/tcoudrat/Scripts/vs_scripts/mapsPocket.icm"
    ligandVlsci = "/vlsci/VR0024/tcoudrat/Scripts/vs_scripts/mapsLigand.icm"

    pocketDesktop = "/home/thomas/Copy/Tools/vs_scripts/mapsPocket.icm"
    ligandDesktop = "/home/thomas/Copy/Tools/vs_scripts/mapsLigand.icm"

    pocketMcc = "/nfs/home/hpcpharm/tcoudrat/Scripts/vs_scripts/mapsPocket.icm"
    ligandMcc = "/nfs/home/hpcpharm/tcoudrat/Scripts/vs_scripts/mapsLigand.icm"

    # Paths to the icm executables
    icmVlsci = "/vlsci/VR0024/tcoudrat/bin/icm-3.7-3b/icm64"
    icmDesktop = "/usr/icm-3.7-3b/icm64"
    icmMcc = "/nfs/home/hpcpharm/tcoudrat/bin/icm-3.7-3b/icm64"

    # Get the hostname
    hostname = socket.gethostname()

    # Select the right path, depending on the platform where this is executed
    if hostname == "barcoo":
        pocketScript = pocketVlsci
        ligScript = ligandVlsci
        icm = icmVlsci
    elif hostname == "linux-T1650":
        pocketScript = pocketDesktop
        ligScript = ligandDesktop
        icm = icmDesktop
    elif hostname == "msgln6.its.monash.edu.au":
        pocketScript = pocketMcc
        ligScript = ligandMcc
        icm = icmMcc
    else:
        print "ICM scripts or ICM executable not found"
        sys.exit()

    return icm, pocketScript, ligScript


def parseArgs():
    """
    Get the arguments, define the script's help
    """

    # Definition of the arguments
    descr = "Creates maps for VS with ICM, requires target in .ob and " \
            "library index in .inx format"
    descr_obPath = "Provide the path to the .ob file containing the target " \
                   "for the VS"
    descr_inxPath = "Provide the path to the .inx index to the ligand " \
                    "to be used for this VS"
    descr_pocketMap = "Use this flag if the maps are to be created using" \
                      " the ICMpocketFinder method"
    descr_ligandMap = "Use this flag if the maps are to be created using" \
                      " the residues around the bound ligand (also deletes " \
                      "the ligand)"
    descr_pocket = "Define the pocket number to use (determined by " \
                   "ICMpocket finder). Default is pocket #1."

    # Parse arguments
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("obPath", help=descr_obPath)
    parser.add_argument("inxPath", help=descr_inxPath)
    parser.add_argument("-pocketMap", action="store_true",
                        help=descr_pocketMap)
    parser.add_argument("-ligandMap", action="store_true",
                        help=descr_ligandMap)
    parser.add_argument("--pocket", help=descr_pocket)
    args = parser.parse_args()

    # Store arguments
    obPath = args.obPath
    inxPath = args.inxPath
    pocketMap = args.pocketMap
    ligandMap = args.ligandMap
    pocket = args.pocket

    # Deal with arguments
    if not pocket:
        pocket = "1"

    if not pocketMap and not ligandMap:
        print "A map mode must be selected, either -pocketMap or -ligandMap"
        sys.exit()
    elif pocketMap and ligandMap:
        print "Only one of -pocketMap or -ligandMap options must be selected"
        sys.exit()
    elif pocketMap:
        mapMode = "pocket"
    elif ligandMap:
        mapMode = "ligand"

    return obPath, inxPath, mapMode, pocket


def runPocketScript(icm, pocketScript, obPath, pocket):
    """
    Copy a temp copy of the .icm script, modify it, and run it
    """
    projName = obPath.replace(".ob", "")

    # Copy
    shutil.copy(pocketScript, "./temp.icm")
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


def runLigScript(icm, ligScript, obPath):
    """
    This modifies and runs the script creating the maps around the bound ligand
    """
    projName = obPath.replace(".ob", "")

    # Copy
    shutil.copy(ligScript, "./temp.icm")
    # Modify
    os.system("sed -e 's|VS_PROJECT|" + projName + "|g' ./temp.icm -i")
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
