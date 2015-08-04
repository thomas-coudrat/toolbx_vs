#!/usr/bin/env python

# ------------------------------------------------------------
#
#   This creates maps for a VS from a .ob object containing
#   the target protein only. The .inx index file pointing
#   to the library to be used must also be present in the
#   directory.
#   It uses an .icm script to create the ICM docking maps
#
#   Thomas Coudrat, February 2014
#
# ------------------------------------------------------------

import os
import sys
import argparse
import socket
from subprocess import check_output, STDOUT, CalledProcessError
import json


def main():

    # Get the paths to the scripts to be executed
    icm = getPath()

    # Get the arguments
    obPath, inxPath, dbType, mapMode, pocket = parseArgs()

    # Generate script (will return the script relevant to the mode chosen)
    script = generateScript(mapMode, obPath, pocket, icm)

    # Run the .icm script
    runScript(icm, script)

    # Set of calls to modify the .dtb file
    modifyDtb("i_maxHdonor", "  15", obPath)
    modifyDtb("i_maxLigSize", "  1000", obPath)
    modifyDtb("i_maxNO", "  20", obPath)
    modifyDtb("i_maxTorsion", "  20", obPath)
    modifyDtb("i_ringFlexLevel", "  1", obPath)
    modifyDtb("r_ScoreThreshold", "  -20", obPath)
    modifyDtb("r_maxPk", "  15", obPath)
    modifyDtb("r_minPk", "  -10", obPath)
    modifyDtb("s_chargeGroups", "  auto", obPath)
    modifyDtb("s_dbIndex", inxPath, obPath)
    if dbType == "3D":
        modifyDtb("s_dbType", "mol 3D", obPath)
        modifyDtb("l_sampleRacemic", "  no", obPath)
    elif dbType == "2Drac":
        modifyDtb("s_dbType", "mol 2D", obPath)
        modifyDtb("l_sampleRacemic", "  yes", obPath)


def parseArgs():
    """
    Get the arguments, define the script's help
    """

    # Definition of the arguments
    descr = "Creates maps for VS with ICM, requires target in .ob as an ICM " \
            "object and library index in .inx format"
    descr_obPath = "Provide the path to the .ob file containing the target " \
                   "for the VS"
    descr_inxPath = "Provide the path to the .inx index to the ligand " \
                    "to be used for this VS"
    descr_dbType = "Database type: '3D' for 3D chiral compounds or '2Drac' " \
                   " for 2D racemic compounds, to be sampled"
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
    parser.add_argument("dbType", help=descr_dbType)
    parser.add_argument("-ligandMap", action="store_true",
                        help=descr_ligandMap)
    parser.add_argument("--pocket", help=descr_pocket)
    args = parser.parse_args()

    # Store arguments
    obPath = args.obPath
    inxPath = args.inxPath
    dbType = args.dbType
    pocketMap = args.pocketMap
    ligandMap = args.ligandMap
    pocket = args.pocket

    # Deal with arguments
    if not pocket:
        pocket = "1"

    if not pocketMap and not ligandMap:
        print("A map mode must be selected, either -pocketMap or -ligandMap")
        sys.exit()
    elif pocketMap and ligandMap:
        print ("Only one of -pocketMap or -ligandMap options must be selected")
        sys.exit()
    elif pocketMap:
        mapMode = "pocket"
    elif ligandMap:
        mapMode = "ligand"

    if dbType not in ("3D", "2Drac"):
        print("For the ligand database type, use either '3D' or '2Drac'")
        sys.exit()

    return obPath, inxPath, dbType, mapMode, pocket


def getPath():
    """
    Devide which platform this is executed on, and select the proper script
    paths
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


def generateScript(mapMode, obPath, pocket, icm):
    """
    The scripts are generated here, and their content is modified to fit the
    tasks they are supposed to carry out
    """

    # Ligand-script base
    scrLig_string = """#!ICM_EXEC

call "_startup"

openFile "VS_PROJ.ob"

# Create sphere around the bound ligand
as_graph = Sphere( a_VS_PROJ.2 , ( ! a_VS_PROJ.2 ) & Obj( a_VS_PROJ.2) 4. )
# Delete the ligand
delete a_VS_PROJ.2

# Setup docking project
currentDockProj.data[8] = "yes"
tempsel = as_graph
dock2SetupReceptor "VS_PROJ" a_ tempsel no "none"
dock5CalcMaps "VS_PROJ" 0.5 4.0 no
currentDockProj.data[1] = "VS_PROJ"

quit
"""

    # Pocket-script base
    scrPok_string = """#!ICM_EXEC

call "_startup"

openFile "VS_PROJ.ob"

# Get pocket
icmPocketFinder Mol(a_*.//DD) & a_*.!H,W 4.6 no no
# Create sphere around fist (best) pocket
as_graph = Sphere( g_pocketPOCKET_NUM a_VS_PROJ. 2.5)

# Setup docking project
currentDockProj.data[8] = "yes"
tempsel = as_graph
dock2SetupReceptor "VS_PROJ" a_ tempsel no "none"
dock5CalcMaps "VS_PROJ" 0.5 4.0 no
currentDockProj.data[1] = "VS_PROJ"

quit
"""

    # Modify the required script and return its path
    projName = obPath.replace(".ob", "")
    if mapMode == "ligand":
        scr_string = scrLig_string
        scr_string = scr_string.replace("ICM_EXEC", icm)
        scr_string = scr_string.replace("VS_PROJ", projName)
    elif mapMode == "pocket":
        scr_string = scrPok_string
        scr_string = scr_string.replace("ICM_EXEC", icm)
        scr_string = scr_string.replace("VS_PROJ", projName)
        scr_string = scr_string.replace("POCKET_NUM", pocket)

    # Write the selected script to a file
    workDir = os.getcwd()
    scr_path = workDir + "/temp.icm"
    with open(scr_path, "w") as scr_file:
        scr_file.write(scr_string)

    return scr_path


def runScript(icm, script):
    """
    This modifies and runs the script creating the maps around the bound ligand
    """

    # Execute
    try:
        check_output(icm + " -s " + script, stderr=STDOUT, shell=True)
    except CalledProcessError as e:
        print(e.output)
        sys.exit()

    # Delete temp script
    os.remove(script)


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

    # print dtbLines[i]

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
