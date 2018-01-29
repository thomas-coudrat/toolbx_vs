#!/usr/bin/env python

# This creates maps for a VS from a .ob object containing
# the target protein only. The .inx index file pointing
# to the library to be used must also be provided.
# It uses an .icm script to create the ICM docking maps.
#
# https://github.com/thomas-coudrat/toolbx_vs
# Thomas Coudrat <thomas.coudrat@gmail.com>

import os
import sys
import argparse
import socket
from subprocess import check_output, STDOUT, CalledProcessError
import json

def main():
    """
    Run script
    """

    # Get the paths to the scripts to be executed
    icm = getPath()

    # Get the arguments
    obPath, inxPath, dbType, mapMode, pocket, residues = parseArgs()

    # Generate script (will return the script relevant to the mode chosen)
    script = generateScript(mapMode, obPath, pocket, icm, residues)

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
    descr_mapMode = "Choose map creation mode: 'ligand' create a map around " \
                    "bound ligand (on object containing ligand). Deletes " \
                    "ligand." \
                    "'pocket' searches for binding pocket using " \
                    "ICMpocketFinder (on object not containing ligand)." \
                    "'residues' selects a list of residues submitted through " \
                    "a text file"
    descr_dbType = "Database type: '3D' for 3D chiral compounds or '2Drac' " \
                   " for 2D racemic compounds, to be sampled"
    descr_pocket = "Define the pocket number to use (determined by " \
                   "ICMpocket finder). Default is pocket #1."
    descr_residues = "Provide a file containing the list of residues numbers" \
                     " to be used for map creation. Format: 1,2,3:10,20"

    # Parse arguments
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("obPath", help=descr_obPath)
    parser.add_argument("inxPath", help=descr_inxPath)
    parser.add_argument("dbType", help=descr_dbType)
    parser.add_argument("mapMode", help=descr_mapMode)
    parser.add_argument("--pocket", help=descr_pocket)
    parser.add_argument("--resPath", help=descr_residues)
    args = parser.parse_args()

    # Store arguments
    obPath = args.obPath
    inxPath = args.inxPath
    dbType = args.dbType
    mapMode = args.mapMode
    pocket = args.pocket
    resPath = args.resPath

    # Deal with arguments
    if not pocket:
        pocket = "1"

    if mapMode not in ("ligand", "pocket", "residues"):
        print ("Either use option 'pocket', 'ligand' or 'residues' for map mode creation")
        sys.exit()

    if dbType not in ("3D", "2Drac"):
        print("For the ligand database type, use either '3D' or '2Drac'")
        sys.exit()

    if mapMode == "residues" and resPath == False:
        print("Option 'residues' requires the user to provide a residues list using --residues flag")
        sys.exit()

    # Read in the residue list
    if mapMode == "residues":
        with open(resPath) as f:
            lines = f.readlines()
            residues = "".join(lines).strip()
    # otherwise initialise residues to empty string
    else:
        residues = ""

    return obPath, inxPath, dbType, mapMode, pocket, residues


def getPath():
    """
    Get the ICM executable from the ICMHOME environment variable
    """

    # Get environment variable. It returns None if not defined on the system.
    icmHome = os.environ.get('ICMHOME')

    # Return path to executable if the environment variable was found
    if icmHome == None:
        "The ICMHOME environment variable must be set for your system. Exiting."
        sys.exit()
    else:
        icm = icmHome + "/icm64"

    return icm


def generateScript(mapMode, obPath, pocket, icm, residues):
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

    # resList-script base
    resList_string = """#!ICM_EXEC

call "_startup"

openFile "VS_PROJ.ob"

# Select residues from the given list
as_graph = a_VS_PROJ./RESIDUE_LIST

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
    elif mapMode == "residues":
        scr_string = resList_string
        scr_string = scr_string.replace("ICM_EXEC", icm)
        scr_string = scr_string.replace("VS_PROJ", projName)
        scr_string = scr_string.replace("RESIDUE_LIST", residues)

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
