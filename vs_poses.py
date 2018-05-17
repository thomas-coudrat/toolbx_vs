#!/usr/bin/env python

# Uses the results.csv previously generated to locate the X docking poses
# following a VS. Loads them using ICM and saves the poses wanted to a single
# .pdb file. Also saves the receptor to that .pdb file. That file can then be
# opened using ICM or an other molecular viewer.

# https://github.com/thomas-coudrat/toolbx_vs
# Thomas Coudrat <thomas.coudrat@gmail.com>

import csv
import sys
import os
import glob
import shutil
import argparse
import socket
from subprocess import check_output, STDOUT, CalledProcessError
import json


def main():
    """
    Run script
    """

    # Get the arguments and paths
    resultsPath, ligs_from_A, ligs_to_B, ligIDs, label = parseArgs()
    icm = getPath()

    # Fix paths
    cwd = os.getcwd()
    vsPath = os.path.dirname(cwd + "/" + resultsPath)
    #projName = os.path.basename(os.path.dirname(cwd + "/" + resultsPath))
    projName = os.path.basename(glob.glob(cwd + "/*.log")[0]).replace(".log", "")

    # Parse the VS results
    resDataAll = parseResultsCsv(resultsPath)

    # Select results, the X docked poses and/or optionally specific ligIDs
    resDataSel = selectResults(resDataAll, ligs_from_A, ligs_to_B, ligIDs)

    # Print the selected results
    printResults(resDataSel)

    # Get the list of poses to pick per repeat directory
    repeatsRes = posesPerRepeat(resDataSel)

    # Load those poses and save them in the /poses directory
    loadAnswersWritePoses(repeatsRes, vsPath, projName, label, icm)

    print("\n")


def parseArgs():
    """
    Create arguments and parse them
    """
    # Parsing description
    descr = "Extract docking poses from a VS in .pdb format"
    descr_resultsPath = "Results file of the VS in .csv format"
    descr_from_A = "Extract ligands from A (indexed at 0)"
    descr_to_B = "Extract ligands to B (indexed at 0, up to but not including)"
    descr_ligIDs = "Optional ligIDs docking poses to be extracted. Format " \
        "e.g. 1-10,133,217-301"
    descr_label = "Add a label suffix to all poses generated"

    # Define arguments
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("resultsPath", help=descr_resultsPath)
    parser.add_argument("from_A", help=descr_from_A)
    parser.add_argument("to_B", help=descr_to_B)
    parser.add_argument("--ligIDs", help=descr_ligIDs)
    parser.add_argument("--label", help=descr_label)

    # Parse arguments
    args = parser.parse_args()
    resultsPath = args.resultsPath
    from_A = int(args.from_A)
    to_B = int(args.to_B)
    ligIDs = args.ligIDs

    if args.label:
        label = args.label
    else:
        label = ""

    # Make ligIDs a list of IDs integers
    if ligIDs:
        ligIDs = makeIDlist(ligIDs)

    return resultsPath, from_A, to_B, ligIDs, label


def makeIDlist(stringID):
    """
    Get a string defining which IDs to be generated into a list
    """

    # This stores the range of IDs into a list
    rangeID = []

    IDportions = stringID.split(",")

    for portion in IDportions:
        # Treat ranges of IDs
        if "-" in portion:
            start, end = portion.split("-")
            start = int(start)
            end = int(end)
            # Do not add the value 0 to the list
            if start == 0 or end == 0:
                pass
            else:
                rangeID = rangeID + range(start, end + 1)
        # Treat single IDs
        else:
            portion = int(portion)
            # Do not add the value 0 to the list
            if portion == 0:
                pass
            else:
                rangeID.append(int(portion))

    return rangeID


def getPath():
    """
    Get the ICM executable from the ICMHOME environment variable
    """

    # Get environment variable. It returns None if not defined on the system.
    icmHome = os.environ.get('ICMHOME')

    # Return path to executable if the environment variable was found
    if icmHome is None:
        print("ICMHOME env. var. must be set for your system. Exiting.")
        sys.exit()
    else:
        icm = icmHome + "/icm64"

    return icm


def parseResultsCsv(resPath):
    """
    Read the resultsPath, and return the location (which repeat) of each ligand
    for which the docking pose should be extracted
    """

    # Read data with the csv reader
    resFile = open(resPath, "r")
    resData = csv.reader(resFile)
    resData = [row for row in resData]
    resFile.close()

    resLen = len(resData)

    # Select the following:
    # overall rank (ID),
    # ligandID row[0],
    # the ligand name row[12],
    # ICM score [9],
    # repeat directory [12]
    # note: skip the first line of the CSV file resData (header)
    return [[ID, eval(row[0]), row[11], eval(row[9]), row[12]] for row, ID in
            zip(resData[1:], range(1, resLen))]


def selectResults(resDataAll, ligs_from_A, ligs_to_B, ligIDs):
    """
    Make a selection of the X ligands from the VS results, also as optional
    ligIDs
    """

    # First select X ligands
    resDataSel = [row for row in resDataAll[ligs_from_A:ligs_to_B]]

    # If the ligID flag was used, select by ligIDs, the ligID is in row[1]
    # and return a combined list: if a selected ligID is within the selected X,
    # there will be a duplicate. However when processed the .pdb docking pose
    # will be written then overwritten, resulting in the result expected
    if ligIDs:
        resDataID = [row for row in resDataAll if row[1] in ligIDs]
        return resDataSel + resDataID
    # Otherwise just return the 'sel' list
    else:
        return resDataSel


def printResults(resData):
    """
    Takes in results and prints them out
    """

    print()
    print("Rank\tID\tName\tScore\tRepeat")
    for res in resData:
        print(str(res[0]) + "\t" +
              str(res[1]) + "\t" +
              str(res[2]) + "\t" +
              str(res[3]) + "\t" +
              str(res[4]))
        if res[3] > -25.:
            print("\t (score above -25., may not have been saved)")
        else:
            print("\t")


def posesPerRepeat(resData):
    """
    Create a list for each repeat directory with the ordered list of ligIDs
    that should be extracted to be part of the results
    """

    # Store the ligand data into a dictionary where the keys are repeat numbers
    repeatsRes = {}

    for row in resData:
        rep = row[4]
        keys = repeatsRes.keys()
        # if the repeat already exists in the keys, add the row to the list
        if rep in keys:
            repeatsRes[rep].append(row)
            repeatsRes[rep].sort
        # otherwise create a new list with that row
        else:
            repeatsRes[rep] = [row]

    return repeatsRes


def loadAnswersWritePoses(repeatsRes, vsPath, projName, label, icm):
    """
    Walk through repeat directories, and load each
    """
    # Create the results directory, delete it if already exists
    resultsPath = vsPath + "/poses/"
    if os.path.exists(resultsPath):
        shutil.rmtree(resultsPath)
    os.makedirs(resultsPath)

    print()

    for key in repeatsRes.keys():
        # Update progress
        # Get the ob file list
        repPath = vsPath + "/" + key + "/"
        # Ligands found in that repeat
        ligsInfo = [[row[1], ""] for row in repeatsRes[key]]
        obFileList = getAnswersList(repPath, ligsInfo)

        # Get the pdb file list
        pdbFileList = []
        for row in repeatsRes[key]:
            # print(row)
            ligRank = row[0]
            ligID = row[1]
            ligName = row[2]
            ligScore = row[3]
            pdbFilePath = resultsPath + str(ligRank) + "_" + str(ligScore) + \
                "_" + str(ligID) + "_" + ligName + "_" + label + ".pdb"
            selectionName = projName.replace("-", "_") + str(ligID)
            icmName = "a_" + selectionName + "."

            pdbFileList.append([selectionName, icmName, pdbFilePath])

        readAndWrite(obFileList, pdbFileList, projName, vsPath, icm)


def getAnswersList(repPath, ligsInfo):
    """
    Given a repeat directory path and a list of ligand IDs, return the list
    of VS answers (.ob files) that contain all the ligand IDs provided.
    """

    # Get a list of all the files, and make is an sorted list
    allObFiles = glob.glob(repPath + "*_answers*.ob")
    allObFiles = [[obFile, int(obFile.split("_answers")[1].replace(".ob", ""))]
                  for obFile in allObFiles]
    sortedObFiles = sorted(allObFiles, key=lambda obFile: obFile[1],
                           reverse=True)
    for obFile in sortedObFiles:
        obFilePath = obFile[0]
        obFileStart = obFile[1]

    # Loop over the ligsInfo list, which contains
    # ligInfo[0] = the ligand ID number
    # ligInfo[1] = an empty string that is used to store the path to the
    # corresponding .ob file to be loaded
    for lig in ligsInfo:
        ligID = lig[0]
        # answerPath = ligInfo[1]

        previousObStart = 999999999999999999999999999999999999999999999999

        # Loop over all answers file, and for each
        for obFile in sortedObFiles:
            obFileStart = obFile[1]
            obFilePath = obFile[0]

            # If the ligand ID number is between between the start of this file
            # and the start of the previous file (part of the current file)
            if ligID >= obFileStart and ligID < previousObStart:
                lig[1] = obFilePath

            previousObStart = obFileStart

    return list(set([lig[1] for lig in ligsInfo]))


def readAndWrite(obFileList, pdbFileList, projName, vsPath, icm):
    """
    Get the information of which *.ob files to read, and which *.pdb files from
    the loaded molecules to write.
    Edit a temporary version of the ICM script with that information, and
    execute it.
    This script is writen and executed for each repeat, selecting only the
    poses from each of those repeats
    """

    # Create temp script
    icmScript = open("./temp.icm", "w")
    icmScript.write('call "_startup"\n')

    # Load the receptor
    recObName = projName + "_rec"
    recObPath = vsPath + "/vs_setup/" + recObName + ".ob"
    print("loading receptor: " + recObPath)
    icmScript.write('openFile "' + recObPath + '"\n')

    # Load each OB file containing the ligands
    icmScript.write("\n# OPEN FILES\n")
    for obFile in obFileList:
        print("loading:" + obFile)
        icmScript.write('openFile "' + obFile + '"\n')

    # Create the receptor-ligand complexes, save complexes
    icmScript.write("\n# MODIFY/WRITE FILES\n")
    for pdbFile in pdbFileList:
        # Renumber ligand 'residue' number to avoid overlap with rec res nums
        icmScript.write('align number ' + pdbFile[1] + 'm 9999\n')
        # Rename receptor to 'a'
        icmScript.write('rename a_' + recObName + '.* Name(Name("a" simple),unique)\n')
        # Create a temporary copy of the receptor
        icmScript.write('copy Obj( a_' + recObName + '.a ) Name( "temp_rec" object ) delete display\n')
        # Move the receptor this ligand ICM object
        icmScript.write('move a_temp_rec.a ' + pdbFile[1] + '\n')
        # Save this new receptor-ligand complex to a pdb file
        icmScript.write('write pdb ' + pdbFile[1] + ' "' + pdbFile[2] + '"\n')

    icmScript.write("\nquit")
    icmScript.close()

    # Execute temp script
    try:
        check_output(icm + " -s ./temp.icm", stderr=STDOUT, shell=True)
    except CalledProcessError as e:
        print(e.output)
        sys.exit()

    # Delete temp script
    os.remove("./temp.icm")


if __name__ == "__main__":
    sys.exit(main())
