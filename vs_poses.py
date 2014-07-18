#!/usr/bin/env python


###############################################################################
#
#   Uses the results.csv previously generated to locate the top X docking poses
#   following a VS. Loads them using ICM and saves the poses wanted to a single
#   .pdb file. Also saves the receptor to that .pdb file. That file can then be
#   opened using ICM or an other molecular viewer.
#
#   Thomas Coudrat, July 2014
#
###############################################################################


import csv
import sys
import os
import glob
import shutil
from subprocess import check_output, STDOUT, CalledProcessError


def main():
    """
    Run the scripts
    """

    X = 10
    icmBin = "/usr/icm-3.7-3b/icm64"

    # Get paths
    resultsPath = sys.argv[1]
    cwd = os.getcwd()
    vsPath = os.path.dirname(cwd + "/" + resultsPath)
    projName = os.path.basename(os.path.dirname(cwd + "/" + resultsPath))

    # Parse the VS results
    resData = parseResultsCsv(resultsPath, X)

    # Get the list of poses to pick per repeat directory
    repeatsRes = posesPerRepeat(resData)

    # Load those poses and save them in the /poses directory
    loadPoses(repeatsRes, vsPath, projName, icmBin)

    # Now load and write the receptor (binding pocket)
    recObName = projName + "_rec"
    recObPath = vsPath + "/vs_setup/" + recObName + ".ob"
    recPdbPath = vsPath + "/poses/" + recObName + ".pdb"
    print "Extracting", recObName
    print
    readAndWrite([recObPath], [["a_" + recObName + ".", recPdbPath]], icmBin)


def parseResultsCsv(resPath, X):
    """
    Read the resultsPath, and extract the location (which repeat) of each ligand
    for which the docking pose should be extracted
    """

    # Read data with the csv reader
    resFile = open(resPath, "r")
    resData = csv.reader(resFile)
    resData = [row for row in resData]
    resFile.close()

    # Generate a selection of the top X docked ligands, selecting only the
    # ligandID [0], the ICM score [9], and the repeat directory [-1]
    resData = [[ID, eval(row[0]), eval(row[9]), row[-1]] for row, ID in
               zip(resData[1:X+1], range(1, X+1))]

    return resData


def posesPerRepeat(resData):
    """
    Create a list for each repeat directory with the ordered list of ligIDs that
    should be extracted to be part of the results
    """

    # Store the ligand data into a dictionary where the keys are repeat numbers
    repeatsRes = {}

    for row in resData:
        rep = row[3]
        keys = repeatsRes.keys()
        # if the repeat already exists in the keys, add the row to the list
        if rep in keys:
            repeatsRes[rep].append(row)
            repeatsRes[rep].sort
        # otherwise create a new list with that row
        else:
            repeatsRes[rep] = [row]

    return repeatsRes


def loadPoses(repeatsRes, vsPath, projName, icmBin):
    """
    Walk through repeat directories, and load each
    """
    resultsPath = vsPath + "/poses/"
    # Create the results directory, delete it if already exists
    if os.path.exists(resultsPath):
        shutil.rmtree(resultsPath)
    os.makedirs(resultsPath)

    print

    for key in repeatsRes.keys():
        # Update progress
        print "Extracting", len(repeatsRes[key]), "poses from repeat #", key

        # Get the ob file list
        repPath = vsPath + "/" + key + "/"
        obFileList = glob.glob(repPath + "*_answers*.ob")

        # Get the pdb file list
        pdbFileList = []
        for row in repeatsRes[key]:
            # print row
            ligPos = row[0]
            ligScore = row[2]
            ligID = row[1]
            pdbFilePath = resultsPath + str(ligPos) + "_" + str(ligScore) + \
                "_" + str(ligID) + ".pdb"
            icmName = "a_" + projName + str(ligID) + "."

            pdbFileList.append([icmName, pdbFilePath])

        readAndWrite(obFileList, pdbFileList, icmBin)


def readAndWrite(obFileList, pdbFileList, icmBin):
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
    # Add the ob file loading part
    icmScript.write("\n# OPENING FILES\n")
    for obFile in obFileList:
        icmScript.write('openFile "' + obFile + '"\n')
    # Add the pdb file saving part
    icmScript.write("\n# WRITING FILES\n")
    for pdbFile in pdbFileList:
        icmScript.write('write pdb ' + pdbFile[0] + ' "' + pdbFile[1] + '"\n')
    icmScript.write("\nquit")
    icmScript.close()

    # Execute temp script
    try:
        check_output(icmBin + " -s ./temp.icm", stderr=STDOUT, shell=True)
    except CalledProcessError, e:
        print e.output
        sys.exit()

    # Delete temp script
    os.remove("./temp.icm")


if __name__ == "__main__":
    main()
