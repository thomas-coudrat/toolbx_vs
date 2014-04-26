#!/usr/bin/env python

# ------------------------------------------------
#
#   Extracts the results from all .ou files contained
#   in the repeats of the current VS directory
#   Regroups the repeats together and extracts only
#   the best score for each ligand.
#
#   Thomas Coudrat, February 2014
#
# -------------------------------------------------

import glob
import os
import argparse


def main():
    """
    Execute VS results script
    """

    # Get arguments
    vsDir = parseArguments()

    # Get the project name out of the vsDir
    projName = os.path.basename(os.path.normpath(vsDir))

    # Create the dictionary storing ligand info
    # based on ligandID: for each ligandID key there
    # is a number of ligangInfo lists equal to the
    # number of repeats
    ligDict = {}

    print "\nPARSING:\n"

    maxRepeatNum = -1

    # Get all .ou files in each repeat directory
    ouFiles = glob.glob(vsDir + "/*/*.ou")
    # Loop through them and look for the 'SCORES' line
    for ouFilePath in ouFiles:
        # Open file containing text result of the VLS
        file = open(ouFilePath, "r")
        lines = file.readlines()
        file.close()

        print "\t", ouFilePath
        vs_dir = os.path.dirname(os.path.dirname(ouFilePath))
        repeatNum = os.path.dirname(ouFilePath).replace(vs_dir + "/", "")
        # print ouFilePath
        # print repeatNum

        # Loop through each line of the file
        for line in lines:
            # We take only the lines that contain "SCORE>"
            if "SCORES>" in line:
                parseScoreLine(ligDict, line, repeatNum)

        # Update the repeat number in order to grab the max repeat number
        if maxRepeatNum < int(repeatNum):
            maxRepeatNum = int(repeatNum)

    # Getting rid of the ligands that were not docking in all repeats attempted
    removeFailed(ligDict, maxRepeatNum)

    # Sort each ligand docking amongst repeats
    sortRepeats(ligDict)

    # Write the results in a .csv file
    writeResultFile(ligDict, projName, vsDir)


def parseArguments():

    # Parsing description of arguments
    descr = "Extract VS results, write results and ROC data to file"
    descr_vsDir = "Directory of the VS to be analysed"

    # Defining the arguments
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("vsDir", help=descr_vsDir)

    # Parsing arguments
    args = parser.parse_args()
    vsDir = args.vsDir

    return vsDir


def parseScoreLine(ligDict, line, repeatNum):
    """
    Populate the ligDict dictionary in the following manner:
    ligDict{ligandID, [[ligInfo_rep1], [ligInfo_rep2], ...]}
    """

    ll = line.split()
    # Store ligID unique identifyer
    ligID = int(ll[2])

    # Will contain all the info for 1 ligand
    ligInfo = []

    # The fist info is the ligID
    ligInfo.append(ligID)

    # Give a generic name for when the ligand does
    # not have one
    ligName = "none"

    # The rest of the info relates to the scoring
    for i, split in enumerate(ll):
        if "Name=" in split:
            ligName = ll[i + 1]
            break
        if "completed" in split or "FINISHED" in split:
            break
        # Store the values following each tag
        # (determined by the presnce of a '=')
        if "=" in split:
            val = ll[i + 1].rstrip("%FINISHED")
            # The score has to be stored as a float,
            # because it is used for sorting
            if split.strip() == "Score=":
                val = float(val)
            ligInfo.append(val)

    # Add the ligand name, which can be none when it is
    # not provided in the original .sdf library
    ligInfo.append(ligName)
    # Lastly adding the repeat number info
    ligInfo.append(repeatNum)

    # Add that ligInfo to the ligDict, if it already exists
    # just append to the list, otherwise create a new list
    keys = ligDict.keys()
    if ligID not in keys:
        ligDict[ligID] = [ligInfo]
    else:
        ligDict[ligID].append(ligInfo)


def removeFailed(ligDict, maxRepeatNum):
    """
    Loop over all results and remove those not successful for all repeats
    attempted. Print the information about the failed dockings.
    """

    print "\nINCOMPLETE DOCKINGS:\n"

    keys = ligDict.keys()
    for key in keys:
        if len(ligDict[key]) != maxRepeatNum:
            print "\tid:", key, "# of sucessful repeats:", len(ligDict[key])
            del ligDict[key]


def sortRepeats(ligDict):
    """
    For each ligandID, get the repeat that got the best score, this will
    represent that ligand in this VS scoring
    """

    # For each ligID, sort each repeat based on score (lig[9])
    # The result is a ligDict for which the first of each ligID is
    # the one with the best score
    keys = ligDict.keys()
    for key in keys:
        repeatsLigInfo = ligDict[key]
        repeatsLigInfo = sorted(repeatsLigInfo, key=lambda lig: lig[9])
        ligDict[key] = repeatsLigInfo


def writeResultFile(ligDict, projName, vsDir):
    """
    Write out the results of this VS
    """
    # Write the ligand info
    keys = ligDict.keys()
    vsResult = []
    for key in keys:
        # Get only the first in the list of repeats information
        # for this ligand
        # for ligInfo in ligDict[key]:
        ligInfo = ligDict[key][0]
        vsResult.append(ligInfo)

    # Sort the vsResult based on score, for the sorted full VS result
    vsResult = sorted(vsResult, key=lambda lig: lig[9])

    print "\nWRITING:\n"

    # Write result file
    if projName == ".":
        projName = os.path.basename(os.getcwd())
    print "\tresults_" + projName + ".csv"

    fileResult = open(vsDir + "/results_" + projName + ".csv", "w")
    fileResult.write("No,Nat,Nva,dEhb,dEgrid,dEin,dEsurf" +
                     ",dEel,dEhp,Score,mfScore,Name,Run#\n")

    for ligInfo in vsResult:
        ligInfoStr = []
        for val in ligInfo:
            ligInfoStr.append(str(val))

        fileResult.write(",".join(ligInfoStr))
        fileResult.write("\n")
    fileResult.close()


if __name__ == "__main__":
    main()
