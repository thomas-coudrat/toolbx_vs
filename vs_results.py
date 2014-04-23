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
    vsDir, knownIDfirst, knownIDlast, \
        ommitIDfirst, ommitIDlast = parseArguments()

    # Get the project name out of the vsDir
    projName = os.path.basename(os.path.normpath(vsDir))

    # Create the dictionary storing ligand info
    # based on ligandID: for each ligandID key there
    # is a number of ligangInfo lists equal to the
    # number of repeats
    #
    ligDict = {}

    print "\nPARSING:\n"

    # Get all .ou files in each repeat directory
    #
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

    # Sort each ligand docking amongst repeats
    #
    sortRepeats(ligDict)

    # Write the results in a .csv file
    #
    vsResult = writeResultFile(ligDict, projName, vsDir)
    # Write ROC curve data to .roc file
    #
    writeROCfile(vsResult, projName, vsDir,
                 knownIDfirst, knownIDlast,
                 ommitIDfirst, ommitIDlast)


def parseArguments():

    # Parsing description of arguments
    descr = "Extract VS results, write results and ROC data to file"
    descr_vsDir = "Directory of the VS to be analysed"
    descr_knownIDrange = "Provide the ID range of known actives lig" \
                         "lib (format: 1-514)"
    descr_ommitIDrange = "Provide the ID range of ligands to ommit " \
                         "from the ROC curve data"

    # Defining the arguments
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("vsDir", help=descr_vsDir)
    parser.add_argument("knownIDrange", help=descr_knownIDrange)
    parser.add_argument("ommitIDrange", help=descr_ommitIDrange)

    # Parsing arguments
    args = parser.parse_args()
    vsDir = args.vsDir
    knownIDrange = args.knownIDrange
    knownIDfirst, knownIDlast = knownIDrange.split("-")
    ommitIDrange = args.ommitIDrange
    ommitIDfirst, ommitIDlast = ommitIDrange.split("-")

    return vsDir, int(knownIDfirst), int(knownIDlast), \
        int(ommitIDfirst), int(ommitIDlast)


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
    print "\tresults_" + projName + ".csv"

    fileResult = open(vsDir + "/results_" + projName + ".csv", "w")
    fileResult.write("No,Nat,Nva,dEhb,dEgrid,dEin,dEsurf" +
                     ",dEel,dEhp,Score,mfScore,Name,Run#\n")

    for ligInfo in vsResult:
        for val in ligInfo:
            fileResult.write(str(val) + ",")
        fileResult.write("\n")
    fileResult.close()

    return vsResult


def writeROCfile(vsResult, projName, vsDir,
                 knownIDfirst, knownIDlast,
                 ommitIDfirst, ommitIDlast):
    """
    Given this VS result, and information about the ID of known actives
    in the library, write in a file the information to plot a ROC curve
    """

    X = 0
    Y = 0
    knowns = "knowns_" + str(knownIDfirst) + "-" + str(knownIDlast)
    ommits = "ommits_" + str(ommitIDfirst) + "-" + str(ommitIDlast)

    # Create filename
    rocFileName = "roc_" + knowns + "_" + ommits + "_" + projName + ".csv"
    print "\t", rocFileName
    rocDataFile = open(vsDir + "/" + rocFileName, "w")

    # Get the total knowns and total library to calculate percentages
    totalKnowns = knownIDlast - knownIDfirst + 1
    totalLibrary = len(vsResult) - (ommitIDlast - ommitIDfirst + 1)

    print totalKnowns
    print totalLibrary - totalKnowns

    for ligInfo in vsResult:
        ligID = int(ligInfo[0])

        # Skip if ligID is part of the range that needs to be ommited
        if ligID in range(ommitIDfirst, ommitIDlast + 1):
            continue
        # Otherwise proceed normally
        else:
            # When the sorted ligID corresponds to a known, increase
            # the value of Y by 1
            if ligID in range(knownIDfirst, knownIDlast + 1):
                Y += 1

            # For each ligand in the full VS, increase X and write
            # the X,Y pair to the data file
            X += 1

            # Calculate percentage X and Y
            Xpercent = (X * 100.0) / totalLibrary
            Ypercent = (Y * 100.0) / totalKnowns
            rocDataFile.write(str(X) + "," + str(Y) + "," +
                              str(Xpercent) + "," + str(Ypercent) + "\n")

    rocDataFile.close()

    print


if __name__ == "__main__":
    main()
