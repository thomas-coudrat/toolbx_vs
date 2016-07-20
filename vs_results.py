#!/usr/bin/env python

# ------------------------------------------------
#
#   Extracts the results from all .ou files contained
#   in the repeats of the current VS directory
#   Regroups the repeats together and extracts either only
#   the best score for each ligand, or all repeats.
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
    vsDir, minRep, allRep = parseArguments()

    # Get the project name out of the vsDir
    projName = os.path.basename(os.path.normpath(vsDir))
    # Fix project name if script was called within its directory
    if projName == ".":
        projName = os.path.basename(os.getcwd())

    # Create the dictionary storing ligand info
    # based on ligandID: for each ligandID key there
    # is a number of ligangInfo lists equal to the
    # number of repeats
    ligDict = {}

    # Goes through repeat directories to gather the score data
    # Returns ligDict (VS results) total number of repeats
    ligDict, totalRepeatNum = collectScoreData(vsDir, ligDict)

    # Getting rid of the ligands that were not docking in all repeats attempted
    ligDict = removeFailed(ligDict, totalRepeatNum, minRep)

    # Sort each ligand docking amongst repeats
    ligDict = sortRepeats(ligDict)

    # Write the results in a .csv file
    writeResultFiles(ligDict, projName, vsDir)

    # Write out individual results files for each repeat, if requested
    if allRep:
        # For each repeat,
        for repeat in range(1, totalRepeatNum + 1):
            # Initialise a results text file
            repFileName = "repeat{}_results_{}.csv".format(repeat, projName)
            print("\t" + repFileName)
            repFile = open(vsDir + "/" + repFileName, "w")
            repFile.write("No,Nat,Nva,dEhb,dEgrid,dEin,dEsurf" +
                          ",dEel,dEhp,Score,mfScore,Name,Run#\n")

            # Then, extract ligDict data corresponding to that repeat, sort,
            # and write to text file
            writeRepeatFile(repFile, repeat, ligDict)


def parseArguments():

    # Parsing description of arguments
    descr = "Extract VS results, write results and ROC data to file"
    descr_vsDir = "Directory of the VS to be analysed"
    descr_minRep = "Minimum number of repeats required to be included in" \
        " the results. Default is max number of repeats"
    descr_allRep = "Print out all results from each repeat in a different text" \
        " file"

    # Defining the arguments
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("vsDir", help=descr_vsDir)
    parser.add_argument("--minRep", help=descr_minRep)
    parser.add_argument("-allRep", action="store_true", help=descr_allRep)

    # Parsing arguments
    args = parser.parse_args()
    vsDir = args.vsDir
    minRep = args.minRep
    allRep = args.allRep

    # Deal with minRep in case the option was not used in which case use a very
    # large int number. Otherwise make the minRep an int.
    if minRep:
        minRep = int(minRep)
    else:
        # This could be improved, but it does work well this way (never will the
        # repeat number be that high)
        minRep = 999999999999999999999

    return vsDir, minRep, allRep


def collectScoreData(vsDir, ligDict):
    """
    Go through the repeat directories and collect the score data
    """

    print("\nPARSING:\n")

    maxRepeatNum = -1

    # Get all .ou files in each repeat directory
    ouFiles = glob.glob(vsDir + "/*/*.ou")
    # Loop through them and look for the 'SCORES' line
    for ouFilePath in ouFiles:
        # Open file containing text result of the VLS
        file = open(ouFilePath, "r")
        lines = file.readlines()
        file.close()

        vs_dir = os.path.dirname(os.path.dirname(ouFilePath))
        repeatNum = os.path.dirname(ouFilePath).replace(vs_dir + "/", "")
        # print ouFilePath
        # print repeatNum

        # Loop through each line of the file
        ligDockedNum = 0
        for line in lines:
            # We take only the lines that contain "SCORE>"
            if "SCORES>" in line:
                ligDockedNum += 1
                ligDict = parseScoreLine(ligDict, line, repeatNum)

        print("\t" + ouFilePath + "\t" + str(ligDockedNum) + " ligands")

        # Update the repeat number in order to grab the max repeat number
        if maxRepeatNum < int(repeatNum):
            maxRepeatNum = int(repeatNum)

    return ligDict, maxRepeatNum


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

    return ligDict


def removeFailed(ligDict, totalRepeatNum, minRepeatNum):
    """
    Loop over all results and remove those not successful for all repeats
    attempted. Print the information about the failed dockings.
    """

    # Get the max ligID of the docked ligands
    ligIDs = sorted(ligDict.keys())
    minLigID = ligIDs[0]
    maxLigID = ligIDs[-1]
    # Create a range list of IDs, stopping at the max
    rangeIDs = range(minLigID, maxLigID + 1)
    rangeFlag = dict([(ligID, False) for ligID in rangeIDs])
    # print rangeFlag

    print("\nINCOMPLETE DOCKINGS:\n")

    #keys = ligDict.keys()
    for key in ligIDs:
        currRepeatNum = len(ligDict[key])

        # When the number of repeats found is not equal to the max number of
        # repeats expected
        if currRepeatNum != totalRepeatNum:
            print("\tid:" + str(key) + "# of sucessful repeats:" +
                  str(currRepeatNum))
            # For cases where a ligand was docked more than the defined repeat
            # number (when there was mistake in the VS setup)
            if currRepeatNum > totalRepeatNum:
                print("\t\t(included)")
            # For cases where the repeat number of a given ligand is above or
            # equal to the user defined minimum repeat number
            elif currRepeatNum >= minRepeatNum:
                print("\t\t(included)")
            # Otherwise delete the ligand's information from the list
            else:
                print("\t\t(deleted)")
                del ligDict[key]

        # Flag the current ligID when it is found
        if key in rangeFlag.keys():
            rangeFlag[key] = True

    for key in rangeFlag.keys():
        if rangeFlag[key] is False:
            print("\tid:", key, "# of successful repeats: 0 (not included)")

    print("\nSUMMARY:\n")

    print("\tTotal ligands docked:" + str(len(ligDict.keys())))

    return ligDict


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

    return ligDict


def writeResultFiles(ligDict, projName, vsDir):
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

    print("\nWRITING:\n")

    # Create results file
    print("\tresults_" + projName + ".csv")
    fileResult = open(vsDir + "/results_" + projName + ".csv", "w")
    fileResult.write("No,Nat,Nva,dEhb,dEgrid,dEin,dEsurf" +
                     ",dEel,dEhp,Score,mfScore,Name,Run#\n")

    # Loop over the sorted results, and write to the result(s) file(s)
    for ligInfo in vsResult:
        # Write single repeat result (the best repeat)
        writeResultLine(ligInfo, fileResult)

    fileResult.close()

def writeRepeatFile(repFile, repeat, ligDict):
    """
    Get a repeat result file and repeat number. Extract VS data corresponding
    to that repeat from ligDict, sort the results according to score, and
    write those to the corresponding repeat results text file.
    """
    # Store VS repeat results in a list
    repeatResults = []

    # Get the repeat results from the ligDict, looping over ligand ID
    for ligID in ligDict.keys():
        # For each ligID, loop over all repeats
        for repeatLigInfo in ligDict[ligID]:
            # If this is the repeat that matches the one we want (last item)
            if repeat == int(repeatLigInfo[-1]):
                # Append to the list
                repeatResults.append(repeatLigInfo)

    # Sort the repeat VS result list, based on score
    repeatResults = sorted(repeatResults, key=lambda lig: lig[9])

    # Write results to file
    for repeatResultsLine in repeatResults:
        writeResultLine(repeatResultsLine, repFile)
    repFile.close()


def writeResultLine(ligInfo, fileResult):
    """
    Write single results line to file
    """
    ligInfoStr = []
    for val in ligInfo:
        ligInfoStr.append(str(val))

    fileResult.write(",".join(ligInfoStr))
    fileResult.write("\n")


if __name__ == "__main__":
    main()
