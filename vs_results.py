#!/vlsci/VR0024/tcoudrat/bin/python275/bin/python

#------------------------------------------------
#
#   Extracts the results from all .ou files contained
#   in the repeats of the current VS directory
#   Regroups the repeats together and extracts only
#   the best score for each ligand.
#
#   Thomas Coudrat, February 2014
#
#-------------------------------------------------

import glob
import os
import argparse

def main():

    # Get arguments
    knownIDfirst, knownIDlast = parseArguments()

    # Get all .ou files in each repeat directory
    ouFiles = glob.glob("*/*.ou")

    print
    print "PARSING:"
    print

    # Create the dictionary storing ligand info
    # based on ligandID: for each ligandID key there
    # is a number of ligangInfo list equal to the
    # number of repeats
    ligDict = {}

    # Get the results from those .ou files
    parseResults(ligDict, ouFiles)

    # Sort each ligand docking amongst repeats
    sortRepeats(ligDict)

    print
    print "WRITING:"
    print

    # Write the results in a .csv file
    vsResult = writeResultFile(ligDict)

    # Write ROC curve data to .roc file
    writeROCdata(vsResult, knownIDfirst, knownIDlast)

    print


def parseArguments():

    # Parsing description of arguments
    descr = "Extract VS results, write results and ROC data to file"
    descr_knownIDrange = "Provide the ID range of known actives lig lib (format: 1-514)"

    # Defining the arguments
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("knownIDrange", help=descr_knownIDrange)

    # Parsing arguments
    args = parser.parse_args()
    knownIDrange = args.knownIDrange
    knownIDfirst, knownIDlast = knownIDrange.split("-")

    return int(knownIDfirst), int(knownIDlast)


def parseResults(ligDict, ouFiles):
    """
    Populate the ligDict dictionary in the following manner:
    ligDict{ligandID, [[ligInfo_rep1], [ligInfo_rep2], ...]}
    """

    # Loop over all .ou files, store ligand docking info
    for ouFilePath in ouFiles:

        # Open file containing text result of the VLS
        file = open(ouFilePath, "r")
        lines = file.readlines()
        file.close()

        print "\t", ouFilePath
        repeatNum = os.path.dirname(ouFilePath)

        #Loop through each line of the file
        for line in lines:

            # We take only the lines that contain "SCORE>"
            if "SCORES>" in line:

                ll = line.split()
                # Store ligID unique identifyer
                ligID = int(ll[2])

                # Will contain all the info for 1 ligand
                ligInfo = []

                # The fist info is the ligID
                ligInfo.append(ligID)

                # The rest of the info relates to the scoring
                for i, split in enumerate(ll):
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

                # Lastly adding the repeat number info
                ligInfo.append(repeatNum)

                # Add that ligInfo to the ligDict, if it already exists
                # just append to the list, otherwise create a new list
                keys = ligDict.keys()
                if ligID not in keys:
                    ligDict[ligID] = [ligInfo]
                else:
                    ligDict[ligID].append(ligInfo)
                ## Adding the information of 1 single ligand to the list
                #ligList.append(ligInfo)


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


def writeResultFile(ligDict):


    # Write the ligand info
    keys = ligDict.keys()
    vsResult = []
    for key in keys:
        # Get only the first in the list of repeats information
        # for this ligand
        #for ligInfo in ligDict[key]:
        ligInfo = ligDict[key][0]
        vsResult.append(ligInfo)

    # Sort the vsResult based on score, for the sorted full VS result
    vsResult = sorted(vsResult, key=lambda lig: lig[9])

    # Write result file
    print "\tranked_results.csv"

    fileResult = open("ranked_results.csv", "w")
    fileResult.write("No,Nat,Nva,dEhb,dEgrid,dEin,dEsurf,dEel,dEhp,Score,mfScore,Name,Run#\n")

    for ligInfo in vsResult:
        for val in ligInfo:
            fileResult.write(str(val) + ",")
        fileResult.write("\n")
    fileResult.close()

    return vsResult


def writeROCdata(vsResult, knownIDfirst, knownIDlast):
    """
    Given this VS result, and information about the ID of known actives in the library,
    write in a file the information to plot a ROC curve
    """

    X = 0
    Y = 0
    #libSize = 0

    rocFileName = "rocKnown_" + str(knownIDfirst) + "-" + str(knownIDlast) + ".csv"
    print "\t", rocFileName
    rocDataFile = open(rocFileName, "w")

    for ligInfo in vsResult:
        ligID = int(ligInfo[0])
        # Update the libSize in order to get the full libSize
        #if libSize < ligID:
        #    libSize = ligID
        # When the sorted ligID corresponds to a known, increase
        # the value of Y by 1
        if ligID in range(knownIDfirst, knownIDlast + 1):
            Y += 1
        # For each ligand in the full VS, increase X and write
        # the X,Y pair to the data file
        X += 1
        rocDataFile.write(str(X) + "," + str(Y) + "\n")

    rocDataFile.close()


if __name__ == "__main__":
    main()
