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

def main():

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
    writeResultFile(ligDict)


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

    # Write result file
    print "\tranked_results.csv"
    print
    fileResult = open("ranked_results.csv", "w")
    fileResult.write("No,Nat,Nva,dEhb,dEgrid,dEin,dEsurf,dEel,dEhp,Score,mfScore,Name,Run#\n")


    # After all information has been gathered,
    # sorting the list based on the score of each ligand
    #ligList = sorted(ligList, key=lambda lig: lig[9])

    # Write the ligand info
    keys = ligDict.keys()
    for key in keys:
        # Get only the first in the list of repeats information
        # for this ligand
        #for ligInfo in ligDict[key]:
        ligInfo = ligDict[key][0]
        for val in ligInfo:
            fileResult.write(str(val) + ",")
        fileResult.write("\n")
    fileResult.close()


if __name__ == "__main__":
    main()
