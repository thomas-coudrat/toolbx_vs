#!/vlsci/VR0024/tcoudrat/bin/python275/bin/python
#   Browse VS results contained in .ou files of several
#   repeats, and write a results list
#
#   Thomas Coudrat, February 2014
#
#--------------------------------------------------------

import os
import glob


def main():

    # Get all .ou files in each repeat directory
    ouFiles = glob.glob("*/*.ou")

    print
    print "PARSING:"
    print

    # Get the results from those .ou files
    ligList = gatherResults(ouFiles)

    print
    print "WRITING:"
    print

    # Write the results in a .csv file
    writeResultFile(ligList)


def gatherResults(ouFiles):

    # Storing ligand info
    ligList = []

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
                # Will contain all the info for 1 ligand
                ligInfo = []

                # Store "Lig number"
                ligInfo.append(int(ll[2]))

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
                # Adding the run number information to the ligInfo
                ligInfo.append(repeatNum)
                # Adding the information of 1 single ligand to the list
                ligList.append(ligInfo)

    return ligList


def writeResultFile(ligList):

    # Write result file
    fileResult = open("ranked_result.csv", "w")
    fileResult.write("No,Nat,Nva,dEhb,dEgrid,dEin,dEsurf,dEel,dEhp,Score,mfScore,Name,Run#\n")


    # After all information has been gathered,
    # sorting the list based on the score of each ligand
    ligList = sorted(ligList, key=lambda lig: lig[9])

    for lig in ligList:
        for info in lig:
            fileResult.write(str(info) + ",")
        fileResult.write("\n")

    fileResult.close()


if __name__ == "__main__":
    main()
