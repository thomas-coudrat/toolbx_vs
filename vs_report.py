#!/usr/bin/env python

# Run in a VS repeat directory, checks all .ou files and
# compiles the number of occurence of the 'SCORE' word in
# order to inform about the status of that repeat
#
# https://github.com/thomas-coudrat/toolbx_vs
# Thomas Coudrat <thomas.coudrat@gmail.com>

import glob
import os
import argparse

def main():
    """
    Run script
    """

    # Setting up variables
    workDir = os.getcwd()

    print("\n************************\n")

    specs, skipCount = loopOverRepeats(workDir)

    # Print the skipped ligands data
    printSkipped(specs, skipCount)

    # Print out the content of the slurm files, which
    # will have ERROR information, or will be blank
    # if no Error
    printSlurmOuts(workDir)


def loopOverRepeats(workDir):
    """
    Loop over the repeats in this VS directory, gathering
    information about the VS status
    """

    # Dictionary containing specs as keys, and [MIN,MAX,COUNT] as values
    specs = {}
    skipCount = 0

    # Loop through the directories in this VS
    for subDir in os.listdir(workDir):
        # Check only repeat directories
        if subDir.isdigit():
            dirPath = os.path.join(workDir, subDir)

            # Get all the *.ou files in this dir
            ouFiles = glob.glob(dirPath + "/*.ou")

            scoreCount = 0
            # Looping over .ou files in the current dir
            for file in ouFiles:

                # Read all lines of that .ou file
                f = open(file, "r")
                lines = f.readlines()
                f.close()

                # Loop through the lines in that file,
                # and collect both the SCORE count and
                # the Skipped information
                for line in lines:
                    # Update "SCORE" count
                    if "SCORE" in line:
                        scoreCount +=1
                    # Update "Skipping" count
                    specs, skipCount = countSkipped(line, specs, skipCount)

            printCompleted(subDir, scoreCount)

    return specs, skipCount


def printCompleted(subDir, scoreCount):
    """
    Print the score count for the current VS repeat
    """

    print("SCORE COUNT FOR " + subDir + ": " + str(scoreCount))


def countSkipped(line, specs, skipCount):
    """
    Keep track of the number of skipped ligands, and the criteria that
    got them rejected
    """

    # If the word 'Skipping' is found,
    if "Skipping" in line:
        skipCount += 1
        endLine = line.split(",")[1].strip()
        endll = endLine.split()

        # Special case for LogP, the line is not written the same way
        if endll[0] == "LogP":
            currSpec = endll[0]
            currVal = float(endll[1])
            #print endll[0]
        else:
            currSpec = ''.join(endll[0:2])
            currVal = int(endll[-3].rstrip("."))
            #print ''.join(endll[0:2])


        # If the key did not exist before, create it (not seen before spec)
        keys = specs.keys()
        if currSpec not in keys:
            specs[currSpec] = [currVal, currVal, 0]
        # If it had already been created
        else:
            # Increment its count by one
            specs[currSpec][2] += 1
            # And update the MIN, MAX values
            if specs[currSpec][0] > currVal:
                specs[currSpec][0] = currVal
            if specs[currSpec][1] < currVal:
                specs[currSpec][1] = currVal

        # This is to pring the full line, for debugging
        #print endLine

    return specs, skipCount


def printSkipped(specs, skipCount):
    """
    Print out the results gathered related to the ligands skipped
    """
    # Print out the result
    keys = specs.keys()

    print("\n************************")
    print("THERE WERE", skipCount, "LIGANDS SKIPPED:\n")

    for key in keys:
        print(key)
        print("MIN:{:>10} \t MAX:{:>10} \t COUNT:{:>10}".format(specs[key][0],
                                                                specs[key][1],
                                                                specs[key][2]))
        #print key, "| MIN:", [0], ", MAX:", specs[key][1], ", COUNT:", specs[key][2]
    print("\n")

def printSlurmOuts(workDir):
    """
    Go through the repeats directories and print out their content,
    they contain errors that might have occured
    """

    print("\n************************")
    print("ERRORS?\n")

    # Loop through the directories in this VS
    for subDir in os.listdir(workDir):
        # Check only repeat directories
        if subDir.isdigit():
            dirPath = os.path.join(workDir, subDir)

            # Get all the *.ou files in this dir
            slurmOutPaths = glob.glob(dirPath + "/*.out")

            for slurmOutPath in slurmOutPaths:
                slurmOutFile = open(slurmOutPath, "r")
                slurmLines = slurmOutFile.readlines()
                slurmOutFile.close()

                if len(slurmLines) > 0:
                    print(slurmOutPath)
                    for line in slurmLines:
                        print(line)
    print("\n")


if __name__ == "__main__":
    main()
