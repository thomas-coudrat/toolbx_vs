#!/usr/bin/env python

#--------------------------------------------------------
#
# Create .slurm files to slice a VS into several equal
# portions, for parallelisation
#
# Thomas Coudrat, February 2014
#
#--------------------------------------------------------

import os
import argparse
import glob
import shutil
import re


def main():
    # Getting all the args
    args = parsing()

    # Library params
    libSize = int(args.libSize)
    sliceSize = int(args.sliceSize)
    repeatNum = int(args.repeatNum)
    # VS params
    walltime = args.walltime
    thor = args.thor
    # Project info
    setupDir = args.setupDir
    projName = glob.glob(setupDir + "/*.dtb")[0].replace(".dtb", "").split("/")[1]

    # Store all the report lines in this file, will be used for printout
    # and to write to file
    reportLines = []
    # Get current working directory
    workDir = os.getcwd()

    reportLines.append("\nPARAMETERS:\n")
    reportLines.append("\t libSize: " + str(libSize))
    reportLines.append("\t sliceSize: " + str(sliceSize))
    reportLines.append("\t repeatNum: " + str(repeatNum))
    reportLines.append("\t walltime: " + walltime)
    reportLines.append("\t thoroughness: " + thor)
    reportLines.append("\t setupDir: " + setupDir)
    reportLines.append("\t projName: " + projName)
    reportLines.append("\n")

    # grep the parameters to lookout for in the .dtb file, and print them out
    reportLines = printParams(setupDir, reportLines)

    reportLines.append("\n***********************\n")

    # Creating the repeats directories, which are copies of the setupDir
    reportLines = createRepeats(repeatNum, setupDir, reportLines)

    reportLines.append("\n***********************\n")

    # Create the .slurm slices
    reportLines = createSlices(libSize, sliceSize, repeatNum, walltime, thor, projName, reportLines)

    reportLines.append("\n")

    # Print and write report
    printWriteReport(reportLines, workDir, projName)


def parsing():
    # Parsing descriptions for arguments
    descr = "Slices a VS scripts into portions of library for slurm submission"
    descr_libSize = "Size of the full library"
    descr_sliceSize = "Size of the slices"
    descr_repeatNum = "Number of repeats"
    descr_walltime = "Walltime for a single slice (format: 1-24:00:00)"
    descr_thor = "Thoroughness of the docking (format: 5.)"
    descr_setupDir = "Name of the directory containing setup files"

    # Defining the arguments
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("libSize", help=descr_libSize)
    parser.add_argument("sliceSize", help=descr_sliceSize)
    parser.add_argument("repeatNum", help=descr_repeatNum)
    parser.add_argument("walltime", help=descr_walltime)
    parser.add_argument("thor", help=descr_thor)
    parser.add_argument("setupDir", help=descr_setupDir)

    # Parsing and storing into variables
    args = parser.parse_args()

    return args


def printParams(setupDir, reportLines):
    """
    Use the grep command to print out common parameters to check
    when running a VS
    """
    dtbPath = glob.glob(setupDir + "/*.dtb")[0]

    dtbFile = open(dtbPath, "r")
    dtbLines = dtbFile.readlines()
    dtbFile.close()

    regEx = "maxHdonors|maxLigSize|maxNO|maxTorsion|ringFlexLevel|sampleRacemic"\
            "|scoreThreshold|maxPk|minPk|chargeGroups|dbIndex|dbType"

    lineNumber=0
    for line in dtbLines:
        if re.search(regEx, line):
            param = dtbLines[lineNumber].strip()
            val = dtbLines[lineNumber + 1].strip()
            reportLines.append("\t" + param + " : " + val)
        lineNumber +=1
    reportLines.append("\n")

    return reportLines


def createRepeats(repeatNum, setupDir, reportLines):

    # Get files in setupDir
    filePaths = glob.glob(setupDir + "/*")

    # Create each repeat dirs, and populate them with files
    for repeatDir in range(1, repeatNum + 1):

        repeatDir = str(repeatDir)
        # Create the directory if it doesn't already exist
        if not os.path.exists(repeatDir):
            os.makedirs(repeatDir)

        reportLines.append("\n")
        reportLines.append("REPEAT:" + repeatDir + "\n")

        # Copy each file into this new directory
        for filePath in filePaths:
            fileName = os.path.basename(filePath)
            shutil.copy(filePath, repeatDir + "/" + fileName)
            reportLines.append("\t COPYING:" + fileName)

    return reportLines


def createSlices(libSize, sliceSize, repeatNum, walltime, thor, projName, reportLines):
    repeat = 1
    # Loop over repeat directories
    while repeat <= repeatNum:

        # Initialize variables for this repeat

        repeatDir = os.getcwd() + "/" + str(repeat) + "/"
        upperLimit = 0
        sliceCount = 0
        exit = False

        reportLines.append("\n")
        reportLines.append("REPEAT:" + repeatDir + "\n")

        # Loop over the slices
        while True:

            # Update upperLimit and sliceCount
            lowerLimit = upperLimit + 1
            upperLimit += sliceSize
            sliceCount += 1

            # Exit statement of the loop, when
            if upperLimit > libSize:
                upperLimit = libSize
                exit = True

            # Create sliceName for job name and slurm file name
            sliceName = projName + "_rep" + str(repeat) + "_sl" + str(upperLimit)

            # Content of the slurm file to be created
            lines = []
            lines.append("#!/bin/bash")
            lines.append("#SBATCH -p main")
            lines.append("#SBATCH --ntasks=1")
            lines.append("#SBATCH --mem-per-cpu=1024")
            lines.append("#SBATCH --time=" + walltime)
            lines.append("#SBATCH --job-name=" + sliceName)
            lines.append("")
            lines.append("ICMHOME=/vlsci/VR0024/tcoudrat/bin/icm-3.7-3b")
            lines.append("$ICMHOME/icm64 -vlscluster $ICMHOME/_dockScan " + projName +
                " thorough=" + thor +
                " from=" + str(lowerLimit) +
                " to=" + str(upperLimit) +
                " >& " + projName + "_" + str(upperLimit) + ".ou")

            # Create the slurm file in current repeat directory
            slurmFile = open(repeatDir + sliceName + ".slurm", "w")
            for line in lines:
                slurmFile.write(line + "\n")
            slurmFile.close()

            reportLines.append("\t SLICE:" + sliceName + ".slurm")

            # Exit statement when end of the library is reached
            if exit:
                break

        # Update the repeat number
        repeat += 1

    return reportLines


def printWriteReport(reportLines, workDir, projName):
    """
    Go through the report lines and print them to standard output and write them to
    a text file for future ref
    """

    reportFile = open(workDir + "/" + projName + ".log", "w")

    for line in reportLines:
        print line
        reportFile.write(line + "\n")

    reportFile.close()

if  __name__ == "__main__":
    main()
