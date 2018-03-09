#!/usr/bin/env python

# Builds the files to split a VS into separate slices to be ran in parallel on
# an HPC cluster using the SLURM or PBS queuing system.
#
# https://github.com/thomas-coudrat/toolbx_vs
# Thomas Coudrat <thomas.coudrat@gmail.com>

import os
import argparse
import glob
import shutil
import re
import sys
import socket
import json
import datetime
import time

def main():
    """
    Run the following script
    """

    # Getting all the args
    libStart, libEnd, sliceSize, repeatNum, thor, \
        walltime, setupDir, projName, queue = parsing()

    # Get the path from the Json file
    icmHome = getPath()

    # Store all the report lines in this file, will be used for printout
    # and to write to file
    reportLines = []

    # Get current working directory
    workDir = os.getcwd()

    # Clean files present in the current repeat directories, if any
    for repeatDir in glob.glob(workDir + "/[0-9]*"):
        cleanRepeatDir(repeatDir)

    reportLines.append("\nPARAMETERS:\n")
    reportLines.append("\t libStart: " + str(libStart))
    reportLines.append("\t libEnd: " + str(libEnd))
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
    reportLines = createSlices(libStart, libEnd, sliceSize, walltime, thor,
                               projName, repeatNum, queue, reportLines, icmHome)

    reportLines.append("\n")

    # Print and write report
    printWriteReport(reportLines, workDir, projName)


def parsing():
    """
    Defining the arguments and parsing them
    """

    # Parsing descriptions for arguments
    descr = "Slices a VS scripts into portions of library for slurm submission"
    descr_libStart = "Ligand library ID where to START the VS"
    descr_libEnd = "Ligand library ID where to END the VS"
    descr_sliceSize = "Size of the slices"
    descr_repeatNum = "Number of repeats"
    descr_thor = "Thoroughness of the docking (format: 5.)"
    descr_walltime = "Walltime for a single slice (format: 1-24:00:00)"
    descr_setupDir = "Name of the directory containing setup files"
    descr_queue = "Queuing system to be used (sge/slurm/slurm-srun)"

    # Defining the arguments
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("libStart", help=descr_libStart)
    parser.add_argument("libEnd", help=descr_libEnd)
    parser.add_argument("sliceSize", help=descr_sliceSize)
    parser.add_argument("repeatNum", help=descr_repeatNum)
    parser.add_argument("thor", help=descr_thor)
    parser.add_argument("walltime", help=descr_walltime)
    parser.add_argument("setupDir", help=descr_setupDir)
    parser.add_argument("queue", help=descr_queue)

    # Parsing and storing into variables
    args = parser.parse_args()

    # Library params
    libStart = int(args.libStart)
    libEnd = int(args.libEnd)
    sliceSize = int(args.sliceSize)
    repeatNum = int(args.repeatNum)
    thor = args.thor
    # VS params
    walltime = args.walltime
    queue = args.queue
    # Project info
    setupDir = args.setupDir
    dtbFileName = glob.glob(setupDir + "/*.dtb")[0]
    projName = dtbFileName.replace(".dtb", "").split("/")[1]

    if queue not in ("sge", "slurm", "slurm-srun"):
        print("'sge', 'slurm' and 'slurm-srun' are the queuing system options")
        sys.exit()

    return libStart, libEnd, sliceSize, repeatNum, thor, walltime, setupDir, \
        projName, queue


def getPath():
    """
    Get the ICM executable from the ICMHOME environment variable
    """

    # Get environment variable. It returns None if not defined on the system.
    icmHome = os.environ.get('ICMHOME')

    # Return path to executable if the environment variable was found
    if icmHome == None:
        "The ICMHOME environment variable must be set for your system. Exiting."
        sys.exit()
    else:
        return icmHome


def cleanRepeatDir(repeatDir):
    """
    Checks for file already present in this VS directory, and prompts
    the user to confirm removal, before creating a new set of files
    """

    # Get files in workDir
    allPaths = glob.glob(repeatDir + "/*")
    filePaths = [f for f in allPaths if not os.path.basename(f).startswith("backup_")]

    # If files were found in this directory, display what they are and prompt
    # for deletion
    if len(filePaths) > 0:

        # Displaying existing files
        print("\n#### VS files were found in this directory ####\n")
        print(repeatDir)
        for filePath in filePaths:
            print("\t" + os.path.basename(filePath))
        print("\n")

        # Prompting for removal
        answer = input('Three options:\n' +
                       '(abort) operation,\n' +
                       '(delete) all existing files and continue,\n' +
                       '(keep) existing files and continue (note that ' +
                       'you risk overwriting some files)\n' +
                       '(backup) existing file in a timestamped folder ' +
                       'per repeat dir\n\n'
                       'abort/delete/keep/backup? (default=abort) ')
        if answer == "delete":
            print("DELETING PREVIOUS FILES...")
            for filePath in filePaths:
                os.remove(filePath)
        elif answer == "keep":
            print("CONTINUE WITHOUT DELETING FILES...")
        elif answer == 'backup':
            print("BACKING UP FILES...")
            t = time.time()
            humanTime = datetime.datetime.fromtimestamp(int(t)).strftime('%Y-%m-%d_%H:%M:%S')
            backupDir = repeatDir + "/backup_" + humanTime
            os.makedirs(backupDir)
            for filePath in filePaths:
                shutil.move(filePath, backupDir)
        else:
            print("NOT DELETING, QUIT...")
            sys.exit()


def printParams(setupDir, reportLines):
    """
    Use the grep command to print out common parameters to check
    when running a VS
    """
    dtbPath = glob.glob(setupDir + "/*.dtb")[0]

    dtbFile = open(dtbPath, "r")
    dtbLines = dtbFile.readlines()
    dtbFile.close()

    regEx = "maxHdonors|maxLigSize|maxNO|maxTorsion|ringFlexLevel|" \
            "sampleRacemic|scoreThreshold|maxPk|minPk|chargeGroups|" \
            "dbIndex|dbType"

    lineNumber = 0
    for line in dtbLines:
        if re.search(regEx, line):
            param = dtbLines[lineNumber].strip()
            val = dtbLines[lineNumber + 1].strip()
            reportLines.append("\t" + param + " : " + val)
        lineNumber += 1
    reportLines.append("\n")

    return reportLines


def createRepeats(repeatNum, setupDir, reportLines):
    """
    Copy the content of the setup directory to however many
    repeat directories wanted by the user
    """

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


def createSlices(libStart, libEnd, sliceSize, walltime, thor, projName,
                 repeatNum, queue, reportLines, icmHome):
    """
    Create the .slurm slices to split the VS job into portions for submission
    to the cluster
    """

    repeat = 1

    # Loop over repeat directories
    while repeat <= repeatNum:

        # Initialize variables for this repeat
        cwd = os.getcwd()
        repeatDir = cwd + "/" + str(repeat) + "/"
        # Update the report
        reportLines.append("\n")
        reportLines.append("REPEAT:" + repeatDir + "\n")

        # Initialize variables for the first slice
        lowerLimit = libStart
        upperLimit = libStart + sliceSize - 1
        sliceCount = 1
        keepLooping = True

        # Loop over the slices
        while keepLooping:

            # Exit statement of the loop, when upperLimit has reached the
            # size of the ligand library
            if upperLimit >= libEnd:
                upperLimit = libEnd
                # Once the end of the libSize has been reached, stop the next
                # loop
                keepLooping = False

            # Create sliceName for job name and slurm file name
            sliceName = projName + "_rep" + str(repeat) + \
                "_sl" + str(upperLimit)

            # Create a slice, check for submission system, run the appropriate
            # command
            if queue == "slurm-srun":
                reportLines = slurmSrunSlice(sliceCount, projName, thor,
                                             lowerLimit, upperLimit,
                                             libStart, libEnd,
                                             repeatDir, reportLines, icmHome)
            elif queue == "sge":
                reportLines = sgeSlice(walltime, sliceName, projName, thor,
                                       lowerLimit, upperLimit, repeatDir,
                                       reportLines, icmHome)
            elif queue == "slurm":
                reportLines = slurmSlice(walltime, sliceName, projName, thor,
                                         lowerLimit, upperLimit, repeatDir,
                                         reportLines, icmHome)

            # Update upperLimit and sliceCount
            lowerLimit += sliceSize
            upperLimit += sliceSize
            sliceCount += 1

        # Update the repeat number
        repeat += 1

        # Combine these slices in a call srun
        if queue == "slurm-srun":
            slurmSrun(projName, libStart, libEnd, walltime,
                      repeatDir, repeat, sliceCount - 1)

    return reportLines


def slurmSrun(projName, libStart, libEnd,  walltime, repeatDir, repeat, sliceCount):
    """
    Create the srun SLURM script which will group all SLURM submissions together
    """

    slurmName = projName + "_" + str(repeat)

    libRange = str(libStart) + "-" + str(libEnd)

    lines = []
    lines.append("#!/bin/bash")
    lines.append("#SBATCH -p main")
    lines.append("#SBATCH --ntasks=" + str(sliceCount))
    lines.append("#SBATCH --mem-per-cpu=1024")
    lines.append("#SBATCH --time=" + walltime)
    lines.append("#SBATCH --job-name=" + slurmName)
    lines.append("")
    lines.append("for i in `seq 1 $SLURM_NTASKS`")
    lines.append("do")
    lines.append('\tsrun --nodes=1 --ntasks=1 --cpus-per-task=1 ' +
                 'sh -c "(sh slice_' + libRange + '_$i.sh)" &')
    lines.append("done")
    lines.append("wait")

    with open(repeatDir + "srun_" + libRange  + ".slurm", "w") as f:
        f.write("\n".join(lines))


def slurmSrunSlice(sliceCount, projName, thor, lowerLimit, upperLimit,
                   libStart, libEnd, repeatDir, reportLines, icmHome):
    """
    Create a slurm slice that will be used as part of a bundled SRUN command
    and write to a file with the info provided
    """

    lines = []
    lines.append("#!/bin/bash")
    lines.append("")
    lines.append("ICMHOME=" + icmHome)
    lines.append("$ICMHOME/icm64 -vlscluster $ICMHOME/_dockScan " + projName +
                 " thorough=" + thor +
                 " from=" + str(lowerLimit) +
                 " to=" + str(upperLimit) +
                 " >& " + projName + "_" + str(upperLimit) + ".ou")

    # WRITE SLURM LINES TO FILE
    sliceName = str(libStart) + "-" + str(libEnd) + "_" + str(sliceCount)
    with open(repeatDir + "slice_" + sliceName + ".sh", "w") as f:
        f.write("\n".join(lines))

    # Update report
    reportLines.append("\tproject: " + projName +
                       ", repeat:" + os.path.relpath(repeatDir) +
                       ", slice:" + str(sliceCount))

    return reportLines


def slurmSlice(walltime, sliceName, projName, thor, lowerLimit, upperLimit,
               repeatDir, reportLines, icmHome):
    """
    Create a slurm slice and write to a file with the info provided
    """

    lines = []
    lines.append("#!/bin/bash")
    lines.append("#SBATCH --mem=1024")
    lines.append("#SBATCH --time=" + walltime)
    lines.append("#SBATCH --job-name=" + sliceName)
    lines.append("")
    lines.append("ICMHOME=" + icmHome)
    lines.append("$ICMHOME/icm64 -vlscluster $ICMHOME/_dockScan " + projName +
                 " thorough=" + thor +
                 " from=" + str(lowerLimit) +
                 " to=" + str(upperLimit) +
                 " >& " + projName + "_" + str(upperLimit) + ".ou")

    # WRITE SLURM LINES TO FILE
    with open(repeatDir + sliceName + ".slurm", "w") as f:
        f.write('\n'.join(lines))

    # Update report
    reportLines.append("\tproject: " + projName +
                       ", repeat:" + os.path.relpath(repeatDir) +
                       ", slice:" + sliceName)

    return reportLines


def sgeSlice(walltime, sliceName, projName, thor, lowerLimit, upperLimit,
             repeatDir, reportLines, icmHome):
    """
    Create a SGE slice given the info provided
    """
    lines = []
    lines.append("#!/bin/sh")
    lines.append("#$ -S /bin/sh")
    lines.append("#$ -l h_rt=" + walltime)
    lines.append("#$ -l h_vmem=1G")
    lines.append("#$ -q hqu9")
    lines.append("#$ -l dpod=1")
    lines.append("#$ -cwd")
    lines.append("#$ -N " + str(sliceName))
    lines.append("")
    lines.append("ICMHOME=" + icmHome)
    lines.append("$ICMHOME/icm64 -vlscluster $ICMHOME/_dockScan " + projName +
                 " thorough=" + thor +
                 " from=" + str(lowerLimit) +
                 " to=" + str(upperLimit) +
                 " >& " + projName + "_" + str(upperLimit) + ".ou")

    # WRITE SLURM LINES TO FILE
    with open(repeatDir + sliceName + ".sge", "w") as f:
        f.write("\n".join(lines))

    # Update report
    reportLines.append("\t SLICE:" + sliceName + ".sge")

    return reportLines


def printWriteReport(reportLines, workDir, projName):
    """
    Go through the report lines and print them to standard output and
    write them to a text file for future ref
    """

    reportFile = open(workDir + "/" + projName + ".log", "w")

    for line in reportLines:
        print(line)
        reportFile.write(line + "\n")

    reportFile.close()


if __name__ == "__main__":
    main()
