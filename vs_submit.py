#!/usr/bin/env python

# Execute within a VS directory, will crawl through
# all its subdirs and submit all .slurm or .sge
# files found there, while pausing for 1 second
# between each submission
#
# https://github.com/thomas-coudrat/toolbx_vs
# Thomas Coudrat <thomas.coudrat@gmail.com>

import os
import time
import argparse
import sys
import socket
import json

def main():
    """
    Run script
    """

    # Return the queuing system chosen
    vsDir, queue = parsing()

    # Get the current working directory
    cwd = os.getcwd()

    # Store all queueing scripts to be submitted in this directory
    queuePaths = getQueueScripts(vsDir, queue)

    # Ask for confirmation to submit run
    confirmSubmit(queuePaths)

    # Submit all those scripts (using the proper queueing system)
    submitQueueScripts(queuePaths, cwd, queue)

    print("")


def parsing():
    """
    Define arguments, parse and return them
    """

    # Define and collect arguments
    descr = "Submits a VS using either -slurm or -sge queuing system"
    descr_vsDir = "VS directory to be submitted to the queue"
    descr_queue = "Queuing system to be used (sge/slurm)"

    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("vsDir", help=descr_vsDir)
    parser.add_argument("queue", help=descr_queue)

    args = parser.parse_args()

    vsDir = args.vsDir
    queue = args.queue

    if queue not in ("sge", "slurm"):
        print("Only 'sge' and 'slurm' are accepted queuing system options")
        sys.exit()

    return vsDir, queue


def confirmSubmit(queuePaths):
    """
    Ask for user input confirmation to submit the jobs, listing information
    on what is about to be submitted
    """

    print("\nYou are about to submit " + str(len(queuePaths)) + " jobs.")

    answer = input("Do you want to proceed? (yes/no) ")

    if answer == "yes":
        print("\nSubmitting jobs...\n")
    else:
        print("Quit...\n")
        sys.exit()


def getQueueScripts(vsDir, queue):
    """
    Make a list of the scripts to be submited
    """

    queuePaths = []
    # Listing direct subdirectories to the dir where this was executed
    for subDir in os.listdir(vsDir):
        if subDir.isdigit():
            path = os.path.join(vsDir, subDir)
            files = os.listdir(path)

            # For each of these, save every file that ends with .slurm or
            # .sge in a list, by saving its full path
            for file in files:
                if file.endswith("." + queue):
                    queuePaths.append(os.path.join(path, file))

    return queuePaths


def submitQueueScripts(queuePaths, cwd, queue):
    """
    Submit all the queueing scripts
    """

    # Loop over the saved .slurm or .sge paths
    for queuePath in queuePaths:

        # Get the full path relative to the root
        queueFullPath = os.path.join(cwd, queuePath)

        # Get the directory name and the file name separately
        queueDir = os.path.dirname(queueFullPath)
        queueFile = os.path.basename(queuePath)
        # print queuePath, queueDir

        # Change to that repeat directory
        # queuePath = os.path.joindir(cwd, queueDir)
        os.chdir(queueDir)
        # And submit using either SLURM or SGE queueing
        # system depending on what was chosen
        if queue == "slurm":
            # print "sbatch " + queueFile
            os.system("sbatch " + queueFile)
        elif queue == "sge":
            # print "qsub " + queueFile
            os.system("qsub " + queueFile)
        time.sleep(1)


if __name__ == "__main__":
    main()
