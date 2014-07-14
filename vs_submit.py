#!/usr/bin/env python

#----------------------------------------------------
#
#   Execute within a VS directory, will crawl through
#   all its subdirs and submit all .slurm or .pbs
#   files found there, while pausing for 1 second
#   between each submission
#
#   Thomas Coudrat, February 2014
#
#----------------------------------------------------

import os
import time
import argparse
import sys
import socket


def main():
    """
    Run the VS submission script
    """

    # Return the queuing system chosen
    vsDir, queue = parsing()

    # Get the current working directory
    cwd = os.getcwd()

    # Store all queueing scripts to be submitted in this directory
    queuePaths = getQueueScripts(vsDir, queue)

    print "\n SUBMITTING", str(len(queuePaths)), "JOBS: \n"

    # Submit all those scripts (using the proper queueing system)
    submitQueueScripts(queuePaths, cwd, queue)

    print


def parsing():
    """
    Define arguments, parse and return them
    """

    # Define and collect arguments
    descr = "Submits a VS using either -slurm or -pbs queuing system"
    descr_vsDir = "VS directory to be submitted to the queue"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("vsDir", help=descr_vsDir)
    args = parser.parse_args()
    vsDir = args.vsDir

    # Get queuing system from hostname
    hostname = socket.gethostname()
    if hostname == "msgln6.its.monash.edu.au":
        queue = ".pbs"
    elif hostname == "barcoo":
        queue = ".slurm"
    else:
        print "Queuing system could not be indentified on your system", hostname
        sys.exit()

    return vsDir, queue


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
            # .pbs in a list, by saving its full path
            for file in files:
                if file.endswith(queue):
                    queuePaths.append(os.path.join(path, file))

    return queuePaths


def submitQueueScripts(queuePaths, cwd, queue):
    """
    Submit all the queueing scripts
    """

    # Loop over the saved .slurm or .pbs paths
    for queuePath in queuePaths:

        # Get the full path relative to the root
        queueFullPath = os.path.join(cwd, queuePath)

        # Get the directory name and the file name separately
        queueDir = os.path.dirname(queueFullPath)
        queueFile = os.path.basename(queuePath)
        #print queuePath, queueDir

        # Change to that repeat directory
        #queuePath = os.path.joindir(cwd, queueDir)
        os.chdir(queueDir)
        # And submit using either SLURM or PBS queueing
        # system depending on what was chosen
        if queue == ".slurm":
            #print "sbatch " + queueFile
            os.system("sbatch " + queueFile)
        elif queue == ".pbs":
            #print "qsub " + queueFile
            os.system("qsub " + queueFile)
        time.sleep(1)


if __name__ == "__main__":
    main()
