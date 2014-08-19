#!/usr/bin/env python

# ----------------------------------------------------
#
#   Execute within a VS directory, will crawl through
#   all its subdirs and submit all .slurm or .pbs
#   files found there, while pausing for 1 second
#   between each submission
#
#   Thomas Coudrat, February 2014
#
# ----------------------------------------------------

import os
import time
import argparse
import sys
import socket
import json


def main():
    """
    Run the VS submission script
    """

    # Return the queuing system chosen
    vsDir = parsing()

    # Get the queuing system on this cluster
    queue = getQueuingSys()

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

    return vsDir


def getQueuingSys():
    """
    Figure out which queuing system to use depending on the platform this script
    is executed on
    """

    # This Json file stores the ICM executable locations for each platform
    queuesJson = os.path.dirname(os.path.realpath(__file__)) + "/queue_sys.json"

    # Read content of .json file
    with open(queuesJson, "r") as jsonFile:
        queues = json.load(jsonFile)

    # Get the hostname to know which computer this is executed on
    hostname = socket.gethostname()

    # Assign the ICM executable path corresponding to the hostname, if it is not
    # defined then stop the execution
    if hostname in queues.keys():
        queue = queues[hostname]
    else:
        print("The queuing system could not be assigned", hostname)
        sys.exit()

    return queue


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
                if file.endswith("." + queue):
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
        # print queuePath, queueDir

        # Change to that repeat directory
        # queuePath = os.path.joindir(cwd, queueDir)
        os.chdir(queueDir)
        # And submit using either SLURM or PBS queueing
        # system depending on what was chosen
        if queue == "slurm":
            # print "sbatch " + queueFile
            os.system("sbatch " + queueFile)
        elif queue == "pbs":
            # print "qsub " + queueFile
            os.system("qsub " + queueFile)
        time.sleep(1)


if __name__ == "__main__":
    main()
