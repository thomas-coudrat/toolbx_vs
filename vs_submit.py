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


def main():
    """
    Run the VS submission script
    """

    # Return the queuing system chosen
    queue = parsing()

    workDir = os.getcwd()
    queuePaths = []

    # Listing direct subdirectories to the dir where this was executed
    for subDir in os.listdir(workDir):
        if subDir.isdigit():
            path = os.path.join(workDir, subDir)

            # If this is a directory, then list the files/directories in it
            if os.path.isdir(path):
                #print path
                files = os.listdir(path)

                # For each of these, save every file that ends with .slurm or
                # .pbs in a list, by saving its full path
                for file in files:
                    if file.endswith(queue):
                        queuePaths.append(os.path.join(path, file))

    print
    print "SUBMITTING", str(len(queuePaths)), "JOBS:"
    print

    # Loop over the saved .slurm or .pbs paths
    for queuePath in queuePaths:

        queueDir = os.path.dirname(queuePath)
        queueFile = os.path.basename(queuePath)
        #print queuePath, queueDir

        # Change to that repeat directory
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

    print


def parsing():
    """
    Define arguments, parse and return them
    """
    descr = "Submits a VS using either -slurm or -pbs queuing system"
    descr_pbs = "Optional argument, use the PBS queuing system"
    descr_slurm = "Optional argument, use the SLURM queuing system"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("-pbs", action="store_true", help=descr_pbs)
    parser.add_argument("-slurm", action="store_true", help=descr_slurm)
    args = parser.parse_args()

    pbs = args.pbs
    slurm = args.slurm

    # Check if a queuing system was chosen
    if not pbs and not slurm:
        print "A queuing system -pbs or -slurm must be chosen"
        sys.exit()
    elif pbs and slurm:
        print "Only one of -pbs or -slurm option must be chosen"
        sys.exit()
    elif pbs:
        queue = ".pbs"
    elif slurm:
        queue = ".slurm"

    return queue


if __name__ == "__main__":
    main()
