#!/usr/bin/env python

#----------------------------------------------------
#
#   Execute within a VS directory, will crawl through
#   all its subdirs and submit all .slurm files found
#   there, while pausing for 1 second between each
#   submission
#
#   Thomas Coudrat, February 2014
#
#----------------------------------------------------

import os
import time

workDir = os.getcwd()
slurmPaths = []

# Listing direct subdirectories to the dir where this was executed
for subDir in os.listdir(workDir):
    if subDir.isdigit():
        path = os.path.join(workDir, subDir)

        # If this is a directory, then list the files/directories in it
        if os.path.isdir(path):
            #print path
            files = os.listdir(path)

            # For each of these, save every file that ends with .slurm in
            # a list, by saving its full path
            for file in files:
                if file.endswith(".slurm"):
                    slurmPaths.append(os.path.join(path, file))

print
print "SUBMITTING", str(len(slurmPaths)), "JOBS:"
print

# Loop over the saved .slurm paths
for slurmPath in slurmPaths:

    slurmDir = os.path.dirname(slurmPath)
    slurmFile = os.path.basename(slurmPath)
    #print slurmPath, slurmDir

    # Change to that repeat directory
    os.chdir(slurmDir)
    os.system("sbatch " + slurmFile)
    time.sleep(1)

print
