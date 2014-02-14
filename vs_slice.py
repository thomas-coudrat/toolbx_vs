#!/vlsci/VR0024/tcoudrat/bin/python275/bin/python

#--------------------------------------------------------
#
# Create .slurm files to slice a VS into several equal
# portions, for parallelisation
#
# Thomas Coudrat, February 2014
#
#--------------------------------------------------------

import os

# Ligand library, slice and repeats
libSize = 17720
sliceSize = 3000
repeatNum = 3

# VS parameters
walltime = "0-24:00:00"
thor = "5."

# Project info
projName = "ADORA2A_2YDV_allLoop"

repeat = 1
# Loop over repeat directories
while repeat <= repeatNum:

    # Initialize variables for this repeat

    repeatDir = os.getcwd() + "/" + str(repeat) + "/"
    if not os.path.exists(repeatDir):
        os.makedirs(repeatDir)
    upperLimit = 0
    sliceCount = 0
    exit = False

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

        # Exit statement when end of the library is reached
        if exit:
            break

    # Update the repeat number
    repeat += 1
