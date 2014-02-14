#!/vlsci/VR0024/tcoudrat/bin/python275/bin/python

#--------------------------------------------------------------
#
#   Run on output (.ou) from VS in order to obtain the number
#   and type of skipped ligands during the VS. Can be ran while
#   VS is still running
#
#   Thomas Coudrat, February 2014
#
#--------------------------------------------------------------

import glob
import os

# Setting up variables
workDir = os.getcwd()
ouFiles = glob.glob(workDir + "/*.ou")

# Loop over all .ou files and store them all in allLines
allLines = []
print
print "CHECKING FILES:\n"
for ouFile in ouFiles:
    print ouFile
    file = open(ouFile, "r")
    lines = file.readlines()
    file.close()

    allLines = allLines + lines
print

# Dictionary containing specs as keys, and [MIN,MAX,COUNT] as values
specs = {}
skipCount = 0
# Looping over all stored lines
for line in allLines:

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

# Print out the result
keys = specs.keys()
print "THERE WERE", skipCount, "LIGANDS SKIPPED:"
print
for key in keys:
    print key, "\t\tMIN:", specs[key][0], "\t\tMAX:", specs[key][1], "\t\tCOUNT:", specs[key][2]
print
