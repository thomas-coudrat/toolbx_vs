#!/usr/bin/env python

#-------------------------------------------------------------
#
#   Run in a VS repeat directory, checks all .ou files and
#   compiles the number of occurence of the 'SCORE' word in
#   order to inform about the status of that repeat
#
#   Thomas Coudrat, February 2014
#
#-------------------------------------------------------------

import glob
import os

# Setting up variables
workDir = os.getcwd()
ouFiles = glob.glob(workDir + "/*.ou")
scoreCount = 0

# Looping over .ou files in the current dir
for file in ouFiles:

    # Read all lines of that .ou file
    f = open(file, "r")
    lines = f.readlines()
    f.close()

    # Count the number of occurence of the word SCORE in it
    for line in lines:
        if "SCORE" in line:
            scoreCount +=1

print
print "COMPLETED DOCKS:", scoreCount
print "DIRECTORY:", workDir
print
