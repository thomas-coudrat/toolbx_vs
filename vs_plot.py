#!/usr/bin/env python

from matplotlib import pyplot as plt
import argparse
import os

descr = "Feed rocData (however many files), plots ROC curves"
descr_rocPaths = "Provide rocDataFiles: data1.csv,data2.csv,data4.csv"
parser = argparse.ArgumentParser(description=descr)
parser.add_argument("rocPaths", help=descr_rocPaths)
args = parser.parse_args()
rocPaths = args.rocPaths.split(",")

rocList = []

for rocPath in rocPaths:
    rocFile = open(rocPath, "r")
    rocLines = rocFile.readlines()
    rocFile.close()

    X = []
    Y = []
    rocName = os.path.basename(rocPath).replace(".csv", "")
    for line in rocLines:
        ll = line.split(",")
        X.append(int(ll[0]))
        Y.append(int(ll[1]))

    plt.plot(X, Y, label=rocName)

plt.legend()
plt.show()
