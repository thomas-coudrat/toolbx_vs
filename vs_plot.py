#!/usr/bin/env python

from matplotlib import pyplot as plt
import argparse
import os


def main():

    rocPaths, gui = parseArgs()

    plot(rocPaths, gui)


def parseArgs():
    """
    Parsing and returning arguments
    """

    descr = "Feed rocData (however many files), plots ROC curves"
    descr_gui = "Use this flag to display plot: saves to .png by the default"
    descr_rocPaths = "Provide rocDataFiles data1.csv data2.csv data4.csv"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("-gui", action="store_true", help=descr_gui)
    parser.add_argument("rocPaths", help=descr_rocPaths, nargs="+")

    args = parser.parse_args()
    rocPaths = args.rocPaths
    gui = args.gui

    return rocPaths, gui


def plot(rocPaths, gui):
    """
    Plot the data provided as argument, to draw ROC curves
    """

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

    plt.xlabel("Total count (" + str(X[-1]) + ")", fontsize=16)
    plt.ylabel("True positive count (" + str(Y[-1]) + ")", fontsize=16)
    plt.title("Recovery of known actives (ROC curve)")
    plt.legend(loc="lower right", prop={'size': 12})

    if gui:
        plt.show()
    else:
        plt.savefig('rocPlot.png', bbox_inches="tight")


if __name__ == "__main__":
    main()
