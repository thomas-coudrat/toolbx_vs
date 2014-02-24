#!/usr/bin/env python

from matplotlib import pyplot as plt
import argparse
import os


def main():

    gui, title, rocPaths = parseArgs()

    # Plot the values 0 and 1, which correspond to X and Y numbers
    plot(title, rocPaths, gui, "number")

    # Plot the values 2 and 3, which correspond to the percentage X and Y
    plot(title, rocPaths, gui, "percent")


def parseArgs():
    """
    Parsing and returning arguments
    """

    descr = "Feed rocData (however many files), plots ROC curves"
    descr_gui = "Use this flag to display plot: saves to .png by the default"
    descr_title = "Provide a title for the graph, also used as filename"
    descr_rocPaths = "Provide rocDataFiles data1.csv data2.csv data4.csv"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("-gui", action="store_true", help=descr_gui)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("rocPaths", help=descr_rocPaths, nargs="+")

    args = parser.parse_args()
    gui = args.gui
    title = args.title
    rocPaths = args.rocPaths

    return gui, title, rocPaths


def plot(title, rocPaths, gui, mode):
    """
    Plot the data provided as argument, to draw ROC curves
    """

    plt.figure()


    for rocPath in rocPaths:
        rocFile = open(rocPath, "r")
        rocLines = rocFile.readlines()
        rocFile.close()

        X = []
        Y = []
        rocName = os.path.basename(rocPath).replace(".csv", "")
        for line in rocLines:
            ll = line.split(",")
            # Pick the numbers corresponding to the mode selected
            if mode == "number":
                X.append(int(ll[0]))
                Y.append(int(ll[1]))
            elif mode == "percent":
                X.append(float(ll[2]))
                Y.append(float(ll[3]))

        plt.plot(X, Y, label=rocName)

    plt.xlabel("Total count (" + str(X[-1]) + ")", fontsize=16)
    plt.ylabel("True positive count (" + str(Y[-1]) + ")", fontsize=16)
    plt.title(title)
    plt.legend(loc="lower right", prop={'size': 12})

    if gui:
        plt.show()
    else:
        fileName = title.replace(" ", "_") + "_" + mode + ".png"
        plt.savefig(fileName, bbox_inches="tight")


if __name__ == "__main__":
    main()
