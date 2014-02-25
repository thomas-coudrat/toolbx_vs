#!/usr/bin/env python

from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import argparse
import os


def main():

    gui, log, title, rocPaths = parseArgs()

    # Plot the values 0 and 1, which correspond to X and Y numbers
    #plot(title, rocPaths, gui, log, "number")

    # Plot the values 2 and 3, which correspond to the percentage X and Y
    plot(title, rocPaths, gui, log, "percent")


def parseArgs():
    """
    Parsing and returning arguments
    """

    descr = "Feed rocData (however many files), plots ROC curves"
    descr_gui = "Use this flag to display plot: saves to .png by the default"
    descr_log = "Draw this plot on a log scale for the X axis"
    descr_title = "Provide a title for the graph, also used as filename"
    descr_rocPaths = "Provide rocDataFiles data1.csv data2.csv data4.csv"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("-gui", action="store_true", help=descr_gui)
    parser.add_argument("-log", action="store_true", help=descr_log)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("rocPaths", help=descr_rocPaths, nargs="+")

    args = parser.parse_args()
    gui = args.gui
    log = args.log
    title = args.title
    rocPaths = args.rocPaths

    return gui, log, title, rocPaths


def plot(title, rocPaths, gui, log, mode):
    """
    Plot the data provided as argument, to draw ROC curves
    """

    fig = plt.figure(figsize=(13, 12), dpi=100)
    ax = fig.add_subplot(111)
    colors = ["#E82F3B", "#3340FF", "#2E2E33"]

    for i, rocPath in enumerate(rocPaths):
        # Read the ROC data file
        rocFile = open(rocPath, "r")
        rocLines = rocFile.readlines()
        rocFile.close()

        # Get the values total ligands in the library, and total known ligands
        # from the last line of the ROC data
        totalLib = int(rocLines[-1].split(",")[0])
        totalKnown = int(rocLines[-1].split(",")[1])
        # Define increment to plot an average curve
        randomIncr = totalKnown / float(totalLib)

        # Write this curve based on the data contained in the ROC data file
        X = []
        Y = []
        random = []
        randomVal = 0
        rocName = os.path.basename(rocPath).replace(".csv", "")
        for line in rocLines:
            ll = line.split(",")
            # Pick the numbers corresponding to the mode selected
            if mode == "number":
                # Create the data curve
                X.append(int(ll[0]))
                Y.append(int(ll[1]))
            elif mode == "percent":
                # Create the data curve
                X.append(float(ll[2]))
                Y.append(float(ll[3]))
                # Create the random curve
                #randomVal += 1
                #random.append(randomVal)
                #print randomVal

            # Create the random curve
            randomVal += randomIncr
            random.append(randomVal)

        plt.plot(X, Y, label=rocName.replace("_", " "),
                 linewidth=2, color=colors[i])
        plt.plot(X, X, "--", color="grey")

    plt.xlabel("Total count (" + str(totalLib) + ")", fontsize=16)
    plt.ylabel("True positive count (" + str(totalKnown) + ")", fontsize=16)
    plt.minorticks_on()
    if log:
        plt.xscale("log")
        ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
    #plt.xticks([0.1, 1, 10, 100])
    plt.title(title, fontsize=18)
    plt.legend(loc="upper left", prop={'size': 12})
    plt.axis('tight')

    if gui:
        plt.show()
    else:
        fileName = title.replace(" ", "_") + "_" + mode + ".png"
        plt.savefig(fileName, bbox_inches="tight")


if __name__ == "__main__":
    main()
