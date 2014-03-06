#!/usr/bin/env python

from matplotlib import pyplot as plt
import argparse


def main():

    gui, log, title, rocData = parseArgs()

    # Extrac the ROC paths and legends from the rocData
    rocPaths = []
    rocLegends = []
    i = 0
    while i < len(rocData):
        rocLegends.append(rocData[i])
        rocPaths.append(rocData[i + 1])
        i += 2

    # Plot the values 0 and 1, which correspond to X and Y numbers
    #plot(title, rocPaths, rocLegends, gui, log, "number")

    # Plot the values 2 and 3, which correspond to the percentage X and Y
    plot(title, rocPaths, rocLegends, gui, log, "percent")


def parseArgs():
    """
    Parsing and returning arguments
    """

    descr = "Feed rocData (however many files), plots ROC curves"
    descr_gui = "Use this flag to display plot: saves to .png by the default"
    descr_log = "Draw this plot on a log scale for the X axis"
    descr_title = "Provide a title for the graph, also used as filename"
    descr_rocData = "Provide rocDataFiles and legend titles for each curve:" \
        " 'legend1!' data1.csv 'legend2?' data2.csv 'legend4!!' data4.csv"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("-gui", action="store_true", help=descr_gui)
    parser.add_argument("-log", action="store_true", help=descr_log)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("rocData", help=descr_rocData, nargs="+")

    args = parser.parse_args()
    gui = args.gui
    log = args.log
    title = args.title
    rocData = args.rocData

    return gui, log, title, rocData


def plot(title, rocPaths, rocLegends, gui, log, mode):
    """
    Plot the data provided as argument, to draw ROC curves
    """

    fig = plt.figure(figsize=(13, 12), dpi=100)
    ax = fig.add_subplot(111)
    colors = ["#E82F3B", "#3340FF", "#2E2E33"]

    for i, (rocPath, rocLegend) in enumerate(zip(rocPaths, rocLegends)):
        # Read the ROC data file
        rocFile = open(rocPath, "r")
        rocLines = rocFile.readlines()
        rocFile.close()

        # Get the values total ligands in the library, and total known ligands
        # from the last line of the ROC data
        totalLib = int(rocLines[-1].split(",")[0])
        totalKnown = int(rocLines[-1].split(",")[1])
        # Define increment to plot an average curve
        #perfectIncr = totalKnown / float(totalLib)

        # Write this curve based on the data contained in the ROC data file
        X = []
        Y = []
        perfect = []
        val = 0
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

            # Create the perfect curve
            if val < totalKnown:
                val += 1
            perfect.append((val * 100.0) / totalKnown)

        # Plot this curve
        plt.plot(X, Y, label=rocLegend,
                 linewidth=2, color=colors[i])

    # Now plot random and perfect curves, common for all plotted curves
    plt.plot(X, X, "--", color="grey")
    plt.plot(X, perfect, color="grey")

    plt.xlabel("% of ranked database (total=" + str(totalLib) + ")",
               fontsize=16)
    plt.ylabel("% of known ligands found (total=" + str(totalKnown) + ")",
               fontsize=16)
    plt.minorticks_on()
    if log:
        plt.xscale("log")
        plt.xticks([0.1, 1, 10, 100])
        ax.set_xticklabels([0.1, 1, 10, 100])
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
