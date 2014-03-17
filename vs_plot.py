#!/usr/bin/env python

from matplotlib import pyplot as plt
import matplotlib
import argparse


def main():

    title, rocLegends, rocPaths, zoom, log, gui = parseArgs()

    # Extract the data from the ROC files
    rocData, perfect, totalLib, totalKnown, xLim, yLim = getData(rocPaths,
                                                                 rocLegends,
                                                                 zoom)
    # Plot the values 2 and 3, which correspond to the percentage X and Y
    plot(title, rocData, perfect, xLim, yLim,
         totalLib, totalKnown, gui, log, zoom)


def parseArgs():
    """
    Parsing and returning arguments
    """

    # Definition of arguments
    descr = "Feed rocData (however many files), plots ROC curves"
    descr_title = "Provide a title for the graph, also used as filename"
    descr_roc = "Provide rocDataFiles and legend titles for each curve:" \
        " 'legend1!' data1.csv 'legend2?' data2.csv 'legend4!!' data4.csv"
    descr_zoom = "Give the percent of ranked database to be displayed in the" \
        " zoomed subplot"
    descr_log = "Draw this plot on a log scale for the X axis"
    descr_gui = "Use this flag to display plot: saves to .png by the default"

    # adding arguments to the parser
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("-log", action="store_true", help=descr_log)
    parser.add_argument("-gui", action="store_true", help=descr_gui)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("roc", help=descr_roc, nargs="+")
    parser.add_argument("zoom", help=descr_zoom)

    # parsing args
    args = parser.parse_args()
    title = args.title
    roc = args.roc
    zoom = float(args.zoom)
    log = args.log
    gui = args.gui

    # Extrac the ROC paths and legends from the roc variable
    rocPaths = []
    rocLegends = []
    i = 0
    while i < len(roc):
        rocLegends.append(roc[i])
        rocPaths.append(roc[i + 1])
        i += 2

    return title, rocLegends, rocPaths, zoom, log, gui


def getData(rocPaths, rocLegends, zoom):
    """
    Read the ROC data files, return the data for plotting
    """

    # Variables that define the x and y limites for the zoomed in subplot
    xLim = 0.0
    yLim = 0.0

    rocData = []

    for rocPath, rocLegend in zip(rocPaths, rocLegends):
        # Read the ROC data file
        rocFile = open(rocPath, "r")
        rocLines = rocFile.readlines()
        rocFile.close()

        # Get the values total ligands in the library, and total known ligands
        # from the last line of the ROC data
        totalLib = int(rocLines[-1].split(",")[0])
        totalKnown = int(rocLines[-1].split(",")[1])

        X = []
        Y = []
        perfect = []
        val = 0
        for line in rocLines:
            # Get the data from the file
            ll = line.split(",")
            xPercent = float(ll[2])
            yPercent = float(ll[3])

            # Create the data curve
            X.append(xPercent)
            Y.append(yPercent)

            # Create the perfect curve
            if val < totalKnown:
                val += 1
            perfect.append((val * 100.0) / totalKnown)

            # Create the subplot limits
            if xPercent <= zoom:
                if yLim < yPercent:
                    xLim = xPercent
                    yLim = yPercent

        rocData.append((X, Y, rocLegend))

    return rocData, perfect, totalLib, totalKnown, xLim, yLim


def plot(title, rocData, perfect, xLim, yLim,
         totalLib, totalKnown, gui, log, zoom):
    """
    Plot the data provided as argument, to draw ROC curves
    """

    # Setting up the figure
    fig = plt.figure(figsize=(13, 12), dpi=100)
    ax = fig.add_subplot(111)

    # Setting up color scheme
    cm = plt.get_cmap("Set1")
    cNorm = matplotlib.colors.Normalize(vmin=0, vmax=len(rocData) - 1)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    #ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(len(rocData))])

    # Drawing data on the figure
    for i, rocDatum in enumerate(rocData):
        color = scalarMap.to_rgba(i)
        X = rocDatum[0]
        Y = rocDatum[1]
        rocLegend = rocDatum[2]

        # Plot this curve
        ax.plot(X, Y, label=rocLegend, linewidth=2, color=color)

        # Plot a blow up of the first X%
        if zoom != 0.0:
            ax2 = plt.axes([.17, .25, .2, .2])
            ax2.semilogx(X, Y, color=color)
            ax2.semilogx(X, perfect, color="grey")
            ax2.semilogx(X, X, "--", color="grey")
            xLimRound = int(xLim * 100) / 100.0
            yLimRound = int(yLim * 100) / 100.0
            plt.setp(ax2, xlim=(0, xLim), ylim=(0, yLim),
                     xticks=[0, xLimRound], yticks=[0, yLimRound])
            ax2.tick_params(axis="both", which="major", labelsize=8)

    # Now plot random and perfect curves, common for all plotted curves
    ax.plot(X, X, "--", color="grey")
    ax.plot(X, perfect, color="grey")

    # Here axis and ticks are improved
    ax.set_xlabel("% of ranked database (total=" + str(totalLib) + ")",
                  fontsize=16)
    ax.set_ylabel("% of known ligands found (total=" + str(totalKnown) + ")",
                  fontsize=16)
    ax.minorticks_on()
    if log:
        ax.set_xscale("log")
        ax.set_xticks([0.1, 1, 10, 100])
        ax.set_xticklabels([0.1, 1, 10, 100])
    ax.set_title(title, fontsize=18)
    ax.legend(loc="upper left", prop={'size': 12})
    ax.axis('tight')

    if gui:
        plt.show()
    else:
        fileName = title.replace(" ", "_") + ".png"
        plt.savefig(fileName, bbox_inches="tight")


if __name__ == "__main__":
    main()
