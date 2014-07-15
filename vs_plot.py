#!/usr/bin/env python

from matplotlib import pyplot as plt
import matplotlib
import argparse
import scipy.integrate
import sys
import os
import math


def main():

    title, rocLegends, resultPaths, zoom, \
        knownIDstr, ommitIDstr, ref, log, gui = parseArgs()

    #
    # Get the knownID range in list format
    #
    knownIDlist = makeIDlist(knownIDstr)
    ommitIDlist = makeIDlist(ommitIDstr)
    print "\n", "Known ID string", knownIDstr
    print "Known ID list", knownIDlist
    print "Ommit ID string", ommitIDstr
    print "Ommit ID list", ommitIDlist, "\n"

    #
    # Generate a dictionary containing the refinement ligands, if any
    # refinement ligand was submitted
    #
    if ref:
        # print ref
        refDict = makeRefDict(ref)
        # print refDict
    else:
        refDict = {}

    #
    # Read the results of each VS and keep only the ligIDs that are common
    # to all of them
    #
    allVsResultsIntersect = intersectResults(resultPaths)
    # print resultPaths
    rocPaths = []
    allTotalLibs = []
    allTotalKnowns = []
    # Calculate ROC curves for each of these (write file + return data)
    for resultPath, vsResult in zip(resultPaths, allVsResultsIntersect):
        vsDir = os.path.dirname(resultPath)
        # print knownIDfirst, knownIDlast, ommitIDfirst, ommitIDlast
        rocPath, totalLib, totalKnown = writeRocFile(vsResult, vsDir,
                                                     knownIDstr,
                                                     knownIDlist,
                                                     ommitIDstr,
                                                     ommitIDlist,
                                                     refDict)
        print rocPath, totalLib, totalKnown
        rocPaths.append(rocPath)
        allTotalLibs.append(totalLib)
        allTotalKnowns.append(totalKnown)

    #
    # Make sure the total library size and the total number of knowns is the
    # same between all vsResults. Exit and print statement if it isn't
    #
    for totalL in allTotalLibs:
        if totalLib != totalL:
            print "Total library size not matching between VS experiments"
            sys.exit()
    for totalK in allTotalKnowns:
        if totalKnown != totalK:
            print "Total number of knowns not matching between VS experiments"
            sys.exit()

    #
    # Extract the data from the ROC files
    #
    rocData, perfect, xLim, yLim = extractRocData(rocPaths,
                                                  rocLegends,
                                                  totalKnown,
                                                  zoom)
    getAUC_NSQ(rocData, perfect)

    #
    # Plot the values 2 and 3, which correspond to the percentage X and Y
    #
    plot(title, rocData, perfect, xLim, yLim, totalLib, totalKnown,
         gui, log, zoom)

    # Write the command used to execute this script into a log file
    writeCommand()


def parseArgs():
    """
    Parsing and returning arguments
    """

    # Definition of arguments
    descr = "Feed VS result data (however many files), plots ROC curves"
    descr_title = "Provide a title for the graph, also used as filename"
    descr_results = "Provide resultDataFiles.csv and 'legend titles' for" \
        " each curve: 'legend1!' data1.csv 'legend2?' data2.csv" \
        " 'legend4!!' data4.csv"
    descr_zoom = "Give the percent of ranked database to be displayed in the" \
        " zoomed subplot"
    descr_knownIDstr = "Provide the IDs of known actives ligands" \
        " lib (format: 1-514,6001,6700-6702)"
    descr_ommitIDstr = "Provide the IDs of ligands to ommit" \
        " from the ROC curve data, same format at knownIDs"
    descr_ref = "Refinement ligand(s) used on this GPCR binding pocket" \
        " refinement. Provide ligand name and ID in the following format:" \
        " lig1:328,lig2:535"
    descr_log = "Draw this plot on a log scale for the X axis"
    descr_gui = "Use this flag to display plot: saves to .png by the default"

    # adding arguments to the parser
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("results", help=descr_results, nargs="+")
    parser.add_argument("zoom", help=descr_zoom)
    parser.add_argument("knownIDstr", help=descr_knownIDstr)
    parser.add_argument("ommitIDstr", help=descr_ommitIDstr)
    parser.add_argument("--ref", help=descr_ref)
    parser.add_argument("-log", action="store_true", help=descr_log)
    parser.add_argument("-gui", action="store_true", help=descr_gui)

    # parsing args
    args = parser.parse_args()
    title = args.title
    results = args.results
    zoom = float(args.zoom)
    knownIDstr = args.knownIDstr
    ommitIDstr = args.ommitIDstr
    ref = args.ref
    log = args.log
    gui = args.gui

    # Extrac the ROC paths and legends from the roc variable
    resultPaths = []
    rocLegends = []
    i = 0
    while i < len(results):
        rocLegends.append(results[i])
        resultPaths.append(results[i + 1])
        i += 2

    return title, rocLegends, resultPaths, zoom, \
        knownIDstr, ommitIDstr, ref, log, gui


def makeIDlist(stringID):
    """
    Get a string defining which IDs to be generated into a list
    """

    # This stores the range of IDs into a list
    rangeID = []

    IDportions = stringID.split(",")

    for portion in IDportions:
        # Treat ranges of IDs
        if "-" in portion:
            start, end = portion.split("-")
            start = int(start)
            end = int(end)
            # Do not add the value 0 to the list
            if start == 0 or end == 0:
                pass
            else:
                rangeID = rangeID + range(start, end + 1)
        # Treat single IDs
        else:
            portion = int(portion)
            # Do not add the value 0 to the list
            if portion == 0:
                pass
            else:
                rangeID.append(int(portion))

    return rangeID


def makeRefDict(refStr):
    """
    Get a string describing refinement ligands and their ID, and generate a
    dictionary to store that information, which is used when plotting the ROC
    curve
    """
    # This dictionary will be used for plotting
    refDict = {}

    refList = refStr.split(",")

    print refList

    for ref in refList:
        refName, refID = ref.split(":")
        refDict[int(refID)] = refName

    return refDict


def intersectResults(resultPaths):
    """
    Read in the results provided in .csv format, and figure out the intersect
    between each of those results set based on the ligID. Then return the
    results set containing only the intersect results for each set.
    """

    allVsResults = []
    allLigIDs = []

    # Read results and populate vsResultsRaw
    for resultPath in resultPaths:
        resultFile = open(resultPath, 'r')
        resultLines = resultFile.readlines()
        resultFile.close()

        vsResult = []
        ligIDs = []
        # loop over the vs result lines omitting the first row
        print resultPath, len(resultLines)
        for line in resultLines[1:]:
            ligInfo = line.strip().split(",")
            vsResult.append(ligInfo)
            ligIDs.append(ligInfo[0])

        allVsResults.append(vsResult)
        allLigIDs.append(set(ligIDs))

    # Get the intersection set
    intersectLigID = set.intersection(*allLigIDs)

    allVsResultsIntersect = []
    # Loop over vsResults and keep only the ones present in the intersection
    for resultPath, vsResult in zip(resultPaths, allVsResults):
        print resultPath, len(vsResult)
        vsResultIntersect = []
        for i, ligInfo in enumerate(vsResult):
            if ligInfo[0] in intersectLigID:
                vsResultIntersect.append(ligInfo)
        allVsResultsIntersect.append(vsResultIntersect)

    return allVsResultsIntersect


def writeRocFile(vsResult, vsDir,
                 knownIDstr, knownIDlist,
                 ommitIDstr, ommitIDlist,
                 refDict):
    """
    Given this VS result, and information about the ID of known actives
    in the library, write in a file the information to plot a ROC curve
    Also collect for each
    """

    known = "knowns_" + knownIDstr
    ommit = "ommits_" + ommitIDstr

    # Create filename
    rocPath = vsDir + "/roc_" + known + "_" + ommit + "_" + vsDir + ".csv"
    print "\t", rocPath
    rocDataFile = open(rocPath, "w")

    # Loop over the results once, in order to check for the actual total number
    # of knowns present in them
    totalKnowns = 0
    for ligInfo in vsResult:
        ligID = int(ligInfo[0])
        if ligID in knownIDlist:
            totalKnowns += 1

    # Get the total library size
    if len(ommitIDlist):
        totalLibrary = len(vsResult)
    else:
        totalLibrary = len(vsResult) - len(ommitIDlist)

    print "\nTotal knowns:", totalKnowns
    print "Total library - knowns:", totalLibrary, totalKnowns

    print vsDir, len(vsResult)

    X = 0
    Y = 0
    for ligInfo in vsResult:
        # print ligInfo
        ligID = int(ligInfo[0])

        # Skip if ligID is part of the range that needs to be ommited
        # If the ommit values are '0', then there is no ligand to ommit
        if ligID in ommitIDlist:
            # print "ligand skipped", ligInfo
            continue
        # Otherwise proceed normally
        else:
            # When the sorted ligID corresponds to a known, increase
            # the value of Y by 1
            if ligID in knownIDlist:
                # print "known ligand", ligInfo
                Y += 1

            # For each ligand in the full VS, increase X and write
            # the X,Y pair to the data file
            X += 1

            # Calculate percentage X and Y
            Xpercent = (X * 100.0) / totalLibrary
            Ypercent = (Y * 100.0) / totalKnowns

            #
            # Find if the current ligand is one of the refinement ligands, if
            # so save its Xpercent value, in order to add a marker at the
            # position where it was recovered in the VS screen
            #
            if ligID in refDict.keys():
                ligName = refDict[ligID]
                rocDataFile.write(str(X) + "," + str(Y) + "," +
                                  str(Xpercent) + "," + str(Ypercent) + "," +
                                  ligName + "\n")
            else:
                rocDataFile.write(str(X) + "," + str(Y) + "," +
                                  str(Xpercent) + "," + str(Ypercent) + "\n")

    rocDataFile.close()

    return rocPath, totalLibrary, totalKnowns


def extractRocData(rocPaths, rocLegends, totalKnown, zoom):
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

        refPlot = {}

        print "ROC PATH:", rocPath

        perfect = []
        X = []
        Y = []
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

            # Collect the X position of the refinement ligand(s)
            if len(ll) == 5:
                ligName = ll[4]
                ligXY = [xPercent, yPercent]
                refPlot[ligName] = ligXY

        rocData.append((X, Y, rocLegend, refPlot))

    return rocData, perfect, xLim, yLim


def getAUC_NSQ(rocData, perfect):
    """
    Calculate AUC and NSQ_AUC for each curve, and return a list with those
    values (corresponds to the order of rocData)
    """

    # aucData = []

    print "perfect=", len(perfect)
    # perfectSq = [math.sqrt(i) for i in perfect]

    for rocDatum in rocData:
        X = rocDatum[0]
        Y = rocDatum[1]
        legend = rocDatum[2]

        Xsq = [math.sqrt(i) for i in X]
        # Ysq = [math.sqrt(i) for i in Y]

        print "X=", len(X)
        print "Y=", len(Y)

        auc = scipy.integrate.trapz(Y, X)
        aucSq = scipy.integrate.trapz(Y, Xsq)
        # auc2 = scipy.integrate.simps(Y, X)

        perf = scipy.integrate.trapz(perfect, X)
        perfSq = scipy.integrate.trapz(perfect, Xsq)
        # perf2 = scipy.integrate.simps(perfect, X)

        rand = scipy.integrate.trapz(X, X)
        randSq = scipy.integrate.trapz(X, Xsq)
        # rand2 = scipy.integrate.simps(X, X)

        print "**************"

        print legend
        print "auc", auc        # , auc2
        print "aucSq", aucSq
        print "perfect", perf   # , perf2
        print "perfectSq", perfSq
        print "rand", rand      # , rand2
        print "randSq", randSq
        print

        nsq_auc = (aucSq - randSq) / (perfSq / randSq)
        nsq_auc_perf = (perfSq - randSq) / (perfSq / randSq)
        nsq_auc_rand = (randSq - randSq) / (perfSq / randSq)

        print "NSQ_AUC:", nsq_auc
        print "NSQ_AUC - perf:", nsq_auc_perf
        print "NSQ_AUC - rand:", nsq_auc_rand

        print "**************"


def plot(title, rocData, perfect, xLim, yLim,
         totalLib, totalKnown, gui, log, zoom):
    """
    Plot the data provided as argument, to draw ROC curves
    """

    # Setting up the figure
    fig = plt.figure(figsize=(13, 12), dpi=100)
    ax = fig.add_subplot(111)
    # Create the ZOOMED graph, if requested
    if zoom != 0.0:
        ax2 = plt.axes([.17, .35, .2, .2])

    # Setting up color scheme
    cm = plt.get_cmap("spectral")
    cNorm = matplotlib.colors.Normalize(vmin=0, vmax=len(rocData))
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    # ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(len(rocData))])

    # Add a value 0 to the perfect curve
    perfect = [0.0] + perfect

    # Drawing data on the figure
    for i, rocDatum in enumerate(rocData):
        # Set color for crystal structures, and the LDM results have their
        # colors defined by the colormap
        if i == 0:
            color = 'black'
            lw = 4
        elif i == 1:
            color = 'grey'
            lw = 4
        else:
            color = scalarMap.to_rgba(i)
            lw = 2
        X = rocDatum[0]
        Y = rocDatum[1]
        X = [0.1] + X
        Y = [0.1] + Y
        rocLegend = rocDatum[2]
        refPlot = rocDatum[3]

        # Plot this curve
        ax.plot(X, Y, label=rocLegend, linewidth=lw, color=color)

        # Plot a blow up of the first X%
        if zoom != 0.0:
            ax2.plot(X, Y, color=color)

        # Plot a vertical line for each refinement ligand
        for ligName in refPlot.keys():
            xPos, yPos = refPlot[ligName]
            ax.axvline(x=xPos, ymax=yPos/100., color=color, linewidth=3)
            print ligName, xPos, yPos
            # ax.text(xPos, -2, ligName, rotation=-90,
            #        color=color, transform=ax.transData)

    # Plot the RANDOM and PERFECT curves on the zoomed and main graph
    if zoom != 0.0:
        ax2.plot(X, perfect, color="grey")
        ax2.plot(X, X, "--", color="grey")
        ax2.tick_params(axis="both", which="major", labelsize=15)
        ax2.set_title("Zoom of the first " + str(zoom) + "%", fontsize=15)

    # Now plot random and perfect curves, common for all plotted curves
    ax.plot(X, X, "--", color="grey")
    ax.plot(X, perfect, color="grey")
    # print X
    # print perfect

    # Here axis and ticks are improved
    ax.set_xlabel("% of ranked database (total=" + str(totalLib) + ")",
                  fontsize=30)
    ax.set_ylabel("% of known ligands found (total=" + str(totalKnown) + ")",
                  fontsize=30)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=30)
    ax.set_title(title, fontsize=30)
    ax.legend(loc="upper left", prop={'size': 30})
    ax.axis('tight')

    if log:
        ax.set_xscale("symlog", linthreshx=0.0001)
        # Quick and dirty fix for the axis
        ax.set_xticks([0.1, 1, 10, 100])
        ax.set_xticklabels([0.1, 1, 10, 100])

        # Setting ZOOMED ax2 graph
        if zoom != 0.0:
            ax2.set_xscale("log")
            ax2.set_xlim([0, zoom])
            ax2.set_ylim([0, yLim])
            # xLimRound = int(xLim * 100) / 100.0
            yLimRound = int(yLim * 100) / 100.0
            ax2.set_yticks([0, yLimRound])
            # ax2.set_xticklabels([])
            # print xLimRound
            # plt.setp(ax2, xlim=(0, zoom), ylim=(0, yLim),
            #         xticks=[0, zoom], yticks=[0, yLimRound])

    if gui:
        plt.show()
    else:
        fileName = title.replace(" ", "_") + ".png"
        plt.savefig(fileName, bbox_inches="tight")


def writeCommand():
    """
    Write down the command that was used to exectute this script in a log
    file, at the location where the script is executed. Also write the
    current working directory at the time of execution
    """

    cwd = os.getcwd()
    logFile = open("plot.log", "w")
    # Write the directory location: this is not executed upong sh call of the
    # plot.log, but serves as information
    logFile.write(cwd + "\n")
    logFile.write(sys.argv[0].split("/")[-1] + " ")
    for arg in sys.argv[1:]:
        if len(arg) > 0:
            # Deal with argument options (starting with '-')
            if arg[0] == "-":
                logFile.write(arg + " ")
            # Do not add "'" on argument if it already has them
            elif arg[0] == "'" and arg[-1] == "'":
                logFile.write(arg + " ")
            # Add the "'" around each other argument
            else:
                logFile.write("'" + arg + "' ")
    logFile.close()


if __name__ == "__main__":
    main()
