#!/usr/bin/env python

from matplotlib import pyplot as plt
import matplotlib
import argparse
import scipy.integrate
import sys
import os
import math


class col:
    head = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    end = '\033[0m'
    BOLD = '\033[1m'
    red = '\033[31m'


def main():
    """
    Exectute the vs_plot script
    """

    title, vsLegends, vsPaths, zoom, \
        libraryIDstr, truePosIDstr, trueNegIDstr, ommitIDstr, \
        ref, log, gui, roc, enrich = parseArgs()

    # Get the truePosID range in list format
    libraryIDlist = makeIDlist(libraryIDstr, "Library IDs (not displayed): ",
                               False)
    truePosIDlist = makeIDlist(truePosIDstr, "True positive ID list: ", True)
    trueNegIDlist = makeIDlist(trueNegIDstr, "True negative ID list: ", True)
    ommitIDlist = makeIDlist(ommitIDstr, "Ommit ID list: ", True)

    # Generate a dictionary containing the refinement ligands, if any
    # refinement ligand was submitted
    if ref:
        refDict = makeRefDict(ref)
    else:
        refDict = {}

    # Read the results of each VS and keep only the ligIDs that are common
    # to all of them
    vsIntersects, libraryCount, truePosCount, trueNegCount, ommitCount \
        = intersectResults(vsPaths, truePosIDlist, trueNegIDlist, ommitIDlist)

    # Calculate % of total curves for each of these (write file + return data)
    percPaths = []
    for vsPath, vsIntersect in zip(vsPaths, vsIntersects):
        vsDir = os.path.dirname(vsPath)
        # print knownIDfirst, knownIDlast, ommitIDfirst, ommitIDlast
        if enrich:
            percPath = writePercFile(vsIntersect, vsDir, "enrich", refDict,
                                     "library", libraryIDstr,
                                     libraryIDlist, libraryCount,
                                     "true_pos", truePosIDstr,
                                     truePosIDlist, truePosCount,
                                     ommitIDstr, ommitIDlist)
        elif roc:
            percPath = writePercFile(vsIntersect, vsDir, "ROC", refDict,
                                     "true_neg", trueNegIDstr,
                                     trueNegIDlist, trueNegCount,
                                     "true_pos", truePosIDstr,
                                     truePosIDlist, truePosCount,
                                     ommitIDstr, ommitIDlist)

        percPaths.append(percPath)
        # writeROCfile()

    # Extract the data from the vs percent data (in both enrichment curves and
    # ROC curves, the truePositive count would be used to draw the perfect curve
    plotData, perfect, random, xLim, yLim = extractPlotData(percPaths,
                                                            vsLegends,
                                                            truePosCount,
                                                            zoom)

    # FIX AND COMPUTE ON ONE CURVE AT A TIME, on percent vs data?
    # getAUC_NSQ(plotData, perfect)

    # Define title and axis names based on mode
    if enrich:
        xAxisName = "% of ranked database (total=" + str(libraryCount) + ")"
        yAxisName = "% of known ligands found (total=" + str(truePosCount) + ")"
    elif roc:
        xAxisName = "% true negatives (total=" + str(trueNegCount) + ")"
        yAxisName = "% true positives (total=" + str(truePosCount) + ")"

    # Plot the data calculated by writePercFile, and read in by extracPlotData
    plot(title, plotData, perfect, random, xLim, yLim,
         xAxisName, yAxisName, gui, log, zoom)

    # Write the command used to execute this script into a log file
    writeCommand()

    print("\n")


def parseArgs():
    """
    Parsing and returning arguments
    """

    # Definition of arguments
    descr = "Feed VS result data (however many files), plots ROC curves or" \
        " Enrichment curves"
    descr_title = "Provide a title for the graph, also used as filename"
    descr_results = "Provide resultDataFiles.csv and 'legend titles' for" \
        " each curve: 'legend1!' data1.csv 'legend2?' data2.csv" \
        " 'legend4!!' data4.csv"
    descr_zoom = "X-axis percentage to be displayed in the zoomed subplot"
    descr_libraryIDstr = "Provide the ID range of the full library screened"
    descr_truePosIDstr = "Provide the IDs of true positive ligands" \
        " lib (format: 1-514,6001,6700-6702)"
    descr_trueNegIDstr = "Provide the IDs of true negative ligands" \
        " lib (format: 1-514,6001,6700-6702)"
    descr_ommitIDstr = "Provide the IDs of ligands to ommit" \
        " from the VS data, same format at knownIDs"
    descr_ref = "Refinement ligand(s) used on this GPCR binding pocket" \
        " refinement. Provide ligand name and ID in the following format:" \
        " lig1:328,lig2:535"
    descr_log = "Draw this plot on a log scale for the X axis"
    descr_gui = "Use this flag to display plot: saves to .png by the default"
    descr_modROC = "This mode will plot a ROC curve, which compares the true" \
        " positive rate to the true negative rate"
    descr_modEnrich = "This mode will plot an Enrichment curve, which" \
        " evaluates the rate of recovery of true positives as a function of" \
        " full library"

    # adding arguments to the parser
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("results", help=descr_results, nargs="+")
    parser.add_argument("zoom", help=descr_zoom)
    parser.add_argument("libraryIDstr", help=descr_libraryIDstr)
    parser.add_argument("truePosIDstr", help=descr_truePosIDstr)
    parser.add_argument("trueNegIDstr", help=descr_trueNegIDstr)
    parser.add_argument("ommitIDstr", help=descr_ommitIDstr)
    parser.add_argument("--ref", help=descr_ref)
    parser.add_argument("-log", action="store_true", help=descr_log)
    parser.add_argument("-gui", action="store_true", help=descr_gui)
    parser.add_argument("-modROC", action="store_true", help=descr_modROC)
    parser.add_argument("-modEnrich", action="store_true", help=descr_modEnrich)

    # parsing args
    args = parser.parse_args()
    title = args.title
    results = args.results
    zoom = float(args.zoom)
    libraryIDstr = args.libraryIDstr
    truePosIDstr = args.truePosIDstr
    trueNegIDstr = args.trueNegIDstr
    ommitIDstr = args.ommitIDstr
    ref = args.ref
    log = args.log
    gui = args.gui
    roc = args.modROC
    enrich = args.modEnrich

    if not roc and not enrich:
        print("\nOne the modes (ROC or Enrichment curve) has to be chosen\n")
        sys.exit()
    elif roc and enrich:
        print("\nOnly one the modes (ROC or Enrichment curve) can be chosen\n")
        sys.exit()

    # Extrac the VS results paths and legends
    vsPaths = []
    vsLegends = []
    i = 0
    while i < len(results):
        vsLegends.append(results[i])
        vsPaths.append(results[i + 1])
        i += 2

    return title, vsLegends, vsPaths, zoom, \
        libraryIDstr, truePosIDstr, trueNegIDstr, ommitIDstr, \
        ref, log, gui, roc, enrich


def makeIDlist(stringID, blurb, printOut):
    """
    Get a string defining which IDs to be generated into a list
    """

    print(col.head + "\n\t*MAKING LIGAND ID LIST*" + col.end)

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

    # Output this rangeID to the prompt
    if printOut:
        print("\n" + blurb + str(rangeID))
    else:
        print("\n" + blurb)

    return rangeID


def makeRefDict(refStr):
    """
    Get a string describing refinement ligands and their ID, and generate a
    dictionary to store that information, which is used when plotting the curves
    """

    print(col.head + "\n\t*MAKING REFERENCE LIGANDS*" + col.end)

    # This dictionary will be used for plotting
    refDict = {}

    refList = refStr.split(",")

    for ref in refList:
        refName, refID = ref.split(":")
        refDict[int(refID)] = refName

    # Output the reference ligands
    print("\nReference ligand list: " + str(refDict))

    return refDict


def intersectResults(vsPaths, truePosIDlist, trueNegIDlist, ommitIDlist):
    """
    Read in the results provided in .csv format, and figure out the intersect
    between each of those results set based on the ligIDs.
    Also check that the truePositive and trueNegative sets (if provided) are
    fully present in the intersect set: send a WARNING if they are not
    Then return the results set containing only the intersect results for each
    set.
    """

    print(col.head + "\n\t*DEFINING INTERSECT*" + col.end)

    allVsResults = []
    allLigIDs = []

    # Read results and populate vsResult
    for resultPath in vsPaths:
        resultFile = open(resultPath, 'r')
        resultLines = resultFile.readlines()
        resultFile.close()

        vsResult = []
        ligIDs = []
        # loop over the vs result lines omitting the first row
        # print resultPath, len(resultLines)
        for line in resultLines[1:]:
            ligInfo = line.strip().split(",")
            vsResult.append(ligInfo)
            ligIDs.append(int(ligInfo[0]))

        allVsResults.append(vsResult)
        allLigIDs.append(set(ligIDs))

    # Get the intersection set
    intersectLigID = set.intersection(*allLigIDs)
    libraryCount = len(intersectLigID)

    # Verify that truePositives are contained within the intersect of ligIDs
    truePosCount = getDockedLigandSet(truePosIDlist,
                                      "true positive IDs: ", intersectLigID)

    # Verify that trueNegatives are contained within the intersect of ligIDs
    trueNegCount = getDockedLigandSet(trueNegIDlist,
                                      "true negative IDs: ", intersectLigID)

    # Verify that ommits are contained within the intersect of ligIDs
    ommitCount = getDockedLigandSet(ommitIDlist,
                                    "ommit IDs: ", intersectLigID)

    allVsResultsIntersect = []
    # Loop over vsResults and keep only the ones present in the intersection
    for resultPath, vsResult in zip(vsPaths, allVsResults):
        # print resultPath, len(vsResult)
        vsResultIntersect = []
        for i, ligInfo in enumerate(vsResult):
            if int(ligInfo[0]) in intersectLigID:
                vsResultIntersect.append(ligInfo)
        allVsResultsIntersect.append(vsResultIntersect)

    return allVsResultsIntersect, libraryCount, \
        truePosCount, trueNegCount, ommitCount


def getDockedLigandSet(ligIDlist, ligType, intersectLigID):
    """
    From a list of ligands expected to be docked, and the actual list of all
    ligands docked, return an updated ligand set of actual docked ligands (and
    their count)
    """

    ligIDset = set(ligIDlist)
    intersect_ligID = set.intersection(ligIDset, intersectLigID)
    missingLigs = ligIDset - intersect_ligID
    if len(missingLigs) > 0:
        print(col.red + "\nWARNING: " + col.end +
              "missing IDs " + ligType + str(missingLigs))
    ligCount = len(intersect_ligID)

    return ligCount


def writePercFile(vsIntersect, vsDir, mode, refDict,
                  xAxisName, xAxisIDstr, xAxisIDlist, xCount,
                  yAxisName, yAxisIDstr, yAxisIDlist, yCount,
                  ommitIDstr, ommitIDlist):
    """
    Given this VS result, and information about the ID of known actives
    in the library, write in a file the information to plot an enrichment curve
    """

    print(col.head + "\n\t*WRITING ENRICHMENT DATA*" + col.end)

    m = mode + "_"
    x = xAxisName + xAxisIDstr + "_"
    y = yAxisName + yAxisIDstr + "_"
    ommit = "ommits_" + ommitIDstr + "_"

    # Create filename
    percentPath = vsDir + "/" + m + x + y + ommit + vsDir + ".csv"
    print("\n" + percentPath)
    percentDataFile = open(percentPath, "w")

    # Initialize file with values of 0,0,0.0,0.0,0.0
    # Those are x,y,xPercent,yPercent,yPerfect,yRandom values
    percentDataFile.write("0,0,0.0,0.0,0.0,0.0\n")

    # Build the % data
    X = 0
    Y = 0
    val = 0
    for ligInfo in vsIntersect:
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
            if ligID in yAxisIDlist:
                # print "known ligand", ligInfo
                Y += 1

            # For each ligand in the full VS, increase X and write
            # the X,Y pair to the data file
            if ligID in xAxisIDlist:
                X += 1

            # Calculate percentage X and Y
            Xpercent = (X * 100.0) / xCount
            Ypercent = (Y * 100.0) / yCount

            # Calculate what the perfect line should be
            if val < yCount:
                val += 1
            perfect = (val * 100.0) / yCount

            # Find if the current ligand is one of the refinement ligands, if
            # so save its Xpercent value, in order to add a marker at the
            # position where it was recovered in the VS screen
            if ligID in refDict.keys():
                ligName = refDict[ligID]
                percentDataFile.write(str(X) + "," + str(Y) + "," +
                                      str(Xpercent) + "," + str(Ypercent) +
                                      "," + str(perfect) + "," + str(Xpercent) +
                                      "," + ligName + "\n")
            else:
                percentDataFile.write(str(X) + "," + str(Y) + "," +
                                      str(Xpercent) + "," + str(Ypercent) +
                                      "," + str(perfect) + "," + str(Xpercent) +
                                      "\n")

    percentDataFile.close()

    return percentPath


def extractPlotData(percentPaths, vsLegends, truePosCount, zoom):
    """
    Read the % result data files, return the data for plotting
    """

    print(col.head + "\n\t*GETTING PLOTTING DATA*" + col.end)

    # Variables that define the x and y limits for the zoomed in subplot
    xLim = 0.0
    yLim = 0.0

    plotData = []

    for percentPath, vsLegend in zip(percentPaths, vsLegends):
        # Read the % data file
        dataFile = open(percentPath, "r")
        dataLines = dataFile.readlines()
        dataFile.close()

        refPlot = {}

        print "\nData path:", percentPath

        perfect = []
        random = []
        X = []
        Y = []
        for line in dataLines:
            # Get the data from the file
            ll = line.split(",")
            xPercent = float(ll[2])
            yPercent = float(ll[3])
            perfVal = float(ll[4])
            randVal = float(ll[5])

            # Create the data curve
            X.append(xPercent)
            Y.append(yPercent)
            # Create the perfect curve
            perfect.append(perfVal)
            # Create the random curve
            random.append(randVal)

            # Create the subplot limits
            if xPercent <= zoom:
                if yLim < yPercent:
                    xLim = xPercent
                    yLim = yPercent

            # Collect the X position of the refinement ligand(s)
            if len(ll) == 7:
                ligName = ll[6]
                ligXY = [xPercent, yPercent]
                refPlot[ligName] = ligXY

        plotData.append((X, Y, vsLegend, refPlot))

    return plotData, perfect, random, xLim, yLim


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


def plot(title, plotData, perfect, random, xLim, yLim,
         xAxis, yAxis, gui, log, zoom):
    """
    Plot the data provided as argument, to draw curves
    """

    print(col.head + "\n\t*PLOTTING DATA*" + col.end)

    # Setting up the figure
    fig = plt.figure(figsize=(13, 12), dpi=100)
    ax = fig.add_subplot(111)
    # Create the ZOOMED graph, if requested
    if zoom != 0.0:
        ax2 = plt.axes([.17, .35, .2, .2])
    else:
        ax2 = None

    # Setting up color scheme
    cm = plt.get_cmap("spectral")
    cNorm = matplotlib.colors.Normalize(vmin=0, vmax=len(plotData))
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    # ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(len(plotData))])

    # Drawing data on the figure
    for i, plotDatum in enumerate(plotData):
        X, Y = drawLine(ax, ax2, plotDatum, i, zoom, scalarMap)

    # Plot the RANDOM and PERFECT curves on the zoomed and main graph
    if zoom != 0.0:
        ax2.plot(X, perfect, color="grey")
        ax2.plot(X, random, "--", color="grey")
        ax2.tick_params(axis="both", which="major", labelsize=15)
        ax2.set_title("Zoom of the first " + str(zoom) + "%", fontsize=15)

    # Now plot random and perfect curves, common for all plotted curves
    ax.plot(X, perfect, color="grey")
    ax.plot(X, random, "--", color="grey")
    # print X
    # print perfect

    # Here axis and ticks are improved
    ax.set_xlabel(xAxis, fontsize=30)
    ax.set_ylabel(yAxis, fontsize=30)

    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=30)
    ax.set_title(title, fontsize=30)
    ax.legend(loc="upper left", prop={'size': 30})
    ax.axis('tight')

    if log:
        ax.set_xscale("symlog", linthreshx=0.01)
        # Quick and dirty fix for the axis
        ax.set_xticks([0.1, 1, 10, 100])
        ax.set_xticklabels([0.1, 1, 10, 100])

        # Setting ZOOMED ax2 graph
        if zoom != 0.0:
            ax2.set_xscale("log")
            ax2.set_xlim([0, zoom])
            ax2.set_ylim([0, yLim])
            # xLimRound = int(xLim * 100) / 100.0
            # yLimRound = int(yLim * 100) / 100.0
            # ax2.set_yticks([0, yLimRound])
            # ax2.set_xticklabels([])
            # print xLimRound
            # plt.setp(ax2, xlim=(0, zoom), ylim=(0, yLim),
            #         xticks=[0, zoom], yticks=[0, yLimRound])

    if gui:
        plt.show()
    else:
        fileName = title.replace(" ", "_") + ".png"
        plt.savefig(fileName, bbox_inches="tight")


def drawLine(ax, ax2, plotDatum, i, zoom, scalarMap):
    """
    Draw the line corresponding to the set of data passed in arguments
    """
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

    X = plotDatum[0]
    Y = plotDatum[1]
    # Add the value 0 in order to have curves that start at the origin
    X = X
    Y = Y
    plotLegend = plotDatum[2]
    refPlot = plotDatum[3]

    # Plot this curve
    ax.plot(X, Y, label=plotLegend, linewidth=lw, color=color)

    # Plot a blow up of the first X%
    if zoom != 0.0:
        ax2.plot(X, Y, color=color)

    # Plot a vertical line for each refinement ligand
    for ligName in refPlot.keys():
        xPos, yPos = refPlot[ligName]
        ax.axvline(x=xPos, ymax=yPos/100., color=color,
                   linewidth=3, linestyle='--')
        # print ligName, xPos, yPos
        # ax.text(xPos, -2, ligName, rotation=-90,
        #        color=color, transform=ax.transData)

    return X, Y


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
