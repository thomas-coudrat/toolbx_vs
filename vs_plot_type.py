#!/usr/bin/env python

from matplotlib import pyplot as plt
import matplotlib
import argparse
import scipy.integrate
import sys
import os
import math


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
    writeCommand(title)

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
