#!/usr/bin/env python

import argparse
import sys
import os
import plotting


def main():
    """
    Exectute the vs_plot_enrich script
    """

    title, vsLegends, vsPaths, \
        truePosIDstr, trueNegIDstr, ommitIDstr, \
        ref, gui = parseArgs()

    # Define mode
    mode = "ROC"
    # Define zoom
    zoom = 0.0
    # Define log
    log = False

    # Creating a plotting instance for access to all methods
    p = plotting.plotting()

    # Get the truePosID range in list format
    libraryIDlist = p.makeIDlist("0-0", "Library IDs (not displayed): ",
                                 False)

    truePosIDlist = p.makeIDlist(truePosIDstr, "True positive ID list: ", True)
    trueNegIDlist = p.makeIDlist(trueNegIDstr, "True negative ID list: ", True)
    ommitIDlist = p.makeIDlist(ommitIDstr, "Ommit ID list: ", True)

    # Generate a dictionary containing the refinement ligands, if any
    # refinement ligand was submitted
    if ref:
        refDict = p.makeRefDict(ref)
    else:
        refDict = {}

    # Read the results of each VS and keep only the ligIDs that are common
    # to all of them
    vsIntersects, libraryCount, truePosCount, trueNegCount, ommitCount \
        = p.intersectResults(vsPaths, truePosIDlist, trueNegIDlist, ommitIDlist)

    # Calculate % of total curves for each of these (write file + return data)
    percPaths = []
    for vsPath, vsIntersect in zip(vsPaths, vsIntersects):
        vsDir = os.path.dirname(vsPath)

        percPath = p.writePercFile(vsIntersect, vsDir, mode, refDict,
                                   "true_neg", trueNegIDstr,
                                   trueNegIDlist, trueNegCount,
                                   "true_pos", truePosIDstr,
                                   truePosIDlist, truePosCount,
                                   ommitIDstr, ommitIDlist)

        percPaths.append(percPath)

    # Extract the data from the vs percent data (in both enrichment curves and
    # ROC curves, the truePositive count would be used to draw the perfect curve
    plotData, xLim, yLim, scatterData = p.extractPlotData(percPaths, vsLegends,
                                                          zoom)

    # FIX AND COMPUTE ON ONE CURVE AT A TIME, on percent vs data?
    # p.getAUC_NSQ(plotData, perfect)

    # Define title and axis names based on mode
    xAxisName = "% true negatives (total=" + str(trueNegCount) + ")"
    yAxisName = "% true positives (total=" + str(truePosCount) + ")"

    # Plot the data calculated by writePercFile, and read in by extracPlotData
    p.plot(title, plotData, libraryCount, truePosCount,
           xLim, yLim, xAxisName, yAxisName, gui, log, zoom, mode, scatterData)

    # Write the command used to execute this script into a log file
    p.writeCommand(title)

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
    descr_truePosIDstr = "Provide the IDs of true positive ligands" \
        " lib (format: 1-514,6001,6700-6702)"
    descr_trueNegIDstr = "Provide the IDs of true negative ligands" \
        " lib (format: 1-514,6001,6700-6702)"
    descr_ommitIDstr = "Provide the IDs of ligands to ommit" \
        " from the VS data, same format at knownIDs"
    descr_ref = "Refinement ligand(s) used on this GPCR binding pocket" \
        " refinement. Provide ligand name and ID in the following format:" \
        " lig1:328,lig2:535"
    descr_gui = "Use this flag to display plot: saves to .png by the default"

    # adding arguments to the parser
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("results", help=descr_results, nargs="+")
    parser.add_argument("truePosIDstr", help=descr_truePosIDstr)
    parser.add_argument("trueNegIDstr", help=descr_trueNegIDstr)
    parser.add_argument("ommitIDstr", help=descr_ommitIDstr)
    parser.add_argument("-gui", action="store_true", help=descr_gui)
    parser.add_argument("--ref", help=descr_ref)

    # parsing args
    args = parser.parse_args()
    title = args.title
    results = args.results
    truePosIDstr = args.truePosIDstr
    trueNegIDstr = args.trueNegIDstr
    ommitIDstr = args.ommitIDstr
    ref = args.ref
    gui = args.gui

    # Extrac the VS results paths and legends
    vsPaths = []
    vsLegends = []
    i = 0
    while i < len(results):
        vsLegends.append(results[i])
        vsPaths.append(results[i + 1])
        i += 2

    return title, vsLegends, vsPaths, \
        truePosIDstr, trueNegIDstr, ommitIDstr, \
        ref, gui

if __name__ == "__main__":
    main()
