#!/usr/bin/env python

# Plot ROC curves describing the true positives vs. false positive rate.

# https://github.com/thomas-coudrat/toolbx_vs
# Thomas Coudrat <thomas.coudrat@gmail.com>

import argparse
import sys
import os
import plotting


def main():
    """
    Run script
    """

    title, vsLegends, vsPaths, vsColors, vsLines, \
        truePosIDstr, falsePosIDstr, xAxisName, yAxisName, \
        ref, gui, log, showAUC = parseArgs()

    # Define mode
    mode = "ROC"
    # Define zoom
    zoom = 0.0
    # Define scatterData
    scatterData = False

    # Creating a plotting instance for access to all methods
    p = plotting.plotting(title)

    # Get the truePosID range in list format
    truePosIDlist = p.makeIDlist(truePosIDstr,
                                 "True positive ID list: ",
                                 True)
    falsePosIDlist = p.makeIDlist(falsePosIDstr,
                                  "False positive ID list: ",
                                  True)
    libraryIDlist = truePosIDlist + falsePosIDlist

    # Generate a dictionary containing the refinement ligands, if any
    # refinement ligand was submitted
    if ref:
        refDict = p.makeRefDict(ref)
    else:
        refDict = {}

    # Read the results of each VS and keep only the ligIDs that are common
    # to all of them
    vsIntersects, ligIDintersectSet = p.intersectResults(vsPaths,
                                                         libraryIDlist)

    # Get updated true positive, true negative and library counts given the
    # intersect results
    truePosCount = p.updatedLigCounts(ligIDintersectSet,
                                      truePosIDlist,
                                      "true positives")
    falsePosCount = p.updatedLigCounts(ligIDintersectSet,
                                       falsePosIDlist,
                                       "false positives")
    # This value is actually not used, but it complies with the plot() function
    # libraryCount = truePosCount + falsePosCount
    libraryCount = p.updatedLigCounts(ligIDintersectSet,
                                      libraryIDlist,
                                      "whole library")

    # Calculate % of total curves for each of these (write file + return data)
    vsPockets = []
    for vsPath, vsIntersect in zip(vsPaths, vsIntersects):
        vsPocket = p.writePercFile(vsIntersect, vsPath, mode, refDict,
                                   "true_neg", falsePosIDstr,
                                   falsePosIDlist, falsePosCount,
                                   "true_pos", truePosIDstr,
                                   truePosIDlist, truePosCount)

        vsPockets.append(vsPocket)

    # Extract the data from the vs percent data (in both enrichment curves and
    # ROC curves)
    plotData, xLim, yLim = p.extractPlotData(vsPockets, vsLegends, zoom)

    # Calculate NSQ_AUC values for each curve and append data to plotData list
    p.getAUC_NSQ(plotData)

    # Define title and axis names based on mode
    yAxisName = yAxisName + " (total=" + str(truePosCount) + ")"
    xAxisName = xAxisName + " (total=" + str(falsePosCount) + ")"

    # Plot the data calculated by writePercFile, and read in by extracPlotData
    p.plotROC(title, plotData, vsColors, vsLines, libraryCount, truePosCount,
              xLim, yLim, xAxisName, yAxisName, gui, log,
              zoom, mode, showAUC)

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
        " each curve: 'legend1!' data1.csv 'blue' 'cont'" \
        " 'legend2?' data2.csv 'red' 'hyph'"
    descr_truePosIDstr = "Provide the IDs of true positive ligands" \
        " lib (format: 1-514,6001,6700-6702)"
    descr_falsePosIDstr = "Provide the IDs of true negative ligands" \
        " lib (format: 1-514,6001,6700-6702)"
    descr_yAxisName = "Name of the Y-axis in the ROC curve"
    descr_xAxisName = "Name of the X-axis in the ROC curve"
    descr_ref = "Refinement ligand(s) used on this GPCR binding pocket" \
        " refinement. Provide ligand name and ID in the following format:" \
        " lig1:328,lig2:535"
    descr_gui = "Use this flag to display plot: saves to .png by the default"
    descr_log = "Display with X-axis in log scale"
    descr_showAUC = "Use this flag to display calculated NSQ-AUC in the inset"

    # adding arguments to the parser
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("results", help=descr_results, nargs="+")
    parser.add_argument("truePosIDstr", help=descr_truePosIDstr)
    parser.add_argument("falsePosIDstr", help=descr_falsePosIDstr)
    parser.add_argument("yAxisName", help=descr_yAxisName)
    parser.add_argument("xAxisName", help=descr_xAxisName)
    parser.add_argument("-gui", action="store_true", help=descr_gui)
    parser.add_argument("-log", action="store_true", help=descr_log)
    parser.add_argument("--ref", help=descr_ref)
    parser.add_argument("-showAUC", action="store_true", help=descr_showAUC)

    # parsing args
    args = parser.parse_args()
    title = args.title
    results = args.results
    truePosIDstr = args.truePosIDstr
    falsePosIDstr = args.falsePosIDstr
    yAxisName = args.yAxisName
    xAxisName = args.xAxisName
    ref = args.ref
    gui = args.gui
    log = args.log
    showAUC = args.showAUC

    # Extrac the VS results paths and legends
    vsPaths = []
    vsLegends = []
    vsColors = []
    vsLines = []
    i = 0
    while i < len(results):
        vsLegends.append(results[i])
        vsPaths.append(results[i + 1])
        vsColors.append(results[i + 2])
        vsLines.append(results[i + 3])
        i += 4

    return title, vsLegends, vsPaths, vsColors, vsLines, \
        truePosIDstr, falsePosIDstr, xAxisName, yAxisName, \
        ref, gui, log, showAUC


if __name__ == "__main__":
    main()
