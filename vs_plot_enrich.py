#!/usr/bin/env python


# Plot an enrichment curve
#
# DEPRECIATED
#
# https://github.com/thomas-coudrat/toolbx_vs
# Thomas Coudrat <thomas.coudrat@gmail.com>

import argparse
import sys
import os
import plotting


def main():
    """
    Exectute the vs_plot_enrich script
    """

    title, vsLegends, vsPaths, \
        libraryIDstr, truePosIDstr, ref, zoom, gui, showAUC = parseArgs()

    # Define mode
    mode = "enrich"
    # Define log
    log = True

    # Creating a plotting instance for access to all methods
    p = plotting.plotting(title)

    # Get the truePosID range in list format
    truePosIDlist = p.makeIDlist(truePosIDstr,
                                 "True positive ID list",
                                 printOut=True)
    libraryIDlist = p.makeIDlist(libraryIDstr,
                                 "Library IDs (not displayed)",
                                 printOut=False)

    # Generate a dictionary containing the refinement ligands, if any
    # refinement ligand was submitted
    if ref:
        refDict = p.makeRefDict(ref)
    else:
        refDict = {}

    # Make zoom a float if it was passed as an argument, otherwise make it "0.0"
    # to have no zoomed window on the plot
    if zoom:
        zoom = float(zoom)
    elif zoom == None:
        zoom = 0.0

    # Read the results of each VS and keep only the ligIDs that are common
    # to all of them (create an interesect result list)
    vsIntersects, ligIDintersectSet = p.intersectResults(vsPaths, libraryIDlist)

    # Get updated true positive, true negative and library counts given the
    # intersect results
    truePosCount = p.updatedLigCounts(ligIDintersectSet,
                                      truePosIDlist,
                                      "true positives")
    #trueNegCount = p.updatedLigCounts(ligIDintersectSet,
    #                                  trueNegIDlist,
    #                                  "true negatives")
    libraryCount = p.updatedLigCounts(ligIDintersectSet,
                                      libraryIDlist,
                                      "whole library")

    # Calculate % of total curves for each of these (write file + return data)
    vsPockets = []
    for vsPath, vsIntersect in zip(vsPaths, vsIntersects):
        #vsDir = os.path.dirname(vsPath)
        vsPocket = p.writePercFile(vsIntersect, vsPath, mode, refDict,
                                   "library", libraryIDstr,
                                   libraryIDlist, libraryCount,
                                   "true_pos", truePosIDstr,
                                   truePosIDlist, truePosCount)

        vsPockets.append(vsPocket)

    # Extract the data from the vs percent data (in both enrichment curves and
    # ROC curves, the truePositive count would be used to draw the perfect curve
    plotData, xLim, yLim = p.extractPlotData(vsPockets, vsLegends, zoom)

    # FIX AND COMPUTE ON ONE CURVE AT A TIME, on percent vs data?
    # p.getAUC_NSQ(plotData, perfect)

    # Define title and axis names based on mode
    xAxisName = "% of ranked database (total=" + str(libraryCount) + ")"
    yAxisName = "% of known ligands found (total=" + str(truePosCount) + ")"

    # Plot the data calculated by writePercFile, and read in by extracPlotData
    p.plot(title, plotData, libraryCount, truePosCount, xLim, yLim,
           xAxisName, yAxisName, gui, log, zoom,
           mode, showAUC, scatterData=False)

    # Write the command used to execute this script into a txt file
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
    descr_libraryIDstr = "Provide the ID range of the full library screened"
    descr_truePosIDstr = "Provide the IDs of true positive ligands" \
        " lib (format: 1-514,6001,6700-6702)"
    descr_ref = "Refinement ligand(s) used on this GPCR binding pocket" \
        " refinement. Provide ligand name and ID in the following format:" \
        " lig1:328,lig2:535"
    descr_zoom = "X-axis percentage to be displayed in the zoomed subplot"
    descr_gui = "Use this flag to display plot: saves to .png by the default"
    descr_showAUC = "Use this flag to display calculated NSQ-AUC in the inset"

    # adding arguments to the parser
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("results", help=descr_results, nargs="+")
    parser.add_argument("libraryIDstr", help=descr_libraryIDstr)
    parser.add_argument("truePosIDstr", help=descr_truePosIDstr)
    parser.add_argument("--ref", help=descr_ref)
    parser.add_argument("--zoom", help=descr_zoom)
    parser.add_argument("-gui", action="store_true", help=descr_gui)
    parser.add_argument("-showAUC", action="store_true", help=descr_showAUC)

    # parsing args
    args = parser.parse_args()
    title = args.title
    results = args.results
    libraryIDstr = args.libraryIDstr
    truePosIDstr = args.truePosIDstr
    ref = args.ref
    zoom = args.zoom
    gui = args.gui
    showAUC = args.showAUC

    # Extrac the VS results paths and legends
    vsPaths = []
    vsLegends = []
    i = 0
    while i < len(results):
        vsLegends.append(results[i])
        vsPaths.append(results[i + 1])
        i += 2

    return title, vsLegends, vsPaths, \
        libraryIDstr, truePosIDstr, ref, zoom, gui, showAUC

if __name__ == "__main__":
    main()
