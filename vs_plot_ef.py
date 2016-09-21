#!/usr/bin/env python

import argparse
import sys
import os
import plotting


def main():
    """
    Exectute the vs_plot_enrich script
    """

    title, vsLegends, vsPaths, vsColors, \
        truePosIDstr, falsePosIDstr, ligLibsJson, \
        ref, gui = parseArgs()

    # Define mode
    mode = "EF"
    # Define zoom
    zoom = 0.0
    # Define log
    log = True
    # Define ef_cutoffs
    ef_cutoffs = [2, 5, 10]

    # Creating a plotting instance for access to all methods
    p = plotting.plotting(title)

    # Create a library ID string, combining true positive and true negative
    # strings
    libraryIDstr = truePosIDstr + "," + falsePosIDstr

    # Get the truePosID range in list format
    truePosIDlist = p.makeIDlist(truePosIDstr, "True positive ID list: ",
                                 printOut=True)
    falsePosIDlist = p.makeIDlist(falsePosIDstr, "False positive ID list: ",
                                 printOut=True)
    libraryIDlist = truePosIDlist + falsePosIDlist

    #print(len(truePosIDlist), len(falsePosIDlist), len(libraryIDlist))


    # Generate a dictionary containing the refinement ligands, if any
    # refinement ligand was submitted
    if ref:
        refDict = p.makeRefDict(ref)
    else:
        refDict = {}

    # Get ligand ID list from sdf file(s)
    # Each of the SDF file represents a ligand "type" can be chemotype,
    # pharmacology, molecularWeight, or interaction pattern.
    # The information of this ligand "type" is stored in the dictionary key,
    # which points to a .sdf path containing ligands of that "type"
    lig_types = p.getLigandListFromJson(ligLibsJson)

    # Read the results of each VS and keep only the ligIDs that are common
    # to all of them
    vsIntersects, ligIDintersectSet = p.intersectResults(vsPaths, libraryIDlist)

    # Get updated true positive, true negative and library counts given the
    # intersect results
    truePosCount = p.updatedLigCounts(ligIDintersectSet,
                                      truePosIDlist,
                                      "true positives")
    falsePosCount = p.updatedLigCounts(ligIDintersectSet,
                                      falsePosIDlist,
                                      "false positives")
    libraryCount = p.updatedLigCounts(ligIDintersectSet,
                                      libraryIDlist,
                                      "full library")

    # Calculate % of total curves for each of these (write file + return data)
    vsPockets = []
    for vsPath, vsIntersect in zip(vsPaths, vsIntersects):
        vsDir = os.path.dirname(vsPath)
        vsPocket = p.writePercFile(vsIntersect, vsDir, mode, refDict,
                                   "full_lib", libraryIDstr,
                                   libraryIDlist, libraryCount,
                                   "true_pos", truePosIDstr,
                                   truePosIDlist, truePosCount)
        vsPockets.append(vsPocket)

    # Extract the data from the vs percent data (in both enrichment curves and
    # ROC curves, the truePositive count would be used to draw the perfect curve
    #plotData, xLim, yLim = p.extractPlotData(vsPockets, vsLegends, zoom)

    # Extract data related to ligand type (plotting and barplot data)
    enrichFactorData = p.extractLigTypeData(vsPockets,
                                            vsLegends,
                                            lig_types,
                                            libraryCount,
                                            ef_cutoffs)

    #import pprint
    #pprint.pprint(enrichFactorData)
    #pprint.pprint(lig_types)

    # Plot the barplot represeting the enrochment factors (EFs) in known ligands
    # at ef_cutoffs of the screened library
    p.barPlot(title, enrichFactorData, vsLegends, ef_cutoffs,
              vsColors, lig_types, gui)

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
    descr_falsePosIDstr = "Provide the IDs of true negative ligands" \
        " lib (format: 1-514,6001,6700-6702)"
    descr_ligLibs = "JSON file containing a dictionary structure of ligand" \
        " library names (keys) pointing to library paths in .sdf"
    descr_ref = "Refinement ligand(s) used on this GPCR binding pocket" \
        " refinement. Provide ligand name and ID in the following format:" \
        " lig1:328,lig2:535"
    descr_gui = "Use this flag to display plot: saves to .png by the default"

    # adding arguments to the parser
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("title", help=descr_title)
    parser.add_argument("results", help=descr_results, nargs="+")
    parser.add_argument("truePosIDstr", help=descr_truePosIDstr)
    parser.add_argument("falsePosIDstr", help=descr_falsePosIDstr)
    parser.add_argument("ligLibs", help=descr_ligLibs)
    parser.add_argument("--ref", help=descr_ref)
    parser.add_argument("-gui", action="store_true", help=descr_gui)

    # parsing args
    args = parser.parse_args()
    title = args.title
    results = args.results
    truePosIDstr = args.truePosIDstr
    falsePosIDstr = args.falsePosIDstr
    ligLibsJson = args.ligLibs
    ref = args.ref
    gui = args.gui

    # Extrac the VS results paths and legends
    vsPaths = []
    vsLegends = []
    vsColors = []
    i = 0
    while i < len(results):
        vsLegends.append(results[i])
        vsPaths.append(results[i + 1])
        vsColors.append(results[i + 2])
        i += 3

    return title, vsLegends, vsPaths, vsColors, \
        truePosIDstr, falsePosIDstr, ligLibsJson, \
        ref, gui

if __name__ == "__main__":
    main()
