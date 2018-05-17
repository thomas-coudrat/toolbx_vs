#!/usr/bin/env python

# Classes used for plotting of ROC curves and bargraphs of VS results
#
# https://github.com/thomas-coudrat/toolbx_vs
# Thomas Coudrat <thomas.coudrat@gmail.com>

from matplotlib import pyplot as plt
import matplotlib as mpl
import math
import os
import sys
import numpy as np
import json
from sklearn.metrics import roc_curve, auc
import random

# Get matplotlib to save SVG text as text, not paths
mpl.rcParams['svg.fonttype'] = 'none'


class col:
    """
    Adding some colours to stdout
    """
    head = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    end = '\033[0m'
    BOLD = '\033[1m'
    red = '\033[31m'


class plotting:
    """
    Contains all methods for plotting VS results
    """

    # Defines the path to a log file used to write output. Initialised as False
    # it is set as a file path upon initialisation of this class. After that
    # all subsequent calls to log information use the same file path and append
    # to that file.
    log_file = False

    def __init__(self, title):
        """
        Initialises the plotting class by storing the title given to the
        executed script. Also generates a log file to which output is
        appended during script execution
        """

        # Name of the log file based on the executed script's title
        self.log_file = title.replace(" ", "_") + ".log"
        # Write string to file
        with open(self.log_file, "w") as f:
            f.write("**************************")
            f.write("\n*** LOG FILE for: " + title)
            f.write("\n**************************")

    def makeIDlist(self, stringID, blurb, printOut):
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
                    rangeID = rangeID + list(range(start, end + 1))
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
            print("\n" + blurb + " : " + str(rangeID))
        else:
            print("\n" + blurb)

        return rangeID

    def makeRefDict(self, refStr):
        """
        Get a string describing refinement ligands and their ID, and generate a
        dictionary to store that information, which is used when plotting the
        curves
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

    def intersectResults(self, vsPaths, libraryIDlist):
        """
        Read in the results provided in .csv format, and figure out the
        intersect between each of those results set based on the ligIDs.
        Also check that the truePositive and falsePositive sets (if provided)
        are fully present in the intersect set: send a WARNING if they are not
        Then return the results set containing only the intersect results for
        each set.
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
        ligIDintersectSet = set.intersection(*allLigIDs)

        # Loop over vsResults and keep only the ones present in the
        # intersection. Also identify which binding pocket does not have which
        # docked ligand information
        vsIntersects = []
        for resultPath, vsResult, ligIDs in zip(vsPaths,
                                                allVsResults,
                                                allLigIDs):
            vsResultIntersect = []
            for ligInfo in vsResult:
                ligID = int(ligInfo[0])
                # Get only ligands that were docked in all binding pockets, and
                # also part of the library ID list provided
                if ligID in ligIDintersectSet and ligID in libraryIDlist:
                    vsResultIntersect.append(ligInfo)
            # Complete docking data for each pocket
            vsIntersects.append(vsResultIntersect)

            # Identify which of the binding pocket is missing which ligand
            # print resultPath
            libraryIDset = set(libraryIDlist)
            # Get the set of ligIDs found in the total library ID set, but not
            # in the current binding pocket docked data
            notDockedLigIDlist = sorted(list(libraryIDset.difference(ligIDs)))
            # Print the binding pocket name, along with the IDs not docked
            bpName = os.path.basename(resultPath)
            if len(notDockedLigIDlist) > 0:
                print("\n" + bpName + col.blue +
                      " missing the following ligands:" + col.end + "\t" +
                      str(notDockedLigIDlist))

        print("\nNumber of intersecting docked ligands: " +
              str(len(ligIDintersectSet)))

        # Return the intersect VS data and intersect ligand ID set
        return vsIntersects, ligIDintersectSet

    def updatedLigCounts(self, ligIDintersectSet, ligIDlist, lig_type):
        """
        Get list of docked ligands kept for analysis (intersection of all
        binding pockets' results). Check which of the true positive, true
        negative and library ligands is present in that intersection set: print
        warnings when ligands are not present from these three original sets,
        and return actual counts for each.
        """

        print(col.head + "\n\t*UPDATE LIGAND COUNTS*" + col.end)

        ligIDset = set(ligIDlist)
        intersectLig_IDs = set.intersection(ligIDset, ligIDintersectSet)
        missingLigs = sorted(list(ligIDset - intersectLig_IDs))
        if len(missingLigs) > 0:
            print(col.red + "\nWARNING: " + col.end +
                  "missing IDs " + lig_type + " : " + str(missingLigs))
        ligCount = len(intersectLig_IDs)

        print("\nUpdated " + lig_type + ": " + str(ligCount))

        return ligCount

    def writePercFile(self, vsIntersect, vsPath, mode, refDict,
                      xAxisName, xAxisIDstr, xAxisIDlist, xCount,
                      yAxisName, yAxisIDstr, yAxisIDlist, yCount):
        """
        Given this VS result, and information about the ID of known actives
        in the library, write in a file the information to plot a
        ROC/enrichment curve (depending on what is supplied as xAxisIDlist,
        falsePositives or full library)
        """

        print(col.head + "\n\t*WRITING ENRICHMENT DATA*" + col.end)

        m = mode + "_"
        x = xAxisName + xAxisIDstr + "_"
        y = yAxisName + yAxisIDstr + "_"
        vsDataName = os.path.basename(vsPath)
        vsDir = os.path.dirname(vsPath)

        # Create filename
        percentPath = vsDir + "/" + m + x + y + vsDataName
        print("\n" + percentPath)
        percentDataFile = open(percentPath, "w")

        # Build the % data
        X = 0
        Y = 0
        val = 0
        for i, ligInfo in enumerate(vsIntersect):
            # print ligInfo
            ligID = int(ligInfo[0])

            # When the sorted ligID corresponds to a known, increase
            # the value of Y by 1
            if ligID in yAxisIDlist:
                # print "known ligand", ligInfo
                Y += 1

            # For each ligand in the full VS (or in the case of a ROC curve,
            # in the truePositives/Negatives), increase X
            if ligID in xAxisIDlist:
                X += 1

            # Calculate percentage X and Y
            Xpercent = (X * 100.0) / xCount
            Ypercent = (Y * 100.0) / yCount

            # Line to be saved to file
            percLine = ",".join([str(ligID), str(X), str(Y), str(Xpercent),
                                str(Ypercent)])
            # If one of the refinement ligand corresponds to that line, add
            # its name on the line
            if ligID in refDict.keys():
                ligName = refDict[ligID]
                percLine = percLine + "," + ligName + "\n"
            else:
                percLine = percLine + "\n"

            # Saving behaviour of enrichment vs. ROC curve vs. EF scatter
            # All points are saved for an enrichment curve,
            # whereas only Y and X increments are saved in a ROC curve and EF.
            """
            if mode in ("EF"):
                # Check again for presence of the current ligID in the y
                # axis list. If present then write line, otherwise don't
                if ligID in yAxisIDlist or len(vsIntersect) == i + 1:
                    percentDataFile.write(percLine)
            """
            if mode in ("enrich"):
                percentDataFile.write(percLine)
            elif mode in ("ROC", "EF"):
                if ligID in yAxisIDlist or ligID in xAxisIDlist:
                    percentDataFile.write(percLine)

        percentDataFile.close()

        return percentPath

    def extractPlotData(self, vsPockets, vsLegends, zoom):
        """
        Read the % result data files, return the data for plotting
        Also collect scattering data for plot_type: lig_types is a dictionary
        with keys representing ligand types, asssociated with lists
        representing ligand IDs.
        """

        print(col.head + "\n\t*GETTING PLOTTING DATA*" + col.end)

        # Variables that define the x and y limits for the zoomed in subplot
        xLim = 0.0
        yLim = 0.0
        plotData = []

        # List of X coordinates for reference ligands in this dataset. This is
        # used to display ligands recovered at the same X coordinate at
        # different angles on the figure and thus avoid overlap
        ref_coords = []

        for vsPocket, vsLegend in zip(vsPockets, vsLegends):

            print("\nOpening:", vsLegend, "from path:", vsPocket)

            # initialising variables
            refPlot = {}
            X = []
            Y = []

            # Read the % data file
            with open(vsPocket) as dataFile:
                for line in dataFile:
                    # Get the data from the file
                    ll = line.split(",")
                    ligID = int(ll[0])
                    xPercent = float(ll[3])
                    yPercent = float(ll[4])
                    X.append(xPercent)
                    Y.append(yPercent)

                    # Create the subplot limits
                    if xPercent <= zoom:
                        if yLim < yPercent:
                            xLim = xPercent
                            yLim = yPercent

                    # Collect the X position of the refinement ligand(s)
                    if len(ll) == 6:
                        ligName = ll[5].strip()
                        # Store the reference ligand's coordinates in a tuple
                        ligXY = (xPercent, yPercent)

                        # Reference ligand information includes coordinates
                        # and a boolean that identifies if there is more than
                        # one reference ligand on the figure that share there
                        # same coordinates

                        # Count the number of reference ligands already
                        # sharing that X coordinate, and store that number
                        # This number will be used as a multiplier to change
                        # the display parameters of the ligand name on the
                        # figure
                        refPlot[ligName] = (ligXY,
                                            ref_coords.count(xPercent))

                        # Update the reference ligand X coordinate list
                        ref_coords.append(xPercent)

            # Position 4 is initialised to None as it can hold the NSQ_AUC
            # value in the future (if calculated)
            plotData.append([X, Y, vsLegend, refPlot, None])

        return plotData, xLim, yLim

    def getAUC_NSQ(self, rocData):
        """
        Calculate the normalised square root area under the curve (NSQ_AUC) for
        each ROC curve (Y data). The Xsq are the square root of the X axis
        values. AUC is calculated using Xsq for the ROC curve (Y), the perfect
        curve (perf), and the random curve (rand).
        """

        self.log_and_print(col.head + "\n\t*CALCULATING NSQ_AUC*" + col.end)

        # Get values of first ROC curve in order to calculate random and
        # perfect curve values
        maxValY = len(rocData[0][1])
        maxValX = len(rocData[0][0])
        X = np.array(rocData[0][0])
        # print("len X", len(X))
        # Calculate random and perfect curve values
        perf = np.array([100.] * len(X))
        # print("len perf", len(perf))
        # print(perf)
        # rand = X
        rand = np.arange(0., 99.9999, 100.00/len(X))
        # print("len rand", len(rand))
        # print(rand)
        # Calculate AUCs for perfect and random curves
        # aucRand = np.trapz(rand)
        # aucPerf = np.trapz(perf)
        aucRand = auc(rand, rand)
        aucPerf = auc(rand, perf)
        # self.log_and_print("\n")
        # self.log_and_print("Random curve AUC: {:.3f}".format(aucRand))
        # self.log_and_print("Perfect curve AUC: {:.3f}".format(aucPerf))
        # self.log_and_print("\n")

        # print(rand)
        # Get square root of X values
        Xsq = np.sqrt(rand)
        # print(Xsq)

        # Calculate AUCs for random and perfect curves
        # aucRandSq = np.trapz(rand) # / 100.
        # aucPerfSq = np.trapz(perf) # / 100.
        aucRandSq = auc(Xsq, rand)
        aucPerfSq = auc(Xsq, perf)
        # Calculate and print random NSQ_AUC
        # self.log_and_print("Random curve AUC_sq: {:.3f}".format(aucRandSq))
        # self.log_and_print("Perfect curve AUC_sq: {:.3f}".format(aucPerfSq))
        # self.log_and_print("\n")

        # Calculate NSQ_AUCs
        nsqauc_rand = (100 * (aucRandSq - aucRandSq) /
                       (aucPerfSq - aucRandSq))
        randStr = "Random curve NSQ_AUC: {:.3f}".format(round(nsqauc_rand, 3))
        # self.log_and_print(randStr)
        nsqauc_perf = (100 * (aucPerfSq - aucRandSq) /
                       (aucPerfSq - aucRandSq))
        perfStr = "Perfect curve NSQ_AUC: {:.3f}".format(round(nsqauc_perf, 3))
        # self.log_and_print(perfStr)
        # self.log_and_print("\n")

        self.log_and_print("\nPocket,NSQ_AUC,AUC")

        for rocDatum in rocData:
            X = np.array(rocDatum[0])
            Y = np.array(rocDatum[1])
            legend = rocDatum[2]

            # Calculate AUC for the current curve
            # print(Y, len(Y))

            current_auc = auc(X, Y) / 100

            Xsq = np.sqrt(X)
            aucSq = auc(Xsq, Y)  # / 100.

            # Normalised square root AUC: 100 is perfect, 0 is random, negative
            # values are below random
            nsq_auc = (100 * (aucSq - aucRandSq) / (aucPerfSq - aucRandSq))

            # Append the nsq_auc value to the current rocDatum. Position 4
            # is initialised with to None
            rocDatum[4] = "{:.1f}".format(round(nsq_auc, 3))

            # print("AUC_sq:", aucSq)
            self.log_and_print("{},{},{}".format(legend.replace(" ", "_"),
                                                 nsq_auc,
                                                 current_auc))

    def plot(self, title, plotData, libraryCount, truePosCount, xLim, yLim,
             xAxis, yAxis, gui, log, zoom, mode, showAUC, scatterData=False):
        """
        Plot the data provided as argument, to draw curves
        """

        print(col.head + "\n\t*PLOTTING DATA*" + col.end)

        dpiVal = 800
        lineWidth = 4
        alphaVal = 0.8

        # Setting up the figure
        fig = plt.figure(figsize=(13, 12), dpi=dpiVal)
        ax = fig.add_subplot(111)

        # Create the ZOOMED graph, if requested
        if zoom != 0.0:
            ax2 = plt.axes([.17, .35, .2, .2])
        else:
            ax2 = None

        # Create a scalar map matching the data to be plotted
        # nonXrayPockets = [p for p in plotData if "X-ray" not in p]
        scalMapPlot = self.getColorMap("jet", plotData)

        # Drawing data on the figure
        for i, plotDatum in enumerate(plotData):
            X, Y = self.drawLine(ax, ax2, plotDatum, "black",
                                 "-", i, zoom, mode, lineWidth,
                                 alphaVal, False, showAUC)

        #                         ax, ax2, plotDatum, color, lineStyle, i,
        #                         zoom, mode, lineWidth, alphaVal,
        #                         ligVerticStyleFollow, showAUC)

        # Plot the scatter data if it was provided (only for plot_type)
        if scatterData:
            scalMapScat = self.getColorMap("jet", scatterData)
            for i, data in enumerate(scatterData):
                color = scalMapScat.to_rgba(i)
                scat_X = data[0]
                scat_Y = data[1]
                lib_name = data[2]
                ax.scatter(scat_X, scat_Y, label=lib_name, color=color)

        # Now plot random and perfect curves, get a range of X values from
        # 0 to 100, with 0.1 increments. These values are submitted to the
        # equations to get corresponding Y values
        xValues = np.arange(0, 100, 0.001)
        # Do not plot a perfect curve for ROC plots, as it won't show, being
        # along the x axis
        if mode in ("enrich", "type"):
            yPerfect = self.formulaPerfect(xValues, libraryCount, truePosCount)
            ax.plot(xValues, yPerfect, ":", color="grey",
                    alpha=alphaVal, linewidth=lineWidth)

        yRandom = self.formulaRandom(xValues)
        ax.plot(xValues, yRandom, "--", color="grey",
                alpha=alphaVal, linewidth=lineWidth)

        # Plot the RANDOM and PERFECT curves on the zoomed and main graph
        if zoom != 0.0:
            ax2.plot(X, perfect, color="grey")
            ax2.plot(X, random, ":", color="grey")
            ax2.tick_params(axis="both", which="major", labelsize=15)
            ax2.set_title("Zoom of the first " + str(zoom) + "%", fontsize=15)

        # Here axis and ticks are improved
        ax.set_xlabel(xAxis, fontsize=30)
        ax.set_ylabel(yAxis, fontsize=30)

        ax.minorticks_on()
        ax.tick_params(axis="both", which="major", labelsize=30)
        ax.set_title(title, fontsize=35, y=1.08)
        ax.legend(loc="upper left", prop={'size': 30})
        ax.axis('tight')
        # Needed when plotting scatterplots (redundant for plotting lines)
        plt.ylim(0, 100)
        plt.xlim(0, 100)

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
            fileName = title.replace(" ", "_")
            # Save svg version
            plt.savefig(fileName + ".svg", bbox_inches="tight", dpi=dpiVal)
            # Save png version
            plt.savefig(fileName + ".png", bbox_inches="tight", dpi=dpiVal)

    def plotROC(self, title, plotData, vsColors, vsLines,
                libraryCount, truePosCount, xLim, yLim, xAxis, yAxis,
                gui, log, zoom, mode, showAUC):
        """
        Plot the data provided as argument, to draw curves
        """

        print(col.head + "\n\t*PLOTTING DATA*" + col.end)

        dpiVal = 800
        lineWidth = 6
        alphaVal = 1

        # Setting up the figure
        fig = plt.figure(figsize=(13, 12), dpi=dpiVal)
        ax = fig.add_subplot(111)

        # Create the ZOOMED graph, if requested
        if zoom != 0.0:
            ax2 = plt.axes([.17, .35, .2, .2])
        else:
            ax2 = None

        # If all lines are two be drawn continuous, then vertical lines that
        # identify the ligand should all be hyphenated. Otherwise, follow line
        # styles of the VS curve.
        if all(x == "cont" for x in vsLines):
            ligVerticStyleFollow = False
        else:
            ligVerticStyleFollow = True

        # Drawing data on the figure
        for i, (plotDatum, color, line) in enumerate(zip(plotData,
                                                         vsColors,
                                                         vsLines)):
            if line == "cont":
                lineStyle = "-"
            elif line == "hyph":
                lineStyle = "--"
            elif line == "dots":
                lineStyle = ":"

            X, Y = self.drawLine(ax, ax2, plotDatum, color, lineStyle, i,
                                 zoom, mode, lineWidth, alphaVal,
                                 ligVerticStyleFollow, showAUC)

        # Now plot random and perfect curves, get a range of X values from
        # 0 to 100, with 0.1 increments. These values are submitted to the
        # equations to get corresponding Y values
        xValues = np.arange(0, 100, 0.001)

        yRandom = self.formulaRandom(xValues)
        ax.plot(xValues, yRandom, linestyle="-", color="black",
                alpha=alphaVal, linewidth=4)

        # Plot the RANDOM and PERFECT curves on the zoomed and main graph
        if zoom != 0.0:
            ax2.plot(X, perfect, color="grey")
            ax2.plot(X, random, ":", color="grey")
            ax2.tick_params(axis="both", which="major", labelsize=15)
            # ax2.set_title("Zoom of the first " + str(zoom) + "%",
            #               fontsize=15)

        # Here axis and ticks are improved
        ax.set_xlabel(xAxis, fontsize=30)
        ax.set_ylabel(yAxis, fontsize=30)

        ax.minorticks_on()
        ax.tick_params(axis="both", which="major", labelsize=30)
        # ax.set_title(title, fontsize=35, y=1.08)
        ax.legend(loc="best", prop={'size': 30},
                  borderpad=0.1, labelspacing=0.3, handletextpad=0.4)
        ax.axis('tight')
        # Needed when plotting scatterplots (redundant for plotting lines)
        plt.ylim(0, 100)
        plt.xlim(0, 100)

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
            fileName = title.replace(" ", "_")
            # Save svg version
            plt.savefig(fileName + ".svg", bbox_inches="tight", dpi=dpiVal)
            # Save png version
            plt.savefig(fileName + ".png", bbox_inches="tight", dpi=dpiVal)

    def drawLine(self, ax, ax2, plotDatum, color, lineStyle, i, zoom,
                 mode, lineWidth, alphaVal,
                 ligVerticStyleFollow, showAUC):
        """
        Draw the line corresponding to the set of data passed in arguments
        """
        # Extract plotting data
        X = plotDatum[0]
        Y = plotDatum[1]

        plotLegend = plotDatum[2]
        refPlot = plotDatum[3]

        # If NSQ_AUC value is found, then add it to the legend
        nsq_auc = plotDatum[4]
        # If showAUC is True, then add the calculated NSQ_AUC value to the
        # legend. Otherwise, show the legend as is.
        if nsq_auc is not None and showAUC:
            plotLegend = plotLegend + ": " + nsq_auc

        """
        # Set color for crystal structures, and the LDM results have their
        # colors defined by the colormap
        if i == 0 and "X-ray" in plotLegend:
            color = 'black'
        elif i == 1 and "X-ray" in plotLegend:
            color = 'grey'
        else:
            color = colorDatum
        """

        # Plot this curve: scatter plot if plotting type, curves otherwise
        if mode in ("type"):
            # first plot individual ligands with scatter plot
            # ax.scatter(X, Y, linewidth=1, color='black')
            # Then draw line that starts from the origin and goes through
            # each of the scatter dots plotted above
            # X = [0.0] + X
            # Y = [0.0] + Y
            ax.plot(X, Y, label=plotLegend,
                    linewidth=lineWidth, color=color, alpha=alphaVal)
        elif mode in ("enrich", "ROC"):
            # X = [0.0] + X
            # Y = [0.0] + Y
            ax.plot(X, Y, label=plotLegend,
                    linewidth=lineWidth, color=color, alpha=alphaVal,
                    linestyle=lineStyle)

        # Plot a blow up of the first X%
        if zoom != 0.0:
            ax2.plot(X, Y, color=color)

        # Decide whether to have vertical lines the same style as VS curves
        # or all hyphenated
        if ligVerticStyleFollow:
            pass
        else:
            lineStyle = "--"

        # Plot a vertical line for each refinement ligand
        for ligName in refPlot.keys():
            (xPos, yPos), rot_multipl = refPlot[ligName]
            ax.axvline(x=xPos, ymax=yPos/100., color=color, alpha=alphaVal,
                       linewidth=lineWidth, linestyle=lineStyle)

            # print("Reference lig:", ligName, color,
            #       xPos, 5, 70 + rot_multipl * 7)

            # Print the ligand name at its X coordinate, along the X axis
            # Ligands that share the same X coordinate are offset by a
            # different angle
            ax.text(xPos, 5, ligName, rotation=70 + rot_multipl * 7,
                    alpha=alphaVal, fontsize=25, color=color,
                    transform=ax.transData)

        return X, Y

    def writeCommand(self, title):
        """
        Write down the command that was used to exectute this script in a log
        file, at the location where the script is executed. Also write the
        current working directory at the time of execution
        """

        filename = title.replace(" ", "_") + ".sh"
        cwd = os.getcwd()
        logFile = open(filename, "w")
        # Write the directory location: this is not executed upong sh call of
        # the plot.log, but serves as information
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

    def formulaRandom(self, x):
        """
        Return the y value corresponding to random enrichment
        """
        return x

    def formulaPerfect(self, x, libraryCount, truePosCount):
        """
        Return the y values corresponding to a perfect enrichment
        """
        # Calculate the percentage value of the
        percentageXmax = ((truePosCount * 1.0 / libraryCount) * 100)

        # Calculate the slope of the line going from the origin to
        # x = (truePosCount / libraryCount) / 100, and y = 100
        slope = 100 / percentageXmax

        return np.multiply(x, slope)

    def getLigandListFromJson(self, jsonFilePath):
        """
        Get a json file as argument that contains the name and path to .sdf
        ligand libraries. Go through each of these libraries and collect all
        ligand ID present in the file. Return a dictionary of library names
        refering to lists ligand IDs.
        """

        jsonDirPath = os.path.dirname(jsonFilePath)

        # Extract information from the json file to a python dictionary
        with open(jsonFilePath, "r") as jsonRead:
            ligand_libraries_paths = json.load(jsonRead)

        # This dictionary will store ligand name (as keys) associated with
        # lists of ligand IDs
        lig_libraries_content = {}

        # For each library, open the file and collect ligand IDs
        for lib_name in ligand_libraries_paths.keys():
            ligFilePath = ligand_libraries_paths[lib_name]
            # Using the relative path to the json file, and the relative path
            # to the sdf file, get access to that sdf file.
            ligFileRelPath = os.path.join(jsonDirPath, ligFilePath)
            # Read the .sdf file and store ligand ID information
            ligand_IDs = []
            with open(ligFileRelPath) as f:
                for line in f:
                    # Get the line directly after the identifier "<lig_ID>"
                    if "<lig_ID>" in line:
                        ligID = int(next(f).strip())
                        ligand_IDs.append(ligID)

            # Storing the information collected: adding empty lists which
            # will be used later to store X and Y plotting information
            lig_count = str(len(ligand_IDs))
            lib_name = lib_name + " (" + lig_count + ")"
            lig_libraries_content[lib_name] = (ligand_IDs, [], [])

        return lig_libraries_content

    def getColorMap(self, color_range, data):
        """
        Get a color map matching the data to be plotted
        """

        # Setting up color scheme
        cm = plt.get_cmap(color_range)
        cNorm = mpl.colors.Normalize(vmin=0, vmax=len(data))
        scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=cm)
        # ax.set_color_cycle([scalarMap.to_rgba(i)
        #                    for i in range(len(plotData))])

        return scalarMap

    def barPlot(self, title, enrichFactorData, pocketNames, ef_cutoffs,
                vsColors, lig_types, gui, labelBars):
        """
        Plot a bar graph showing the EF0.1, EF1 and EF10 (enrichment factors)
        values for each binding pocket compared.
        Plot this enrichment factor for each of the library types compared
        """
        print(col.head + "\n\t*PLOTTING BAR GRAPH DATA*" + col.end)

        # Default graphic values
        dpiVal = 800
        alphaVal = 1

        # Values of EF_a, EF_b and EF_c
        EF_a = ef_cutoffs[0]
        EF_b = ef_cutoffs[1]
        EF_c = ef_cutoffs[2]

        # Setting up the barchart figure
        fig_bar = plt.figure(figsize=(30, 12), dpi=dpiVal)
        ax_bar = fig_bar.add_subplot(111)

        # Default values
        groups = 3
        ind = np.arange(groups*2, step=2.4)
        # print(ind)
        # ind = np.array([0, 2.7, 5.4])
        # ind = np.array([0, 2.4, 4.8])
        width = 0.07
        efNumber = len(enrichFactorData.keys())
        # Create a list of pockets and ligand library that matches the order
        # submitted in arguments
        efPockets = []
        for pocket in pocketNames:
            for lib in sorted(lig_types.keys()):
                efPockets.append([pocket, lib])

        # Get the count of the largest library
        # max_count = 0
        # Create the list of hatches for each ligand type
        # Hatches will be stored in this dictionary
        ligLibHatches = {}
        # Finite set of patterns to be associated with ligand types
        patterns = ('/', '*', 'x', 'O', '|', '-', "\\\\", '\\', '.', 'o')

        # print lig_types
        for i, lig_lib in enumerate(sorted(lig_types.keys())):
            # Create the name as "ligType (lig_count)"
            # current_count = len(lig_types[lig_lib][0])
            # lig_lib = lig_lib + " (" + str(current_count)+ ")"
            # Add a dictionary value, associate it to a new pattern
            if i <= len(patterns) - 1:
                ligLibHatches[lig_lib] = patterns[i]
            else:
                ligLibHatches[lig_lib] = ""
            # Update max_count to know the size of the largest library
            # if current_count > max_count:
            #    max_count = current_count

        # Go through a sorted list of the pocket-ligType combinations
        allBars = []
        # Store largest EF value
        # max_ef_val = 0
        # Loop over binding pocket/ligand library combinations in the order
        # defined by the user in arguments
        for i, efName in enumerate(efPockets):
            for efKey in enrichFactorData.keys():
                if enrichFactorData[efKey][4] == efName:
                    # Get the data to be plotted
                    efData = enrichFactorData[efKey][0]
                    # Get ligand count to be written above each bar
                    ligCountData = enrichFactorData[efKey][1]
                    # Get total count for that ligand library
                    libTotalCount = enrichFactorData[efKey][3]
                    # Choose the bar color: match the pocket
                    curr_pocket_name = enrichFactorData[efKey][4][0]
                    # Loop over pocket names to get the pocket index
                    # (use 0 if pocket was not found)
                    num = 0
                    for j, pocketName in enumerate(pocketNames):
                        if curr_pocket_name == pocketName:
                            num = j

                    # Get the color using the pocket index
                    if num < len(vsColors):
                        color = vsColors[num]
                    else:
                        color = "grey"

                    # Choose bar hatch: match the ligand types
                    curr_lib_name = enrichFactorData[efKey][4][1]
                    # print lib_name
                    # print ligLibHatches.keys()
                    if curr_lib_name in ligLibHatches.keys():
                        hatch = ligLibHatches[curr_lib_name]
                        # print lib_name
                    else:
                        hatch = ""

                    # Plotting totals in white bars
                    # barTots = ax_bar.bar(ind + i*(width), efTotals, width,
                    #                      alpha=alphaVal, color="white",
                    #                      align="center", linewidth=0)

                    # Plotting bar (with matching color and hatch)
                    bars = ax_bar.bar(ind + i*(width), efData, width,
                                      alpha=alphaVal, color=color,
                                      align="center", hatch=hatch)

                    # Keep the largest value to update figure size
                    # for ef_val in efData:
                    #    if max_ef_val < ef_val:
                    #        max_ef_val = ef_val

                    # Display values above bars, only if there was at least
                    # one bar above 0
                    # if max_ef_val != 0:
                    for j, (value,
                            ligCount,
                            bar) in enumerate(zip(efData,
                                                  ligCountData,
                                                  bars)):
                        # If no ligand of that type was found, avoid division
                        # by 0 ind + i*width + j*width +
                        x_position = bar.get_x()
                        if labelBars and ligCount > 0:
                            # Write information on ligands found over total
                            # above bar
                            ax_bar.text(x=x_position + 0.073,
                                        y=value + 0.5,
                                        s="{}/{}".format(ligCount,
                                                         libTotalCount),
                                        fontsize=17, color="black",
                                        verticalalignment="bottom",
                                        horizontalalignment="right",
                                        rotation=90)
                        # Write information about ligand library below bar
                        # Blended transform: x in data untis, y in axes
                        # fraction
                        transf = ax_bar.get_xaxis_transform()
                        ax_bar.annotate(s=efName[1].split()[0],
                                        horizontalalignment="right",
                                        fontsize=17,
                                        xy=(x_position + 0.065, -0.025),
                                        color=color,
                                        fontweight="bold",
                                        xycoords=transf)
                    allBars.append(bars)

        # Setting ticks and limits
        # ax_bar.set_title(title, fontsize=35, y=1.08)
        # Defining the text and positions for enrichment factor (EF) labels
        ticksLabels = ('EF' + str(EF_a), 'EF' + str(EF_b), 'EF' + str(EF_c))
        ticksXpos = ind + (efNumber * width)/2
        ticksYpos = np.array([-4] * 3)
        ax_bar.set_xticks(ticksXpos)
        ax_bar.set_xticklabels(labels=ticksLabels, y=-0.04)
        ax_bar.tick_params(axis="both", which="both",
                           top=False, right=False,
                           labelsize=30)
        ax_bar.spines['right'].set_visible(False)
        ax_bar.spines['top'].set_visible(False)
        ax_bar.set_ylabel(r"EF$_x$ (%)", fontsize=30)
        # Set the upperlimit at 20% more than the maximum value of the graph
        # if max_ef_val != 0:
        #    ax_bar.set_ylim(0, max_ef_val * 1.20)
        # else:
        #    ax_bar.set_ylim(bottom=0)
        # Set margins left and right of the bar groups
        ax_bar.margins(x=.1)

        # Setting legend
        # fig_leg = plt.figure(dpi=dpiVal)

        # Generate the two lists that store the custom legend information
        # This could have been done in the for loop above, but is abstracted
        # here for code simplicity
        legRects = []
        legNames = []
        # Binding-pockets legend
        for i, (pocketName, color) in enumerate(zip(pocketNames, vsColors)):
            legRects.append(plt.Rectangle((0, 0), 10, 10, facecolor=color,
                                          alpha=alphaVal))
            legNames.append(pocketName)

        # Ligand-types legend
        for hatchKey in sorted(ligLibHatches.keys()):
            legRects.append(plt.Rectangle((0, 0), 10, 10,
                            hatch=ligLibHatches[hatchKey],
                            facecolor="white"))
            legNames.append(hatchKey)
        # Create the custom figure legend
        # plt.figlegend(legRects, legNames, prop={'size': 30},
        #              loc="center")
        ax_bar.legend(legRects, legNames,  # ncol=2,"
                      loc="best", prop={"size": 30}, frameon=False)

        # Display or save barchart and legend
        if gui:
            plt.show()
        else:
            barFile = title.replace(" ", "_")
            # legFile = title.replace(" ", "_") + "_barLeg"
            # Save pdf version
            # fig_bar.savefig(barFile + ".pdf", bbox_inches="tight",
            #                 dpi=dpiVal)
            # fig_leg.savefig(legFile + ".png", bbox_inches="tight",
            #                 format="png", dpi=dpiVal)
            # Save PNG version. Neither SVG nor PDF backends produce the right
            # figure.
            fig_bar.savefig(barFile + ".png", bbox_inches="tight", dpi=dpiVal)
            # fig_leg.savefig(legFile + ".pdf", bbox_inches="tight",
            #                 format="pdf", dpi=dpiVal)

    def extractLigTypeData(self, percentPaths, vsLegends,
                           lig_types, libraryCount, ef_cutoffs):
        """
        Extract enrichment factor data from each binding pocket VS data file
        File format
        """

        print(col.head + "\n\t*GETTING LIGAND TYPE DATA*" + col.end)

        # Values of EF_a, EF_b and EF_c
        EF_a = ef_cutoffs[0]
        EF_b = ef_cutoffs[1]
        EF_c = ef_cutoffs[2]

        # Dictionary containing enrichment factor information
        enrichFactorData = {}

        # Read each binding pocket VS data file
        for percentPath, vsLegend in zip(percentPaths, vsLegends):

            # Flags used for EF
            EF_a_notReached = True
            EF_b_notReached = True
            EF_c_notReached = True

            print("\nOpening:", vsLegend, "from path:", percentPath)

            with open(percentPath) as dataFile:
                # Open the binding pocket VS data line by line
                for line in dataFile:
                    # Get the data from the file
                    ll = line.split(",")
                    ligID = int(ll[0])
                    X = int(ll[1])
                    xPercent = float(ll[3])
                    yPercent = float(ll[4])

                    # Read each ligand type ligand ID data
                    for i, lib_name in enumerate(sorted(lig_types.keys())):

                        # Collect the list of ligand IDs that correspond to
                        # that type
                        lib_IDs = lig_types[lib_name][0]
                        # Count number of ligands in that library
                        ligCount = len(lib_IDs)
                        # Get number of ligand libraries
                        libTypesCount = len(lig_types)

                        # Initialise enrichFactorData, dictionary with keys
                        # combining bindingPocket name and ligLibrary name.
                        # It stores in position 0 and 1 two lists of three
                        # elements. Each of these elements correspond to
                        # cutoffs EF_a, EF_b and EF_c, respectively. The
                        # commonly used cutoffs are 2, 5 and 10 (defined in
                        # vs_plot_ef.py). These cutoffs are refered below as
                        # X %.
                        # Description of what is in the dictionary:
                        # [0] is initiated to False, and will contain the EF
                        # (enrichment factor) values which are calculated after
                        # this loop, and used for plotting the barchart.
                        # [1] contains the true positive ligand counters at
                        # X % of the screened database (TPx)
                        # [2] contains the number of ligands from the library
                        # at X % of screened database (Nx)
                        # [3] contains the total count of ligand
                        # of that library, it is used to calculate the EF.
                        # [4] is a list storing the
                        # strings for legend plotting, binding pocket name
                        # and ligand library name.

                        # libNameNum = lib_name + " (" + str(ligCount) + ")"
                        efName = vsLegend + " - " + lib_name
                        if efName not in enrichFactorData.keys():
                            enrichFactorData[efName] = [False,
                                                        [0, 0, 0],
                                                        [0, 0, 0],
                                                        ligCount,
                                                        [vsLegend, lib_name]]

                        # Populate the library counts at EF a, b and c.
                        # If the current ligand ID is in that list, then store
                        # its X and Y data in the lig_types dictionary
                        if xPercent - EF_a <= 0 and EF_a_notReached:
                            # update library ligand count
                            enrichFactorData[efName][2][0] = X
                            # update ligand type count
                            if ligID in lib_IDs:
                                enrichFactorData[efName][1][0] += 1
                            # update found flag
                            if int(xPercent) == EF_a and \
                                    i + 1 == libTypesCount:
                                EF_a_notReached = False

                        if xPercent - EF_b <= 0 and EF_b_notReached:
                            # update library ligand count
                            enrichFactorData[efName][2][1] = X
                            # update ligand type count
                            if ligID in lib_IDs:
                                enrichFactorData[efName][1][1] += 1
                            # update found flag
                            if int(xPercent) == EF_b and \
                                    i + 1 == libTypesCount:
                                EF_b_notReached = False

                        if xPercent - EF_c <= 0 and EF_c_notReached:
                            # update library ligand count
                            enrichFactorData[efName][2][2] = X
                            # update ligand type count
                            if ligID in lib_IDs:
                                enrichFactorData[efName][1][2] += 1
                            # update found flag
                            if int(xPercent) == EF_c and \
                                    i + 1 == libTypesCount:
                                EF_c_notReached = False

        # Calculate EF values now that ligand type counts at 0.1, 1 and 10
        # have been gathered. The equation for EF used is the following:
        # EF = (TPx / Nx) / (TP / N) where TPx is the count of true positives
        # at x % of the screened library, Nx is the count of ligands at x %,
        # and TP is the total true positives in the full library, and N is the
        # count of ligands in the library.
        for efName in sorted(enrichFactorData.keys()):

            # Get the total ligand count for this ligand type, on that
            # binding pocket
            ligCount = enrichFactorData[efName][3]

            # Calculating TP / N
            ligTypeTotalPerc = float(ligCount) / float(libraryCount)

            # Get the true positive and library count at X (0.1, 1 and 10)
            ligCountsX = enrichFactorData[efName][1]
            libraryCountsX = enrichFactorData[efName][2]

            efValues = []
            for ligCountX, libraryCountX in zip(ligCountsX, libraryCountsX):
                # Calculating the numerator of EF equation: TPx / Nx
                ligTypeXPerc = float(ligCountX) / float(libraryCountX)
                # Calculate EF and append to efValues list
                efValues.append(ligTypeXPerc / ligTypeTotalPerc)

            # Adding the new efValues list to the dictionary
            enrichFactorData[efName][0] = efValues

        return enrichFactorData

    def log_and_print(self, string):
        """
        Receives a string to be printed to stdout, as well as logged in a file.
        Initialises the log file the first time it's called, then uses a class
        variable to append to that file
        """

        # Append the string to the logging file
        with open(self.log_file, "a") as f:
            f.write("\n" + string)

        # Print to stdout the string passed as argument
        print(string)


if __name__ == "__main__":
    print("Plotting class: create a instance of the plotting object \
          to access functions ")
