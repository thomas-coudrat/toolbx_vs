#!/usr/bin/env python

from matplotlib import pyplot as plt
import matplotlib
import scipy.integrate
import math
import os, sys
import numpy as np
import json


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


    def intersectResults(self, vsPaths, truePosIDlist, trueNegIDlist,
                         ommitIDlist):
        """
        Read in the results provided in .csv format, and figure out the
        intersect between each of those results set based on the ligIDs.
        Also check that the truePositive and trueNegative sets (if provided)
        are fully present in the intersect set: send a WARNING if they are not
        Then return the results set containing only the intersect results for
        each set.
        Remove from the set returned the ligands that correspond to the ommited
        set
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
        truePosCount = self.getDockedLigandSet(truePosIDlist,
                                        "true positive IDs: ", intersectLigID)

        # Verify that trueNegatives are contained within the intersect of ligIDs
        trueNegCount = self.getDockedLigandSet(trueNegIDlist,
                                        "true negative IDs: ", intersectLigID)

        # Verify that ommits are contained within the intersect of ligIDs
        ommitCount = self.getDockedLigandSet(ommitIDlist,
                                        "ommit IDs: ", intersectLigID)

        allVsResultsIntersect = []
        # Loop over vsResults and keep only the ones present in the intersection
        for resultPath, vsResult in zip(vsPaths, allVsResults):
            # print resultPath, len(vsResult)
            vsResultIntersect = []
            for i, ligInfo in enumerate(vsResult):
                # Get only ligands that were docked in all binding pockets, and
                # also do not store the ommitted ligand information
                if int(ligInfo[0]) in intersectLigID and \
                       int(ligInfo[0]) not in ommitIDlist:
                    vsResultIntersect.append(ligInfo)
            allVsResultsIntersect.append(vsResultIntersect)

        # Return the library count as the total count - ommit count
        return allVsResultsIntersect, libraryCount - ommitCount, \
            truePosCount, trueNegCount, ommitCount


    def getDockedLigandSet(self, ligIDlist, ligType, intersectLigID):
        """
        From a list of ligands expected to be docked, and the actual list of
        all ligands docked, return an updated ligand set of actual docked
        ligands (and their count)
        """

        ligIDset = set(ligIDlist)
        intersect_ligID = set.intersection(ligIDset, intersectLigID)
        missingLigs = ligIDset - intersect_ligID
        if len(missingLigs) > 0:
            print(col.red + "\nWARNING: " + col.end +
                "missing IDs " + ligType + str(missingLigs))
        ligCount = len(intersect_ligID)

        return ligCount


    def writePercFile(self, vsIntersect, vsDir, mode, refDict,
                      xAxisName, xAxisIDstr, xAxisIDlist, xCount,
                      yAxisName, yAxisIDstr, yAxisIDlist, yCount,
                      ommitIDstr, ommitIDlist):
        """
        Given this VS result, and information about the ID of known actives
        in the library, write in a file the information to plot an enrichment
        curve
        """

        print(col.head + "\n\t*WRITING ENRICHMENT DATA*" + col.end)

        m = mode + "_"
        x = xAxisName + xAxisIDstr + "_"
        y = yAxisName + yAxisIDstr + "_"
        ommit = "ommits_" + ommitIDstr + "_"
        fname = os.path.basename(vsDir)

        # Create filename
        percentPath = vsDir + "/" + m + x + y + ommit + fname + ".csv"
        print("\n" + percentPath)
        percentDataFile = open(percentPath, "w")

        # Build the % data
        X = 0
        Y = 0
        val = 0
        for i, ligInfo in enumerate(vsIntersect):
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

                # For each ligand in the full VS (or in the case of a ROC curve,
                # in the truePositives/Negatives), increase X
                if ligID in xAxisIDlist:
                    X += 1

                # Calculate percentage X and Y
                Xpercent = (X * 100.0) / xCount
                Ypercent = (Y * 100.0) / yCount

                # Line to be saved to file
                percLine = ",".join([str(X), str(Y), str(Xpercent),
                                    str(Ypercent)])
                # If one of the refinement ligand corresponds to that line, add
                # its name on the line
                if ligID in refDict.keys():
                    ligName = refDict[ligID]
                    percLine = percLine + ligName + "\n"
                else:
                    percLine = percLine + "\n"

                # Saving behaviour of enrichment vs. ROC curve is different:
                # only the Y increments are saved in an enrichment curve,
                # whereas both Y and X increments are saved in a ROC curve.
                if mode in ("type"):
                    # Check again for presence of the current ligID in the y
                    # axis list. If present then write line, otherwise don't
                    if ligID in yAxisIDlist or len(vsIntersect) == i + 1:
                        percentDataFile.write(percLine)
                elif mode in ("enrich"):
                    percentDataFile.write(percLine)
                elif mode in ("ROC"):
                    if ligID in yAxisIDlist or ligID in xAxisIDlist:
                        percentDataFile.write(percLine)

        percentDataFile.close()

        return percentPath


    def extractPlotData(self, percentPaths, vsLegends, truePosCount, zoom,
                        mode):
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

            X = []
            Y = []

            for line in dataLines:
                # Get the data from the file
                ll = line.split(",")
                xPercent = float(ll[2])
                yPercent = float(ll[3])

                X.append(xPercent)
                Y.append(yPercent)

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

        return plotData, xLim, yLim


    def getAUC_NSQ(self, rocData, perfect):
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


    def plot(self, title, plotData, libraryCount, truePosCount, xLim, yLim,
             xAxis, yAxis, gui, log, zoom, mode):
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
            X, Y = self.drawLine(ax, ax2, plotDatum, i, zoom, scalarMap, mode)

        # Now plot random and perfect curves, get a range of X values from
        # 0 to 100, with 0.1 increments. These values are submitted to the
        # equations to get corresponding Y values
        xValues = np.arange(0, 100, 0.1)
        # Do not plot a perfect curve for ROC plots, as it won't show, being
        # along the x axis
        if mode in ("enrich", "type"):
            yPerfect = self.formulaPerfect(xValues, libraryCount, truePosCount)
            ax.plot(xValues, yPerfect, "--", color="grey")

        yRandom = self.formulaRandom(xValues)
        ax.plot(xValues, yRandom, ":", color="grey")

        # Plot the RANDOM and PERFECT curves on the zoomed and main graph
        if zoom != 0.0:
            ax2.plot(X, perfect, color="grey")
            ax2.plot(X, random, "--", color="grey")
            ax2.tick_params(axis="both", which="major", labelsize=15)
            ax2.set_title("Zoom of the first " + str(zoom) + "%", fontsize=15)

        # Here axis and ticks are improved
        ax.set_xlabel(xAxis, fontsize=30)
        ax.set_ylabel(yAxis, fontsize=30)

        ax.minorticks_on()
        ax.tick_params(axis="both", which="major", labelsize=30)
        ax.set_title(title, fontsize=30)
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
            fileName = title.replace(" ", "_") + ".png"
            plt.savefig(fileName, bbox_inches="tight")


    def drawLine(self, ax, ax2, plotDatum, i, zoom, scalarMap, mode):
        """
        Draw the line corresponding to the set of data passed in arguments
        """
        # Set color for crystal structures, and the LDM results have their
        # colors defined by the colormap
        if i == 0:
            color = 'black'
            lw = 2
        elif i == 1:
            color = 'grey'
            lw = 2
        else:
            color = scalarMap.to_rgba(i)
            lw = 2

        X = plotDatum[0]
        Y = plotDatum[1]

        plotLegend = plotDatum[2]
        refPlot = plotDatum[3]

        # Plot this curve: scatter plot if plotting type, curves otherwise
        if mode in ("type"):
            # first plot individual ligands with scatter plot
            ax.scatter(X, Y, label=plotLegend, linewidth=lw, color=color)
            # Then draw line that starts from the origin and goes through
            # each of the scatter dots plotted above
            X = [0.0] + X
            Y = [0.0] + Y
            ax.plot(X, Y, linewidth=0.2, color="grey")
        elif mode in ("enrich", "ROC"):
            X = [0.0] + X
            Y = [0.0] + Y
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


    def writeCommand(self, title):
        """
        Write down the command that was used to exectute this script in a log
        file, at the location where the script is executed. Also write the
        current working directory at the time of execution
        """

        filename = title.replace(" ", "_") + ".txt"
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


    def ligID_from_sdf(self, sdf_path):
        """
        This function takes in the path to an sdf file and returns a list of
        ligandIDs that were found in this file
        """

        # Read sdf file line by line, and store ID information
        lig_IDs = []
        with open(sdf_path) as f:
            for line in f:
                # Get the line directly after the identifier lig_ID
                if "<lig_ID>" in line:
                    lig_IDs.append(next(f).strip())

        return lig_IDs


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

        # Extrac information from the json file to a python dictionary

        print jsonFilePath
        with open(jsonFilePath, "r") as jsonRead:
            ligand_libraries_paths = json.load(jsonRead)

        # This dictionary will store ligand name (as keys) associated with
        # lists of ligand IDs
        lig_libraries_content = {}

        # For each library, open the file and collect ligand IDs
        for lib_name in ligand_libraries_paths.keys():
            ligFilePath = ligand_libraries_paths[lib_name]

            ligand_IDs = []
            with open(ligFilePath) as f:
                for line in f:
                    if "<lig_ID>" in line:
                        ligID = int(next(f).strip())
                        ligand_IDs.append(ligID)

            # Storing the information collected
            lig_libraries_content[lib_name] = ligand_IDs

        return lig_libraries_content


if __name__ == "__main__":
    p = plotting()
    list = p.ligID_from_sdf("./ADORA2A_inhib_cluster_zm.sdf")
    print(list)
