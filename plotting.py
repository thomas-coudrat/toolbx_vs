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
        Also check that the truePositive and trueNegative sets (if provided)
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

        # Loop over vsResults and keep only the ones present in the intersection
        vsIntersects = []
        for resultPath, vsResult in zip(vsPaths, allVsResults):
            # print resultPath, len(vsResult)
            vsResultIntersect = []
            for i, ligInfo in enumerate(vsResult):
                # Get only ligands that were docked in all binding pockets, and
                # also part of the library ID list provided
                if int(ligInfo[0]) in ligIDintersectSet and \
                       int(ligInfo[0]) in libraryIDlist:
                    vsResultIntersect.append(ligInfo)
            vsIntersects.append(vsResultIntersect)

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

        ligIDset = set(ligIDlist)
        intersectLig_IDs = set.intersection(ligIDset, ligIDintersectSet)
        missingLigs = sorted(list(ligIDset - intersectLig_IDs))
        if len(missingLigs) > 0:
            print(col.red + "\nWARNING: " + col.end +
                "missing IDs " + lig_type + " : " + str(missingLigs))
        ligCount = len(intersectLig_IDs)

        return ligCount


    def writePercFile(self, vsIntersect, vsDir, mode, refDict,
                      xAxisName, xAxisIDstr, xAxisIDlist, xCount,
                      yAxisName, yAxisIDstr, yAxisIDlist, yCount):
        """
        Given this VS result, and information about the ID of known actives
        in the library, write in a file the information to plot an enrichment
        curve
        """

        print(col.head + "\n\t*WRITING ENRICHMENT DATA*" + col.end)

        m = mode + "_"
        x = xAxisName + xAxisIDstr + "_"
        y = yAxisName + yAxisIDstr + "_"
        fname = os.path.basename(vsDir)

        # Create filename
        percentPath = vsDir + "/" + m + x + y + fname + ".csv"
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

            # Saving behaviour of enrichment vs. ROC curve vs. type scatter
            # All points are saved for an enrichment curve,
            # whereas only Y and X increments are saved in a ROC curve.
            # The "type" scatter collects only Y increments plus the last
            # ligand of the VS (even if it's not a Y increment)
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

        for vsPocket, vsLegend in zip(vsPockets, vsLegends):

            print "\nOpening:", vsLegend, "from path:", vsPocket

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
             xAxis, yAxis, gui, log, zoom, mode, scatterData=False):
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

        # Create a scalar map matching the data to be plotted
        scalMapPlot = self.getColorMap("spectral", plotData)

        # Drawing data on the figure
        for i, plotDatum in enumerate(plotData):
            X, Y = self.drawLine(ax, ax2, plotDatum, i, zoom, scalMapPlot, mode)

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
        # Extract plotting data
        X = plotDatum[0]
        Y = plotDatum[1]

        plotLegend = plotDatum[2]
        refPlot = plotDatum[3]

        # Set color for crystal structures, and the LDM results have their
        # colors defined by the colormap
        if i == 0 and "X-ray" in plotLegend:
            color = 'black'
            lw = 2
        elif i == 1 and "X-ray" in plotLegend:
            color = 'grey'
            lw = 2
        else:
            color = scalarMap.to_rgba(i)
            lw = 2

        # Plot this curve: scatter plot if plotting type, curves otherwise
        if mode in ("type"):

            # first plot individual ligands with scatter plot
            # ax.scatter(X, Y, linewidth=1, color='black')
            # Then draw line that starts from the origin and goes through
            # each of the scatter dots plotted above
            X = [0.0] + X
            Y = [0.0] + Y
            ax.plot(X, Y, label=plotLegend, linewidth=1, color=color)
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
            ax.text(xPos, -2, ligName, rotation=-70,
                    color=color, transform=ax.transData)

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

        # Extract information from the json file to a python dictionary
        with open(jsonFilePath, "r") as jsonRead:
            ligand_libraries_paths = json.load(jsonRead)

        # This dictionary will store ligand name (as keys) associated with
        # lists of ligand IDs
        lig_libraries_content = {}

        # For each library, open the file and collect ligand IDs
        for lib_name in ligand_libraries_paths.keys():
            ligFilePath = ligand_libraries_paths[lib_name]

            # Read the .sdf file and store ligand ID information
            ligand_IDs = []
            with open(ligFilePath) as f:
                for line in f:
                    # Get the line directly after the identifier "<lig_ID>"
                    if "<lig_ID>" in line:
                        ligID = int(next(f).strip())
                        ligand_IDs.append(ligID)

            # Storing the information collected: adding empty lists which
            # will be used later to store X and Y plotting information
            lig_libraries_content[lib_name] = (ligand_IDs, [], [])

        return lig_libraries_content


    def getColorMap(self, color_range, data):
        """
        Get a color map matching the data to be plotted
        """

        # Setting up color scheme
        cm = plt.get_cmap(color_range)
        cNorm = matplotlib.colors.Normalize(vmin=0, vmax=len(data))
        scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
        # ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(len(plotData))])

        return scalarMap


    def barPlot(self, title, enrichFactorData, pocketNames, lig_types, gui):
        """
        Plot a bar graph showing the EF0.1, EF1 and EF10 (enrichment factors)
        values for each binding pocket compared.
        Plot this enrichment factor for each of the library types compared
        """
        print(col.head + "\n\t*PLOTTING BAR GRAPH DATA*" + col.end)

        # Setting up the barchart figure
        fig_bar = plt.figure(figsize=(13, 12), dpi=100)
        ax_bar = fig_bar.add_subplot(111)

        # Default values
        groups = 3
        ind = np.arange(groups)
        width = 0.05
        separator = 0.02
        efNames = sorted(enrichFactorData.keys())
        efNumber = len(efNames)

        # Create the list of hatches for each ligand type
        # Hatches will be stored in this dictionary
        ligLibHatches = {}
        # Finite set of patterns to be associated with ligand types
        patterns = ('/', '.', 'x', '\\', '|', '-', "\\\\", 'o', '*', 'O')
        # print lig_types
        for i, lig_lib in enumerate(lig_types.keys()):
            # Create the name as "ligType (lig_count)"
            lig_lib_name = lig_lib + " (" + str(len(lig_types[lig_lib][0]))+ ")"
            # Add a dictionary value, associate it to a new pattern
            if i <= len(patterns) - 1:
                ligLibHatches[lig_lib_name] = patterns[i]
            else:
                ligLibHatches[lig_lib_name] = ""

        # Create a scalar map matching the binding pockets to be plotted
        scalMapEF = self.getColorMap("spectral", pocketNames)

        # Go through a sorted list of the pocket-ligType combinations
        allBars = []
        maxTotal = 0
        for i, efName in enumerate(sorted(enrichFactorData.keys())):
            # Get the data to be plotted
            efData = enrichFactorData[efName][0]
            efTotals = enrichFactorData[efName][1]

            # Choose the bar color: match the pocket
            curr_pocket_name = enrichFactorData[efName][2][0]
            # Loop over pocket names to get the pocket index (use 0 if pocket
            # nawas not found)
            num = 0
            for j, pocketName in enumerate(pocketNames):
                if curr_pocket_name == pocketName:
                    num = j
            # X-ray pockets use black and grey, remaining pockets use colors
            # matching the scalMapEF defined above
            if num == 0 and "X-ray" in curr_pocket_name:
                color = 'black'
            elif num == 1 and "X-ray" in curr_pocket_name:
                color = 'grey'
            else:
                color = scalMapEF.to_rgba(num)

            # Choose bar hatch: match the ligand types
            curr_lib_name = enrichFactorData[efName][2][1]
            # print lib_name
            # print ligLibHatches.keys()
            if curr_lib_name in ligLibHatches.keys():
                hatch = ligLibHatches[curr_lib_name]
                # print lib_name
            else:
                hatch = ""

            # Plotting totals in white bars
            barTots = ax_bar.bar(ind + i*(width), efTotals, width,
                                 alpha=0.5, color="white", align="center")
            # Plotting bar (with matching color and hatch)
            bars = ax_bar.bar(ind + i*(width), efData, width,
                              alpha=0.5, color=color, align="center",
                              hatch=hatch)

            # Keep the max total information to set the y limit of the final
            # plot
            maxTotalCurrent = max(efTotals[0], efTotals[1], efTotals[2])
            if maxTotalCurrent > maxTotal:
                maxTotal = maxTotalCurrent

            allBars.append(bars)

        # Setting ticks and limits
        ax_bar.set_xticks(ind + ((groups * width) + width / 2.0))
        y_axis = ax_bar.get_yaxis()
        y_axis.set_major_locator(plt.MaxNLocator(integer=True))
        ax_bar.set_xticklabels( ('EF 0.1 %', 'EF 1 %', 'EF 10 %') )
        ax_bar.tick_params(axis="both", which="major", labelsize=30)
        # Set the upperlimit at 10% more than the maximum value of the graph
        ax_bar.set_ylim(0, maxTotal * 1.10)

        # Setting legend
        fig_leg = plt.figure(figsize=(13, 12), dpi=100)
        plt.figlegend([bar[0] for bar in allBars], efNames,
                      loc="upper left", prop={'size': 30})

        # Display or save barchart and legend
        if gui:
            plt.show()
        else:
            barFile = title.replace(" ", "_") + "_bar.png"
            legFile = title.replace(" ", "_") + "_barLeg.png"
            fig_bar.savefig(barFile, bbox_inches="tight")
            fig_leg.savefig(legFile, bbox_inches="tight")


    def extractLigTypeData(self, percentPaths, vsLegends,
                           lig_types, libraryCount):
        """
        Extract enrichment factor data from each binding pocket VS data file
        File format
        """

        print(col.head + "\n\t*GETTING LIGAND TYPE DATA*" + col.end)

        # Dictionary containing enrichment factor information
        enrichFactorData = {}

        # Calculate the ligand count at each enrichment factor (0.1, 1, and
        # 10 %). This is evaluated to plot values against totals at each EF.
        totalLibTenthPercent = int(0.001 * libraryCount)
        totalLibOnePercent = int(0.01 * libraryCount)
        totalLibTenPercent = int(0.1 * libraryCount)

        # Read each binding pocket VS data file
        for percentPath, vsLegend in zip(percentPaths, vsLegends):

            print "\nOpening:", vsLegend, "from path:", percentPath

            with open(percentPath) as dataFile:
                # Open the binding pocket VS data line by line
                for line in dataFile:
                    # Get the data from the file
                    ll = line.split(",")
                    ligID = int(ll[0])
                    xPercent = float(ll[3])
                    yPercent = float(ll[4])

                    # Read each ligand type ligand ID data
                    for lib_name in lig_types.keys():

                        # Collect the list of ligand IDs that correspond to that
                        # type
                        lib_IDs = lig_types[lib_name][0]
                        # Count number of ligands in that library
                        ligCount = len(lib_IDs)

                        # Create the enrichmentFactorData data structure
                        # Dictionary with keys combining bindingPocket name and
                        # ligLibrary name. Point to two lists, each initialised
                        # to 0,0,0. For each list, values at [0], [1], [2]
                        # correspond to EF0.1, EF1 and EF10 respectively.
                        # The first list stores values at each EF level, the
                        # second list stores total hypothetical values at each
                        # EF (total ligand library numbers normalised to EF 0.1,
                        # 1 and 10).
                        libNameNum = lib_name + " (" + str(ligCount) + ")"
                        efName = vsLegend + " - " + libNameNum
                        if efName not in enrichFactorData.keys():
                            efTenth = min(ligCount, totalLibTenthPercent)
                            efOne = min(ligCount, totalLibOnePercent)
                            efTen = min(ligCount, totalLibTenPercent)
                            enrichFactorData[efName] = [[0,0,0],
                                                        [efTenth,efOne,efTen],
                                                        [vsLegend, libNameNum]]

                        # If the current ligand ID is in that list, then store
                        # its X and Y data in the lig_types dictionary
                        if ligID in lib_IDs:
                            print lib_name
                            # Store the X value
                            lig_types[lib_name][1].append(xPercent)
                            # Store the Y value
                            lig_types[lib_name][2].append(yPercent)

                            # Create a dictionary entry for the pocket-ligType
                            # combination if it doesn't exist yet. Otherwise
                            # populate the enrichment factor data (EF0.1, EF1
                            # and EF10) if the current line corresponds to a
                            # known ligand found in the current ligand type list
                            if xPercent <= 0.1:
                                enrichFactorData[efName][0][0] += 1
                                print "EF0.1", xPercent, yPercent, enrichFactorData[efName][0][0]
                            if xPercent <= 1:
                                enrichFactorData[efName][0][1] += 1
                                print "EF1", xPercent, yPercent, enrichFactorData[efName][0][1]
                            if xPercent <= 10:
                                enrichFactorData[efName][0][2] += 1
                                print "EF10", xPercent, yPercent, enrichFactorData[efName][0][2]

        scatterData = []
        for lib_name in lig_types.keys():
            scat_X = lig_types[lib_name][1]
            scat_Y = lig_types[lib_name][2]
            scatterData.append((scat_X, scat_Y, lib_name))

        return scatterData, enrichFactorData


if __name__ == "__main__":
    print("Plotting class: create a instance of the plotting object \
          to access functions ")
