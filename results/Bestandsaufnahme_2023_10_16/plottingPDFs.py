# IMPORTANT Note: only works with single mu-values!

###################
# imports
###################
import numpy as np
import matplotlib.pyplot as plt
import os
from DataSets import *


###################
# user settings
###################
dirThisFile = os.path.dirname(__file__) + "/"

# plot settings
saveFig = True
savePDF = True
showPlot = False
scale= 1.75
sizeHochkant = 15
dpi = 200

plotHochkant = True


# NOCH GRV AT SQRT2 UND GRV AT SQRT20 BERECHNEN!!!!!!


##########################
# data import & handling
##########################
rangeParticles = 6 

for DataSetsKey in DataSetsDict.keys():
    DataSets      = DataSetsDict[DataSetsKey]["Data"]
    ErrorDataSets = DataSetsDict[DataSetsKey]["Errors"]
    pltName       = dirBestandsaufnahme + DataSetsDict[DataSetsKey]["pltName"]

    xVals = np.array([])
    muVal = 0



    for Set in DataSets.keys():
        DataSets[Set]["MuVal"]    = float(np.loadtxt(open(DataSets[Set]["FilePath"]), delimiter=",", unpack=True, skiprows=1, max_rows=1))
        DataSets[Set]["NumXVals"] =   int(np.loadtxt(open(DataSets[Set]["FilePath"]), delimiter=",", unpack=True, skiprows=3, max_rows=1))

        if not "XColumn" in DataSets[Set]:
            DataSets[Set]["XColumn"] = 0 # set 0 as default

        if not "DataColumn" in DataSets[Set]:
            DataSets[Set]["DataColumn"] = 1 # set 1 as default

        if not "Label" in DataSets[Set]:
            #DataSets[Set]["Label"] = "Data Set: " + str(Set) # make default Label
            DataSets[Set]["Label"] = str(Set) # make default Label

        numXVals = DataSets[Set]["NumXVals"]

        DataSets[Set]["X"]    = np.zeros(numXVals)
        DataSets[Set]["Data"] = np.zeros((rangeParticles+1, numXVals)) # rangeParticles+1 because of Singlet

        for particle in range(rangeParticles):
            dataArrayTemp = np.loadtxt(open(DataSets[Set]["FilePath"]), delimiter=",", unpack=True, skiprows=(particle * numXVals + 5), max_rows=numXVals)
            DataSets[Set]["Data"][particle] = dataArrayTemp[DataSets[Set]["DataColumn"]]

            if particle == 0: # get x values (only once needed)
                DataSets[Set]["X"] = dataArrayTemp[DataSets[Set]["XColumn"]]
            else: # make Singlet data (gluon not included)
                DataSets[Set]["Data"][rangeParticles] += 2 * DataSets[Set]["Data"][particle]

        # check if all sets have same x vals
        if len(xVals) == 0:
            xVals = DataSets[Set]["X"]
        elif xVals.all() != DataSets[Set]["X"].all():
            sameXVals = False

        # check if all sets have same mu val
        if muVal == 0:
            muVal = DataSets[Set]["MuVal"]
        elif muVal != DataSets[Set]["MuVal"]:
            sameMuVal = False

    for Set in ErrorDataSets.keys():
        numXVals = DataSets[Set]["NumXVals"]
        ErrorDataSets[Set]["Data"] = np.zeros((rangeParticles+1, numXVals)) # rangeParticles+1 because of Singlet

        for particle in range(rangeParticles):
            dataArrayTemp = np.loadtxt(open(DataSets[Set]["FilePath"]), delimiter=",", unpack=True, skiprows=(particle * numXVals + 5), max_rows=numXVals)
            ErrorDataSets[Set]["Data"][particle] = dataArrayTemp[ErrorDataSets[Set]["DataColumn"]]

            if particle != 0: # make Singlet data (gluon not included)
                ErrorDataSets[Set]["Data"][rangeParticles] += 2 * ErrorDataSets[Set]["Data"][particle]




    ########################
    # prepare plot and plot
    ########################
    ## setting up the layout of the plot
    if DataSetsKey in ["DataSets03","DataSets04"]:
        fig, axs = plt.subplots(2, 2, sharex=True, figsize=(sizeHochkant, sizeHochkant))
        subplt = [axs[1,1], axs[0,1], axs[0,0], axs[1,0]]
        noPlotList = [4, 5, 6]
    else:
        fig, axs = plt.subplots(3, 2, sharex=True, figsize=(sizeHochkant, sizeHochkant))
        subplt = [axs[2,1], axs[0,1], axs[0,0], axs[1,0], axs[1,1], axs[2,0]]
        noPlotList = [5]

    namesPlots = ["gluon", "down quark", "up quark", "strange quark\n", "charm quark", "bottom quark", "singlet"]
    yLabelsPlots = ["g", "d", "u", "s", "c", "b", "\Sigma(x)"]

    ## placing the legend
    if DataSetsKey in ["DataSets01","DataSets02"]:
        legendLoc = 7*["upper right"]
    elif DataSetsKey in ["DataSets03","DataSets04","DataSets09"]:
        legendLoc = 7*["upper left"]
        legendLoc[0] = "upper right" #gluon
    elif DataSetsKey in ["DataSets05","DataSets06","DataSets07","DataSets08","DataSets10"]:
        legendLoc = 7*["upper left"]
        legendLoc[0] = "upper right" #gluon
        legendLoc[3] = "upper right" #strange
    else:
        legendLoc = 7*["upper right"]

    # plot the values for the imported pdfs
    i = 0
    for particle in range(7):

        # bottom quark is not plotted
        if particle not in noPlotList:
            maxval = 0

            NumberOfSet=0
            lines = ["-", "--", "-.", ":"]

            for Set in DataSets.keys():
                if maxval < max(DataSets[Set]["Data"][particle]):
                    maxval=max(DataSets[Set]["Data"][particle])

                subplt[i].plot(xVals, DataSets[Set]["Data"][particle],  label=DataSets[Set]["Label"], linewidth=2, linestyle=lines[NumberOfSet])
                # plot ErrorPDFs
                if Set in ErrorDataSets.keys():
                    subplt[i].fill_between(xVals, DataSets[Set]["Data"][particle] - ErrorDataSets[Set]["Data"][particle], DataSets[Set]["Data"][particle] + ErrorDataSets[Set]["Data"][particle], alpha=0.2)
                else: 
                    subplt[i].fill_between(xVals, DataSets[Set]["Data"][particle], DataSets[Set]["Data"][particle], alpha=0.2)
            
                NumberOfSet+=1


            if DataSetsKey in ["DataSets05","DataSets06","DataSets07","DataSets08","DataSets09","DataSets10","DataSets15","DataSets16"]:
                subplt[i].set_ylim(bottom = -0.01, top = 1.3*maxval)
            else:
                subplt[i].set_ylim(bottom = -0.01, top = 1.1*maxval)

            if DataSetsKey in ["DataSets03","DataSets04"]:
                if i == 3:
                    heightTitle = 0.93
                else: 
                    heightTitle = 0.96
            else:
                if i == 3:
                    heightTitle = 0.9
                else: 
                    heightTitle = 0.95

            #subplt[i].set_title(namesPlots[particle]+" PDF")
            #subplt[i].text(0.43, 0.9, namesPlots[particle]+" PDF", transform=subplt[i].transAxes, size=10)
            subplt[i].text(0.5, heightTitle, namesPlots[particle]+" PDF", transform=subplt[i].transAxes, size=10, ha='center')
            subplt[i].tick_params('x', labelbottom=True)
            subplt[i].legend(loc=legendLoc[particle])
            subplt[i].set_xlabel("$x$")
            subplt[i].set_ylabel("$x"+yLabelsPlots[particle]+"/\\alpha_{QED}$")
            i+=1

    plt.xscale('log')
    plt.xlim(left=10**(-4), right=1)
    if saveFig:
        print("Save "+DataSetsKey+" / "+DataSetsDict[DataSetsKey]["pltName"])
        if savePDF:
            plt.savefig(pltName+".pdf", bbox_inches='tight', dpi=dpi)
        plt.savefig(pltName+".png", bbox_inches='tight', dpi=dpi)
    if showPlot:
        plt.show()