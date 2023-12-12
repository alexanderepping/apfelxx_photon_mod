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

# plot ratio firstSet/allOtherSets
# no file output
plotRatio = False

# plot ratios given in the DataSets
plotRatioAlt = False

# output file settings
writeFile = False # only useful if data from multiple files is combined
outputFileName = dirCurrent + "dataEvolvedPDFsAllHO.txt"

# plot settings
saveFig = False
savePDF = True
#pltName = dirThisFile + "plotTest"
#pltName = dirBestandsaufnahme3 + "plotEvolvedPDFs"+name+"only"
pltName = dirBestandsaufnahme3 + "plotInitialPDFs"+order
scale= 1.75
sizeHochkant = 15
dpi = 200

plotHochkant = True

# input data sets settings
# some DataSets can be found in plottingEvolvedPDFsDataSets.py

#DataSets      = {name+", 2 GeV²": data_EvolvedPDFsSqrt2_1, name+", 10 GeV²": data_EvolvedPDFsSqrt10_1}
#ErrorDataSets = {name+", 2 GeV²": data_EvolvedPDFsSqrt2_4, name+", 10 GeV²": data_EvolvedPDFsSqrt10_4}

#DataSets      = {"APFEL++, (1.3)² GeV²": data_EvolvedPDFs13_1, 
#                 "APFEL++, 10 GeV²": data_EvolvedPDFsSqrt10_1, 
#                 "GRV, (1.3)² GeV²": data_EvolvedPDFs13_2, 
#                 "GRV, 10 GeV²": data_EvolvedPDFsSqrt10_2}
#ErrorDataSets = {}

DataSets      = {"SAL5"+order+", (1.3)² GeV²": data_EvolvedPDFs13_SAL5, 
                 "SAL4"+order+", (1.3)² GeV²": data_EvolvedPDFs13_SAL4, 
                 "SAL3"+order+", (1.3)² GeV²": data_EvolvedPDFs13_SAL3}
ErrorDataSets = {}


###################
# definitions
###################
# only include the first 6 particles: gluon, down, up, strange, charm, bottom
# the next particle is usually top or singlet
rangeParticles = 6 

xVals = np.array([])
muVal = 0

sameXVals = True
sameMuVal = True



##########################
# data import & handling
##########################
for Set in DataSets.keys():
    DataSets[Set]["MuVal"]    = float(np.loadtxt(open(DataSets[Set]["FilePath"]), delimiter=",", unpack=True, skiprows=1, max_rows=1))
    DataSets[Set]["NumXVals"] = int(np.loadtxt(open(DataSets[Set]["FilePath"]), delimiter=",", unpack=True, skiprows=3, max_rows=1))

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



######################
# write data to file
######################
# only write data to file, if all sets have same x vals and mu val
if writeFile and sameXVals and sameMuVal and not plotRatio:

    OutputFile = open(outputFileName, "w")

    OutputFile.write("# mu value:\n")
    OutputFile.write(str(muVal)+"\n")
    OutputFile.write("# numXVals\n")
    OutputFile.write(str(len(xVals))+"\n")

    OutputFile.write("# x")
    for Set in DataSets.keys():
        OutputFile.write(", "+str(Set))
    OutputFile.write("\n")

    for particle in range(rangeParticles):
        for i in range(len(xVals)):
            OutputFile.write(str(xVals[i]))
            for Set in DataSets.keys():
                OutputFile.write(", "+str(DataSets[Set]["Data"][particle][i]))
            OutputFile.write("\n")

    OutputFile.close()



########################
# prepare plot and plot
########################
# setting up the layout of the plot
if plotHochkant:
    fig, axs = plt.subplots(3, 2, sharex=True, figsize=(sizeHochkant, sizeHochkant))
    subplt = [axs[2,1], axs[0,1], axs[0,0], axs[1,0], axs[1,1], axs[2,0]]
else:
    fig, axs = plt.subplots(2, 3, sharex=True, figsize=(12*scale,7*scale))
    subplt = [axs[0,2], axs[0,0], axs[1,0], axs[0,1], axs[1,1], axs[1,2]]

namesPlots = ["gluon", "down quark", "up quark", "strange quark", "charm quark", "bottom quark", "singlet"]
yLabelsPlots = ["g", "d", "u", "s", "c", "b", "\Sigma(x)"]
#legendLoc = 7*["upper right"]
legendLoc = 7*["upper left"]
legendLoc[0] = "upper right" #gluon
#legendLoc[3] = "upper right" #strange

# plot the values for the imported pdfs/ratios
i = 0
for particle in range(rangeParticles+1):

    # bottom quark is not plotted
    if (particle != 5):
        # plot the ratio:
        if plotRatio:
            CompareSet = list(DataSets.keys())[0]
            for Set in DataSets.keys():
                if Set == CompareSet:
                    subplt[i].plot(xVals, DataSets[CompareSet]["Data"][particle]/DataSets[CompareSet]["Data"][particle])
                else:
                    subplt[i].plot(xVals, DataSets[Set]["Data"][particle]/DataSets[CompareSet]["Data"][particle],  label=DataSets[Set]["Label"])
        # plot the data:
        else:
            maxval = 0
            for Set in DataSets.keys():
                if maxval < max(DataSets[Set]["Data"][particle]):
                    maxval=max(DataSets[Set]["Data"][particle])

                subplt[i].plot(xVals, DataSets[Set]["Data"][particle],  label=DataSets[Set]["Label"])
                # plot ErrorPDFs
                if Set in ErrorDataSets.keys():
                    subplt[i].fill_between(xVals, DataSets[Set]["Data"][particle] - ErrorDataSets[Set]["Data"][particle], DataSets[Set]["Data"][particle] + ErrorDataSets[Set]["Data"][particle], alpha=0.2)
                else: 
                    subplt[i].fill_between(xVals, DataSets[Set]["Data"][particle], DataSets[Set]["Data"][particle], alpha=0.2)


        if i == 4:
            subplt[i].set_ylim(bottom = -0.01, top = 0.01)
        else:
            subplt[i].set_ylim(bottom = 0, top = 1.1*maxval)


#        subplt[i].set_ylim(bottom = 0, top = 1.1*maxval)
        subplt[i].set_title(namesPlots[particle]+" PDF")
        subplt[i].tick_params('x', labelbottom=True)
        subplt[i].legend(loc=legendLoc[particle])
        subplt[i].set_xlabel("$x$")
        if plotRatio:
            subplt[i].set_ylabel("ratio to "+DataSets[CompareSet]["Label"])
        elif plotRatioAlt:
            subplt[i].set_ylabel("")
        else:
            subplt[i].set_ylabel("$x"+yLabelsPlots[particle]+"/\\alpha_{QED}$")
        i+=1

plt.xscale('log')
plt.xlim(left=10**(-4), right=1)
if saveFig:
    if plotRatio or plotRatioAlt:
        pltName += "_ratio"
    if savePDF:
        plt.savefig(pltName+".pdf", bbox_inches='tight', dpi=dpi)
    plt.savefig(pltName+".png", bbox_inches='tight', dpi=dpi)
plt.show()