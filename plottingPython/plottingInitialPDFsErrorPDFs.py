#   found out: the ErrorParametersMinus[4] is the one responsible for the strange behaviour!
#               also: it is just in the pointlike part!
#
#       PROBLEM: If B_Q_PL is negative, the pointlike part will go to infinity for some x!!!
#                   therefore we need to require B_Q_PL to be always positive!!!
#                   already done for initial params, but not for the error params!
#       This means, that we need to require, that the delta aio is not too big. Since we cant change the Hessian matrix,
#       we need to require a maximum value for DeltaChi2: ~1.23
#

#########################
# Imports & Definitions
#########################
from plottingPythonHelperFunctions import *
#from pickle import TRUE
import numpy as np
import matplotlib.pyplot as plt
import os

dirThisFile = os.path.dirname(__file__) + "/"
dirApfel = "/home/alexander/Uni/apfelxx_photon_mod/"



###################
# Change Options
###################
dirName = dirApfel+"results/Bestandsaufnahme_2022_11_13/"

#input_file = dirThisFile+"dataInitialPDFs.txt"
input_file = dirName+"dataInitialPDFsSAL5HO.txt"
# input_file = dirName+"dataInitialPDFsSAL4HO.txt"

plotWhichPart = "PL"


showSALInitialPDFs = False #SALInitialPDFs are at sqrt2 GeV, whereas out InitialPDFs are most likely at 1.3 GeV

startingLine = 1 # line in which # INITIALPDFS_... is written, usually 1
scale = 1.75
dpi = 200 # default is 100



###################
# Definitions
###################
x = np.sort(np.outer(np.array([1e-1, 1e-2, 1e-3, 1e-4]), np.array([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5])).flatten())

# order and name of particles
array_particles = ["gluon", "d", "u", "s"]

# dictionary to save which function etc belongs to which particle
BigDataDictionary = {0: {"Name": "Gluon",
                         "legend loc": "upper right",
                         "InitialPDF": InitialPDFgluonSAL},
                     1: {"Name": "Strange Quark",
                         "legend loc": "upper right",
                         "InitialPDF": InitialPDFstrangeSAL },
                     2: {"Name": "Down Quark",
                         "legend loc": "upper right",
                         "InitialPDF": InitialPDFdownSAL },
                     3: {"Name": "Up Quark",
                         "legend loc": "upper right",
                         "InitialPDF": InitialPDFupSAL }}

# parameters if the SAL input PDFs 
# taken from photon_pdfs_v3, eq. (5.3), same as from SAL Table1 Zeus-TR, but B+1 because there they have f and not x*f
#K_S, B_G_HAD, C_G_HAD, A_Q_HAD, B_Q_HAD, C_Q_HAD, A_Q_PL, B_Q_PL, A_G_HAD
SALParameters = [0.3, -0.57, 3, 0.065, -0.16, 1, 4.45, 1.9, 0.027] 



###################
# import data
###################

usedInitialPDFs      = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+1, max_rows=1)

InitialPDFsType = "SAL"

usedExperimentalData = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+3, max_rows=1)
ParametersNames      = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+5, max_rows=1)
Parameters           = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+7,    max_rows=1)

# maybe only plot parts of the PDFs
if plotWhichPart == "Had":
    Parameters[6] = 0
    Parameters[7] = 0
if plotWhichPart == "PL":
    Parameters[3] = 0
    Parameters[8] = 0

# get the number of free parameters and the number of overall parameters
a = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+9,    max_rows=1)
# remove the fixed variables (finalErrorParameters are the same as the finalParameters if the variables are fixed)
a -= Parameters
# the AGHad doesn't count as free parameter, but is changed. Therefore ignore it
aTemp = a[0:-1]
# get the number of free parameters that can contribute to error pdfs
numErrorParams   = len(aTemp[aTemp != 0])
# prepare the arrays for the error parameters
ErrorParametersPlus  = np.zeros((numErrorParams, len(a)))
ErrorParametersMinus = np.zeros((numErrorParams, len(a)))
for i in range(numErrorParams):
    ErrorParametersPlus[i]  = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+9+i,    max_rows=1)
    ErrorParametersMinus[i] = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+10+numErrorParams+i,    max_rows=1)

    # maybe only plot parts of the PDFs
    if plotWhichPart == "Had":
        ErrorParametersPlus[i][6] = 0
        ErrorParametersPlus[i][7] = 0
        ErrorParametersMinus[i][6] = 0
        ErrorParametersMinus[i][7] = 0
    if plotWhichPart == "PL":
        ErrorParametersPlus[i][3] = 0
        ErrorParametersPlus[i][8] = 0
        ErrorParametersMinus[i][3] = 0
        ErrorParametersMinus[i][8] = 0

# check if the file includes chi2 for each experiment
if (str(np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+12+2*numErrorParams, max_rows=1, comments="§#§")) != "## chi2/NumberOfDataPoints:"):
    additionalLines = 2
else:
    additionalLines = 0

chi2             = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+11+2*numErrorParams,    max_rows=1)
chi2PerDP        = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+additionalLines+13+2*numErrorParams,   max_rows=1)
deltaChi2        = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+2*additionalLines+15+2*numErrorParams,   max_rows=1)


###################
# plotting
###################
for n in range(numErrorParams+1):
    # setting up the layout of the plot
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(12*scale,7*scale))

    #         gluon,    strange,  down,     up
    subplt = [axs[0,0], axs[0,1], axs[1,0], axs[1,1]]

    # making the Title
    title = makeTitle(usedInitialPDFs, chi2, chi2PerDP, usedExperimentalData, Parameters, ParametersNames, True, deltaChi2)


    # plotting the SAL-type InitialPDFs
    for i in range(len(subplt)):
        InitialPDF = BigDataDictionary[i]["InitialPDF"]
    
        subplt[i].plot(x, InitialPDF(x, Parameters), label="Apfel++ result")

        if n == max(range(numErrorParams+1)):
            for s in range(len(ErrorParametersMinus)):
                labelMinus = "-" + str(s)
                subplt[i].plot(x, InitialPDF(x, ErrorParametersMinus[s]), label=labelMinus)

            for s in range(len(ErrorParametersPlus)):
                labelPlus = "+" + str(s)
                subplt[i].plot(x, InitialPDF(x, ErrorParametersPlus[s]), label=labelPlus)
        else:
            labelMinus = "-" + str(n)
            labelPlus = "+" + str(n)
            subplt[i].plot(x, InitialPDF(x, ErrorParametersMinus[n]), label=labelMinus)
            subplt[i].plot(x, InitialPDF(x, ErrorParametersPlus[n]), label=labelPlus)


        if showSALInitialPDFs: #plot SAL InitialPDFs
            subplt[i].plot(x, InitialPDF(x, SALParameters), label="SAL result")


        subplt[i].set_title(BigDataDictionary[i]["Name"] + " InitialPDF")
        subplt[i].tick_params('x', labelbottom=True)
        subplt[i].legend(loc=BigDataDictionary[i]["legend loc"])


    plt.suptitle(title)
    plt.xscale('log')
    plt.xlim(left=10**(-4), right=1)
    plt.show()