#########################
# Imports & Definitions
#########################
from plottingPythonHelperFunctions import *
from pickle import TRUE
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
#input_file = dirName+"dataInitialPDFsSAL4VADIMHO.txt"
startingLine = 1 # line in which # INITIALPDFS_... is written, usually 1

ratioError = False
showSALInitialPDFs = False #SALInitialPDFs are at sqrt2 GeV, whereas out InitialPDFs are most likely at 1.3 GeV

save_fig = False
pltname=dirThisFile+"../plots/plotInitialPDFsSAL3LO_woutOPAL2"
#pltname=dirName+"plotInitialPDFsSAL5HOErrors"
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

if (str(usedInitialPDFs)[0:15] == "INITIALPDFS_SAL"):
    InitialPDFsType = "SAL"
else: 
    InitialPDFsType = "0"

usedExperimentalData = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+3, max_rows=1)
ParametersNames      = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+5, max_rows=1)
Parameters           = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+7,    max_rows=1)

if (str(np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+8, max_rows=1, comments="§#§")) != "## chi2:"):
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

    # check if the file includes chi2 for each experiment
    if (str(np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+12+2*numErrorParams, max_rows=1, comments="§#§")) != "## chi2/NumberOfDataPoints:"):
        additionalLines = 2
    else:
        additionalLines = 0

    chi2             = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+11+2*numErrorParams,    max_rows=1)
    chi2PerDP        = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+additionalLines+13+2*numErrorParams,   max_rows=1)
    deltaChi2        = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+2*additionalLines+15+2*numErrorParams,   max_rows=1)
    ErrorPDFs = True
else:
    # check if the file includes chi2 for each experiment
    if (str(np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+10, max_rows=1, comments="§#§")) != "## chi2/NumberOfDataPoints:"):
        additionalLines = 2
    else:
        additionalLines = 0

    chi2             = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+9,    max_rows=1)
    chi2PerDP        = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+additionalLines+11,   max_rows=1)
    ErrorPDFs = False


###################
# plotting
###################
# setting up the layout of the plot
fig, axs = plt.subplots(2, 2, sharex=True, figsize=(12*scale,7*scale))

#         gluon,    strange,  down,     up
subplt = [axs[0,0], axs[0,1], axs[1,0], axs[1,1]]

# making the Title
if ErrorPDFs:
    title = makeTitle(usedInitialPDFs, chi2, chi2PerDP, usedExperimentalData, Parameters, ParametersNames, ErrorPDFs, deltaChi2)
else:
    title = makeTitle(usedInitialPDFs, chi2, chi2PerDP, usedExperimentalData, Parameters, ParametersNames)


# plotting the SAL-type InitialPDFs
if (InitialPDFsType == "SAL"):
    for i in range(len(subplt)):
        if not ratioError: #plot data

            InitialPDF = BigDataDictionary[i]["InitialPDF"]
        
            subplt[i].plot(x, InitialPDF(x, Parameters), label="Apfel++ result")

            if ErrorPDFs: #plot ErrorPDFs
                subplt[i].fill_between(x, LowerErrorPDF(InitialPDF, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), UpperErrorPDF(InitialPDF, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), alpha=0.2)

            if showSALInitialPDFs: #plot SAL InitialPDFs
                subplt[i].plot(x, InitialPDF(x, SALParameters), label="SAL result")

        else: #plot Errro-Ratio
            subplt[i].plot(x, DeltaErrorPDF(InitialPDF, x, ErrorParametersPlus, ErrorParametersMinus)/InitialPDF(x, Parameters), label="ratio DeltaError/Apfel")

        subplt[i].set_title(BigDataDictionary[i]["Name"] + " InitialPDF")
        subplt[i].tick_params('x', labelbottom=True)
        subplt[i].legend(loc=BigDataDictionary[i]["legend loc"])

# plotting other InitialPDFs
else:
    print("Not supported InitialPDF!")
    exit()
     


plt.suptitle(title)
plt.xscale('log')
plt.xlim(left=10**(-4), right=1)
if save_fig:
    plt.savefig(pltname+".pdf", bbox_inches='tight', dpi=dpi)
    plt.savefig(pltname+".png", bbox_inches='tight', dpi=dpi)
plt.show()