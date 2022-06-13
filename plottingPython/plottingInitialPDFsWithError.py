###################
# imports
###################
from pickle import TRUE
import numpy as np
import matplotlib.pyplot as plt
import math
import sys



###################
# definitions
###################
# variables changeable by user
input_file = "/home/alexander/Uni/apfelxx_photon_mod/plottingPython/data_InitialPDFs.txt"
# input_file = "/home/alexander/Uni/apfelxx_photon_mod/plottingPython/data_InitialPDFs_Sample.txt"
# input_file = "/home/alexander/Uni/apfelxx_photon_mod/minimizationMinuit/outputDataErrorPDFsHO.md"
# input_file = "/home/alexander/Uni/apfelxx_photon_mod/minimizationMinuit/outputDataErrorPDFsLO.md"
startingLine = 1 # line in which # INITIALPDFS_... is written

ratioError = False

c = 1.75
showSALInitialPDFs = True

# save figure options
save_fig = False
pltdir = "/home/alexander/Uni/apfelxx_photon_mod/plots/mainPlots/" 
pltname = pltdir + "PlotInitialPdfs_LO_SAL3_Ks03_pto0_pl0.pdf"
dpi = 200 # default is 100

# x values to plot
#x = np.array([0.0001, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.105, 0.110, 0.115, 0.120, 0.125, 0.130, 0.135, 0.140, 0.145, 0.150, 0.155, 0.160, 0.165, 0.170, 0.175, 0.180, 0.185, 0.190, 0.195, 0.200, 0.205, 0.210, 0.215, 0.220, 0.225, 0.230, 0.235, 0.240, 0.245, 0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.285, 0.290, 0.295, 0.300, 0.305, 0.310, 0.315, 0.320, 0.325, 0.330, 0.335, 0.340, 0.345, 0.350, 0.355, 0.360, 0.365, 0.370, 0.375, 0.380, 0.385, 0.390, 0.395, 0.400, 0.405, 0.410, 0.415, 0.420, 0.425, 0.430, 0.435, 0.440, 0.445, 0.450, 0.455, 0.460, 0.465, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495, 0.500, 0.505, 0.510, 0.515, 0.520, 0.525, 0.530, 0.535, 0.540, 0.545, 0.550, 0.555, 0.560, 0.565, 0.570, 0.575, 0.580, 0.585, 0.590, 0.595, 0.600, 0.605, 0.610, 0.615, 0.620, 0.625, 0.630, 0.635, 0.640, 0.645, 0.650, 0.655, 0.660, 0.665, 0.670, 0.675, 0.680, 0.685, 0.690, 0.695, 0.700, 0.705, 0.710, 0.715, 0.720, 0.725, 0.730, 0.735, 0.740, 0.745, 0.750, 0.755, 0.760, 0.765, 0.770, 0.775, 0.780, 0.785, 0.790, 0.795, 0.800, 0.805, 0.810, 0.815, 0.820, 0.825, 0.830, 0.835, 0.840, 0.845, 0.850, 0.855, 0.860, 0.865, 0.870, 0.875, 0.880, 0.885, 0.890, 0.895, 0.900])
#x = np.exp(np.linspace(0.0001, 0.9, 100))
#print(x)

a1 = np.array([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5])
a2 = np.array([1e-1, 1e-2, 1e-3, 1e-4])#, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
x = np.sort(np.outer(a2, a1).flatten())

# order and name of particles
array_particles = ["gluon", "d", "u", "s"]

# parameters if the SAL input PDFs 
# taken from photon_pdfs_v3, eq. (5.3), same as from SAL Table1 Zeus-TR, but B+1 because there they have f and not x*f
#K_S, B_G_HAD, C_G_HAD, A_Q_HAD, B_Q_HAD, C_Q_HAD, A_Q_PL, B_Q_PL, A_G_HAD
SALParameters = [0.3, -0.57, 3, 0.065, -0.16, 1, 4.45, 1.9, 0.027] 


###################
# functions
###################
def InitialPDFgluon0(x, An, A, B):
    return An * x**A * (1-x)**B

def InitialPDFdown0(x, An, A, B):
    return An * x**A * (1-x)**B

def InitialPDFup0(x, An, A, B):
    return An * x**A * (1-x)**B

def InitialPDFstrange0(x, K, An_d, A_d, B_d, An_u, A_u, B_u):
    return K/2 * (InitialPDFdown0(x, An_d, A_d, B_d) + InitialPDFup0(x, An_u, A_u, B_u))


def InitialPDFgluonSAL(x, params):
    A = params[8]
    B = params[1]
    C = params[2]
    return A * x**B * (1-x)**C

def InitialPDFdownSAL(x, params):
    A_had = params[3]
    B_had = params[4]
    C_had = params[5]
    A_pl  = params[6]
    B_pl  = params[7]
    return A_had * x**B_had * (1-x)**C_had + 1/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))

def InitialPDFupSAL(x, params):
    A_had = params[3]
    B_had = params[4]
    C_had = params[5]
    A_pl  = params[6]
    B_pl  = params[7]
    return A_had * x**B_had * (1-x)**C_had + 4/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))

def InitialPDFstrangeSAL(x, params):
    K     = params[0]
    A_had = params[3]
    B_had = params[4]
    C_had = params[5]
    A_pl  = params[6]
    B_pl  = params[7]
    return K * A_had * x**B_had * (1-x)**C_had + 1/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))

def DeltaErrorPDF(InitialPDF, x, errorParamsPlus, errorParamsMinus):
    result = 0.0
    for k in range(len(errorParamsPlus)):
        result += (InitialPDF(x, errorParamsPlus[k]) - InitialPDF(x, errorParamsMinus[k]))**2
    result = 0.5 * np.sqrt(result)
    return result

def LowerErrorPDF(InitialPDF, x, params, errorParamsPlus, errorParamsMinus):
    return InitialPDF(x, params) - DeltaErrorPDF(InitialPDF, x, errorParamsPlus, errorParamsMinus)

def UpperErrorPDF(InitialPDF, x, params, errorParamsPlus, errorParamsMinus):
    return InitialPDF(x, params) + DeltaErrorPDF(InitialPDF, x, errorParamsPlus, errorParamsMinus)



###################
# import data
###################
if input_file == "/home/alexander/Uni/apfelxx_photon_mod/plottingPython/data_InitialPDFs.txt":
    startingLine = 1
if input_file == "/home/alexander/Uni/apfelxx_photon_mod/plottingPython/data_InitialPDFs_Sample.txt":
    startingLine = 1

usedInitialPDFs      = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+1, max_rows=1)

if (str(usedInitialPDFs)[0:15] == "INITIALPDFS_SAL"):
    InitialPDFsType = "SAL"
else: 
    InitialPDFsType = "0"

usedExperimentalData = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+3, max_rows=1)
ParametersNames      = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+5, max_rows=1)
Parameters           = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+7,    max_rows=1)

if (str(np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+8, max_rows=1)) != "## chi2:"):
    # get the number of free parameters and the number of overall parameters
    a = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+9,    max_rows=1)
    # remove the fixed variables
    a -= Parameters
    # get the number of free parameters that can contribute to error pdfs
    numErrorParams   = len(a[a != 0])
    # prepare the arrays for the error parameters
    ErrorParametersPlus  = np.zeros((numErrorParams, len(a)))
    ErrorParametersMinus = np.zeros((numErrorParams, len(a)))
    for i in range(numErrorParams):
        ErrorParametersPlus[i]  = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+9+i,    max_rows=1)
        ErrorParametersMinus[i] = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+10+numErrorParams+i,    max_rows=1)
    chi2             = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+11+2*numErrorParams,    max_rows=1)
    chi2PerDP        = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+13+2*numErrorParams,   max_rows=1)
    deltaChi2        = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+15+2*numErrorParams,   max_rows=1)
    ErrorPDFs = True
else:
    chi2             = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+9,    max_rows=1)
    chi2PerDP        = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+11,   max_rows=1)
    ErrorPDFs = False



###################
# plotting
###################
# setting up the layout of the plot
fig, axs = plt.subplots(2, 2, sharex=True, figsize=(12*c,7*c))
subplt = [axs[0,0], axs[0,1], axs[1,0], axs[1,1]]

if (InitialPDFsType == "0"):
    # plot gluon
    subplt[0].plot(x, InitialPDFgluon0(x, Parameters[9], Parameters[1], Parameters[2]))
    subplt[0].set_title("gluon InitialPDF")
    subplt[0].tick_params('x', labelbottom=True)

    # plot down quark
    subplt[2].plot(x, InitialPDFdown0(x, Parameters[3], Parameters[4], Parameters[5]))
    subplt[2].set_title("down quark InitialPDF")
    subplt[2].tick_params('x', labelbottom=True)

    # plot up quark
    subplt[3].plot(x, InitialPDFup0(x, Parameters[6], Parameters[7], Parameters[8]))
    subplt[3].set_title("up quark InitialPDF")
    subplt[3].tick_params('x', labelbottom=True)

    # plot strange quark
    subplt[1].plot(x, InitialPDFstrange0(x, Parameters[0], Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7], Parameters[8]))
    subplt[1].set_title("strange quark InitialPDF")
    subplt[1].tick_params('x', labelbottom=True)

elif (InitialPDFsType == "SAL"):
    # plot gluon
    if not ratioError:
        subplt[0].plot(x, InitialPDFgluonSAL(x, Parameters), label="Apfel++ result")
        if ErrorPDFs:
            subplt[0].plot(x, LowerErrorPDF(InitialPDFgluonSAL, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), label="Apfel++ lower error PDF", color="red")
            subplt[0].plot(x, UpperErrorPDF(InitialPDFgluonSAL, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), label="Apfel++ upper error PDF", color="red")
        if showSALInitialPDFs:
            subplt[0].plot(x, InitialPDFgluonSAL(x, SALParameters), label="SAL result")
    else:
        subplt[0].plot(x, DeltaErrorPDF(InitialPDFgluonSAL, x, ErrorParametersPlus, ErrorParametersMinus)/InitialPDFgluonSAL(x, Parameters), label="ratio DeltaError/Apfel")
    subplt[0].set_title("gluon InitialPDF")
    subplt[0].tick_params('x', labelbottom=True)
    subplt[0].legend(loc="upper right")

    # plot down quark
    if not ratioError:
        subplt[2].plot(x, InitialPDFdownSAL(x, Parameters), label="Apfel++ result")
        if ErrorPDFs:
            subplt[2].plot(x, LowerErrorPDF(InitialPDFdownSAL, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), label="Apfel++ lower error PDF", color="red")
            subplt[2].plot(x, UpperErrorPDF(InitialPDFdownSAL, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), label="Apfel++ upper error PDF", color="red")
        if showSALInitialPDFs:
            subplt[2].plot(x, InitialPDFdownSAL(x, SALParameters), label="SAL result")
    else:
        subplt[2].plot(x, DeltaErrorPDF(InitialPDFdownSAL, x, ErrorParametersPlus, ErrorParametersMinus)/InitialPDFdownSAL(x, Parameters), label="ratio DeltaError/Apfel")
    subplt[2].set_title("down quark InitialPDF")
    subplt[2].tick_params('x', labelbottom=True)
    subplt[2].legend(loc="upper right")

    # plot up quark
    if not ratioError:
        subplt[3].plot(x, InitialPDFupSAL(x, Parameters), label="Apfel++ result")
        if ErrorPDFs:
            subplt[3].plot(x, LowerErrorPDF(InitialPDFupSAL, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), label="Apfel++ lower error PDF", color="red")
            subplt[3].plot(x, UpperErrorPDF(InitialPDFupSAL, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), label="Apfel++ upper error PDF", color="red")
        if showSALInitialPDFs:
            subplt[3].plot(x, InitialPDFupSAL(x, SALParameters), label="SAL result")
    else:
        subplt[3].plot(x, DeltaErrorPDF(InitialPDFupSAL, x, ErrorParametersPlus, ErrorParametersMinus)/InitialPDFupSAL(x, Parameters), label="ratio DeltaError/Apfel")
    subplt[3].set_title("up quark InitialPDF")
    subplt[3].tick_params('x', labelbottom=True)
    subplt[3].legend(loc="upper right")

    # plot strange quark
    if not ratioError:
        subplt[1].plot(x, InitialPDFstrangeSAL(x, Parameters), label="Apfel++ result")
        if ErrorPDFs:
            subplt[1].plot(x, LowerErrorPDF(InitialPDFstrangeSAL, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), label="Apfel++ lower error PDF", color="red")
            subplt[1].plot(x, UpperErrorPDF(InitialPDFstrangeSAL, x, Parameters, ErrorParametersPlus, ErrorParametersMinus), label="Apfel++ upper error PDF", color="red")
        if showSALInitialPDFs:
            subplt[1].plot(x, InitialPDFstrangeSAL(x, SALParameters), label="SAL result")
    else:
        subplt[1].plot(x, DeltaErrorPDF(InitialPDFstrangeSAL, x, ErrorParametersPlus, ErrorParametersMinus)/InitialPDFstrangeSAL(x, Parameters), label="ratio DeltaError/Apfel")
    subplt[1].set_title("strange quark InitialPDF")
    subplt[1].tick_params('x', labelbottom=True)
    subplt[1].legend(loc="upper right")

title = "InitialPDFs: " + str(usedInitialPDFs)
title += ",\nchi2: " + str(chi2) + ", chi2/NumberOfDataPoints: " + str(chi2PerDP) 
if ErrorPDFs:
    title += ", delta chi2: " + str(deltaChi2)
title += ",\nexperimentalData: "
for text in usedExperimentalData:
    title += text + "," 
title = title[:len(title)-1] + ",\nParameters: "
for i in range(len(Parameters)):
    title += ParametersNames[i] + "=" + str(Parameters[i]) + ", " 
title = title[:len(title)-2]

plt.suptitle(title)
plt.xscale('log')
plt.xlim(left=10**(-4), right=1)
#plt.legend(loc="upper right")
if save_fig:
    plt.savefig(pltname, bbox_inches='tight', dpi=dpi)
plt.show()