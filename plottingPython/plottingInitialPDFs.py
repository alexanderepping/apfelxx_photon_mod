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
input_file = "/home/alexander/Documents/apfelxx_photon_mod/plottingPython/data_InitialPDFs.txt"
#input_file = "/home/alexander/Documents/apfelxx_photon_mod/minimizationMinuit/outputData.md"
c = 1.75
InitialPDFsType = "SAL"
showSALInitialPDFs = True

# x values to plot
#x = np.array([0.0001, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.105, 0.110, 0.115, 0.120, 0.125, 0.130, 0.135, 0.140, 0.145, 0.150, 0.155, 0.160, 0.165, 0.170, 0.175, 0.180, 0.185, 0.190, 0.195, 0.200, 0.205, 0.210, 0.215, 0.220, 0.225, 0.230, 0.235, 0.240, 0.245, 0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.285, 0.290, 0.295, 0.300, 0.305, 0.310, 0.315, 0.320, 0.325, 0.330, 0.335, 0.340, 0.345, 0.350, 0.355, 0.360, 0.365, 0.370, 0.375, 0.380, 0.385, 0.390, 0.395, 0.400, 0.405, 0.410, 0.415, 0.420, 0.425, 0.430, 0.435, 0.440, 0.445, 0.450, 0.455, 0.460, 0.465, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495, 0.500, 0.505, 0.510, 0.515, 0.520, 0.525, 0.530, 0.535, 0.540, 0.545, 0.550, 0.555, 0.560, 0.565, 0.570, 0.575, 0.580, 0.585, 0.590, 0.595, 0.600, 0.605, 0.610, 0.615, 0.620, 0.625, 0.630, 0.635, 0.640, 0.645, 0.650, 0.655, 0.660, 0.665, 0.670, 0.675, 0.680, 0.685, 0.690, 0.695, 0.700, 0.705, 0.710, 0.715, 0.720, 0.725, 0.730, 0.735, 0.740, 0.745, 0.750, 0.755, 0.760, 0.765, 0.770, 0.775, 0.780, 0.785, 0.790, 0.795, 0.800, 0.805, 0.810, 0.815, 0.820, 0.825, 0.830, 0.835, 0.840, 0.845, 0.850, 0.855, 0.860, 0.865, 0.870, 0.875, 0.880, 0.885, 0.890, 0.895, 0.900])
#x = np.exp(np.linspace(0.0001, 0.9, 100))
#print(x)

a1 = np.array([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5])
a2 = np.array([1e-1, 1e-2, 1e-3, 1e-4])#, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
x = np.sort(np.outer(a2, a1).flatten())

# order and name of particles
array_particles = ["gluon", "d", "u", "s"]

# parameters if the SAL input PDFs (taken from photon_pdfs_v3, eq. (5.3))
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


def InitialPDFgluonSAL(x, A, B, C):
    return A * x**B * (1-x)**C

def InitialPDFdownSAL(x, A_had, B_had, C_had, A_pl, B_pl):
    return A_had * x**B_had * (1-x)**C_had + 1/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))

def InitialPDFupSAL(x, A_had, B_had, C_had, A_pl, B_pl):
    return A_had * x**B_had * (1-x)**C_had + 4/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))

def InitialPDFstrangeSAL(x, K, A_had, B_had, C_had, A_pl, B_pl):
    return K * A_had * x**B_had * (1-x)**C_had + 1/9 * A_pl * x * (x**2 + (1-x)**2)/(1-B_pl*np.log(1-x))



###################
# import data
###################
startingLine = 1 # line in which # INITIALPDFS_... is written
usedInitialPDFs      = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+1, max_rows=1)
usedExperimentalData = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+3, max_rows=1)
ParametersNames      = np.genfromtxt(input_file,    delimiter=",", unpack=False, dtype='str',   skip_header=startingLine+5, max_rows=1)
Parameters           = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+7,    max_rows=1)
chi2                 = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+9,    max_rows=1)
chi2PerDP            = np.loadtxt(open(input_file), delimiter=",", unpack=False, dtype='float', skiprows=startingLine+11,   max_rows=1)



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
    subplt[3].plot(x, InitialPDFdown0(x, Parameters[6], Parameters[7], Parameters[8]))
    subplt[3].set_title("up quark InitialPDF")
    subplt[3].tick_params('x', labelbottom=True)

    # plot strange quark
    subplt[1].plot(x, InitialPDFstrange0(x, Parameters[0], Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7], Parameters[8]))
    subplt[1].set_title("strange quark InitialPDF")
    subplt[1].tick_params('x', labelbottom=True)

elif (InitialPDFsType == "SAL"):
    # plot gluon
    subplt[0].plot(x, InitialPDFgluonSAL(x, Parameters[8], Parameters[1], Parameters[2]))
    if showSALInitialPDFs:
        subplt[0].plot(x, InitialPDFgluonSAL(x, SALParameters[8], SALParameters[1], SALParameters[2]))
    subplt[0].set_title("gluon InitialPDF")
    subplt[0].tick_params('x', labelbottom=True)

    # plot down quark
    subplt[2].plot(x, InitialPDFdownSAL(x, Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]))
    if showSALInitialPDFs:
        subplt[2].plot(x, InitialPDFdownSAL(x, SALParameters[3], SALParameters[4], SALParameters[5], SALParameters[6], SALParameters[7]))
    subplt[2].set_title("down quark InitialPDF")
    subplt[2].tick_params('x', labelbottom=True)

    # plot up quark
    subplt[3].plot(x, InitialPDFupSAL(x, Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]))
    if showSALInitialPDFs:
        subplt[3].plot(x, InitialPDFupSAL(x, SALParameters[3], SALParameters[4], SALParameters[5], SALParameters[6], SALParameters[7]))
    subplt[3].set_title("up quark InitialPDF")
    subplt[3].tick_params('x', labelbottom=True)

    # plot strange quark
    subplt[1].plot(x, InitialPDFstrangeSAL(x, Parameters[0], Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]))
    if showSALInitialPDFs:
        subplt[1].plot(x, InitialPDFstrangeSAL(x, SALParameters[0], SALParameters[3], SALParameters[4], SALParameters[5], SALParameters[6], SALParameters[7]))
    subplt[1].set_title("strange quark InitialPDF")
    subplt[1].tick_params('x', labelbottom=True)

title = "InitialPDFs: " + str(usedInitialPDFs)
title += ",\nchi2: " + str(chi2) + ", chi2/NumberOfDataPoints: " + str(chi2PerDP) + ",\nexperimentalData: "
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
#plt.savefig('plotAllFlavors_LO.png', bbox_inches='tight')
plt.show()