###########################
# imports
###########################
import numpy as np
import matplotlib.pyplot as plt
import os
from experimentalData import *



###########################
# definitions, changeable
###########################
# variables changeable by user
dirThisFile = os.path.dirname(__file__) + "/"
input_file = dirThisFile+"dataStructureFunctions.txt"
#input_file = dirThisFile+"../results/Bestandsaufnahme_2022_10_06/dataStructureFunctionsSAL5HOErrors.txt"

# All DataSets and Energies:
#       DataSets = [ALEPH1_DATA, ALEPH2_DATA, AMY_DATA, DELPHI_DATA, JADE_DATA, L3_DATA, OPAL1_DATA, OPAL2_DATA, PLUTO_DATA, TASSO_DATA, TOPAZ_DATA, TPC_DATA]
#       NamesSets= ["ALEPH1", "ALEPH2", "AMY", "DELPHI", "JADE", "L3", "OPAL1", "OPAL2", "PLUTO", "TASSO", "TOPAZ", "TPC"]
#       Energies = [1.86, 1.9, 2.4, 2.8, 3.7, 3.76, 4.3, 5.0, 5.1, 5.3, 5.9, 6.8, 7.5, 8.9, 9.0, 9.2, 9.9, 10.7, 10.8, 12.0, 14.5, 14.7, 15.3, 16.0, 17.3, 17.5, 17.8, 20.7, 23.0, 23.1, 24.0, 30.0, 45.0, 59.0, 67.2, 73.0, 80.0, 100.0, 135.0, 284.0, 390.0, 780.0]

# Used DataSets and Energies:
#DataSets = [ALEPH1_DATA, ALEPH2_DATA, AMY_DATA, DELPHI_DATA, JADE_DATA, L3_DATA, OPAL1_DATA, PLUTO_DATA, TASSO_DATA, TOPAZ_DATA]
#NamesSets= ["ALEPH1", "ALEPH2", "AMY", "DELPHI", "JADE", "L3", "OPAL1", "PLUTO", "TASSO", "TOPAZ"]
DataSets = [ALEPH1_DATA, ALEPH2_DATA]
NamesSets= ["ALEPH1", "ALEPH2"]

# plot SAL Structure Functions?
plotSAL = True

# set ylim?
setYLim = False
ylim = [-0.1, 5]

# scaling factor for plot
c = 1.75
               
               
               



###########################
# functions
###########################

def makeErrorbars(arr):
    x = np.array([])
    errbar = np.array([])
    for i in range(len(arr)-1):
        errbar = np.append(errbar, (arr[i+1] - arr[i])/2. )
        x = np.append(x, arr[i]+errbar[i])
    return x, errbar

def makeDataSet(DataSet):
    for Q in range(len(DataSet["intervals"])):
        DataSet["x_data"][Q], DataSet["x_error"][Q] = makeErrorbars(DataSet["intervals"][Q])
    return DataSet



###########################
# definitions
###########################

# loading basic information on the following data
mu2_vals = np.loadtxt(open(input_file), delimiter=",", unpack=False, skiprows=1, max_rows=1)
num_x_vals = int(np.loadtxt(open(input_file), delimiter=",", unpack=True, skiprows=3, max_rows=1))

mu_vals = np.sqrt(mu2_vals)


# setting up arrays to save the data, val[index of mu] gives the array of values
x       = np.zeros((len(mu_vals), num_x_vals))
SAL     = np.zeros((len(mu_vals), num_x_vals))
SF      = np.zeros((len(mu_vals), num_x_vals))
DeltaSF = np.zeros((len(mu_vals), num_x_vals))



###########################
# main program
###########################

# import the data
for i_mu in range(len(mu_vals)):

    x[i_mu], SAL[i_mu], SF[i_mu], DeltaSF[i_mu]= np.loadtxt(open(input_file), delimiter=",", unpack=True, skiprows=(i_mu * num_x_vals + 5), max_rows=num_x_vals)
    #SAL[i_mu] *= 1/137


#    if (DataSetName == "GRV"):
#
#        subplt[i_mu].plot(x[i_mu], SF[i_mu] + DeltaSF[i_mu], label="StructureFunction + Error")
#        subplt[i_mu].plot(x[i_mu], SF[i_mu], label="StructureFunction")
#        subplt[i_mu].plot(x[i_mu], SF[i_mu] - DeltaSF[i_mu], label="StructureFunction - Error")
#
#        subplt[i_mu].set_title("F_2^gamma/x at µ²="+str(mu2_vals[i_mu])+"GeV")
#        subplt[i_mu].tick_params('x', labelbottom=True)
#
#
#    else:
#
#        lbl="Structure Function values from "+DataSetName
#        subplt[i_mu].plot(x[i_mu], SF[i_mu] + DeltaSF[i_mu], label="StructureFunction + Error")
#        subplt[i_mu].plot(x[i_mu], SF[i_mu], label="StructureFunction")
#        subplt[i_mu].plot(x[i_mu], SF[i_mu] - DeltaSF[i_mu], label="StructureFunction - Error")
#
#        subplt[i_mu].errorbar(DataSet["x_data"][i_mu], DataSet["F2Gamma"][i_mu], xerr=DataSet["x_error"][i_mu], yerr=DataSet["y_error"][i_mu], linestyle=" ", color="red", capsize=6, linewidth=1.65, label=lbl)
#
#        subplt[i_mu].set_title("F_2^gamma/x at µ²="+str(mu2_vals[i_mu])+"GeV")
#        subplt[i_mu].tick_params('x', labelbottom=True)



for i in range(len(DataSets)):
    DataSet = DataSets[i]
    EnergiesTemp = []

    for Q in DataSet["Q2Data"]:
        if not Q in EnergiesTemp:
            EnergiesTemp.append(Q)
    
    xData      = []
    F2Gamma    = []
    F2GammaErr = []
    
    n = 0
    xDataTemp = []
    F2GammaTemp = []
    F2GammaErrTemp = []
    for j in range(len(DataSet["Q2Data"])):
        if n != DataSet["Q2Data"][j]:
            n = DataSet["Q2Data"][j]
            if len(xDataTemp)>0:
                xData.append(xDataTemp)
                F2Gamma.append(F2GammaTemp)
                F2GammaErr.append(F2GammaErrTemp)
            xDataTemp = []
            F2GammaTemp = []
            F2GammaErrTemp = []

        xDataTemp.append(DataSet["xData"][j])
        F2GammaTemp.append(DataSet["F2Gamma"][j])
        F2GammaErrTemp.append(DataSet["F2GammaErr"][j])
    xData.append(xDataTemp)
    F2Gamma.append(F2GammaTemp)
    F2GammaErr.append(F2GammaErrTemp)

    m = 0
    while (len(EnergiesTemp)-m > 3):
        fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12*c,7*c))
        subplt = [axs[0], axs[1], axs[2]] 

        for l in range(3):
            i_mu = np.where(mu2_vals==EnergiesTemp[l+m])[0][0]
            lbl="StructureFunction values from " + NamesSets[i]

            subplt[l].plot(x[i_mu], SF[i_mu], label="StructureFunction, calculated")
            subplt[l].fill_between(x[i_mu], SF[i_mu] - DeltaSF[i_mu], SF[i_mu] + DeltaSF[i_mu], alpha=0.2)
            if plotSAL:
                subplt[l].plot(x[i_mu], SAL[i_mu], label="StructureFunction, SAL")

            #subplt[l].errorbar(DataSet["x_data"][i_mu], DataSet["F2Gamma"][i_mu], xerr=DataSet["x_error"][i_mu], yerr=DataSet["y_error"][i_mu], linestyle=" ", color="red", capsize=6, linewidth=1.65, label=lbl)
            subplt[l].errorbar(xData[l+m], F2Gamma[l+m], yerr=F2GammaErr[l+m], linestyle=" ", color="red", capsize=6, linewidth=1.65, label=lbl)

            subplt[l].set_title("F_2^gamma/x at µ²="+str(mu2_vals[i_mu])+"GeV")
            subplt[l].tick_params('x', labelbottom=True)
            if setYLim:
                subplt[l].set_ylim(bottom = ylim[0], top = ylim[1])

        plt.legend(loc="upper right")
        plt.show()
        
        m += 3
    
    if (len(EnergiesTemp)-m == 1):
        fig, axs = plt.subplots(1, 1, sharex=True, figsize=(12*c,7*c))
        subplt = [axs]
    elif (len(EnergiesTemp)-m == 2):
        fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12*c,7*c))
        subplt = [axs[0], axs[1]] 
    elif (len(EnergiesTemp)-m == 3):
        fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12*c,7*c))
        subplt = [axs[0], axs[1], axs[2]] 
        
    for l in range(len(EnergiesTemp)-m):
        i_mu = np.where(mu2_vals==EnergiesTemp[l+m])[0][0]
        lbl="StructureFunction values from " + NamesSets[i]

        subplt[l].plot(x[i_mu], SF[i_mu], label="StructureFunction, calculated")
        subplt[l].fill_between(x[i_mu], SF[i_mu] - DeltaSF[i_mu], SF[i_mu] + DeltaSF[i_mu], alpha=0.2)
        if plotSAL:
            subplt[l].plot(x[i_mu], SAL[i_mu], label="StructureFunction, SAL")

        #subplt[l].errorbar(DataSet["x_data"][i_mu], DataSet["F2Gamma"][i_mu], xerr=DataSet["x_error"][i_mu], yerr=DataSet["y_error"][i_mu], linestyle=" ", color="red", capsize=6, linewidth=1.65, label=lbl)
        subplt[l].errorbar(xData[l+m], F2Gamma[l+m], yerr=F2GammaErr[l+m], linestyle=" ", color="red", capsize=6, linewidth=1.65, label=lbl)

        subplt[l].set_title("F_2^gamma/x at µ²="+str(mu2_vals[i_mu])+"GeV")
        subplt[l].tick_params('x', labelbottom=True)
        if setYLim:
            subplt[l].set_ylim(bottom = ylim[0], top = ylim[1])

    plt.legend(loc="upper right")
    plt.show()