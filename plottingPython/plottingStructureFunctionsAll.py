##########################
# imports
###########################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from experimentalData import *



###########################
# definitions, changeable
###########################
# variables changeable by user
dirThisFile = os.path.dirname(__file__) + "/"
#input_dir = dirThisFile + "../results/Bestandsaufnahme_2022_11_13/"
input_dir = dirThisFile + "../results/Bestandsaufnahme_2023_08_31/"

ErrorPDFs           = ["SAL4"]
input_fileSAL3      = input_dir + "dataStructureFunctionsSAL3HO.txt"
input_fileSAL4      = input_dir + "dataStructureFunctionsSAL4HO.txt"
#input_fileSAL5      = input_dir + "dataStructureFunctionsSAL5HO.txt"

# All DataSets and Energies:
#       DataSets = [ALEPH1_DATA, ALEPH2_DATA, AMY_DATA, DELPHI_DATA, JADE_DATA, L3_DATA, OPAL1_DATA, OPAL2_DATA, PLUTO_DATA, TASSO_DATA, TOPAZ_DATA, TPC_DATA]
#       NamesSets= ["ALEPH1", "ALEPH2", "AMY", "DELPHI", "JADE", "L3", "OPAL1", "OPAL2", "PLUTO", "TASSO", "TOPAZ", "TPC"]
#       Energies = [1.86, 1.9, 2.4, 2.8, 3.7, 3.76, 4.3, 5.0, 5.1, 5.3, 5.9, 6.8, 7.5, 8.9, 9.0, 9.2, 9.9, 10.7, 10.8, 12.0, 14.5, 14.7, 15.3, 16.0, 17.3, 17.5, 17.8, 20.7, 23.0, 23.1, 24.0, 30.0, 45.0, 59.0, 67.2, 73.0, 80.0, 100.0, 135.0, 284.0, 390.0, 780.0]

# Used DataSets and Energies:
DataSets = [ALEPH1_DATA, ALEPH2_DATA, AMY_DATA, DELPHI_DATA, JADE_DATA, L3_DATA, OPAL1_DATA, OPAL2_less_DATA, PLUTO_DATA, TASSO_DATA, TOPAZ_DATA]
NamesSets= ["ALEPH1", "ALEPH2", "AMY", "DELPHI", "JADE", "L3", "OPAL1", "OPAL2_less", "PLUTO", "TASSO", "TOPAZ"]
#DataSets = [ALEPH1_DATA, ALEPH2_DATA]
#NamesSets= ["ALEPH1", "ALEPH2"]

# set x-axis scale to log?
setLogScale = True

# set ylim?
setYLim = True
ylim = [-0.1, 1.5]

# scaling factor for plot
c = 1.75



###########################
# definitions
###########################
input_files = {"SAL3" : input_fileSAL3, "SAL4" : input_fileSAL4}

# loading basic information on the following data, should be the same for all input files
mu2_vals   =     np.loadtxt(open(input_files[list(input_files.keys())[0]]), delimiter=",", unpack=False, skiprows=1, max_rows=1)
num_x_vals = int(np.loadtxt(open(input_files[list(input_files.keys())[0]]), delimiter=",", unpack=True,  skiprows=3, max_rows=1))

mu_vals = np.sqrt(mu2_vals)

# setting up arrays to save the data
x       = {key : np.zeros((len(mu_vals), num_x_vals)) for key in input_files.keys()}
SAL     = {key : np.zeros((len(mu_vals), num_x_vals)) for key in input_files.keys()}
SF      = {key : np.zeros((len(mu_vals), num_x_vals)) for key in input_files.keys()}
DeltaSF = {key : np.zeros((len(mu_vals), num_x_vals)) for key in input_files.keys()}



###########################
# main program
###########################

# import the data
for i_mu in range(len(mu_vals)):
    for key in input_files.keys():
        x[key][i_mu], SAL[key][i_mu], SF[key][i_mu], DeltaSF[key][i_mu]= np.loadtxt(open(input_files[key]), delimiter=",", unpack=True, skiprows=(i_mu * num_x_vals + 5), max_rows=num_x_vals)



# loop through each DataSet, meaning different experiment
for i in range(len(DataSets)):
    DataSet = DataSets[i]
    EnergiesTemp = []

    # save all unique Q-Values in EnergiesTemp
    for Q in DataSet["Q2Data"]:
        if not Q in EnergiesTemp:
            EnergiesTemp.append(Q)
    
    xData      = []
    F2Gamma    = []
    F2GammaErr = []
    
    xDataTemp = []
    F2GammaTemp = []
    F2GammaErrTemp = []

    currentQ2Value = 0

    # go through every Q2-Data-Point
    for j in range(len(DataSet["Q2Data"])):
        # if we arrive at new Q2Value
        if currentQ2Value != DataSet["Q2Data"][j]:
            currentQ2Value = DataSet["Q2Data"][j]
             
            # save the temp list from previous Q2 value as list in the non-temp lists
            if len(xDataTemp) > 0:
                xData.append(xDataTemp)
                F2Gamma.append(F2GammaTemp)
                F2GammaErr.append(F2GammaErrTemp)

            # clear the temp lists
            xDataTemp = []
            F2GammaTemp = []
            F2GammaErrTemp = []

        # save the Data in a temp list
        xDataTemp.append(DataSet["xData"][j])
        F2GammaTemp.append(DataSet["F2Gamma"][j])
        F2GammaErrTemp.append(DataSet["F2GammaErr"][j])

    # save the temp list from last Q2 value as list in the non-temp lists
    xData.append(xDataTemp)
    F2Gamma.append(F2GammaTemp)
    F2GammaErr.append(F2GammaErrTemp)

    m = 0
    while (len(EnergiesTemp)-m > 0):
        # make different figures depending on how many different Q2Values there are still left
        if (len(EnergiesTemp)-m == 1):
            fig, axs = plt.subplots(1, 1, sharex=True, figsize=(12*c,7*c))
            subplt = [axs]
        elif (len(EnergiesTemp)-m == 2):
            fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12*c,7*c))
            subplt = [axs[0], axs[1]] 
        elif (len(EnergiesTemp)-m >= 3):
            fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12*c,7*c))
            subplt = [axs[0], axs[1], axs[2]] 
        
        # plot the stuff for the different Q2Values, max 3 at a time
        for l in range(min([len(EnergiesTemp)-m,3])):
            i_mu = np.where(mu2_vals==EnergiesTemp[l+m])[0][0]

            for key in input_files.keys():
                # plot Structure Functions
                subplt[l].plot(x[key][i_mu], SF[key][i_mu], label=key)

                if key in ErrorPDFs:
                    # plot ErrorPDFs
                    subplt[l].fill_between(x[key][i_mu], SF[key][i_mu] - DeltaSF[key][i_mu], SF[key][i_mu] + DeltaSF[key][i_mu], alpha=0.2)
                else:
                    # plot ErrorPDFs (but with 0 opacity); done such that colors of ErrorPDFs always match
                    subplt[l].fill_between(x[key][i_mu], SF[key][i_mu] - DeltaSF[key][i_mu], SF[key][i_mu] + DeltaSF[key][i_mu], alpha=0)

            # Plot SAL Structure functions (are the same values for SAL3, SAL4 etc)
            subplt[l].plot(x[list(input_files.keys())[0]][i_mu], SAL[list(input_files.keys())[0]][i_mu], label="SAL, original")

            subplt[l].errorbar(xData[l+m], F2Gamma[l+m], yerr=F2GammaErr[l+m], linestyle=" ", marker="+", color="red", capsize=6, linewidth=1.65, label=NamesSets[i])

            # some stuff to change appearance of plots
            subplt[l].text(0.43, 0.9, "$F_2^{\gamma}/x$ at µ² = "+str(mu2_vals[i_mu])+" GeV²", transform=subplt[l].transAxes, size=14)

            # ylim
            if setYLim:
                subplt[l].set_ylim(bottom = ylim[0], top = ylim[1])
            
            # logscale, including ticks and xlim
            if setLogScale:
                subplt[l].set_xscale("log")
                #subplt[l].set_xlim(left = 5*10**-3, right = 1)
                subplt[l].set_xlim(left = 0.0030, right = 1)
                xticks = [0.005, 0.01, 0.05, 0.1, 0.5, 1]
                subplt[l].xaxis.set_ticks(xticks)
                #subplt[l].xaxis.set_ticklabels( ['%1.e'  % i for i in xticks] )
                subplt[l].xaxis.set_ticklabels( ['$5 \cdot 10^{-3}$', '$10^{-3}$', '$5 \cdot 10^{-2}$', '$10^{-2}$', '$5 \cdot 10^{-1}$', '$1$'] )
            # no logscale xlim
            else:
                subplt[l].set_xlim(left = 0, right = 1)

            #x-ticks on all axis
            subplt[l].tick_params('x', labelbottom=True)

        ## Find maxval and minval for top and bottom of plot
            maxval = 0
            minval = 0
            # SAL3, 4 and 5
            for key in input_files.keys():
                if maxval < max(SF[key][i_mu][1:]):
                    maxval=max(SF[key][i_mu][1:])
                if minval > min(SF[key][i_mu]):
                    minval=min(SF[key][i_mu])
            # SAL
            if maxval < max(SAL[list(input_files.keys())[0]][i_mu][1:]):
                maxval=max(SAL[list(input_files.keys())[0]][i_mu][1:])
            if minval > min(SAL[list(input_files.keys())[0]][i_mu]):
                minval=min(SAL[list(input_files.keys())[0]][i_mu])
            # Data points
            F2Plus  = np.array(F2Gamma[l+m]) + np.array(F2GammaErr[l+m])
            F2Minus = np.array(F2Gamma[l+m]) - np.array(F2GammaErr[l+m])
            if maxval < max(F2Plus):
                maxval = max(F2Plus)
            if minval > min(F2Minus):
                minval = min(F2Minus)




            subplt[l].set_ylim(bottom = minval-0.1, top = maxval+0.1)

            # labels for axis
            if l == min([len(EnergiesTemp)-m,3])-1:
                subplt[l].set_xlabel("x", size=14)
            subplt[l].set_ylabel("$F_2^{\gamma}/x$", size=14)

            if l == 0:
                subplt[l].legend(loc="upper left", ncol=2)

        #plt.legend(loc="lower right", ncol=2)
        plt.show()
        
        m += 3