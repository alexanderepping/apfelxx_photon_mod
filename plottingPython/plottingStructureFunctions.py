###########################
# imports
###########################
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
from experimentalData import *



###########################
# definitions, changeable
###########################
# variables changeable by user
input_file = "data_StructureFunctions.txt"

# dataset to compare the Apfel++ values to
DataSetName = "GRV"

# plot either the data or the ratio of Apfel and GRV
mode = "ratio"

# save figure options
save_fig = False
pltname = "plotStructureFunctions_AMY_2021_10_28.png"
dpi = 200 # default is 100

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

# defining the DataSet depending on user input
if (DataSetName == "GRV"):
            mu2_vals = mu2_vals[0:3]
            mu_vals = np.sqrt(mu2_vals)
elif (DataSetName == "ALEPH"):
            DataSet = makeDataSet(ALEPH_DATA)
            mu2_vals = mu2_vals[0:3]
            mu_vals = np.sqrt(mu2_vals)
elif (DataSetName == "AMY"):
            DataSet = makeDataSet(AMY_DATA)
            mu2_vals = mu2_vals[3:6]
            mu_vals = np.sqrt(mu2_vals)


# setting up the layout of the plot
if len(mu_vals) == 1:
    fig, axs = plt.subplots(1, 1, sharex=True, figsize=(12*c,7*c))
    subplt = [axs[0]]
elif len(mu_vals) == 2:
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12*c,7*c))
    subplt = [axs[0], axs[1]] 
elif len(mu_vals) == 3:
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12*c,7*c))
    subplt = [axs[0], axs[1], axs[2]] 
elif len(mu_vals) == 4:
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(12*c,7*c))
    subplt = [axs[0,0], axs[0,1], axs[1,0], axs[1,1]] 


# setting up arrays to save the data, val[index of mu] gives the array of values
x = np.zeros((len(mu_vals), num_x_vals))
apfel = np.zeros((len(mu_vals), num_x_vals))
lhapdf = np.zeros((len(mu_vals), num_x_vals))
ratio = np.zeros((len(mu_vals), num_x_vals))



###########################
# main program
###########################
# import the data
for i_mu in range(len(mu_vals)):

    x[i_mu], apfel[i_mu], lhapdf[i_mu], ratio[i_mu]= np.loadtxt(open(input_file), delimiter=",", unpack=True, skiprows=(i_mu * num_x_vals + 5), max_rows=num_x_vals)


    if (DataSetName == "GRV"):

        if (mode == "data"):
            subplt[i_mu].plot(x[i_mu], apfel[i_mu], label="Structure function values, PDFs calculated by Apfel++")
            subplt[i_mu].plot(x[i_mu], lhapdf[i_mu], label="Structure function values, PDFs given by LHAPDF")

        elif (mode == "ratio"):
            subplt[i_mu].plot(x[i_mu], ratio[i_mu], label="Structure Function Ratio Apfel++/GRV")

        subplt[i_mu].set_title("F_2^gamma/x at µ²="+str(mu2_vals[i_mu])+"GeV")
        subplt[i_mu].tick_params('x', labelbottom=True)


    else:

        lbl="Structure Function values from "+DataSetName
        subplt[i_mu].plot(x[i_mu], apfel[i_mu], label="Structure function values, PDFs calculated by Apfel++")
        subplt[i_mu].errorbar(DataSet["x_data"][i_mu], DataSet["F2Gamma"][i_mu], xerr=DataSet["x_error"][i_mu], yerr=DataSet["y_error"][i_mu], linestyle=" ", color="red", capsize=6, elinewidth=1.65, label=lbl)

        subplt[i_mu].set_title("F_2^gamma/x at µ²="+str(mu2_vals[i_mu])+"GeV")
        subplt[i_mu].tick_params('x', labelbottom=True)


#plt.xlim(left=10**(-4), right=1)
plt.legend(loc="upper right")
if save_fig:
    plt.savefig(pltname, bbox_inches='tight', dpi=dpi)
plt.show()