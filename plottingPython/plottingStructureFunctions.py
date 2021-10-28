###################
# imports
###################
import numpy as np
import matplotlib.pyplot as plt
import math
import sys



###################
# definitions, changeable
###################
# variables changeable by user
input_file = "data_StructureFunctions.txt"

# dataset to compare the Apfel++ values to
dataset = "ALEPH"

# plot either the data or the ratio of the two datasets
mode = "data"



###################
# functions
###################

def ALEPH(x, Q2):
    if (Q2 == 9.9):
        if (x >= 0.005 and x < 0.080):
            return 0.30
        elif (x >= 0.080 and x < 0.200):
            return 0.40
        elif (x >= 0.200 and x < 0.400):
            return 0.41
        elif (x >= 0.400 and x < 0.800):
            return 0.27
        else:
            return 0
    elif (Q2 == 20.7):
        if (x >= 0.009 and x < 0.120):
            return 0.36
        elif (x >= 0.120 and x < 0.270):
            return 0.34
        elif (x >= 0.270 and x < 0.500):
            return 0.56
        elif (x >= 0.500 and x < 0.890):
            return 0.45
        else:
            return 0
    elif (Q2 == 284):
        if (x >= 0.003 and x < 0.350):
            return 0.65
        elif (x >= 0.350 and x < 0.650):
            return 0.70
        elif (x >= 0.650 and x < 0.970):
            return 1.28
        else:
            return 0
    else:
        return 0

def makeXAleph(x_arr, Q2):
    aleph_arr = np.array([])
    for x in x_arr:
        aleph_arr = np.append(aleph_arr, ALEPH(x, Q2))
    return aleph_arr

def removeZero(x_arr, a_arr, b_arr):
    # remove the data points from all three arrays, if array a is zero 
    x_arr_2 = np.array([])
    a_arr_2 = np.array([])
    b_arr_2 = np.array([])
    for i in range(len(a_arr)):
        if (a_arr[i] != 0):
            x_arr_2 = np.append(x_arr_2, x_arr[i])
            a_arr_2 = np.append(a_arr_2, a_arr[i])
            b_arr_2 = np.append(b_arr_2, b_arr[i])
    return [x_arr_2, a_arr_2, b_arr_2]


###################
# definitions
###################

c = 1.75

# loading basic information on the following data
mu2_vals = np.loadtxt(open(input_file), delimiter=",", unpack=False, skiprows=1, max_rows=1)
mu_vals = np.sqrt(mu2_vals)
num_x_vals = int(np.loadtxt(open(input_file), delimiter=",", unpack=True, skiprows=3, max_rows=1))


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



###################
# main program
###################
# import the data
for i_mu in range(len(mu_vals)):

    x[i_mu], apfel[i_mu], lhapdf[i_mu], ratio[i_mu]= np.loadtxt(open(input_file), delimiter=",", unpack=True, skiprows=(i_mu * num_x_vals + 5), max_rows=num_x_vals)

    if (mode == "data"):
        subplt[i_mu].plot(x[i_mu], apfel[i_mu], label="Structure function values, PDFs calculated by Apfel++")

        if (dataset == "GRV"):
            subplt[i_mu].plot(x[i_mu], lhapdf[i_mu], label="Structure function values, PDFs given by LHAPDF")
        elif (dataset == "ALEPH"):
            data = removeZero(x[i_mu], makeXAleph(x[i_mu], mu2_vals[i_mu]), apfel[i_mu])
            subplt[i_mu].plot(data[0], data[1], label="Structure function values, given by ALEPH")

        subplt[i_mu].set_title("F_2^gamma/alpha at µ²="+str(mu2_vals[i_mu])+"GeV")
        subplt[i_mu].tick_params('x', labelbottom=True)
    elif (mode == "ratio"):
        if (dataset == "GRV"):
            subplt[i_mu].plot(x[i_mu], ratio[i_mu], label="ratio SF of Apfel++ / SF of LHAPDF")
        elif (dataset == "ALEPH"):
            data = removeZero(x[i_mu], makeXAleph(x[i_mu], mu2_vals[i_mu]), apfel[i_mu])
            subplt[i_mu].plot(data[0], data[2]/data[1], label="ratio SF of Apfel++ / SF of ALEPH")
        
        subplt[i_mu].set_title("F_2^gamma/alpha at µ²="+str(mu2_vals[i_mu])+"GeV")
        subplt[i_mu].tick_params('x', labelbottom=True)
    


#plt.xlim(left=10**(-4), right=1)
plt.legend(loc="upper right")
#plt.savefig('plotStructureFunctions_data_LO_2021_10_27.png', bbox_inches='tight')
plt.show()