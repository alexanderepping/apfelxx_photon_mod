###########################
# imports
###########################
import numpy as np
from experimentalData import *



###########################
# definitions, changeable
###########################
DataSetName = ALEPH_DATA
DataSetString = '"ALEPH"'
Energies = [9.9, 20.7, 284]
accuracy = 4


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
# main program
###########################

DataSet = makeDataSet(DataSetName)
nQ = len(Energies)

## Energies 
x = "x_data"
DSet = DataSet[x]
output = ""
for i in range(nQ):
    for j in DSet[i]:
        output = output + str(Energies[i]) + ", "
print('{' + DataSetString + ', {{"Energies", {' + output[:len(output)-2] + '}},')

## xData
x = "x_data"
DSet = DataSet[x]
output = ""
for i in range(nQ):
    for j in DSet[i]:
        output = output + str(round(j, accuracy)) + ", "
print('{"xData",    {' + output[:len(output)-2] + '}},')

## xError
x = "x_error"
DSet = DataSet[x]
output = ""
for i in range(nQ):
    for j in DSet[i]:
        output = output + str(round(j, accuracy)) + ", "
print('{"xError",   {' + output[:len(output)-2] + '}},')

## F2Gamma
x = "F2Gamma"
DSet = DataSet[x]
output = ""
for i in range(nQ):
    for j in DSet[i]:
        output = output + str(round(j, accuracy)) + ", "
print('{"F2Gamma",  {' + output[:len(output)-2] + '}},')

## yError
x = "y_error"
DSet = DataSet[x]
output = ""
for i in range(nQ):
    for j in DSet[i]:
        output = output + str(round(j, accuracy)) + ", "
print('{"yError",   {' + output[:len(output)-2] + '}}}},')