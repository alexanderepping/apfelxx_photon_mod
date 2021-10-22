# plottingPython
Contains various python files to plot the data collected from the Evolution.cc/EvolutionFlavors.cc program. 

Also contains a program to make an array of logarithmic values and compare Vadims and my data.

## plottingEvolution.py
Plots the values from apfel and lhapdf, taken from data_EvolutionFlavors.txt.

## data_Evolution.txt
Contains the output data by Evolution.cc and is read by plottingEvolution.py. 

## plottingEvolutionFlavors.py
Plots the values from apfel and lhapdf, taken from data_EvolutionFlavors.txt. Other than plottingEvolution.py, this program plots the values for multiple particles at multpile energies.

## data_EvolutionFlavors.txt
Contains the output data by EvolutionFlavors.cc and is read by plottingEvolutionFlavors.py. 

## LogArray.py
Calculates an array of logartihmic values, sorted from smallest to biggest, and saves it to the x_data.txt file.

## x_data.txt
Contains an array of logarithmic values sorted from smallest to biggest, calculated by LogArray.py.

## compareWithVadim.py
Quick program to compare Vadims and my data. Not very useful anymore.

## GRV_param_Q2_....txt
Data compared by the compareWithVadim.py program. The numbers after Q2 are the squared energies at which the values are taken. The values are taken for the up quark.
