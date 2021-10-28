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

## plottingStructureFunctions.py
Plots the Structure Functions, taken from data_StructureFunctions.txt. It also calculates the Structure Functions given from experimental data and compares them to the ones given by Apfel++.

## data_StructureFunctions.txt
Contains the output data by StructureFunctions.cc and is read by plottingStructureFunctions.py. The data contains the x data, the Structure Functions calculated using the PDFs that were evolved by Apfel++, the Structure Functions calculated using the LHAPDF data and the ratio of the two. 

## calculationAlphaS.py
Calculates AlphaS(Q, nf, orderQCD). Formula taken from Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01), Eq. (2.5).
