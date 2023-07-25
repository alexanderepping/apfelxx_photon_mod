# mainPrograms
Modified versions of the examples found in [vbertone/APFEL_EXAMPLES](https://github.com/vbertone/APFEL_Examples).

All references to Apfel are deleted in Evolution.cc, StructureFunctions.cc and the Makefile. Furthermore both .cc files are modified to output the name of the used LHAPDF Set and the perturbative order.

## Evolution.cc
In the file Evolution.cc an option to use Apfel++ or to only show the LHAPDF results is added. It can be changed using the constant boolean _includeApfel_. The program outputs the data to the terminal and also to a file in plottingPython/ where it can be read by the plotting program. These programs can be executed by the bashFiles/run_Evolution.cc program.

## EvolutionFlavors.cc
Derivative of Evolution.cc. It can output the values for mutliple different particles/flavors at multiple energies.The program outputs the data to the terminal and also to a file in plottingPython/ where it can be read by the plotting program. These programs can be executed by the bashFiles/run_EvolutionFlavors.cc program.

## StructureFunctions.cc
It evolves some given LHAPDFSet and then outputs it to the fitting directory in the myLHAPDF/share/ directory in the LHAPDF format. With the evolved and the not-evolved PDFs it calculates the Structure Functions. The program outputs the data to the terminal and also to a file in plottingPython/ where it can be read by the plotting program. These programs can be executed by the bashFiles/run_StructureFunctions.cc program.
