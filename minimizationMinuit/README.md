# minimizationMinuit
Directory contains files related to the fitting of the Initial PDFs to the Structure Function data from several experiments (of course after evolving the Initial PDFs and calculating the Structure Functions from them) by minimization of the chi-Square using the MinuitCpp library. 

## files
### minimizationMinuit.h
File that `#include`s all the used header files.

### minimizationStructureFunctions.cc 
Main program to run for the minimization/fitting. It outputs the resuls to ../plottingPython/data_InitialPDFs.txt and to the terminal.

### configMinuit.h
File that contains configurations for all the other files like which Initial PDFs are used, some parameters etc. Can be changed by user.

### experimentalData.h
File that contains the experimental data. Experimental data can be added by user. To help that, see plottingPython/calculateExperimentalData.py

### StructureFunctionFcn.cc/.h
File that contains the class that calculates the StructureFunctions, takes the experimental data and returns the chi-square. It also contains the InitialPDFs and functions to help calculating them.

### Makefile
Not ready to use yet. See "running the program" below.

## running the program 
### compile, execute and time the runtime:
```time g++ -std=c++17 -I/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/include -I/usr/include/eigen3 -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc experimentalData.h ErrorPDFs.cc -L/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF -L/home/alexander/Uni/minuit-cpp/build/lib/ -lminuit-cpp && time ./minimizationStructureFunctions```

### compile and execute:
```g++ -std=c++17 -I/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/include -I/usr/include/eigen3 -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc experimentalData.h ErrorPDFs.cc -L/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF -L/home/alexander/Uni/minuit-cpp/build/lib/ -lminuit-cpp && ./minimizationStructureFunctions```

### compile:
```g++ -std=c++17 -I/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/include -I/usr/include/eigen3 -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc experimentalData.h ErrorPDFs.cc -L/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF -L/home/alexander/Uni/minuit-cpp/build/lib/ -lminuit-cpp```

### execute:
```./minimizationStructureFunctions```

### delete executable file:
```rm minimizationStructureFunctions```

### execute with just comments as output:
```./minimizationStructureFunctions | grep §```
