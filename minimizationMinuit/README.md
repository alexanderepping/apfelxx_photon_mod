# minimizationMinuit
Directory contains files related to the fitting of the Initial PDFs to the Structure Function data from several experiments (of course after evolving the Initial PDFs and calculating the Structure Functions from them) by minimization of the chi-Square using the MinuitCpp library. 

## files
### minimizationMinuit.h
File that `#include`s all the used header files.

### minimizationStructureFunctions.cc 
Main program to run for the minimization/fitting.

### configMinuit.h
File that contains configurations for all the other files like which Initial PDFs are used, some parameters etc. Can be changed by user.

### experimentalData.h
File that contains the experimental data. Experimental data can be added by user. To help that, see plottingPython/calculateExperimentalData.py

### StructureFunctionFcn.cc/.h
File that contains the class that calculates the StructureFunctions, takes the experimental data and returns the chi-square. It also contains the InitialPDFs.

### Makefile
Not ready to use yet. See "running the program" below.

## running the program
### compile:
```g++ -I/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/include -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc StructureFunctionsFcn.h minimizationMinuit.h configMinuit.h experimentalData.h -L/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF -L/home/alexander/Documents/libraries_cpp/minuit-cpp/build/lib/ -lminuit-cpp```

### execute:
```./minimizationStructureFunctions```

### delete executable file:
```rm minimizationStructureFunctions```

### execute with just comments as output:
```./minimizationStructureFunctions | grep §```

### compile and execute:
```g++ -I/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/include -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc StructureFunctionsFcn.h minimizationMinuit.h configMinuit.h experimentalData.h -L/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF && ./minimizationStructureFunctions -L/home/alexander/Documents/libraries_cpp/minuit-cpp/build/lib/ -lminuit-cpp```

### compile and execute with just comments as output:
```g++ -I/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/include -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc StructureFunctionsFcn.h minimizationMinuit.h configMinuit.h experimentalData.h -L/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF -L/home/alexander/Documents/libraries_cpp/minuit-cpp/build/lib/ -lminuit-cpp && ./minimizationStructureFunctions | grep §```