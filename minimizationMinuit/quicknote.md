
# To-Do:
[ ] write LaTeX file further
[ ] write README.md file for this folder

[ ] mainfile execution
[ ] `#pragma once` errors?


[ ] config file: add parameters like (Output? y/n?)

[x] is it okay to leave the cutoff out? - Yes. The cutoff was only used because the output to the LHAPDF file couldn't be zero. If it doesn't produce any errors, we can just leave the negative numbers
[x] initial values? - Don't really matter
[x] add initial values in the config file


# Things to check, if there is some error:
- maybe introduce cutoff in StructureFunctionsFcn.cc in PDFsEvolved if there are some errors bc of negative numbers
	

# commands to compile etc:
## compile:
g++ -I/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/include -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc StructureFunctionsFcn.h minimizationMinuit.h configMinuit.h experimentalData.h InitialPDFs.h -L/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF

## execute:
./minimizationStructureFunctions

## execute with just comments as output:
./minimizationStructureFunctions | grep §

## compile and execute:
g++ -I/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/include -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc StructureFunctionsFcn.h minimizationMinuit.h configMinuit.h experimentalData.h InitialPDFs.h -L/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF && ./minimizationStructureFunctions

## compile and execute with just comments as output:
g++ -I/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/include -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc StructureFunctionsFcn.h minimizationMinuit.h configMinuit.h experimentalData.h InitialPDFs.h -L/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF && ./minimizationStructureFunctions | grep §


# old commands:
## command to run the AIOStructureFunctions (not up to data)
g++ -I/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/include -o AIOStructureFunctions_ AIOStructureFunctions_.cc AIOStructureFunctions.cc AIOStructureFunctions.h -L/home/alexander/Documents/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF
