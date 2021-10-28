# bashFiles
Some useful bash files I use to run the programs and/or change between different installations of Apfel++ and the LHAPDF library.

## run_Evolution.sh
Program to make and run the Evolution program and possibly reinstall the LHAPDF and/or Apfel++ library.
- The program also outputs the path of the currently used Apfel++ and LHAPDF Installation. This output is aligned with he information output of the Evolution.cc file.
- The variable $CURRENT_APFEL should point to the folder with the apfelxx folder inside.
- The variable $CURRENT_LHAPDF should point to the folder with the LHAPDF-6.4.0 folder inside.
- The variable $$CURRENT_APFEL_TEST should point to the folder with the Evolution.cc and its respective Makefile inside.     

## run_EvolutionFlavors.sh
Program to make and run the EvolutionFlavors program and possibly reinstall the LHAPDF and/or Apfel++ library.
- Other than the Evolution program the EvolutionFlavors program outputs the results for multiple particles and final energies.
- The program also outputs the path of the currently used Apfel++ and LHAPDF Installation. This output is aligned with he information output of the EvolutionFlavors.cc file.
- The variable $CURRENT_APFEL should point to the folder with the apfelxx folder inside.
- The variable $CURRENT_LHAPDF should point to the folder with the LHAPDF-6.4.0 folder inside.
- The variable $$CURRENT_APFEL_TEST should point to the folder with the EvolutionFlavors.cc and its respective Makefile inside.   

## run_StructureFunctions.sh
Program to make and run the EvolutionStructureFunctions program to evolve the PDFs and then make, then run the StructureFunctions program and possibly reinstall the LHAPDF and/or Apfel++ library.
- The EvolutionStructureFunctions.cc program makes an LHAPDF grid out of the evolved PDFs it calculates. This grid is then used to calculate the Structure Functions using the StructureFunctions.cc program. The output is then saved to a file and read and plotted by plottingStructureFunctions.py.
- The program also outputs the path of the currently used Apfel++ and LHAPDF Installation. This output is aligned with he information output of the StructureFunctions.cc file.
- The variable $CURRENT_APFEL should point to the folder with the apfelxx folder inside.
- The variable $CURRENT_LHAPDF should point to the folder with the LHAPDF-6.4.0 folder inside.
- The variable $$CURRENT_APFEL_TEST should point to the folder with the EvolutionStructureFunctions.cc file, the StructureFunctions.cc file and its respective Makefile inside.    

## makeLHAPDFLO.sh
Quick program to run the makeLHAPDFLO.f program and copy the file into the LHAPDF share directory.  

## makeLHAPDFHO.sh
Quick program to run the makeLHAPDFHO.f program and copy the file into the LHAPDF share directory.

## change_Apfel_Installation.sh
Convenient program to change the currently used installation of Apfel++.
- Furthermore it will save the current state of the .bash_aliases_apfel file in a backup file
- The .bash_aliases_apfel file should be constucted in a similar way to the one below to work properly. My current .bash_aliases_apfel file can be also found in the Github.
- The paths obiously should be changed accordingly.

## change_LHAPDF_Installation.sh
Convenient program to change the currently used installation of LHAPDF.
- Furthermore it will save the current state of the .bash_aliases file in a backup file
- The .bash_aliases file should be constucted in a similar way to the one below to work properly. My current .bash_aliases file can be also found in the Github.
- The paths obiously should be changed accordingly.

## .bash_aliases_apfel
My currently used .bash_aliases_apfel file which works with the bash files above. For this file to have any effect it should be included into the .bash_aliases or .bashrc file by writing:
```
if [ -f ~/.bash_aliases_apfel ]; then
    . ~/.bash_aliases_apfel
fi
```
