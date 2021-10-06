# bashFiles
Some useful bash files I use to run the programs and/or change between different installations of Apfel++ and the LHAPDF library.

## run_Evolution.sh
Program to make and run the Evolution program and possibly reinstall the LHAPDF and/or Apfel++ library.
- The program also outputs the path of the currently used Apfel++ and LHAPDF Installation. This output is aligned with he information output of the Evolution.cc file.
- The variable $CURRENT_APFEL should point to the folder with the apfelxx folder inside.
- The variable $CURRENT_LHAPDF should point to the folder with the LHAPDF-6.4.0 folder inside.
- The variable $$CURRENT_APFEL_TEST should point to the folder with the Evolution.cc and its respectiveMakefile inside.                               

## change_Apfel_Installation.sh
Convenient program to change the currently used installation of Apfel++.
- Furthermore it will save the current state of the .bash_aliases file in a backup file
- The .bash_aliases file should be constucted in a similar way to the one below to work properly. My current .bash_aliases file can be also found in the Github.
- The paths obiously should be changed accordingly.

## change_LHAPDF_Installation.sh
Convenient program to change the currently used installation of LHAPDF.
- Furthermore it will save the current state of the .bash_aliases file in a backup file
- The .bash_aliases file should be constucted in a similar way to the one below to work properly. My current .bash_aliases file can be also found in the Github.
- The paths obiously should be changed accordingly.

## .bash_aliases
My currently used .bash_aliases file which works with the bash files above
