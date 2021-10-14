# apfelxx_photon_mod
Modification of the [vbertone/apfelxx](https://github.com/vbertone/apfelxx) code to allow the calculation of photon PDFs.

Done in the course of writing my master thesis.

## compiling
in the main folder use
```
cmake .
make
```
- fresh installation: remove `CMakeCache.txt` and the folder `CMakeFiles`
- executables can be added in the `CMakeLists.txt` 
- to switch the `LHAPDF`-installation just change your `PATH` variable to point to the correct `lhapdf-config`

## folders:

### bashFiles
Some useful bash files I use to run the programs and/or change between different installations of Apfel++ and the LHAPDF library.

### makeLHAPDF
Fortran program to write the GRVCustomSet_0000.dat file in the LHAPDF fromat and related files.

### myApfelModified/apfelxx
Folder containing the changed source code and headers of the apfelxx program.

### myLHAPDFModified/share/LHAPDF/
Folder containing folders with the used LHAPDF data sets.

### plottingPython
Folder containing various python files to plot the data collected from the Evolution.cc program. It also contains a program to make an array of logarithmic values.

### testModifiedApfel
Modified versions of the examples found in [vbertone/APFEL_EXAMPLES](https://github.com/vbertone/APFEL_Examples).
