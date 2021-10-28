# apfelxx_photon_mod
Modification of the [vbertone/apfelxx](https://github.com/vbertone/apfelxx) code to allow the calculation of photon PDFs. 

The modification is done by adding the pointlike contribution directly in the EvolveObject function in [matchedevolution.cc](https://github.com/alexanderepping/apfelxx_photon_mod/blob/main/myApfelModified/apfelxx/src/kernel/matchedevolution.cc).
The pointlike contributions can be seein in [pointlikecontributions.cc](https://github.com/alexanderepping/apfelxx_photon_mod/blob/main/myApfelModified/apfelxx/src/kernel/pointlikecontributions.cc).

Done in the course of writing my master thesis.

## folders:

### bashFiles
Some useful bash files I use to run the programs and/or change between different installations of Apfel++ and the LHAPDF library.

### makeLHAPDF
Fortran program to write GRVCustomSet_0000.dat files in the LHAPDF fromat and related files.

### myApfelModified/apfelxx
Folder containing the changed source code and headers of the apfelxx program.

Changed files: 
- matchedevolution.cc/.h: 
    - Removed Set<Distribution> as template class
    - Manually added MatchedEvolution<Set<Distribution>>
    - Added Evaluate and EvolveObject functions with additional Alphas input parameter. The EvolveObject function ( for <Set<Distribution>) with Alphas input includes the pointlike contribution into the evolution.
- tabulateobjects.cc/.h:
    - added functions that can pass the Alphas through to the functions of matchedevolution

added files: 
- pointlikecontributions.cc/.h: 
    - functions etc. to include the pointlike contribution of the photon into the evolution 


### myLHAPDFModified/share/LHAPDF/
Folder containing folders with the used LHAPDF data sets.

### plottingPython
Folder containing various python files to plot the data collected from the Evolution.cc program. It also contains a program to make an array of logarithmic values.

### testModifiedApfel
Modified versions of the examples found in [vbertone/APFEL_EXAMPLES](https://github.com/vbertone/APFEL_Examples).


## compiling using cmake
in the main folder use
```
cmake .
make
```
- fresh installation: remove `CMakeCache.txt` and the folder `CMakeFiles`
- executables can be added in the `CMakeLists.txt` 
- to switch the `LHAPDF`-installation just change your `PATH` variable to point to the correct `lhapdf-config`
