# makeLHAPDF
Fortran program to write GRVCustomSet_0000.dat files in the LHAPDF fromat and related files.

## makeLHAPDFLO.f
Program to write a GRVCustomSetLO_0000.dat file in the LHAPDF fromat.

The program takes a template file and copies the first three lines. 
Then it reads the x- (4th line) and Q-values (5th line) and saves them. 
They are then used with the GRVGLO subroutine to calculate the Photon PDF values at different x- and Q²-values. 
These values are then written to the output file.

The program can be run by running the bashFiles/makeLHAPDFLO.sh program, which will also copy the GRVCustomSetLO_0000.dat directly into the LHAPDF share directory and delete the executable as well as the .dat file in here.

## makeLHAPDFHO.f
Program to write a GRVCustomSetHO_0000.dat file in the LHAPDF fromat.

The program takes a template file and copies the first three lines.
Then it reads the x- (4th line) and Q-values (5th line) and saves them.
They are then used with the subroutines from the grvphoton.f program to calculate the Photon PDF values at different x- and Q²-values.
These values are then written to the output file.
You can change if you want to add for example the LO terms by changing the integer "mode".

The program can be run by running the bashFiles/makeLHAPDFHO.sh program, which will also copy the GRVCustomSetHO_0000.dat directly into the LHAPDF share directory and delete the executable as well as the .dat file in here.

## grvphoton.f
Program by GRV containing the GRVGLO subroutine which calculates the PDFs at different x- and Q²-values.

## LHAPDF_template.dat
Template for the LHAPDF format.
