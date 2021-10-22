# makeLHAPDF
Fortran program to write the GRVCustomSet_0000.dat file in the LHAPDF fromat and related files.

The program can be run by running the bashFiles/makeLHAPDF.sh program, which will also copy the GRVCustomSet_0000.dat directly into the LHAPDF share directory.

## makeLHAPDF.f
Program to write a GRVCustomSet_0000.dat file in the LHAPDF fromat.

The program takes a template file and copies the first three lines. 
Then it reads the x- (4th line) and Q-values (5th line) and saves them. 
They are then used with the GRVGLO subroutine to calculate the Photon PDF values at different x- and Q²-values. 
These values are then written to the output file.

## grvphoton.f
Program by GRV containing the GRVGLO subroutine which calculates the PDFs at different x- and Q²-values.

## LHAPDF_template.dat
Template for the LHAPDF format.

## GRVCustomSet_0000.dat
GRV Photon PDF at different x- and Q²-values in the LHAPDF format.
It is output by the makeLHAPDF.f file and is readable by the LHAPDF library.
