# makeLHAPDF
Fortran program to write the GRVCustomSet_0000.dat file in the LHAPDF fromat and the other used files.

## makeLHAPDF.f
Program to write a GRVCustomSet_0000.dat file in the LHAPDF fromat.

The program takes a template file and copies the first three lines. 
Then it reads the x- (4th line) and Q-values (5th line) and saves them. 
They are then used with the GRVGLO subroutine to calculate the Photon PDF values at different x- and Q²-values. 
These values are then written to the output file.
