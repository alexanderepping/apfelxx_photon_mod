* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  makeLHAPDF.f                                                                                                     *
*  Author: Alexander Epping: a_eppi01@uni-muenster.de                                                               *
*  6 Oct 2021                                                                                                       *         
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  Program to write a GRVCustomSet_0000.dat file in the LHAPDF fromat.                                              *            
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*  The program takes a template file and copies the first three lines.                                              *
*  Then it reads the x- (4th line) and Q-values (5th line) and saves them.                                          *
*  They are then used with the GRVGLO subroutine to calculate the Photon PDF values at different x- and Q²-values.  *
*  These values are then written to the output file.                                                                *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      program makeLHAPDF 
        implicit none

c       file, which has the same X- and Q-values as the one we want to create    
c       the file should be in a LHAPDF format
        character(len=19) :: template_file = 'LHAPDF_template.dat'

c       file, which we want to create in the LHAPDF format
        character(len=21) :: output_file = 'GRVCustomSet_0000.dat'

        integer :: lenXArr, lenQArr
        character(len=99999) :: lineOne, lineTwo, lineThree
        character(len=99999) :: lineX, lineQ, lineF

        real, allocatable,dimension(:) :: XArr, QArr

c       variables that are allocated inside the do-loops/by the GRVGLO-subroutine
        real :: X, Q2, UL, DL, SL, CL, BL, GL    

c       running index variables for do-loops
        integer :: i, j

c       Structure constant
        real :: alphaQED = 1./137.




c - get the Q- and X-values        
        open(1, file = template_file, status = 'old')
              
        read (1,'(A)') lineOne  
        read (1,'(A)') lineTwo
        read (1,'(A)') lineThree

c       read the next three lines as strings
        read (1,'(A)') lineX
        read (1,'(A)') lineQ
        read (1,'(A)') lineF
        
        close(1)

c       calculate the number of numbers in the respective lines
c       13 is the number of characters per data-point
        lenXArr = ( len(trim(lineX)) + 1 ) / 13 
        lenQArr = ( len(trim(lineQ)) + 1 ) / 13

c       tell the arrays how big they should be  
        allocate(XArr(lenXArr), QArr(lenQArr))


        open(1, file = template_file, status = 'old')
        
c       skip the first three lines        
        read (1,*)
        read (1,*)
        read (1,*) 

c       read the next two lines into the array
        read (1,*) XArr(:) 
        read (1,*) QArr(:)
        
        close(1)

        
c - open the output file and write the first lines

        open(2, file = output_file, status="unknown", action="write")

        write(2,'(A)') trim(lineOne)
        write(2,'(A)') trim(lineTwo)
        write(2,'(A)') trim(lineThree)
        write(2,'(A)') trim(lineX)
        write(2,'(A)') trim(lineQ)
        write(2,'(A)') trim(lineF)

c - run the loop, write all the data and close the file
        do i = 1, lenXArr
          do j = 1, lenQArr

            X = XArr(i)
            Q2 = QArr(j)*QArr(j)

c           call the subroutine from grvphoton.f
            call GRVGLO (X, Q2, UL, DL, SL, CL, BL, GL)

c           grvglo returns 1 / alphaQED * X * PDF ->
            GL = GL * alphaQED / X
            BL = BL * alphaQED / X
            CL = CL * alphaQED / X
            SL = SL * alphaQED / X
            DL = DL * alphaQED / X
            UL = UL * alphaQED / X
c           grvglo returns 1 / alphaQED * X * PDF <-

            write(2,"(10(es12.6, 1x), es12.6)") BL, CL, SL, DL, UL, UL, 
     >       DL, SL, CL, BL, GL
          end do
        end do 

c       write the last line which is the same as the third
        write(2,'(A)') trim(lineThree)

        close(2)


c - deallocate the memory of the arrays
        deallocate(QArr, XArr)

      end program makeLHAPDF
