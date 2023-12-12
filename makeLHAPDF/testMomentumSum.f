      program momentumSumRule 
        implicit none


        real :: pi = 3.14159265

        integer :: numberSteps = 9999999
        real :: Q2 = 2.0
        real :: momentumSum = 0.0
        real :: momentumSumCorrect = 0.0

        integer :: i

        real :: X, UL, DL, SL, CL, BL, GL    

        do i = 1, numberSteps

          X = (i-0.5)/numberSteps

          call GRVGHO (X, Q2, UL, DL, SL, CL, BL, GL)

          momentumSum = momentumSum + 1.0/numberSteps
     >      *(GL+2*(UL+DL+SL+CL+BL))

        end do 

        momentumSumCorrect = 1.0 + 2.0/(3.0*pi)*log(Q2/4)

        write(*,*) "The Momentum Sum for Q^2=", Q2, ":"
        write(*,*) momentumSum

        write(*,*) "It should be:"
        write(*,*) momentumSumCorrect

        write(*,*) "Therefore the difference in percent is:"
        write(*,*) (momentumSum/momentumSumCorrect -1.0) * 100

      end program momentumSumRule
