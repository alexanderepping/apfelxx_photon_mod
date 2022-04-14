      program momentumSumRule 
        implicit none


        integer :: numberSteps = 9999999
        real :: mu = 10.0
        real :: momentumSum = 0.0

        integer :: i

        real :: X, Q2, UL, DL, SL, CL, BL, GL    

        Q2 = mu*mu

        do i = 1, numberSteps

          X = (i-0.5)/numberSteps

          call GRVGLO (X, Q2, UL, DL, SL, CL, BL, GL)

          momentumSum = momentumSum + 1.0/numberSteps
     >      *(GL+2*(UL+DL+SL+CL+BL))

        end do 

        write(*,*) "The Momentum Sum for
        write(*,*) momentumSum

      end program momentumSumRule
