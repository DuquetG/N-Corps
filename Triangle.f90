program Triangle
    integer, parameter :: nbCorps = 3, Nstep = 100000, nbSimulations = 3
    integer :: sim, line
    real(8), parameter :: dt = 1000.
    real(8), dimension(nbCorps) :: M = [1.989e30, 0.33011e24, 4.8675e24]
    real(8), dimension(nbCorps, 4) :: X
    character(100) :: current

    open(1, file="Triangle.csv", action='readwrite', status='replace')
    close(1)
    
    do sim = 1, nbSimulations
        X(1,1) = 0; X(1,2) = 0; X(1,3) = 0; X(1,4) = 0
        X(2,1) = 57.909e9; X(2,2) = sim*17e9; X(2,3) = 0; X(2,4) = 47.36e3
        X(3,1) = 108.209e9; X(3,2) = 0; X(3,3) = 0; X(3,4) = 35.02e3

        call Simulation2D(X, M, nbCorps, Nstep, dt, .true., 'csv', .false.)

        open(2, file="Triangle.csv", action='readwrite', position='append')
        open(3, file="bodies_movement2D.csv", action='read')
        do line = 1, Nstep
            read(3, '(A)') current
            write(2, '(A)') current
        end do
        close(3)
        close(2)
    end do

end program Triangle


include 'Simulation2D.f90'