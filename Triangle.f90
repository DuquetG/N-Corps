program Triangle
    integer, parameter :: nbCorps = 2, Nstep = 100000, nbSimulations = 100
    integer :: sim, line
    real(8) :: dt = 1000.
    real(8), dimension(nbCorps) :: M = [1.989e30, 0.33011e24] !, 4.8675e24, 5.9724e24, 0.64171e24, 1898.19e24, 568.34e24, 86.813e24, 102.413e24]
    real(8), dimension(nbCorps, 4) :: X
    character(100) :: current

    open(1, file="Triangle.csv", action='readwrite', status='replace')
    close(1)
    
    do sim = 1, nbSimulations
        X(1,1) = 0; X(1,2) = 0; X(1,3) = 0; X(1,4) = 0
        X(2,1) = 57.909e9; X(2,2) = (sim-1)*1e8; X(2,3) = 0; X(2,4) = 47.36e3
        ! X(3,1) = 108.209e9; X(3,2) = 0; X(3,3) = 0; X(3,4) = 35.02e3
        ! X(4,1) = 149.596e9; X(4,2) = 0; X(4,3) = 0; X(4,4) = 29.78e3
        ! X(5,1) = 227.923e9; X(5,2) = 0; X(5,3) = 0; X(5,4) = 24.07e3
        ! X(6,1) = 778.570e9; X(6,2) = 0; X(6,3) = 0; X(6,4) = 13e3
        ! X(7,1) = 1433.529e9; X(7,2) = 0; X(7,3) = 0; X(7,4) = 9.68e3
        ! X(8,1) = 2872.463e9; X(8,2) = 0; X(8,3) = 0; X(8,4) = 6.80e3
        ! X(9,1) = 4495.060e9; X(9,2) = 0; X(9,3) = 0; X(9,4) = 5.43e3

        call Simulation2D(X, M, nbCorps, Nstep, dt, .true., 'csv', .false., .false.)

        open(2, file="Triangle.csv", action='readwrite', position='append')
        open(3, file="bodies_movement2D.csv", action='read')
        do line = 1, Nstep
            read(3, '(A)') current
            write(2, '(A)') current
        end do
        close(2)
        close(3)
    end do

end program Triangle


include 'Simulation2D.f90'