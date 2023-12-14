program main2D
    integer, parameter :: nbCorps = 2, Nstep = 100000
    integer :: sim, line, line2, nbSimulations = 5, pasCalcul=10 ! Nombre de pas sautés à chaque calcul de Lyapunov
    integer:: length = 0
    real(8) :: dt = 500.
    real(8), dimension(nbCorps) :: M = [1.989e30, 0.33011e24]!, 4.8675e24, 5.9724e24, 0.64171e24, 1898.19e24]!, 568.34e24, &
    !86.813e24, 102.413e24] !
    real(8), dimension(nbCorps, 4) :: X
    logical, parameter:: wtraj=.true., wenergy=.true., wviriel=.true., wvelocity=.true.
    character(1000) :: current, current_v, command
    character(len=*), parameter::  format='csv'
    integer :: status


    open(1, file="CSVs/all_positions_2D.csv", action='readwrite', status='replace')
    close(1)
    open(6, file="CSVs/all_velocities_2D.csv", action='readwrite', status='replace')
    close(6)

    
    do sim = 1, nbSimulations

        ! Coordonnées initiales du système solaire: premier indice: numéro du corps. Second indice: (1,2,3,4)==(x,y,vx,vy)

        X(1,1) = 0; X(1,2) = 0; X(1,3) = 0; X(1,4) = 0
        X(2,1) = 57.909e9+(sim-1)*5e9; X(2,2) = 0; X(2,3) = 0; X(2,4) = 47.36e3
        ! X(3,1) = 108.209e9; X(3,2) = 0; X(3,3) = 0; X(3,4) = 35.02e3
        ! X(4,1) = 149.596e9; X(4,2) = 0; X(4,3) = 0; X(4,4) = 29.78e3
        ! X(5,1) = 227.923e9; X(5,2) = 0; X(5,3) = 0; X(5,4) = 24.07e3
        ! X(6,1) = 778.570e9; X(6,2) = 0; X(6,3) = 0; X(6,4) = 13e3
        ! X(7,1) = 1433.529e9; X(7,2) = 0; X(7,3) = 0; X(7,4) = 9.68e3
        ! X(8,1) = 2872.463e9; X(8,2) = 0; X(8,3) = 0; X(8,4) = 6.80e3
        ! X(9,1) = 4495.060e9; X(9,2) = 0; X(9,3) = 0; X(9,4) = 5.43e3

        call Simulation2D(X, M, nbCorps, Nstep, dt, wtraj, format, wenergy, wviriel, wvelocity)
        
        open(2, file="CSVs/all_positions_2D.csv", action='readwrite', position='append')
        open(3, file="CSVs/positions_2D.csv", action='read')
        line=0
       
        do line = 0, Nstep-1
            read(3, '(A)') current
            if (mod(line,pasCalcul)==0) then
                length=length+1
                write(2, '(A)') current
            endif
        end do
        close(3)
        close(2)


        open(4, file="CSVs/all_velocities_2D.csv", action='readwrite', position='append')
        open(5, file="CSVs/velocities_2D.csv", action='read')
       
        do line2 = 0, Nstep-1
            read(5, '(A)') current_v
            if (mod(line2,pasCalcul)==0) then
                write(4, '(A)') current_v
            endif
        end do
        close(5)
        close(4)
    enddo


    write(command, '(A,I0)') "python Outils/csv_healer.py ", nbSimulations
    call execute_command_line(command, wait=.true., exitstat=status)
    
    call system("python Outils/Animation2D.py")

    call Lyapunov(nbSimulations, Nstep, nbCorps, M, dt, length, pasCalcul)
    
    call system("python Outils/graphiques.py")

end program main2D


include 'Outils/Simulation2D.f90'
include 'Outils/Lyapunov.f90'
include 'Outils/Integrator.f90'