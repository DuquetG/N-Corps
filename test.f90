program test
    integer, parameter:: nbCorps=4, Nstep=500000
    real(8):: dt=0.01
    real(8), dimension(nbCorps), parameter:: M=[1.989e24, 1.55e24, 2.1e24, 4.7e30]
    real(8), dimension(nbCorps, 4):: X

    X(1,1)=-50e9; X(1,2)=-50e9; X(1,3)=0; X(1,4)=-50e3
    X(2,1)=50e9; X(2,2)=-50e9; X(2,3)=0; X(2,4)=46e3
    X(3,1)=-50e9; X(3,2)=50e9; X(3,3)=0; X(3,4)=50e3
    X(4,1)=0; X(4,2)=0; X(4,3)=0; X(4,4)=0



    call Simulation2D(X, M, nbCorps, Nstep, dt, .true., 'dat', .true.,.false.)
end program test

include 'Simulation2D.f90'