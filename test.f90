program test
    integer, parameter:: nbCorps=2, Nstep=100000
    real(8):: dt=10000
    real(8), dimension(nbCorps), parameter:: M=[1.989e30, 0.33011e24]
    real(8), dimension(nbCorps, 4):: X

    X(1,1)=0; X(1,2)=0; X(1,3)=0; X(1,4)=0
    X(2,1)=57.909e9; X(2,2)=0; X(2,3)=0; X(2,4)=47.36e3


    call Simulation2D(X, M, nbCorps, Nstep, dt, .true., 'dat', .false.)
end program test

include 'Simulation2D.f90'