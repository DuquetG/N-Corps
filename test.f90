program test
    integer, parameter:: nbCorps=3, Nstep=100000
    real(8):: dt=100000
    real(8), dimension(nbCorps), parameter:: M=[1e24, 1e24, 1e24]
    real(8), dimension(nbCorps, 4):: X

    X(1,1)=50e4; X(1,2)=0; X(1,3)=0; X(1,4)=0
    X(2,1)=-50e4; X(2,2)=0; X(2,3)=0; X(2,4)=0
    X(3,1)=20e4; X(3,2)=25e4; X(3,3)=0; X(3,4)=0

    call Simulation2D(X, M, nbCorps, Nstep, dt, .true., 'dat', .false.)
end program test

include 'Simulation2D.f90'