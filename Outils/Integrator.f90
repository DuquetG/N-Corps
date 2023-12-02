subroutine rk4(t,X,dt,N,M,deriv)
    implicit none
    integer , intent (in) :: N
    real (8) , intent (inout) :: t, dt
    real(8), dimension(N), intent(in) :: M
    real (8) , dimension (N,4) , intent ( inout ) :: X
    real (8) :: ddt
    real (8) , dimension (N,4) :: Xp , k1 , k2 , k3 , k4
    ddt = 0.5* dt
    call deriv (t,N,M,X,k1); Xp = X + ddt*k1
    call deriv (t+ddt,N,M,Xp,k2); Xp = X + ddt*k2
    call deriv (t+ddt,N,M,Xp ,k3); Xp = X + dt*k3
    call deriv (t+dt,N,M,Xp ,k4); X = X+dt*(k1+2.0*k2+2.0*k3+k4)/6.0
    t=t+dt
end subroutine rk4

subroutine euler(t,X,dt,N,M,deriv)
    implicit none
    integer , intent (in) :: N
    real (8) , intent (in) :: t, dt
    real(8), dimension(N), intent(in) :: M
    real (8) , dimension (N,4) , intent (inout) :: X
    real (8) , dimension (N,4) :: dX
    call deriv (t,N,M,X,dX)
    X= X+dt*dX
end subroutine euler

subroutine velocity_verlet (t,X,dt,N,M,deriv)
    implicit none
    integer , intent (in) :: n
    real (8) , intent (in) :: t, dt
    real(8), dimension(N), intent(in) :: M
    real (8) , dimension (N,4) , intent ( inout ) :: X
    real (8) :: ddt
    real (8) , dimension (N,4) :: dX
        if (mod (n ,2) .ne. 0) write (*,*) 'WARNING : N should be even for Verlet'
        call deriv (t,n,M,X,dX)
        X(1:N/2,:) = X(1:N/2,:) + dt* dX(1:N/2,:)
        call deriv (t,n,M,X,dX)
        X(N /2+1: N, :) = X(N /2+1: N, :) + 0.5* dt* dX(N/2+1: N, :)
        X(N /2+1: N, :) = X(N /2+1: N, :) + 0.5* dt* dX(N /2+1: N, :)
end subroutine velocity_verlet

subroutine adaptativerk4(t,X,dt,N,M,tolerance)
    implicit none
    integer , intent (in) :: N
    real (8) , intent (inout) :: t, dt
    real(8), dimension(N), intent(in) :: M
    real (8) , dimension (N,4) , intent ( inout ) :: X
    real (8), intent(in) :: tolerance
    real (8), dimension(N,4):: Xnew
    real (8)::  etot1, etot2, epot, ecin
    external:: deriv

    Xnew=X
    call energy(N,M,X,ecin,epot)
    etot1=epot+ecin
    Call rk4(t,Xnew,dt,N,M,deriv)
    call energy(N,M,Xnew,ecin,epot)
    etot2=epot+ecin

    if(abs(etot2-etot1)<=tolerance) then
        dt=dt*1.2
    else 
        dt=dt/1.2
    end if
    !write(*,*) dt
    X=Xnew
    !write(*,*) t
    
end subroutine adaptativerk4

