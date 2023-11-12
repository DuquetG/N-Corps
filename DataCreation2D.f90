module Constant
    implicit none
    real(8):: G=1 !G=6.67430*1e-11
end module Constant

!System's parameters
module System
    implicit none
    integer, parameter:: N=2 !number of body in the system
    integer, parameter:: dim=2 !number of dimension/!!! Currently, the code only supports dim=2
    integer, parameter:: d=2*dim
    real(8), dimension(N), parameter:: M=[100000.0, 1.0] !mass of bodies
end module System

!simulation's parameters
module Simulation
    implicit none
    integer, parameter:: Nstep=100000 !number of steps considered for the simulation
    real(8), parameter:: dt=0.01 !time step 
end module Simulation


program DataCreation2D
    use Constant
    use System
    use Simulation 
    implicit none
    real(8), dimension(N,d):: X
    real(8):: t
    integer:: i,j,k
    external:: deriv

    X(1,1)=0; X(1,2)=0; X(1,3)=0; X(1,4)=0
    X(2,1)=0; X(2,2)=100; X(2,3)=100; X(2,4)=0

    open(1, file='bodies_movement.dat')

    do i=1, Nstep
        t=i*dt

        call rk4(t,X,dt,N,d,deriv)

        write(1,*) X(2,1), X(2,2)

    enddo 

    close(1)

end program DataCreation2D

subroutine deriv(t,X,dx)
    use Constant
    use System
    use Simulation
    implicit none
    real(8), intent(in):: t
    real(8), dimension(N,d), intent(in):: X
    real(8), dimension(N,d), intent(out)::dx
    real(8), dimension(N,N):: xforce, yforce, Xdis
    integer:: i,j
    
    call distance(X,Xdis)
    call force(X,Xdis,xforce,yforce)

    do i=1, N
        dx(i,1)=X(i,3)
        dx(i,2)=X(i,4)
        dx(i,3)=0
        dx(i,4)=0

        do j=1, N
            dx(i,3)=dx(i,3)+1/M(i)*xforce(i,j)
            dx(i,4)=dx(i,4)+1/M(i)*yforce(i,j)
        enddo
    enddo 

end subroutine deriv

subroutine distance(X,Xdis)
    use System
    implicit none
    real(8), dimension(N,4), intent(in):: X
    real(8), dimension(N,N), intent(out):: Xdis
    integer:: i,j

    do i=1, N 
        do j=i, N
            Xdis(i,j)=sqrt((X(i,1)-X(j,1))**2+(X(i,2)-X(j,2))**2) 
            Xdis(j,i)=Xdis(i,j)
        enddo
    enddo

end subroutine distance

subroutine force(X,Xdis,xforce,yforce)
    use Constant
    use System
    implicit none
    real(8), dimension(N,d), intent(in):: X
    real(8), dimension(N,N), intent(in):: Xdis
    real(8), dimension(N,N), intent(out):: xforce
    real(8), dimension(N,N), intent(out):: yforce
    integer:: i,j
    do i=1, N
        do j=i, N
            
            if (i==j) then
                xforce(i,j)=0
                yforce(i,j)=0
            else if (i/=j) then

            xforce(i,j)=G*M(i)*M(j)*(X(j,1)-X(i,1))/Xdis(j,i)**3
            yforce(i,j)=G*M(i)*M(j)*(X(j,2)-X(i,2))/Xdis(j,i)**3

            xforce(j,i)=-xforce(i,j)
            yforce(j,i)=-yforce(i,j)
            end if
        enddo
    enddo

end subroutine

subroutine rk4(t,X,dt,N,d,deriv)
    implicit none
    integer , intent (in) :: d,N
    real (8) , intent (in) :: t, dt
    real (8) , dimension (N,d) , intent ( inout ) :: X
    real (8) :: ddt
    real (8) , dimension (N,d) :: Xp , k1 , k2 , k3 , k4
    ddt = 0.5* dt
    call deriv (t,X,k1); Xp = X + ddt *k1
    call deriv (t+ddt ,Xp ,k2); Xp = X + ddt*k2
    call deriv (t+ddt ,Xp ,k3); Xp = X + dt*k3
    call deriv (t+dt ,Xp ,k4); X = X + dt *( k1 + 2.0* k2 + 2.0* k3 + k4 )/6.0

end subroutine rk4