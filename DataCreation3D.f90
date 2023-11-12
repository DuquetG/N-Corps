module Constant
    implicit none
    real(8):: G=0.00001 !G=6.67430*1e-11
end module Constant

!System's parameters
module System
    implicit none
    integer, parameter:: N=4 !number of body in the system
    integer, parameter:: dim=3 !number of dimension/!!! Currently, the code only supports dim=3
    integer, parameter:: d=2*dim
    real(8), dimension(N), parameter:: M=[5000, 5000, 5000, 5000] !bodies's masses 
end module System

!simulation's parameters
module Simulation
    implicit none
    integer, parameter:: Nstep=1000000 !number of steps considered for the simulation
    real(8), parameter:: dt=0.01 !time step 
end module Simulation

program N_Body
    use Constant
    use System
    use Simulation 
    implicit none
    real(8), dimension(N,d):: X
    real(8):: t
    integer:: i,j,k
    external:: deriv

    X(1,1)=25; X(1,2)=32; X(1,3)=31; X(1,4)=0; X(1,5)=0; X(1,6)=0
    X(2,1)=19.5; X(2,2)=45; X(2,3)=-27; X(2,4)=0; X(2,5)=0; X(2,6)=0
    X(3,1)=17.6; X(3,2)=-30; X(3,3)=-28; X(3,4)=0; X(3,5)=0; X(3,6)=0
    X(4,1)=-10; X(4,2)=25; X(4,3)=-12; X(4,4)=0; X(4,5)=0; X(4,6)=0

    open(1, file='bodies_movement3D.dat')

    do i=1, Nstep
        t=i*dt

        call rk4(t,X,dt,N,d,deriv)

        write(1,*) X(2,1), X(2,2), X(1,1), X(1,2), X(3,1), X(3,2), X(3,1), X(3,2)

    enddo 

    close(1)

end program N_Body

subroutine deriv(t,X,dx)
    use Constant
    use System
    use Simulation
    implicit none
    real(8), intent(in):: t
    real(8), dimension(N,d), intent(in):: X
    real(8), dimension(N,d), intent(out)::dx
    real(8), dimension(N,N):: xforce, yforce, zforce, Xdis
    integer:: i,j
    
    call distance(X,Xdis)
    call force(X,Xdis,xforce,yforce,zforce)

    do i=1, N
        dx(i,1)=X(i,4)
        dx(i,2)=X(i,5)
        dx(i,3)=X(i,6)
        dx(i,4)=0
        dx(i,5)=0
        dx(i,6)=0


        do j=1, N
            dx(i,4)=dx(i,4)+1/M(i)*xforce(i,j)
            dx(i,5)=dx(i,5)+1/M(i)*yforce(i,j)
            dx(i,6)=dx(i,6)+1/M(i)*zforce(i,j)
        enddo
    enddo 

end subroutine deriv

subroutine distance(X,Xdis)
    use System
    implicit none
    real(8), dimension(N,d), intent(in):: X
    real(8), dimension(N,N), intent(out):: Xdis
    integer:: i,j

    do i=1, N 
        do j=i, N
            Xdis(i,j)=sqrt((X(i,1)-X(j,1))**2+(X(i,2)-X(j,2))**2+(X(i,3)-X(j,3))**2) 
            Xdis(j,i)=Xdis(i,j)
        enddo
    enddo

end subroutine distance

subroutine force(X,Xdis,xforce,yforce, zforce)
    use Constant
    use System
    implicit none
    real(8), dimension(N,d), intent(in):: X
    real(8), dimension(N,N), intent(in):: Xdis
    real(8), dimension(N,N), intent(out):: xforce, yforce, zforce
    integer:: i,j
    do i=1, N
        do j=i, N
            
            if (i==j) then
                xforce(i,j)=0
                yforce(i,j)=0
            else if (i/=j) then

            xforce(i,j)=G*M(i)*M(j)*(X(j,1)-X(i,1))/Xdis(j,i)**3
            yforce(i,j)=G*M(i)*M(j)*(X(j,2)-X(i,2))/Xdis(j,i)**3
            zforce(i,j)=G*M(i)*M(j)*(X(j,3)-X(i,3))/Xdis(j,i)**3

            xforce(j,i)=-xforce(i,j)
            yforce(j,i)=-yforce(i,j)
            zforce(j,i)=-zforce(i,j)
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