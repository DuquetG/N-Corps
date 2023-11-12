
!Constantes fondamentales
module Constant
    implicit none
    real(8):: G=6.67430*1e-11
end module Constant

!System's parameters
module System
    implicit none
    integer, parameter:: N=9 !Number of body in the system
    integer, parameter:: dim=2 !Number of dimension/!!! Currently, the code only supports dim=2
    integer, parameter:: d=2*dim
    real(8), dimension(N), parameter:: M=[1.989e30, 0.33011e24, 4.8675e24, 5.9724e24, &
                                          0.64171e24, 1898.19e24, 568.34e24, 86.813e24, 102.413e24] !Mass of bodies
end module System

!Simulation's parameters
module Simulation
    implicit none
    integer, parameter:: Nstep=1000000 !number of steps considered for the simulation
    real(8), parameter:: dt=10000 !time step 
end module Simulation


program DataCreation2D
    use Constant
    use System
    use Simulation 
    implicit none
    real(8), dimension(N,d):: X !position/velocity matrix// X(i,1)=x_i, X(i,2)=y_i, X(i,3)=v_x_i, X(i,4)=v_y_i, i: body's index
    real(8):: t               
    integer:: i,j,k
    external:: deriv  

    !Initialize position and velocity
    X(1,1)=0; X(1,2)=0; X(1,3)=0; X(1,4)=0
    X(2,1)=57.909e9; X(2,2)=0; X(2,3)=0; X(2,4)=47.36e3
    X(3,1)=108.209e9; X(3,2)=0; X(3,3)=0; X(3,4)=35.02e3
    X(4,1)=149.596e9; X(4,2)=0; X(4,3)=0; X(4,4)=29.78e3
    X(5,1)=227.923e9; X(5,2)=0; X(5,3)=0; X(5,4)=24.07e3
    X(6,1)=778.570e9; X(6,2)=0; X(6,3)=0; X(6,4)=13e3
    X(7,1)=1433.529e9; X(7,2)=0; X(7,3)=0; X(7,4)=9.68e3
    X(8,1)=2872.463e9; X(8,2)=0; X(8,3)=0; X(8,4)=6.80e3
    X(9,1)=4495.060e9; X(9,2)=0; X(9,3)=0; X(9,4)=5.43e3

    open(1, file='bodies_movement2D.csv')

    do i=1, Nstep
        t=i*dt

        !Update X for a time step dt
        call rk4(t,X,dt,N,d,deriv)

        write (1, '(*(G0.6,:,";"))') X(1,1), X(1,2), X(2,1), X(2,2), X(3,1), X(3,2), X(4,1), X(4,2), &
                                     X(5,1), X(5,2), X(6,1), X(6,2), X(7,1), X(7,2), X(8,1), X(8,2), X(9,1), X(9,2)

    enddo 

    close(1)

end program DataCreation2D

!Calculate the derivative dX of the X matrix
subroutine deriv(t,X,dX)
    use Constant
    use System
    use Simulation
    implicit none
    real(8), intent(in):: t
    real(8), dimension(N,d), intent(in):: X
    real(8), dimension(N,d), intent(out)::dX !X's derivative
    real(8), dimension(N,N):: xforce, yforce, Xdis
    integer:: i,j
    
    call distance(X,Xdis)
    call force(X,Xdis,xforce,yforce)

    do i=1, N
        dX(i,1)=X(i,3) !dx/dt=v_x
        dX(i,2)=X(i,4) !dy/dt=v_y
        dX(i,3)=0
        dX(i,4)=0
        
        !add all the forces applied for the N bodies
        do j=1, N
            dX(i,3)=dX(i,3)+1/M(i)*xforce(i,j) !dv_x/dt=mi*a_x=G*mi*m1/ri1+G*mi*m2/ri2+...
            dX(i,4)=dX(i,4)+1/M(i)*yforce(i,j) !dv_y/dt=ma_y
        enddo
    enddo 

end subroutine deriv

!Calculate de Xdis matrix where Xdis(i,j) stands for 
!the distance between the i and j bodies.
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

!Calculate the force between each body.
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
                xforce(i,j)=0  !the force of the body applied to itself is null
                yforce(i,j)=0  
            else if (i/=j) then

            xforce(i,j)=G*M(i)*M(j)*(X(j,1)-X(i,1))/Xdis(j,i)**3 !newton's inverse-square law
            yforce(i,j)=G*M(i)*M(j)*(X(j,2)-X(i,2))/Xdis(j,i)**3

            xforce(j,i)=-xforce(i,j) 
            yforce(j,i)=-yforce(i,j)
            end if
        enddo
    enddo

end subroutine

!Implements the fourth-order Runge-Kutta method to update the positions of the bodies over time.
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