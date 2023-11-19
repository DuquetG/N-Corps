
!Constantes fondamentales
module Constant
    implicit none
    real(8):: G=6.67430*1e-11
end module Constant

!System's parameters
module System
    implicit none
    integer, parameter:: dim=2 !Number of dimension
    integer, parameter:: d=2*dim
    
    !integer, parameter:: N=0 !Number of body in the system
    
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
    !real(8), dimension(N,d):: X !position/velocity matrix// X(i,1)=x_i, X(i,2)=y_i, X(i,3)=v_x_i, X(i,4)=v_y_i, i: body's index
    !integer, parameter:: nbCorps=N
    real(8):: t               
    integer:: i,j,k,a,b
    external:: deriv  



    !Initialize position and velocity
    integer :: nbCorps, io_status
    character(len=100) :: line, preset_name, filename
    real(8), allocatable :: X(:,:)
    real(8), allocatable :: M(:) !Mass of bodies
    !integer :: q=0,p=0


    if (command_argument_count() == 1) then
        call get_command_argument(1, filename)
    else
        write(*,*) 'Entrez le nom du fichier de presets : '
        read(*,*) filename
    end if


    open(unit=20, file=filename, status='old', action='read',iostat=io_status)

    if (io_status /= 0) then
        write(*,*) 'Erreur lors de l''ouverture du fichier.'
        stop
    end if

    read(20,*) nbCorps
    allocate(X(nbCorps,d))
    allocate(M(nbCorps))
    read(20,*) M

   
    
    do j = 1, nbCorps
        write(*,*) "prout"
        read(20,*) (X(j,k),k=1,d)
    end do

    close(20)
    !write(*,*) M
    write(*,*) X(9,4)
    open(1, file='bodies_movement2D.csv',iostat=io_status)
    if (io_status /= 0) then
        write(*,*) 'Erreur lors de l''ouverture du fichier.'
        stop
    end if
    write(*,*) "kocessa"
    do i=1, Nstep
        t=i*dt
        
        !Update X for a time step dt
        call rk4(t,X,dt,nbCorps,M,d,deriv)
    
        write(1, '(*(G0.6,:,";"))', advance='no') ((X(b, a), a = 1, 2),b=1,nbCorps)!X(1,1), X(1,2), X(2,1), X(2,2), X(3,1), X(3,2), X(4,1), X(4,2), &
                                     !X(5,1), X(5,2), X(6,1), X(6,2), X(7,1), X(7,2), X(8,1), X(8,2), X(9,1), X(9,2)
        write(1,*)
    enddo 

    close(1)

end program DataCreation2D

!Calculate the derivative dX of the X matrix
subroutine deriv(t,nbCorps,M,X,dX)
    use Constant
    use System
    use Simulation
    implicit none
    integer :: nbCorps
    real(8), intent(in):: t
    real(8), dimension(nbCorps,d), intent(in):: X
    real(8), dimension(nbCorps) :: M
    real(8), dimension(nbCorps,d), intent(out)::dX !X's derivative
    real(8), dimension(nbCorps,nbCorps):: xforce, yforce, Xdis
   
    integer:: i,j
    
    call distance(nbCorps,M,X,Xdis)
    call force(nbCorps,M,X,Xdis,xforce,yforce)

    do i=1, nbCorps
        dX(i,1)=X(i,3) !dx/dt=v_x
        dX(i,2)=X(i,4) !dy/dt=v_y
        dX(i,3)=0
        dX(i,4)=0
        
        !add all the forces applied for the N bodies
        do j=1, nbCorps
            dX(i,3)=dX(i,3)+1/M(i)*xforce(i,j) !dv_x/dt=mi*a_x=G*mi*m1/ri1+G*mi*m2/ri2+...
            dX(i,4)=dX(i,4)+1/M(i)*yforce(i,j) !dv_y/dt=ma_y
        enddo
    enddo 

end subroutine deriv

!Calculate de Xdis matrix where Xdis(i,j) stands for 
!the distance between the i and j bodies.
subroutine distance(nbCorps,M,X,Xdis)
    use System
    implicit none
    integer, intent(in) :: nbCorps
    real(8), dimension(nbCorps), intent(in) :: M
    real(8), dimension(nbCorps,4), intent(in):: X
    real(8), dimension(nbCorps,nbCorps), intent(out):: Xdis
    integer:: i,j

    do i=1, nbCorps 
        do j=i, nbCorps
            Xdis(i,j)=sqrt((X(i,1)-X(j,1))**2+(X(i,2)-X(j,2))**2) 
            Xdis(j,i)=Xdis(i,j)
        enddo
    enddo

end subroutine distance

!Calculate the force between each body.
subroutine force(nbCorps,M,X,Xdis,xforce,yforce)
    use Constant
    use System
    implicit none
    integer, intent(in) :: nbCorps
    real(8), dimension(nbCorps), intent(in) :: M
    real(8), dimension(nbCorps,d), intent(in):: X
    real(8), dimension(nbCorps,nbCorps), intent(in):: Xdis
    real(8), dimension(nbCorps,nbCorps), intent(out):: xforce
    real(8), dimension(nbCorps,nbCorps), intent(out):: yforce
    integer:: i,j
    do i=1, nbCorps
        do j=i, nbCorps
            
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
subroutine rk4(t,X,dt,N,M,d,deriv)
    implicit none
    integer , intent (in) :: d,N
    real (8) , intent (in) :: t, dt
    real(8), dimension(N), intent(in) :: M
    real (8) , dimension (N,d) , intent ( inout ) :: X
    real (8) :: ddt
    real (8) , dimension (N,d) :: Xp , k1 , k2 , k3 , k4
    ddt = 0.5* dt
    call deriv (t,N,M,X,k1); Xp = X + ddt *k1
    call deriv (t+ddt,N,M,Xp ,k2); Xp = X + ddt*k2
    call deriv (t+ddt,N,M,Xp ,k3); Xp = X + dt*k3
    call deriv (t+dt,N,M,Xp ,k4); X = X + dt *( k1 + 2.0* k2 + 2.0* k3 + k4 )/6.0

end subroutine rk4