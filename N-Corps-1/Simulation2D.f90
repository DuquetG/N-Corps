!Constantes fondamentales
module Constant
    implicit none
    real(8):: G=6.67430*1e-11
end module Constant

subroutine simulation2D(X, M, nbCorps, Nstep, dt, wtraj, format, wenergy,integrator)
    implicit none
    integer, intent(in):: nbCorps, Nstep                !Number of bodies/Number of steps for the simulation
    Real(8), intent(in):: dt                            !time step
    Real(8), intent(in), dimension(nbCorps):: M         !M(i)= mass of the body 'i'
    Real(8), intent(inout), dimension(nbCorps,4):: X    !Position/velocity matrix 
    logical, intent(in):: wtraj                         !boolean, if traj='.true.' the program will output the trajectory of the bodies 
    logical, intent(in):: wenergy                       !boolean, if energy='.true.' the program will output the energy fluctuation of the system
    character(len=*), intent(in):: format               !trajectory's output format, 'csv' or 'dat'.
    character(len=*), intent(in):: integrator           !integration methode, 'euler', 'rk4', 

    integer:: i, io_status, a, b
    Real(8):: t, ecin, epot
    Real(8), dimension(nbCorps,nbCorps):: Xdis
    external:: deriv 
    write(*,*) format
    open(1, file='bodies_movement2D.csv',iostat=io_status)
    open(2, file='bodies_movement2D.dat')
    open(3, file='energy.dat')

    if (format/='csv' .and. format/='dat') then
        write(*,*) 'Le format de sortie doit être .dat ou .csv'
        stop
    endif 

    if (io_status /= 0) then
        write(*,*) 'Erreur lors de l''ouverture du fichier.'
        stop
    end if

    do i=0, Nstep
        call rk4(t,X,dt,nbCorps,M,deriv)

        if (wtraj) then

            if (format=='csv') then
                write(1, '(*(G0.6,:,";"))', advance='no') ((X(b, a), a = 1, 2), b=1,nbCorps)
            endif

            if (format=='dat') then
                write(2,*) ((X(b, a), a = 1, 2), b=1,nbCorps)
            endif 

        endif

        if (wenergy) then

            call distance(nbCorps,X,Xdis)
            call energy(nbCorps,M,X,Xdis,ecin,epot)
            write(3,*) ecin, epot, ecin+epot

        end if
    end do

    close(1)
    close(2)

end subroutine simulation2D 

!Calculate the energy of the system
subroutine energy(N,M,X,Xdis,ecin,epot)
    use Constant
    implicit none
    integer :: N
    real(8), dimension(N,4), intent(in):: X
    real(8), dimension(N) :: M
    real(8), dimension(N,N), intent(in):: Xdis
    real(8), intent(out):: epot, ecin
    integer:: i,j

    ecin=0
    epot=0

    do i=1, N
        ecin=ecin+0.5*M(i)*(X(i,3)**2+X(i,4)**2) !kinetic energy
        do j=i+1, N      
            epot=epot-G*M(i)*M(j)/Xdis(i,j)      !potential energy
        enddo
    enddo

end subroutine energy

!Calculate the derivative dX of the X matrix
subroutine deriv(t,nbCorps,M,X,dX)
    use Constant
    implicit none
    integer :: nbCorps
    real(8), intent(in):: t
    real(8), dimension(nbCorps,4), intent(in):: X
    real(8), dimension(nbCorps) :: M
    real(8), dimension(nbCorps,4), intent(out)::dX !X's derivative
    real(8), dimension(nbCorps,nbCorps):: xforce, yforce, Xdis
   
    integer:: i,j
    
    call distance(nbCorps,X,Xdis)
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
subroutine distance(nbCorps,X,Xdis)
    implicit none
    integer, intent(in) :: nbCorps
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
    implicit none
    integer, intent(in) :: nbCorps
    real(8), dimension(nbCorps), intent(in) :: M
    real(8), dimension(nbCorps,4), intent(in):: X
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

include 'integrateur.f90'