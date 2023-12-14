module Constant
    implicit none
    real(8):: G=0.00001 ! G=6.67430*1e-11
end module Constant


! Généralisation 3D pratiquement identique au programme homologue à 2D
subroutine Simulation3D(X, M, nbCorps, Nstep, dt) 
    implicit none
    integer:: nbCorps, Nstep
    integer:: i, io_status, a, b
    Real(8):: dt,t, ecin, epot
    Real(8), dimension(nbCorps):: M
    Real(8), dimension(nbCorps,6):: X
    Real(8), dimension(nbCorps,nbCorps):: Xdis
    external:: deriv 

    open(1, file='bodies_movement2D.csv',iostat=io_status)
    open(2, file='energy.dat')

    if (io_status /= 0) then
        write(*,*) 'Erreur lors de l''ouverture du fichier.'
        stop
    end if

    do i=0, Nstep
        call rk4(t,X,dt,nbCorps,M,deriv)
        write(1, '(*(G0.6,:,";"))', advance='no') ((X(b, a), a = 1, 2, 3), b=1,nbCorps)
        
        if (mod(i,50)==0) then

            call distance(nbCorps,X,Xdis)
            call energy(nbCorps,M,X,Xdis,ecin,epot)
            write(2,*) ecin, epot, ecin+epot

        end if
    end do

    close(1)
    close(2)
end subroutine Simulation3D

subroutine energy(N,M,X,Xdis,ecin,epot)
    use Constant
    implicit none
    integer :: N
    real(8), dimension(N,6), intent(in):: X
    real(8), dimension(N) :: M
    real(8), dimension(N,N), intent(in):: Xdis
    real(8), intent(out):: epot, ecin
    integer:: i,j

    ecin=0
    epot=0

    do i=1, N
        ecin=ecin+0.5*M(i)*(X(i,4)**2+X(i,5)**2+X(i,6)**2) ! Énergie cinétique
        do j=i+1, N      
            epot=epot-G*M(i)*M(j)/Xdis(i,j)      ! Énergie potentielle
        enddo
    enddo

end subroutine energy

subroutine deriv(t,nbCorps,M,X,dX)
    use Constant
    implicit none
    integer, intent(in):: nbCorps
    real(8), intent(in):: t
    real(8), dimension(nbCorps) :: M
    real(8), dimension(Nbcorps,6), intent(in):: X
    real(8), dimension(Nbcorps,6), intent(out)::dX
    real(8), dimension(Nbcorps,Nbcorps):: xforce, yforce, zforce, Xdis
    integer:: i,j
    
    call distance(nbCorps,X,Xdis)
    call force(nbCorps,M,X,Xdis,xforce,yforce,zforce)

    do i=1, Nbcorps
        dX(i,1)=X(i,4)
        dX(i,2)=X(i,5)
        dX(i,3)=X(i,6)
        dX(i,4)=0
        dX(i,5)=0
        dX(i,6)=0


        do j=1, Nbcorps
            dX(i,4)=dX(i,4)+1/M(i)*xforce(i,j)
            dX(i,5)=dX(i,5)+1/M(i)*yforce(i,j)
            dX(i,6)=dX(i,6)+1/M(i)*zforce(i,j)
        enddo
    enddo 

end subroutine deriv

subroutine distance(nbCorps,X,Xdis)
    implicit none
    integer, intent(in):: nbCorps
    real(8), dimension(nbCorps,6), intent(in):: X
    real(8), dimension(nbCorps,nbCorps), intent(out):: Xdis
    integer:: i,j

    do i=1, nbCorps
        do j=i, nbCorps
            Xdis(i,j)=sqrt((X(i,1)-X(j,1))**2+(X(i,2)-X(j,2))**2+(X(i,3)-X(j,3))**2) 
            Xdis(j,i)=Xdis(i,j)
        enddo
    enddo

end subroutine distance

subroutine force(nbCorps,M,X,Xdis,xforce,yforce,zforce)
    use Constant
    implicit none
    integer, intent(in):: nbCorps
    real(8), dimension(nbCorps,6), intent(in):: X
    real(8), dimension(nbCorps) :: M
    real(8), dimension(nbCorps,nbCorps), intent(in):: Xdis
    real(8), dimension(nbCorps,nbCorps), intent(out):: xforce, yforce, zforce
    integer:: i,j
    do i=1, nbCorps
        do j=i, NbCorps
            
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

subroutine rk4(t,X,dt,N,M,deriv)
    implicit none
    integer , intent (in) :: N
    real (8) , intent (in) :: t, dt
    real(8), dimension(N), intent(in) :: M
    real (8) , dimension (N,6) , intent ( inout ) :: X
    real (8) :: ddt
    real (8) , dimension (N,6) :: Xp , k1 , k2 , k3 , k4
    ddt = 0.5* dt
    call deriv (t,N,M,X,k1); Xp = X + ddt *k1
    call deriv (t+ddt,N,M,Xp,k2); Xp = X + ddt*k2
    call deriv (t+ddt,N,M,Xp,k3); Xp = X + dt*k3
    call deriv (t+dt,N,M,Xp ,k4); X = X + dt *( k1 + 2.0* k2 + 2.0* k3 + k4 )/6.0

end subroutine rk4
