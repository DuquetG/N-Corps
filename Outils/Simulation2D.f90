! Ce code fait appel à la méthode d'intégration de Runge-Kutta du 4eme ordre 
! afin de mesurer et inscrire dans un fichier les trajectoires et vitesses des masses d'un système à N-corps.
! Il peut également fournir les variations d'énergie et les énergies moyennes du système.

! Constantes fondamentales
module Constant
    implicit none
    real(8):: G=6.67430*1e-11, eps=1e5, tolerance=1e25
end module Constant

subroutine simulation2D(X, M, nbCorps, Nstep, dt, wtraj, format, wenergy, wviriel, wvelocity)
    use constant

    implicit none
    integer, intent(in):: nbCorps, Nstep                ! Nombre de corps, Nombre de pas pour la simulation
    Real(8), intent(inout):: dt                         ! Pas de temps
    Real(8), intent(in), dimension(nbCorps):: M         ! M(i) = masse du corps 'i'
    Real(8), intent(inout), dimension(nbCorps,4):: X    ! Matrice position/vitesse
    logical, intent(in):: wtraj                         ! Booléen, si traj='.true.' le programme inscrit les positions des corps dans le fichier de sortie
    logical, intent(in):: wenergy                       ! Booléen, si energy='.true.' le programme inscrit les variations d'énergie du système
    logical, intent(in):: wviriel                       ! Booléen, si wviriel=.true., le programme inscrit les énergies moyennes pour vérifier le théoreme du Viriel
    logical, intent(in):: wvelocity                     ! Booléen, si wvelocities=.true., le programme inscrit les vitesses des corps dans le fichier de sortie
    character(len=*), intent(in):: format               ! Format du fichier des positions et vitesses, '.csv' ou '.dat'


    integer:: i, io_status, a, b, u, v
    Real(8):: t=0, ecin, epot, ecinmoy=0, epotmoy=0
    Real(8), dimension(nbCorps,nbCorps):: Xdis
    external:: deriv 
    open(1, file='CSVs/positions_2D.csv',iostat=io_status)
    open(2, file='CSVs/positions_2D.dat',iostat=io_status)
    open(3, file='CSVs/energy.csv',iostat=io_status)
    open(4, file='CSVs/energy.dat',iostat=io_status)
    open(5, file='CSVs/viriel.csv',iostat=io_status)
    open(6, file='CSVs/velocities_2D.csv',iostat=io_status)

    if (format/='csv' .and. format/='dat') then
        write(*,*) 'Le format de sortie doit être .dat ou .csv'
        stop
    endif 

    if (io_status /= 0) then
        write(*,*) 'Erreur lors de l''ouverture du fichier.'
        stop
    end if
    do i=0, Nstep

        ! Choix du type d'intégration
        ! call adaptativerk4(t,X,dt,Nbcorps,M,tolerance)
        call rk4(t,X,dt,Nbcorps,M,deriv)

        if (wtraj) then

            if (format=='csv') then
                write(1, '(*(G0.6,:,";"))', advance='no') ((X(b, a), a = 1, 2), b=1,nbCorps)
                write(1,*)
            endif

            if (format=='dat') then
                write(2,*) ((X(b, a), a = 1, 2), b=1,nbCorps)
            endif 

            if (wvelocity) then
                write(6, '(*(G0.6,:,";"))', advance='no') ((X(v, u), u = 3, 4), v=1,nbCorps)
                write(6,*)
            endif

        endif

        if (wenergy) then

            call energy(nbCorps,M,X,ecin,epot)

            if (format=='csv') then
                write(3, '(*(G0.6,:,";"))', advance='no') ecin, epot, ecin+epot, t
                write(3,*)

            endif

            if (format=='dat') then
                write(4,*) ecin, epot, ecin+epot, t
            endif 


            if (wviriel) then
                ecinmoy=(ecinmoy*(i)+ecin)/(i+1)
                epotmoy=(epotmoy*(i)+epot)/(i+1)
                write(5, '(*(G0.6,:,";"))', advance='no') ecinmoy, epotmoy, i*dt
                write(5,*)

            endif
            

        end if
    end do

    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(6)

end subroutine simulation2D 

! Calcule l'énergie du système
subroutine energy(N,M,X,ecin,epot)
    use Constant
    implicit none
    integer :: N
    real(8), dimension(N,4), intent(in):: X
    real(8), dimension(N) :: M
    real(8), intent(out):: epot, ecin
    integer:: i,j

    ecin=0
    epot=0

    do i=1, N
        ecin=ecin+0.5*M(i)*(X(i,3)**2+X(i,4)**2) !kinetic energy
        do j=i+1, N      
            epot=epot-G*M(i)*M(j)/sqrt((X(i,1)-X(j,1))**2+(X(i,2)-X(j,2))**2)     !potential energy
        enddo
    enddo

end subroutine energy

! Calcule la dérivé dX de la matrice X
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
    
    call force(nbCorps,M,X,xforce,yforce)

    do i=1, nbCorps
        dX(i,1)=X(i,3) !dx/dt=v_x
        dX(i,2)=X(i,4) !dy/dt=v_y
        dX(i,3)=0
        dX(i,4)=0
        
        ! Ajoute toutes les forces appliquées pour les N-corps
        do j=1, nbCorps
            dX(i,3)=dX(i,3)+1/M(i)*xforce(i,j) !dv_x/dt=mi*a_x=G*mi*m1/ri1+G*mi*m2/ri2+...
            dX(i,4)=dX(i,4)+1/M(i)*yforce(i,j) !dv_y/dt=ma_y
        enddo
    enddo 

end subroutine deriv

! Calcule la matrice Xdis dans laquelle Xdis(i,j) représente la distance entre les corps i et j
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

! Calcule la force entre chaque corps
subroutine force(nbCorps,M,X,xforce,yforce)
    use Constant
    implicit none
    integer, intent(in) :: nbCorps
    real(8), dimension(nbCorps), intent(in) :: M
    real(8), dimension(nbCorps,4), intent(in):: X
    real(8), dimension(nbCorps,nbCorps), intent(out):: xforce
    real(8), dimension(nbCorps,nbCorps), intent(out):: yforce
    real(8):: Xdis
    integer:: i,j
    do i=1, nbCorps
        do j=i, nbCorps
            
            if (i==j) then
                xforce(i,j)=0  ! La force du corps appliquée à lui-même est nulle
                yforce(i,j)=0  
            else if (i/=j) then
            Xdis=sqrt((X(i,1)-X(j,1))**2+(X(i,2)-X(j,2))**2)
            xforce(i,j)=G*M(i)*M(j)*(X(j,1)-X(i,1))/((Xdis**2+eps**2)**(1.5)) ! newton's inverse-square law
            yforce(i,j)=G*M(i)*M(j)*(X(j,2)-X(i,2))/((Xdis**2+eps**2)**(1.5))

            xforce(j,i)=-xforce(i,j) 
            yforce(j,i)=-yforce(i,j)
            end if
        enddo
    enddo

end subroutine

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
    t=t+dt
end subroutine adaptativerk4

