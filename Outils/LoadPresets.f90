! Prend un fichier .txt selon un format précis pour permettre l'usage de conditions initiales préréglées 

program DataCreation2D
    implicit none
    integer, parameter:: Nstep=100000     ! number of steps considered for the simulation
    real(8):: dt=10000.                   ! Pas de temps
    integer, parameter:: dim=2            ! Nombre de dimensions du système (un maximum de 2 dimensions sont supportées actuellement)
    integer, parameter:: d=2*dim
    logical, parameter:: wtraj=.true., wenergy=.true., wviriel=.true., wvelocity=.true.
    character(len=*), parameter::  format='csv'

    ! Initialise les positions et vitesses
    integer :: nbCorps, io_status
    character(len=100) :: line, preset_name, filename
    real(8), allocatable :: X(:,:)         ! Matrice des positions et vitesses
    real(8), allocatable :: M(:)           ! Masses des corps
    real(8), allocatable:: Xdis(:,:)       ! Xdis(i,j)= distance entre i et j

    integer :: i,j,k
    

    if (command_argument_count() == 1) then
        call get_command_argument(1, filename)
        
    else
        write(*,*) 'Entrez le nom du fichier de presets : '
        read(*,*) filename
    end if
    

    open(unit=20, file=filename, status='old', action='read', iostat=io_status)
    

    if (io_status /= 0) then
        write(*,*) "Erreur lors de l'ouverture du fichier."
        stop
    end if

   

    read(20,*) nbCorps

    allocate(X(nbCorps,d))
    allocate(M(nbCorps))
    allocate(Xdis(nbCorps,nbCorps))

    read(20,*) M
    do j = 1, nbCorps
        read(20,*) (X(j,k),k=1,d)
    end do

    
    
   
    close(20)

    call simulation2D(X, M, nbCorps, Nstep, dt, wtraj, format, wenergy, wviriel, wvelocity)

end program DataCreation2D

include 'Simulation2D.f90'
