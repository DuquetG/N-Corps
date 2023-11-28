program DataCreation2D
    implicit none
    integer, parameter:: Nstep=1000000 !number of steps considered for the simulation
    real(8), parameter:: dt=10000.      !time step
    integer, parameter:: dim=2         !number of dimension for the problem (currently only support dim=2)
    integer, parameter:: d=2*dim
    logical, parameter:: wtraj=.true., wenergy=.true.
    character(len=*), parameter::  format='csv'

    !Initialize position and velocity
    integer :: nbCorps, io_status
    character(len=100) :: line, preset_name, filename
    real(8), allocatable :: X(:,:)         !Position/velocity matrix
    real(8), allocatable :: M(:)           !Mass of bodies
    real(8), allocatable:: Xdis(:,:)       !Xdis(i,j)= distance between i and j

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

    write(*,*) "prouton 3000"
    
    
   
    close(20)

    call simulation2D(X, M, nbCorps, Nstep, dt, wtraj, format, wenergy)

end program DataCreation2D

include 'Simulation2D.f90'