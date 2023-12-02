

subroutine Lyapunov(nb_simulations, Nstep, nbCorps, chosen_mass, dt, length)
implicit none
    integer, intent(in) :: nb_simulations, Nstep, nbCorps, chosen_mass, length
    integer :: line, k, lol, lmfao, zoinks,l, m, n
    real(8), allocatable :: X(:,:)
    real(8), allocatable :: Xstep(:), Vstep(:), delta_zeros(:), delta_tees(:), lyp(:)
    real(8) :: temps, dt, delZero2, del2

    lol = 1
    lmfao = 1

    open(1, file='CSVs/all_positions_2D.csv')  
    open(3, file='CSVs/all_velocities_2D.csv')      
    open(2, file='CSVs/Lyapunov.csv')

    allocate(X(Nstep, nbCorps * 2 * nb_simulations))
    allocate(Xstep(nbCorps * 2 * nb_simulations))
    allocate(Vstep(nbCorps * 2 * nb_simulations))
    allocate(delta_zeros(nb_simulations))
    allocate(delta_tees(nb_simulations))
    allocate(lyp(nb_simulations))

    read(1,*) (Xstep(k), k = 1, nbCorps * 2 * nb_simulations)
    read(3,*) (Vstep(k), l = 1, nbCorps * 2 * nb_simulations)
    do lol = 1, nb_simulations
        m=0
        delZero2=0
        do m=0,nbCorps-1
            delZero2=delZero2+(Xstep(nbCorps*2*lol+2*m)-Xstep(2*m))**2 &
            + (Xstep(nbCorps*2*lol+2*m+1)-Xstep(2*m+1))**2 &
            + (Vstep(nbCorps*2*lol+2*m)-Vstep(2*m))**2 &
            + (Vstep(nbCorps*2*lol+2*m+1)-Vstep(2*m+1))**2
        enddo

        delta_zeros(lol) = sqrt(delZero2)
    enddo

    do line = 1, length/nb_simulations-1
        read(1,*) (Xstep(k), k = 1, nbCorps * 2 * nb_simulations)
        read(3,*) (Vstep(k), l = 1, nbCorps * 2 * nb_simulations)
        temps = line * dt
        do lmfao = 1, nb_simulations
            n=0
            del2=0
            do n=0,nbCorps-1
                del2=del2+(Xstep(nbCorps*2*lmfao+2*n)-Xstep(2*n))**2 &
                + (Xstep(nbCorps*2*lmfao+2*n+1)-Xstep(2*n+1))**2 &
                + (Vstep(nbCorps*2*lmfao+2*n)-Vstep(2*n))**2 &
                + (Vstep(nbCorps*2*lmfao+2*n+1)-Vstep(2*n+1))**2
            enddo
            delta_tees = sqrt(del2)
            lyp(lmfao) = (1/temps) * log(delta_tees(lmfao)/delta_zeros(lmfao))
        enddo
        write(2,'(*(G0.6,:,";"))', advance='no') temps, (lyp(zoinks), zoinks = 1, nb_simulations)
        write(2,*)
    enddo

    close(1)
    close(2)
    close(3)

end subroutine Lyapunov