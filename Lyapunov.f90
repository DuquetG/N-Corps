

subroutine Lyapunov(nb_simulations, Nstep, nbCorps, chosen_mass, dt, length)
implicit none
    integer, intent(in) :: nb_simulations, Nstep, nbCorps, chosen_mass, length
    integer :: line, k, lol, lmfao, zoinks,l
    real(8), allocatable :: X(:,:)
    real(8), allocatable :: Xstep(:), Vstep(:), delta_zeros(:), delta_tees(:), lyp(:)
    real(8) :: temps, dt

    lol = 1
    lmfao = 1

    open(1, file='Cured_Triangle.csv')  
    open(3, file='Cured_Triangle_v.csv')      
    open(2, file='Lyapunov.csv')

    allocate(X(Nstep, nbCorps * 2 * nb_simulations))
    allocate(Xstep(nbCorps * 2 * nb_simulations))
    allocate(Vstep(nbCorps * 2 * nb_simulations))
    allocate(delta_zeros(nb_simulations))
    allocate(delta_tees(nb_simulations))
    allocate(lyp(nb_simulations))

    read(1,*) (Xstep(k), k = 1, nbCorps * 2 * nb_simulations)
    read(3,*) (Vstep(k), l = 1, nbCorps * 2 * nb_simulations)
    do lol = 1, nb_simulations
        delta_zeros(lol) = sqrt((Xstep(nbCorps*2*lol+chosen_mass)-Xstep(chosen_mass))**2 &
                                 + (Xstep(nbCorps*2*lol+chosen_mass+1)-Xstep(chosen_mass+1))**2&
                                 +(Vstep(nbCorps*2*lol+chosen_mass)-Vstep(chosen_mass))**2 &
                                 + (Vstep(nbCorps*2*lol+chosen_mass+1)-Vstep(chosen_mass+1))**2)
    enddo

    do line = 1, length/nb_simulations-1
        read(1,*) (Xstep(k), k = 1, nbCorps * 2 * nb_simulations)
        read(3,*) (Vstep(k), l = 1, nbCorps * 2 * nb_simulations)
        temps = line * dt
        do lmfao = 1, nb_simulations
            delta_tees = sqrt((Xstep(nbCorps*2*lmfao+chosen_mass)-Xstep(chosen_mass))**2 &
                            + (Xstep(nbCorps*2*lmfao+chosen_mass+1)-Xstep(chosen_mass+1))**2 &
                            + (Vstep(nbCorps*2*lmfao+chosen_mass)-Vstep(chosen_mass))**2 &
                            + (Vstep(nbCorps*2*lmfao+chosen_mass+1)-Vstep(chosen_mass+1))**2)
            lyp(lmfao) = (1/temps) * log(delta_tees(lmfao)/delta_zeros(lmfao))
        enddo
        write(2,'(*(G0.6,:,";"))', advance='no') temps, (lyp(zoinks), zoinks = 1, nb_simulations)
        write(2,*)
    enddo

    close(1)
    close(2)
    close(3)

end subroutine Lyapunov