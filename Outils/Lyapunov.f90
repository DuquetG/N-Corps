

subroutine Lyapunov(nb_simulations, Nstep, nbCorps, mass ,chosen_mass, dt, length, pasCalcul)
implicit none
    integer, intent(in) :: nb_simulations, Nstep, nbCorps, chosen_mass, length, pasCalcul
    real(8), dimension(nbCorps), intent(in) :: mass
    integer :: line, k, i, j, zoinks,l, m, n
    real(8), allocatable :: Xstep(:), Vstep(:), delta_zeros(:), delta_tees(:), lyp(:)
    real(8) :: temps, dt, delZero2, del2

    i = 1
    j = 1
    k=1
    l=1

    open(1, file='CSVs/cured_all_positions_2D.csv')  
    open(3, file='CSVs/cured_all_velocities_2D.csv')      
    open(2, file='CSVs/Lyapunov.csv')

    allocate(Xstep(nbCorps * 2 * nb_simulations))
    allocate(Vstep(nbCorps * 2 * nb_simulations))
    allocate(delta_zeros(nb_simulations-1))
    allocate(delta_tees(nb_simulations-1))
    allocate(lyp(0:nb_simulations-1))

    read(1,*) (Xstep(k), k = 1, nbCorps * 2 * nb_simulations)
    read(3,*) (Vstep(l), l = 1, nbCorps * 2 * nb_simulations)
  
    !Calculating delta zero
    do i = 1, nb_simulations-1
        
        m=0
        delZero2=0
        do m=1,nbCorps
            delZero2=delZero2+(Xstep(nbCorps*2*i+2*m-1)-Xstep(2*m-1))**2 &
            + (Xstep(nbCorps*2*i+2*m)-Xstep(2*m))**2 &
           + ((Vstep(nbCorps*2*i+2*m-1)-Vstep(2*m-1)))**2 &
           + ((Vstep(nbCorps*2*i+2*m)-Vstep(2*m)))**2
        enddo

        delta_zeros(i) = sqrt(delZero2)
    enddo

    

    do line = 1, length/nb_simulations-1
        read(1,*) (Xstep(k), k = 1, nbCorps * 2 * nb_simulations)
        read(3,*) (Vstep(l), l = 1, nbCorps * 2 * nb_simulations)
        temps = line * dt * pasCalcul
        do j = 1, nb_simulations-1
            n=0
            del2=0
            do n=1,nbCorps
                del2=del2+(Xstep(nbCorps*2*j+2*n-1)-Xstep(2*n-1))**2 &
                + (Xstep(nbCorps*2*j+2*n)-Xstep(2*n))**2 &
               + ((Vstep(nbCorps*2*j+2*n-1)-Vstep(2*n-1)))**2 &
               + ((Vstep(nbCorps*2*j+2*n)-Vstep(2*n)))**2
            enddo
            delta_tees(j) = sqrt(del2)
            lyp(j) = (1/temps) * log(delta_tees(j)/delta_zeros(j))
        enddo
        write(2,'(*(G0.6,:,";"))', advance='no') temps, (lyp(zoinks), zoinks = 1, nb_simulations-1)
        write(2,*)
    enddo

    close(1)
    close(2)
    close(3)

end subroutine Lyapunov