module Constant
    implicit none
    real(8):: G=6.67430*1e-11, eps=1e5, tolerance=1e24, SunMass=1.9885e30
end module Constant

subroutine SolarSystemSimulation(X, M, nbCorps, Nstep, dt)
    use Constant
    implicit none
    integer, intent(in):: nbCorps, Nstep                !Number of bodies/Number of steps for the simulation
    Real(8), intent(inout):: dt                         !time step
    Real(8), dimension(nbcorps,6), intent(inout):: X    !Position/velocity matrix 
    Real(8), intent(in), dimension(nbCorps):: M         !M(i)= masses of the bodies
    Real(8)::t,etotorb, epotorb, ecinorb            
    integer(8):: i,a,b
    external:: deriv

    open(1, file='SolarSystemTrajectory.dat')
    open(2, file='OrbEnergy.dat')
    do i=0, Nstep
        write(1,*) ((X(b, a), a = 1,3), b=1,nbCorps)
        !call rk4(t,X,dt,NbCorps,M,deriv,6)
        call OrbitalEnergy(Nbcorps,M,X,etotorb,ecinorb,epotorb,5)
        call adaptativerk4(t,X,dt,NbCorps,M,tolerance,6)
        write(2,*) t, etotorb, epotorb, ecinorb
    end do
    close(1)
    close(2) 

end subroutine SolarSystemSimulation

subroutine deriv(t,nbCorps,M,X,dX)
    use Constant
    implicit none
    integer :: nbCorps
    real(8), intent(in):: t
    real(8):: Xdis
    real(8), dimension(nbCorps,6), intent(in):: X
    real(8), dimension(nbCorps) :: M
    real(8), dimension(nbCorps,6), intent(out)::dX !X's derivative
    real(8), dimension(nbCorps,nbCorps):: xforce, yforce, zforce
   
    integer:: i,j

    call force(nbCorps,M,X,xforce,yforce,zforce)

    do i=1, nbCorps
        dX(i,1)=X(i,4) !dx/dt=v_x
        dX(i,2)=X(i,5) !dy/dt=v_y
        dX(i,3)=X(i,6)
        dX(i,4)=0
        dX(i,5)=0
        dX(i,6)=0
        
        
        !add all the forces applied for the N bodies
        do j=1, nbCorps
            dX(i,4)=dX(i,4)+1/M(i)*xforce(i,j) !dv_x/dt=mi*a_x=G*mi*m1/ri1+G*mi*m2/ri2+...
            dX(i,5)=dX(i,5)+1/M(i)*yforce(i,j) !dv_y/dt=ma_y
            dX(i,6)=dX(i,6)+1/M(i)*zforce(i,j)  
        enddo

         !  Xdis=sqrt(X(i,1)**2+X(i,2)**2+X(i,3)**2)
         !  dX(i,4)=dX(i,4)+SunMass/(((Xdis)**2+eps**2)**(1.5))
         !  dX(i,5)=dX(i,5)+SunMass/(((Xdis)**2+eps**2)**(1.5))
         !  dX(i,6)=dX(i,6)+SunMass/(((Xdis)**2+eps**2)**(1.5))
    enddo 

end subroutine deriv

subroutine force(nbCorps,M,X,xforce,yforce,zforce)
    use Constant
    implicit none
    integer, intent(in) :: nbCorps
    real(8), dimension(nbCorps), intent(in) :: M
    real(8), dimension(nbCorps,6), intent(in):: X
    real(8), dimension(nbCorps,nbCorps), intent(out):: xforce
    real(8), dimension(nbCorps,nbCorps), intent(out):: yforce
    real(8), dimension(nbCorps,nbCorps), intent(out):: zforce
    real(8):: Xdis
    integer:: i,j
    do i=1, nbCorps
        do j=i+1, nbCorps
    
            Xdis=sqrt((X(i,1)-X(j,1))**2+(X(i,2)-X(j,2))**2+(X(i,3)-X(j,3))**2)
            xforce(i,j)=G*M(i)*M(j)*(X(j,1)-X(i,1))/((Xdis**2+eps**2)**(1.5)) !newton's inverse-square law
            yforce(i,j)=G*M(i)*M(j)*(X(j,2)-X(i,2))/((Xdis**2+eps**2)**(1.5))
            zforce(i,j)=G*M(i)*M(j)*(X(j,3)-X(i,3))/((Xdis**2+eps**2)**(1.5))
                
            xforce(j,i)=-xforce(i,j) 
            yforce(j,i)=-yforce(i,j)
            zforce(j,i)=-zforce(i,j)
            
        enddo

    enddo

end subroutine force

subroutine energy(N,M,X,ecin,epot)
    use Constant
    implicit none
    integer :: N
    real(8), dimension(N,6), intent(in):: X
    real(8), dimension(N) :: M
    real(8), intent(out):: epot, ecin
    integer:: i,j

    ecin=0
    epot=0

    do i=1, N
        ecin=ecin+0.5*M(i)*(X(i,4)**2+X(i,5)**2+X(i,6)**2) !kinetic energy
        do j=i+1, N      
            epot=epot-G*M(i)*M(j)/sqrt((X(i,1)-X(j,1))**2+(X(i,2)-X(j,2))**2+(X(i,3)-X(j,3))**2)     !potential energy
        enddo
        !epot=epot-G*M(i)*SunMass/sqrt(X(i,1)**2+X(i,2)**2+X(i,3)**2)
    enddo

end subroutine energy

subroutine OrbitalEnergy(N,M,X,etotorb, ecinorb, epotorb,i)
    use Constant
    implicit none
    integer, intent(in) :: N,i
    real(8), dimension(N,6), intent(in):: X
    real(8), dimension(N) :: M
    real(8), intent(out):: etotorb
    real(8), intent(inout):: ecinorb, epotorb

    ecinorb=0.5*M(i)*(X(i,4)**2+X(i,5)**2+X(i,6)**2)

    epotorb=-G*M(i)*SunMass/sqrt((X(i,1))**2+(X(i,2))**2+(X(i,3))**2)
    etotorb=ecinorb+epotorb





end subroutine OrbitalEnergy

 
include 'integrator.f90'