!===============================================================================
!> \file source_term.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is the subroutine that control the source terms in the equations.
!! \details
!! Contains source_term(), gravity_predictor()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          25-03-2015 
!! \b last \b modified: 25-03-2015
!<
!===============================================================================
!> Compute source term
!===============================================================================
subroutine source_term
  use params
  use variables
  use atmosphere
  use mpi_var
  implicit none

  integer :: i,j,k

  real(dp) :: rhon,rhonp1

  !local variables for balanced scheme
  real(dp) :: r,p,eke
  real(dp) :: rm,rp,pm,pp

  !local variables for cooling function
  real(dp) :: rho,rhovx,rhovy,rhovz,Bxc,Byc,Bzc,E,csq,T,Teq,pressure
  real(dp) :: Q,oneOverTrad,rhoeq,fij

  !local variables for sponge layer
  real(dp) :: Rw,C,zeta,zetaS,Eken

  !local variables for Coriolis force
  real(dp) :: f,lambda

  !local variables for 2D jet
  real(dp) :: Ujet,tauJet,dpJet,dEJet

  if (verbose) write (*,*) 'Entering source_term subroutine...'

!  if ((zposition==nzslice-1).and.(yposition==0).and.(mype==15)) then
!     write (*,*) 'yo',mype,uin(:,1,nz,4)
!  endif

  !Gravity term
  do k=1,nz
     do j=1,ny
        do i=1,nx

           !update momentum
           if (balanced) then

              ! Compute pressure from total energy
              rho=uin_old(i,j,k,1)
              rhovx=uin_old(i,j,k,2) ; rhovy=uin_old(i,j,k,3) ; rhovz=uin_old(i,j,k,4)
              Bxc=half*(uin_old(i,j,k,6)+uin_old(i,j,k,9 ))
              Byc=half*(uin_old(i,j,k,7)+uin_old(i,j,k,10))
              Bzc=half*(uin_old(i,j,k,8)+uin_old(i,j,k,11))
              E=uin_old(i,j,k,5)
              call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,p)

              if (ndim==3) then
                 pm=0.d0 ; pp=0.d0
                 if (z(k)<zbuffer) call getPressureFaces(rho,p,pm,pp,gamma,dz/2.d0)
                 uin(i,j,k,4)=uin(i,j,k,4)+(pm-pp)/dz*dt
              else if (ndim==1) then
                 call getPressureFaces(rho,p,pm,pp,gamma,dx/2.d0)
                 uin(i,j,k,2)=uin(i,j,k,2)+(pm-pp)/dx*dt
              endif
              

           else

              !density
              rhon  =uin_old(i,j,k,1)
              rhonp1=    uin(i,j,k,1)

              uin(i,j,k,2)=uin(i,j,k,2)+half*(rhon+rhonp1)*gravin(i,j,k,1)*dt
              if (ndim>1) uin(i,j,k,3)=uin(i,j,k,3)+half*(rhon+rhonp1)*gravin(i,j,k,2)*dt
              if (ndim>2) uin(i,j,k,4)=uin(i,j,k,4)+half*(rhon+rhonp1)*gravin(i,j,k,3)*dt

           endif

#if ISO==0
           !update energy (1st order in z-direction)
           if (ndim==3) then
              if (z(k)<zbuffer) uin(i,j,k,5) = uin(i,j,k,5) + uin_old(i,j,k,4)*gravin(i,j,k,3)*dt
           else if (ndim==1) then
              uin(i,j,k,5) = uin(i,j,k,5) + uin_old(i,j,k,2)*gravin(i,j,k,1)*dt
           endif
#endif
        end do
     end do
  end do

!  if ((zposition==nzslice-1).and.(yposition==0).and.(mype==15)) then
!     write (*,*) mype,uin(:,1,nz,4)
!  endif
  !stop

  !    ! Coriolis force 
  !    if (ndim==3) then
  !       do k=1,nz
  !          do j=1,ny
  !             do i=1,nx

  !                   f=2.d0*(omegaEq+half*beta*y(j))

  !                   !Euler naive scheme - 16/09/13
  !                   uin(i,j,k,2)=uin(i,j,k,2)+f*uin_old(i,j,k,3)*dt
  !                   uin(i,j,k,3)=uin(i,j,k,3)-f*uin_old(i,j,k,2)*dt

  ! !                  !Crank-Nicholson scheme - 07/10/13
  ! !                  !lambda=f*dt/2.d0
  ! !                  !uin(i,j,k,2)=(1.d0-lambda**2)*uin(i,j,k,2)+2.d0*lambda*uin(i,j,k,3)
  ! !                  !uin(i,j,k,3)=(1.d0-lambda**2)*uin(i,j,k,3)-2.d0*lambda*uin(i,j,k,2)
  ! !                  !uin(i,j,k,2)=uin(i,j,k,2)/(1.d0+lambda**2)
  ! !                  !uin(i,j,k,3)=uin(i,j,k,3)/(1.d0+lambda**2)

  ! !                  !Matrix rotation - 07/10/13
  ! !                  lambda=2.d0*(omegaEq+half*beta*y(j))*dt
  ! !                  uin(i,j,k,2)=+cos(lambda)*uin(i,j,k,2)+sin(lambda)*uin(i,j,k,3)
  ! !                  uin(i,j,k,3)=-sin(lambda)*uin(i,j,k,2)+cos(lambda)*uin(i,j,k,3)

  !             end do
  !          end do
  !       end do
  !    end if

  ! Rayleigh drag
  if (drag) then
     uin(:,:,:,2:4) = uin(:,:,:,2:4) - uin(:,:,:,2:4)/t_drag*dt
  endif

    ! Sponge layer for deep model
!     if (sponge) then
! !       zetaS=7.5d-1 ; C=1.d-1 !ok for deep model with nx=1024 (0.2 in mayne) -> for explicit scheme
!        zetaS=9.d-1 ; C=1.d-1 !ok for deep model with nx=1024 (0.2 in mayne) -> for explicit scheme
!   !     etaS=6.d-1 ; C=0.d-1 !(0.2 in mayne)
!        do k=1,nz
!           zeta=z(k)/(zmax-zmin) ; Rw=0.d0
!           if (zeta>zetaS) Rw=C*(sin(pi/2.d0*(zeta-zetaS)/(1.d0-zetaS)))**2
!           do j=1,ny
!              do i=1,nx
!                 Eken=5.d-1*uin(i,j,k,4)**2/uin(i,j,k,1)
!                 uin(i,j,k,5)=uin(i,j,k,5)-Eken
!                 uin(i,j,k,4)=uin(i,j,k,4)/(1.d0+Rw*dt)
!                 Eken=5.d-1*uin(i,j,k,4)**2/uin(i,j,k,1)
!                 uin(i,j,k,5)=uin(i,j,k,5)+Eken
!              end do
!           end do
!        end do
!     endif

    ! Sponge layer for deep model for all velocity components
    if (sponge) then
       !zetaS=9.d-1 ; C=1.d-4 !ok for deep model with nx=1024 (0.2 in mayne) -> for explicit scheme
       zetaS=8.d-1 ; C=1.d-3 !ok for deep model with nx=1024 (0.2 in mayne) -> for explicit scheme
       do k=1,nz
          zeta=z(k)/(zmax-zmin) ; Rw=0.d0
          if (zeta>zetaS) Rw=C*(sin(pi/2.d0*(zeta-zetaS)/(1.d0-zetaS)))**2
          do j=1,ny
             do i=1,nx
                Eken=5.d-1*(uin(i,j,k,2)**2+uin(i,j,k,3)**2+uin(i,j,k,4)**2)/uin(i,j,k,1)
                uin(i,j,k,5)=uin(i,j,k,5)-Eken
                uin(i,j,k,2)=uin(i,j,k,2)/(1.d0+Rw*dt)
                uin(i,j,k,3)=uin(i,j,k,3)/(1.d0+Rw*dt)
                uin(i,j,k,4)=uin(i,j,k,4)/(1.d0+Rw*dt)
                Eken=5.d-1*(uin(i,j,k,2)**2+uin(i,j,k,3)**2+uin(i,j,k,4)**2)/uin(i,j,k,1)
                uin(i,j,k,5)=uin(i,j,k,5)+Eken
             end do
          end do
       end do
    endif

  ! Newtonian cooling
  if (cooling) then
     do k=1,nz
        do j=1,ny
           do i=1,nx

              ! Compute local temperature
              rho=uin(i,j,k,1) ; E=uin(i,j,k,5)
              rhovx=uin(i,j,k,2) ; rhovy=uin(i,j,k,3) ; rhovz=uin(i,j,k,4)
              Bxc=half*(uin(i,j,k,6)+uin(i+1,j  ,k  ,6 ))
              Byc=half*(uin(i,j,k,7)+uin(i  ,j+1,k  ,7))
              Bzc=half*(uin(i,j,k,8)+uin(i  ,j  ,k+1,8))
              call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,pressure)
              T=pressure/k_mumh/rho

              ! Compute equilibrium temperature
              if (ndim==3) then
                 if (deep) then
                    call getTeq(Teq,x(i),y(j),z(k),pressure)
                    call getTauRad(oneOverTrad,pressure)
                 else
                    call getTeq(Teq,x(i),y(j),z(k),Tsurf,Dtrop,zstra,dTstra,DT_EP,i,j,k)
                    oneOverTrad=1.d0/t_rad
                 endif
              else if (ndim==1) then
                 Teq=Tsurf-Dtrop*(x(i)-xmin)
              endif

              ! Update total energy
!!! Explicit scheme
              !uin(i,j,k,5) = uin(i,j,k,5) + k_mumh*rho*(Teq-T)*oneOverTrad/(gamma-1.d0)*dt
!!! Implicit scheme
              T=(T+oneOverTrad*dt*Teq)/(1.d0+oneOverTrad*dt)
              call get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,gamma*k_mumh*T)
              uin(i,j,k,5)=E

           end do
        end do
     end do
  endif

  ! Shallow water forcing
  if (forcing) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! ! Explicit scheme (fails for large forcing, not sure momentum is correct [Seb, 09/05/14]...)
              ! ! Compute forcing term
              ! Q=(rho0*(1.d0+Q0*exp(-y(j)*y(j)/2.d0)*cos(x(i)-pi))-uin(i,j,k,1))/t_rad
              ! ! Update mass and velocities
              ! uin(i,j,k,2) = uin(i,j,k,2) - (max(Q,0.d0)/uin(i,j,k,1)**2 + 1.d0/t_drag)*uin(i,j,k,2)*dt
              ! uin(i,j,k,3) = uin(i,j,k,3) - (max(Q,0.d0)/uin(i,j,k,1)**2 + 1.d0/t_drag)*uin(i,j,k,3)*dt
              ! uin(i,j,k,1) = uin(i,j,k,1) + Q*dt
              ! ! Compute total energy
              ! uin(i,j,k,5) = 5.d-1*(uin(i,j,k,2)**2+uin(i,j,k,3)**2)/uin(i,j,k,1)+half*uin(i,j,k,1)**2/(gamma-1.d0)
              !Implicit scheme
              !Compute forcing term
              uin(i,j,k,2)=uin(i,j,k,2)/uin(i,j,k,1) !velocity from momentum
              uin(i,j,k,3)=uin(i,j,k,3)/uin(i,j,k,1) !velocity from momentum
              !fij=(1.d0+Q0*exp(-y(j)*y(j)/2.d0/lytherm**2)*cos(x(i)-pi))
              fij=1.d0+Q0*exp(-y(j)*y(j)/2.d0/lytherm**2)*cos(2*pi*x(i)/(xmax-xmin))
              !Update density
              rhoeq=rho0*fij
              uin(i,j,k,1)=(uin(i,j,k,1)+rhoeq*dt/t_rad)/(1.d0+dt/t_rad) ; rhonp1=uin(i,j,k,1)
              ! Update mass and velocities
              Q=(rhoeq-rhonp1)/t_rad
              uin(i,j,k,2)=uin(i,j,k,2)/(1.d0+max(Q,0.d0)/rhonp1*dt+dt/t_drag)
              uin(i,j,k,3)=uin(i,j,k,3)/(1.d0+max(Q,0.d0)/rhonp1*dt+dt/t_drag)
              ! Compute momentum & total energy
              uin(i,j,k,2)=rhonp1*uin(i,j,k,2)
              uin(i,j,k,3)=rhonp1*uin(i,j,k,3)
              uin(i,j,k,5) = 5.d-1*(uin(i,j,k,2)**2+uin(i,j,k,3)**2)/uin(i,j,k,1)+half*uin(i,j,k,1)**2/(gamma-1.d0)
           end do
        end do
     end do
  endif

  if (forcing_jet) then
     tauJet=3.e4
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! Add momentum source term
              call getU(x(i),y(j),Ujet)
              dpJet=(uin(i,j,k,1)*Ujet-uin(i,j,k,2))/tauJet*dt
              dEJet=dpJet*uin(i,j,k,2)/uin(i,j,k,1)
              uin(i,j,k,2)=uin(i,j,k,2)+dpJet
              uin(i,j,k,5)=uin(i,j,k,5)+dEJet
              ! Cooling function
              ! Compute local temperature
              rho=uin(i,j,k,1) ; E=uin(i,j,k,5)
              rhovx=uin(i,j,k,2) ; rhovy=uin(i,j,k,3) ; rhovz=uin(i,j,k,4)
              call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,pressure)
              T=pressure/k_mumh/rho 
              ! Compute equilibrium temperature
              call getTeq(Teq,x(i),y(j),z(k),pressure)
              call getTauRad(oneOverTrad,pressure)
              ! Update total energy
              T=(T+oneOverTrad*dt*Teq)/(1.d0+oneOverTrad*dt)
              call get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,gamma*k_mumh*T)
              uin(i,j,k,5)=E
           end do
        end do
     end do
  endif

  if (verbose) write (*,*) 'End of source_term subroutine...'

  return
end subroutine source_term
!===============================================================================
!> Gravity predictor
!===============================================================================
subroutine gravity_predictor(v, igrav, jgrav, kgrav)
  use params
  use variables
  implicit none

  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,3), intent(inout) :: v
  logical, intent(in) :: igrav, jgrav, kgrav
  integer :: ilo, ihi, jlo, jhi, klo, khi
  integer :: i, j, k

  ilo = min(1,iu1+1); ihi = max(1,iu2-1)
  jlo = min(1,ju1+1); jhi = max(1,ju2-1)
  klo = min(1,ku1+1); khi = max(1,ku2-1)

  ! v = v + 1/2*gravin*dt
  !$acc kernels loop
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           if (igrav) v(i,j,k,1) = v(i,j,k,1)! + half*gravin(i,j,k,1)*dt
           if ((jgrav) .and. (ndim > 1)) then
              v(i,j,k,2) = v(i,j,k,2)! + half*gravin(i,j,k,2)*dt
           endif
           if ((kgrav) .and. (ndim > 2)) then
              v(i,j,k,3) = v(i,j,k,3)! + half*gravin(i,j,k,3)*dt
           endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine gravity_predictor
