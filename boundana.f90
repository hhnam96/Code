!===============================================================================
!> \file boundana.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This subroutine contains the analytical boundary conditions. This version
!! contains dummy subroutines only.
!! \details
!! Contains xinner_ana(), xouter_ana(), yinner_ana(), youter_ana(),
!! zinner_ana(), zouter_ana()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          07-05-2015 
!! \b last \b modified: 
!<
!===============================================================================
!> Compute boundary conditions in the x-direction, inner edge
!===============================================================================
subroutine xinner_ana
  use params
  use variables
  use atmosphere
  implicit none
  integer :: i,j,k,ighost
!
  real(dp) :: rho,rhovx,rhovy,rhovz,vx,vy,vz,Bxc,Byc,Bzc,csq,E,p,rm,rp,pm,pp
!
  Bxc=0.d0 ; Byc=0.d0 ; Bzc=0.d0
!
  do k=ku1,ku2
     do j=ju1,ju2

        ! Compute velocities
!        vx=min(uin(iu1+3,j,k,2)/uin(iu1+3,j,k,1),0.d0)
        vx=0.d0!uin(iu1+3,j,k,2)/uin(iu1+3,j,k,1)
        if (wave>0.d0) vx=wave*sin(8.d0*2.d0*asin(1.d0)*time)
        vy=uin(iu1+3,j,k,3)/uin(iu1+3,j,k,1)
        vz=uin(iu1+3,j,k,4)/uin(iu1+3,j,k,1)
        
        ! Get primitive variables
        rho=uin(iu1+3,j,k,1) ; E=uin(iu1+3,j,k,5)
        rhovx=uin(iu1+3,j,k,2) ; rhovy=uin(iu1+3,j,k,3) ; rhovz=uin(iu1+3,j,k,4)
        call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,p)
        
        ! Extrapolate density and pressure
        ighost=1
        do i=iu1+2,iu1,-1
           call getRhoFaces     (rho,p,rm,rp,gamma,ighost*dx)
           call getPressureFaces(rho,p,pm,pp,gamma,ighost*dx) ; csq=gamma*pp/rp
           call get_E(rp,rp*vx,rp*vy,rp*vz,E,Bxc,Byc,Bzc,gamma,csq)
           uin(i,j,k,1) = rp ; uin(i,j,k,5)=E
           ighost=ighost+1
        end do

        ! Zero gradient velocities
        uin(iu1:iu1+2,j,k,2) = uin(iu1:iu1+2,j,k,1)*vx
        uin(iu1:iu1+2,j,k,3) = uin(iu1:iu1+2,j,k,1)*vy
        uin(iu1:iu1+2,j,k,4) = uin(iu1:iu1+2,j,k,1)*vz

     end do
  end do
!

  return
end subroutine xinner_ana
!===============================================================================
!> Compute boundary conditions in the x-direction, outer edge
!===============================================================================
subroutine xouter_ana
  use params
  use variables
  use atmosphere
  implicit none

  integer ::i,k,j,ighost

  real(dp) :: rho,rhovx,rhovy,rhovz,vx,vy,vz,Bxc,Byc,Bzc,csq,E,p,rm,rp,pm,pp
!
  Bxc=0.d0 ; Byc=0.d0 ; Bzc=0.d0
!
  do k=ku1,ku2
     do j=ju1,ju2

        ! Compute velocities
        vx=0.d0!max(uin(iu2-3,j,k,2)/uin(iu2-3,j,k,1),0.d0)
        vy=uin(iu2-3,j,k,3)/uin(iu2-3,j,k,1)
        vz=uin(iu2-3,j,k,4)/uin(iu2-3,j,k,1)
        
        ! Get primitive variables
        rho=uin(iu2-3,j,k,1) ; E=uin(iu2-3,j,k,5)
        rhovx=uin(iu2-3,j,k,2) ; rhovy=uin(iu2-3,j,k,3) ; rhovz=uin(iu2-3,j,k,4)
        call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,p)
        
        ! Extrapolate density and pressure
        ighost=1
        do i=iu2-2,iu2
           call getRhoFaces     (rho,p,rm,rp,gamma,ighost*dx)
           call getPressureFaces(rho,p,pm,pp,gamma,ighost*dx) ; csq=gamma*pm/rm
           call get_E(rm,rm*vx,rm*vy,rm*vz,E,Bxc,Byc,Bzc,gamma,csq)
           uin(i,j,k,1) = rm ; uin(i,j,k,5)=E
           ighost=ighost+1
        end do

        ! Zero gradient velocities
        uin(iu2-2:iu2,j,k,2) = uin(iu2-2:iu2,j,k,1)*vx
        uin(iu2-2:iu2,j,k,3) = uin(iu2-2:iu2,j,k,1)*vy
        uin(iu2-2:iu2,j,k,4) = uin(iu2-2:iu2,j,k,1)*vz

     end do
  end do
!

  return
end subroutine xouter_ana
!===============================================================================
!> Compute boundary conditions in the y-direction, inner edge
!===============================================================================
subroutine yinner_ana
#if NDIM>1
  use params
  use variables
  use atmosphere , only : B0
  implicit none

  integer ::k,i,j
  real(dp) :: rho,rhovx,rhovy,rhovz,P,Bxc,Byc,Bzc,E,Pmag

  ! Periodic BC
  do k=ku1,ku2-1
     do i=iu1,iu2-1

        ! Set zero gradient for density
        rho=uin(i,ju1+3,k,1) ; E=uin(i,ju1+3,k,5)
        uin(i,ju1:ju1+2,k,1) = rho

        ! Set horizontal velocities to zero
        !uin(i,ju1:ju1+2,k,2:3) = 0.d0
        uin(i,ju1:ju1+2,k,2) = uin(i,ju1+3,k,2)
        uin(i,ju1:ju1+2,k,3) = uin(i,ju1+3,k,3)

        ! Zero gradient for vertial velocity
        rhovz=uin(i,ju1+3,k,4)
        uin(i,ju1:ju1+2,k,4) = rhovz

        ! Compute pressure
        rhovx=uin(i,ju1+3,k,2) ; rhovy=uin(i,ju1+3,k,3)
        ! Bxc=5.d-1*(uin(i,ju1+3,k,6)+uin(i,ju1+3,k,9 ))
        ! Byc=5.d-1*(uin(i,ju1+3,k,7)+uin(i,ju1+3,k,10))
        ! Bzc=5.d-1*(uin(i,ju1+3,k,8)+uin(i,ju1+3,k,11))
        Bxc=5.d-1*(uin(i,ju1+3,k,6)+uin(i+1,ju1+3,k  ,6))
        Byc=5.d-1*(uin(i,ju1+3,k,7)+uin(i  ,ju1+4,k  ,7))
        Bzc=5.d-1*(uin(i,ju1+3,k,8)+uin(i  ,ju1+3,k+1,8))
        call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,P)

        ! BC on magnetic field (perp to boundary)
!        uin(i,ju1:ju1+2,k,7  )=B0
        uin(i,ju1:ju1+2,k,7)=uin(i,ju1+3,k,7)
!        uin(i,ju1:ju1+2,k,6)=0.d0 ; uin(i,ju1:ju1+2,k,8)=0.d0
        uin(i,ju1:ju1+2,k,6)=uin(i,ju1+3,k,6)
        uin(i,ju1:ju1+2,k,8)=uin(i,ju1+3,k,8)

        ! Set zero gradient for pressure
        Pmag=     (5.d-1*(uin(i,ju1+3,k,6)+uin(i+1,ju1+3,k  ,6)))**2
        Pmag=Pmag+(5.d-1*(uin(i,ju1+2,k,7)+uin(i  ,ju1+3,k  ,7)))**2
        Pmag=Pmag+(5.d-1*(uin(i,ju1+3,k,8)+uin(i  ,ju1+3,k+1,8)))**2
        Pmag=Pmag/2.d0
        uin(i,ju1:ju1+2,k,5) = 5.d-1*(rhovx**2+rhovy**2+rhovz**2)/rho + P/(gamma-1.d0) + Pmag

!        if (k==1) write (*,*) i,uin(i,ju1+3,k,7),uin(i,ju1+3,k,4)/uin(i,ju1+3,k,1)

     end do
  end do
!
#endif
  return
end subroutine yinner_ana
!===============================================================================
!> Compute boundary conditions in the y-direction, outer edge
!===============================================================================
subroutine youter_ana
#if NDIM>1
  use params
  use variables
  use atmosphere , only : B0
  implicit none

  integer ::k,i
  real(dp) :: rho,rhovx,rhovy,rhovz,P,Bxc,Byc,Bzc,E,Pmag

  do k=ku1,ku2-1
     do i=iu1,iu2-1

        ! Set zero gradient for density
        rho=uin(i,ju2-3,k,1) ; E=uin(i,ju2-3,k,5)
        uin(i,ju2-2:ju2,k,1) = rho

        ! Set horizontal velocities to zero
        !uin(i,ju2-2:ju2,k,2:3) = 0.d0
        uin(i,ju2-2:ju2,k,2) = uin(i,ju2-3,k,2)
        uin(i,ju2-2:ju2,k,3) = uin(i,ju2-3,k,3)

        ! Zero gradient vertical velocity
        rhovz=uin(i,ju2-3,k,4)
        uin(i,ju2-2:ju2,k,4) = rhovz

        ! Compute pressure
        rhovx=uin(i,ju2-3,k,2) ; rhovy=uin(i,ju2-3,k,3)
        ! Bxc=5.d-1*(uin(i,ju2-3,k,6)+uin(i,ju2-3,k,9 ))
        ! Byc=5.d-1*(uin(i,ju2-3,k,7)+uin(i,ju2-3,k,10))
        ! Bzc=5.d-1*(uin(i,ju2-3,k,8)+uin(i,ju2-3,k,11))
        Bxc=5.d-1*(uin(i,ju2-3,k,6)+uin(i+1,ju2-3,k  ,6))
        Byc=5.d-1*(uin(i,ju2-3,k,7)+uin(i  ,ju2-2,k  ,7))
        Bzc=5.d-1*(uin(i,ju2-3,k,8)+uin(i  ,ju2-3,k+1,8))
        call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,P)

        ! BC on magnetic field (perp to boundary)
!        uin(i,ju2-1:ju2,k,7  )=B0
        uin(i,ju2-1:ju2,k,7  )=uin(i,ju2-2,k,7  )
        uin(i,ju2-2:ju2,k,6  )=uin(i,ju2-3,k,6  )
        uin(i,ju2-2:ju2,k,8  )=uin(i,ju2-3,k,8  )
!        uin(i,ju2-2:ju2,k,6)=0.d0 ; uin(i,ju2-2:ju2,k,8)=0.d0

        ! Set zero gradient for pressure
        Pmag=     (5.d-1*(uin(i,ju2-3,k,6)+uin(i+1,ju2-3,k  ,6)))**2
        Pmag=Pmag+(5.d-1*(uin(i,ju2-3,k,7)+uin(i  ,ju2-2,k  ,7)))**2
        Pmag=Pmag+(5.d-1*(uin(i,ju2-3,k,8)+uin(i  ,ju2-3,k+1,8)))**2
        Pmag=Pmag/2.d0
        uin(i,ju2-2:ju2,k,5) = 5.d-1*(rhovx**2+rhovy**2+rhovz**2)/rho + P/(gamma-1.d0) + Pmag

     end do
  end do
!
#endif
  return
end subroutine youter_ana
!===============================================================================
!> Compute boundary conditions in the z-direction, inner edge
!===============================================================================
subroutine zinner_ana
#if NDIM > 2
  use params
  use variables
  use atmosphere , only : B0, zbuffer
  implicit none

  integer ::i,j,k,kghost
  real(dp) :: rho,rhovx,rhovy,rhovz,vx,vy,vz,vzbp1,vzbp2,vzbp3,Bxc,Byc,Bzc,csq,E
  real(dp) :: p,rm,rp,pm,pp,Bfield,dbxdx,dbydy

  ! Zero gradient tangential B-field
  do j=ju1,ju2
     do i=iu1,iu2
        Bfield=uin(i,j,ku1+3,6) ; uin(i,j,ku1:ku1+2,6)=Bfield
!        uin(i,j,ku1:ku1+2,7)=B0
        Bfield=uin(i,j,ku1+3,7) ; uin(i,j,ku1:ku1+2,7)=Bfield
     end do
  end do

  !Normal B-field
  do j=ju1,ju2-1
     do i=iu1,iu2-1
        dbxdx=(uin(i+1,j  ,ku1+2,6)-uin(i,j,ku1+2,6))/dx
        dbydy=(uin(i  ,j+1,ku1+2,7)-uin(i,j,ku1+2,7))/dy
        uin(i,j,ku1+2,8) = uin(i,j,ku1+3,8) + dz*(dbxdx+dbydy)
        dbxdx=(uin(i+1,j  ,ku1+1,6)-uin(i,j,ku1+1,6))/dx
        dbydy=(uin(i  ,j+1,ku1+1,7)-uin(i,j,ku1+1,7))/dy
        uin(i,j,ku1+1,8) = uin(i,j,ku1+2,8) + dz*(dbxdx+dbydy)
        dbxdx=(uin(i+1,j  ,ku1  ,6)-uin(i,j,ku1  ,6))/dx
        dbydy=(uin(i  ,j+1,ku1  ,7)-uin(i,j,ku1  ,7))/dy
        uin(i,j,ku1  ,8) = uin(i,j,ku1+1,8) + dz*(dbxdx+dbydy)
     end do
  end do

  ! Velocity & Pressure BC  
  do j=ju1,ju2-1
     do i=iu1,iu2-1

        ! Set ghost zone velocities
        vx=uin(i,j,ku1+3,2)/uin(i,j,ku1+3,1) !0.d0
        vy=uin(i,j,ku1+3,3)/uin(i,j,ku1+3,1) !0.d0
        vzbp1=-uin(i,j,ku1+3,4)/uin(i,j,ku1+3,1)
        vzbp2=-uin(i,j,ku1+4,4)/uin(i,j,ku1+4,1)
        vzbp3=-uin(i,j,ku1+5,4)/uin(i,j,ku1+5,1)
!        vz=0.d0! +uin(i,j,ku1+3,4)/uin(i,j,ku1+3,1) !beware...changed on 03/02/14
        
        ! Compute pressure in first active cell
        rho=uin(i,j,ku1+3,1) ; E=uin(i,j,ku1+3,5)
        rhovx=uin(i,j,ku1+3,2) ; rhovy=uin(i,j,ku1+3,3) ; rhovz=uin(i,j,ku1+3,4)
        Bxc=5.d-1*(uin(i,j,ku1+3,6)+uin(i+1,j  ,ku1+3,6))
        Byc=5.d-1*(uin(i,j,ku1+3,7)+uin(i  ,j+1,ku1+3,7))
        Bzc=5.d-1*(uin(i,j,ku1+3,8)+uin(i  ,j  ,ku1+4,8))
        call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,p)

        ! Extrapolate density and pressure
        kghost=1
        do k=ku1+2,ku1,-1

           call getRhoFaces     (rho,p,rm,rp,gamma,kghost*dz)
           call getPressureFaces(rho,p,pm,pp,gamma,kghost*dz) ; csq=gamma*pp/rp
           Bxc=5.d-1*(uin(i,j,k,6)+uin(i+1,j  ,k  ,6))
           Byc=5.d-1*(uin(i,j,k,7)+uin(i  ,j+1,k  ,7))
           Bzc=5.d-1*(uin(i,j,k,8)+uin(i  ,j  ,k+1,8))
           call get_E(rp,rp*vx,rp*vy,rp*vzbp1,E,Bxc,Byc,Bzc,gamma,csq)
           uin(i,j,k,1) = rp ; uin(i,j,k,5)=E
           kghost=kghost+1

        !    call getRhoFaces     (rho,p,rm,rp,gamma,dz)
        !    call getPressureFaces(rho,p,pm,pp,gamma,dz) ; csq=gamma*pp/rp

        !    Bxc=5.d-1*(uin(i,j,k,6)+uin(i+1,j  ,k  ,6))
        !    Byc=5.d-1*(uin(i,j,k,7)+uin(i  ,j+1,k  ,7))
        !    Bzc=5.d-1*(uin(i,j,k,8)+uin(i  ,j  ,k+1,8))
        !    call get_E(rp,rp*vx,rp*vy,-uin(i,j,ku1+2+kghost,4),E,Bxc,Byc,Bzc,gamma,csq)           
        !    uin(i,j,k,1) = rp ; uin(i,j,k,5)=E
           
        !    rho=rp
        !    rhovx=rho*vx ; rhovy=rho*vy ; rhovz=-uin(i,j,ku1+2+kghost,4)
        !    call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,p)

            kghost=kghost+1
         end do

        ! Zero gradient horizontal velocities
        uin(i,j,ku1:ku1+2,2) = uin(i,j,ku1:ku1+2,1)*vx
        uin(i,j,ku1:ku1+2,3) = uin(i,j,ku1:ku1+2,1)*vy

!        ! Wall-symmetry in Vertical velocity
!        uin(i,j,ku1+2,4) = uin(i,j,ku1+2,1)*vzbp1
!        uin(i,j,ku1+1,4) = uin(i,j,ku1+1,1)*vzbp2
!        uin(i,j,ku1  ,4) = uin(i,j,ku1  ,1)*vzbp3
        ! No mass flux in
        uin(i,j,ku1+2,4) = -uin(i,j,ku1+3,4)
        uin(i,j,ku1+1,4) = -uin(i,j,ku1+4,4)
        uin(i,j,ku1  ,4) = -uin(i,j,ku1+5,4)

    end do
  end do
!
#endif 
  return
end subroutine zinner_ana
!===============================================================================
!> Compute boundary conditions in the z-direction, inner edge
!===============================================================================
subroutine zouter_ana
#if NDIM > 2
  use params
  use variables
  use atmosphere , only : B0, zbuffer
  implicit none

  integer ::i,j,k,kghost
  real(dp) :: rho,rhovx,rhovy,rhovz,vx,vy,vz,Bxc,Byc,Bzc,csq,E
  real(dp) :: p,rm,rp,pm,pp,Bfield,dbxdx,dbydy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zero gradient horizontal field, vertical field reconstructed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Zero gradient tangential B-field
  do j=ju1,ju2
     do i=iu1,iu2
        Bfield=uin(i,j,ku2-3,6) ; uin(i,j,ku2-2:ku2,6)=Bfield
!        uin(i,j,ku2-2:ku2,7)=B0
!        uin(i,j,ku2-2:ku2,7)=0.d0
        uin(i,j,ku2-2:ku2,7)=uin(i,j,ku2-3,7)
!        B0=uin(i,j,ku2-3,7) ; uin(i,j,ku2-2:ku2,7)=B0
     end do
  end do

  !Normal B-field
  do j=ju1,ju2-1
     do i=iu1,iu2-1
        dbxdx=(uin(i+1,j  ,ku2-2,6)-uin(i,j,ku2-2,6))/dx
        dbydy=(uin(i  ,j+1,ku2-2,7)-uin(i,j,ku2-2,7))/dy
        uin(i,j,ku2-1,8) = uin(i,j,ku2-2,8)-dz*(dbxdx+dbydy)
        dbxdx=(uin(i+1,j  ,ku2-1,6)-uin(i,j,ku2-1,6))/dx
        dbydy=(uin(i  ,j+1,ku2-1,7)-uin(i,j,ku2-1,7))/dy
        uin(i,j,ku2  ,8) = uin(i,j,ku2-1,8)-dz*(dbxdx+dbydy)
     end do
  end do

  ! Velocity & Pressure BC
  do j=ju1,ju2-1
     do i=iu1,iu2-1

        ! Set ghost zone velocities
        vx=uin(i,j,ku2-3,2)/uin(i,j,ku2-3,1) !0.d0
        vy=uin(i,j,ku2-3,3)/uin(i,j,ku2-3,1) !0.d0
        vz=uin(i,j,ku2-3,4)/uin(i,j,ku2-3,1) !-uin(i,j,ku2-3,4)/uin(i,j,ku2-3,1)

        ! Compute pressure in last active cell
        rho=uin(i,j,ku2-3,1) ; E=uin(i,j,ku2-3,5)
        rhovx=uin(i,j,ku2-3,2) ; rhovy=uin(i,j,ku2-3,3) ; rhovz=uin(i,j,ku2-3,4)
        Bxc=5.d-1*(uin(i,j,ku2-3,6)+uin(i+1,j  ,ku2-3,6))
        Byc=5.d-1*(uin(i,j,ku2-3,7)+uin(i  ,j+1,ku2-3,7))
        Bzc=5.d-1*(uin(i,j,ku2-3,8)+uin(i  ,j  ,ku2-2,8))
        call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,p)
        
        ! Extrapolate density and pressure
        kghost=1
        ! do k=ku2-2,ku2-1
        !    call getRhoFaces     (rho,p,rm,rp,gamma,kghost*dz)
        !    call getPressureFaces(rho,p,pm,pp,gamma,kghost*dz) ; csq=gamma*pm/rm
        !    Bxc=5.d-1*(uin(i,j,k,6)+uin(i+1,j  ,k  ,6))
        !    Byc=5.d-1*(uin(i,j,k,7)+uin(i  ,j+1,k  ,7))
        !    Bzc=5.d-1*(uin(i,j,k,8)+uin(i  ,j  ,k+1,8))
        !    call get_E(rm,rm*vx,rm*vy,rm*vz,E,Bxc,Byc,Bzc,gamma,csq)
        !    uin(i,j,k,1) = rm ; uin(i,j,k,5)=E
        !    kghost=kghost+1
        ! end do

        uin(i,j,ku2-2:ku2,1) = uin(i,j,ku2-3,1)
        do k=ku2-2,ku2-1

           if (z(k)<zbuffer) then

              call getRhoFaces     (rho,p,rm,rp,gamma,dz)
              call getPressureFaces(rho,p,pm,pp,gamma,dz) ; csq=gamma*pm/rm
              
              Bxc=5.d-1*(uin(i,j,k,6)+uin(i+1,j  ,k  ,6))
              Byc=5.d-1*(uin(i,j,k,7)+uin(i  ,j+1,k  ,7))
              Bzc=5.d-1*(uin(i,j,k,8)+uin(i  ,j  ,k+1,8))
              
              !call get_E(rm,rm*vx,rm*vy,-uin(i,j,ku2-2-kghost,4),E,Bxc,Byc,Bzc,gamma,csq)
              call get_E(rm,rm*vx,rm*vy,rm*vz,E,Bxc,Byc,Bzc,gamma,csq)
              uin(i,j,k,1) = rm ; uin(i,j,k,5)=E
              
              rho=rm
              !rhovx=rho*vx ; rhovy=rho*vy ; rhovz=-uin(i,j,ku2-2-kghost,4)
              rhovx=rho*vx ; rhovy=rho*vy ; rhovz=rho*vz
              call get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,p)
              kghost=kghost+1
              
           else

              Bxc=5.d-1*(uin(i,j,k,6)+uin(i+1,j  ,k  ,6))
              Byc=5.d-1*(uin(i,j,k,7)+uin(i  ,j+1,k  ,7))
              Bzc=5.d-1*(uin(i,j,k,8)+uin(i  ,j  ,k+1,8))
              
              call get_E(rho,rho*vx,rho*vy,rho*vz,E,Bxc,Byc,Bzc,gamma,gamma*p/rho)
              uin(i,j,k,5)=E

           endif

        end do

        ! Velocities BC
        uin(i,j,ku2-2:ku2,2) = uin(i,j,ku2-2:ku2,1)*vx
        uin(i,j,ku2-2:ku2,3) = uin(i,j,ku2-2:ku2,1)*vy
        uin(i,j,ku2-2:ku2,4) = uin(i,j,ku2-2:ku2,1)*vz

        ! ! Reflective vertical mass flux
        !uin(i,j,ku2-2,4) = -uin(i,j,ku2-3,4)
        !uin(i,j,ku2-1,4) = -uin(i,j,ku2-4,4)
        !uin(i,j,ku2  ,4) = -uin(i,j,ku2-5,4)

    end do
  end do
!
#endif
  return
end subroutine zouter_ana
