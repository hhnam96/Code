!###########################################################
!###########################################################
!###########################################################
subroutine history
  use params
  use variables
  use mpi_var
#if MPI==1
  use mpi
#endif
  implicit none

  real(dp)::mass,maxwell,reynolds,dtau,magp,divB
  real(dp),dimension(ndim)::mean_B
  real(dp),dimension(:,:,:,:),allocatable :: loc_var
  real(dp),dimension(:,:    ),allocatable :: loc_var_yz_mean,loc_var_xz_mean
  character(LEN=255) :: filename='history.txt'
  integer ::i,j,k,nvar_mean,ipos,ymidm1,ymidp1

  !------ Zonal mean variables -----------------------------------
  real(dp) :: u_mean

#if MPI==1
  !------ MPI variables ------------------------------------------
  integer :: status(MPI_STATUS_SIZE),ierr,iproc
#endif

  if (verbose) write (*,*) 'Entering history...'

#if NDIM==1
  ipos=1
  do while (x(ipos)<0.5)
     ipos=ipos+1
  end do
  if (mype==0) then
     if (time==0.d0) then
        open(unit=2,file=filename,status='unknown')
     else
        open(unit=2,file=filename,status='old',position='append')
     endif
     write(2,110) time,uin(1,ipos,1,1,1),uin(1,ipos,1,1,2),uin(1,ipos,1,1,5)
     close(2)
  endif
110  format  (1x,1pe12.5,1x,3e14.5)  
#endif

#if NDIM==3

  ! Volume averaged quantities
  mass = 0.d0 ; magp = 0.d0 ; maxwell = 0.d0
  mean_B = 0.d0
  do k=1,nz
      do j=1,ny
         do i=1,nx
            mass = mass + uin(i,j,k,1)
            magp = magp + &
                 & 2.5d-1*(uin(i,j,k,6)+uin(i+1,j  ,k  ,6))**2 + &
                 & 2.5d-1*(uin(i,j,k,7)+uin(i  ,j+1,k  ,7))**2 + &
                 & 2.5d-1*(uin(i,j,k,8)+uin(i  ,j  ,k+1,8))**2
            maxwell = maxwell - &
                 & 5.d-1*(uin(i,j,k,6)+uin(i+1,j,k  ,6)) * &
                 & 5.d-1*(uin(i,j,k,8)+uin(i  ,j,k+1,8))
            mean_B(1) = mean_B(1) + uin(i,j,k,6)
            mean_B(2) = mean_B(2) + uin(i,j,k,7)
            mean_B(3) = mean_B(3) + uin(i,j,k,8)
         end do
      end do
  end do
  dtau = dx*dy*dz/(xmax-xmin)/(ymax-ymin)/(zmax-zmin)
  magp    = magp*dtau/2.
  mass = mass*dtau
  maxwell = maxwell*dtau
  mean_B  = mean_B*dtau

  reynolds = 0.d0

  nvar_mean=3
  allocate(loc_var        (iu1:iu2,ju1:ju2,ku1:ku2,nvar_mean))
  allocate(loc_var_yz_mean(iu1:iu2,                nvar_mean))
  loc_var(:,:,:,1)=uin(:,:,:,1)
  loc_var(:,:,:,2)=uin(:,:,:,2)/uin(i,j,k,1)
  loc_var(:,:,:,3)=uin(:,:,:,4)/uin(i,j,k,1)
  call compute_yz_mean(loc_var,loc_var_yz_mean,nvar_mean)
  deallocate(loc_var)

  do k=1,nz
      do j=1,ny
         do i=1,nx
            reynolds=reynolds+uin(i,j,k,1)*dtau*&
                    & (uin(i,j,k,2)/uin(i,j,k,1)-loc_var_yz_mean(i,2))*&
                    & (uin(i,j,k,4)/uin(i,j,k,1)-loc_var_yz_mean(i,3))
         end do
      end do
  end do
  deallocate(loc_var_yz_mean)

  divB=0.d0
  do k=1,nz
      do j=1,ny
         do i=1,nx
            divB=divB + &
              (uin(i+1,j  ,k  ,6)-uin(i,j,k,6))/dx + &
              (uin(i  ,j+1,k  ,7)-uin(i,j,k,7))/dy + &
              (uin(i  ,j  ,k+1,8)-uin(i,j,k,8))/dz
         end do
      end do
  end do

#if MPI==1
  call vol_average(mass,1)
  call vol_average(reynolds,1)
  call vol_average(maxwell,1)
  call vol_average(magp,1)
  call vol_average(divB,1)
  call vol_average(mean_B,3)
#endif      

  ! Zonally averaged quantities at equator
  allocate(loc_var        (iu1:iu2,ju1:ju2,ku1:ku2,3))
  allocate(loc_var_xz_mean(        ju1:ju2        ,3))
  loc_var(:,:,:,1)=uin(:,:,:,2)/uin(:,:,:,1)
  loc_var(:,:,:,2)=uin(:,:,:,3)/uin(:,:,:,1)
  loc_var(:,:,:,3)=uin(:,:,:,4)/uin(:,:,:,1)
  call compute_xz_mean(loc_var,loc_var_xz_mean,3)
#if MPI==1
  if (mod(nyslice,2)==0) then
     ymidm1=ju1+2 ; ymidp1=ymidm1+1
     u_mean=5.d-1*(loc_var_xz_mean(ymidm1,1)+loc_var_xz_mean(ymidp1,1))
  else
     ymidm1=(ju2-ju1-5)/2. ; ymidp1=ymidm1+1
     u_mean=5.d-1*(loc_var_xz_mean(ymidm1,1)+loc_var_xz_mean(ymidp1,1))
  endif
  if (nyslice>1) then
     iproc=(nxslice-1)+nxslice*nzslice*(nyslice/2)
     if (mype==iproc) call MPI_Send(u_mean,1,MPI_DOUBLE_PRECISION,0    ,10,MPI_COMM_WORLD, ierr)
     if (mype==0    ) call MPI_Recv(u_mean,1,MPI_DOUBLE_PRECISION,iproc,10,MPI_COMM_WORLD, status, ierr)
  endif
#else
  ymidm1=(ju2-ju1-5)/2. ; ymidp1=ymidm1+1
  u_mean=5.d-1*(loc_var_xz_mean(ymidm1,1)+loc_var_xz_mean(ymidp1,1))
#endif
  deallocate(loc_var,loc_var_xz_mean)

  if (mype==0) then
     if (time==0.d0) then
        open(unit=2,file=filename,status='unknown')
     else
        open(unit=2,file=filename,status='old',position='append')
     endif
     write(2,120) time,dt,mass,maxwell,reynolds,maxwell+reynolds,magp,&
                & mean_B(1),mean_B(2),mean_B(3),divB,u_mean
     close(2)
  endif
120  format  (1x,1pe12.5,1x,1pe12.5,10e14.5)
#endif

end subroutine history
!===============================================================================
!> This routine computes volume average
!===============================================================================
subroutine vol_average(quant, nquant)
#if MPI == 1
  use params
  use mpi
  implicit none

  integer, intent(in) :: nquant
  real(dp), dimension(nquant), intent(inout) :: quant
  real(dp), dimension(nquant) :: vol_quant
  integer :: ierr

  call MPI_Reduce(quant, vol_quant, nquant, MPI_DOUBLE_PRECISION, MPI_SUM, 0 &
       , MPI_COMM_WORLD, ierr)
  quant = vol_quant

#endif
  return
end subroutine vol_average
!===============================================================================
!> This routine computes mean in the yz-plane
!===============================================================================
subroutine compute_yz_mean(quant, mean_quant, nquant)
#if MPI == 1
  use mpi
  use mpi_var
#endif
  use params

  integer, intent(in) :: nquant
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nquant), intent(in) :: quant
  real(dp), dimension(iu1:iu2,nquant), intent(out) :: mean_quant
  real(dp), dimension(:,:), allocatable :: mean_quant_local
  real(dp) :: mean_quant_loc
  integer :: ivar, i, j, k
#if MPI == 1
  integer :: ierr, commx
#endif

  allocate(mean_quant_local(iu1:iu2, nquant))

  do ivar = 1, nquant
     do i = iu1, iu2
        mean_quant_loc = zero
        !$OMP PARALLEL DO REDUCTION(+: mean_quant_loc) SCHEDULE(RUNTIME)
        do k = 1, nz
           do j = 1, ny
              mean_quant_loc = mean_quant_loc + quant(i,j,k,ivar)
           enddo
        enddo
        !$OMP END PARALLEL DO
        mean_quant_local(i,ivar) = mean_quant_loc
     enddo
  enddo

#if MPI == 1
  call create_comm_plane(commx, 1)
  mean_quant = zero
  call MPI_Reduce(mean_quant_local, mean_quant, (iu2-iu1+1)*nquant &
       , MPI_DOUBLE_PRECISION, MPI_SUM, 0, commx, ierr)
  mean_quant = mean_quant/ny/nz/nyslice/nzslice
  call MPI_Bcast(mean_quant, (iu2-iu1+1)*nquant, MPI_DOUBLE_PRECISION &
       , 0, commx, ierr)
  call MPI_Comm_Free(commx, ierr)
#else
  mean_quant = mean_quant_local/ny/nz
#endif
  
  deallocate(mean_quant_local)

  return
end subroutine compute_yz_mean
!===============================================================================
!> This routine computes mean in the xz-plane
!===============================================================================
subroutine compute_xz_mean(quant, mean_quant, nquant)
#if MPI == 1
  use mpi
  use mpi_var
#endif
  use params

  integer, intent(in) :: nquant
  real(dp), dimension(iu1:iu2,ju1:ju2,ku1:ku2,nquant), intent(in) :: quant
  real(dp), dimension(ju1:ju2,nquant), intent(out) :: mean_quant
  real(dp), dimension(:,:), allocatable :: mean_quant_local
  real(dp) :: mean_quant_loc
  integer :: ivar, i, j, k
#if MPI == 1
  integer :: ierr, commy
#endif

  allocate(mean_quant_local(ju1:ju2,nquant))

  do ivar = 1, nquant
     do j = ju1, ju2
        mean_quant_loc = zero
        !$OMP PARALLEL DO REDUCTION(+: mean_quant_loc) SCHEDULE(RUNTIME)
        do k = 1, nz
           do i = 1, nx
              mean_quant_loc = mean_quant_loc + quant(i,j,k,ivar)
           enddo
        enddo
        !$OMP END PARALLEL DO
        mean_quant_local(j,ivar) = mean_quant_loc
     enddo
  enddo

#if MPI == 1
  call create_comm_plane(commy, 2)
  mean_quant = zero
  call MPI_Reduce(mean_quant_local, mean_quant, (ju2-ju1+1)*nquant &
       , MPI_DOUBLE_PRECISION, MPI_SUM, 0, commy, ierr)
  mean_quant = mean_quant/nx/nz/nxslice/nzslice
  call MPI_Bcast(mean_quant, (ju2-ju1+1)*nquant, MPI_DOUBLE_PRECISION &
       , 0, commy, ierr)
  call MPI_Comm_Free(commy, ierr)
#else
  mean_quant = mean_quant_local/nx/nz
#endif

  deallocate(mean_quant_local)

  return
end subroutine compute_xz_mean
