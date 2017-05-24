!===============================================================================
!> \file user_restart.f90
!! \brief
!! \b DUMSES-Hybrid:
!! This is user's restart subroutines.
!! \details
!! Contains user_restart()
!! \author
!! Marc Joos <marc.joos@cea.fr>, Sébastien Fromang, Romain Teyssier, 
!! Patrick Hennebelle
!! \copyright
!! Copyrights 2013-2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          06-22-2015
!! \b last \b modified: 06-22-2015
!<
!===============================================================================
!> User customization at restart
!===============================================================================
subroutine user_restart
  use precision
  use params
  use variables
  implicit none

  integer :: k, j, i
  real(dp) :: B0  
  logical :: addmag
  namelist /restart_params/addmag,B0

  !Add B-field
  addmag=.false.
  open(unit=1,file='input' ,status='old')
  read(1,restart_params)
  close(1)

  if (addmag) then
     uin(:,:,:,7     )=B0
     uin(:,:,:,nvar+2)=B0
     uin(:,:,:,5)=uin(:,:,:,5)+B0**2/2.d0
  endif

  return
end subroutine user_restart
