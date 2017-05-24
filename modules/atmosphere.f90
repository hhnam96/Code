module atmosphere
  use precision
  implicit none

  real(dp) :: g,Mp,k_mumh,wave,DT_EP,beta,omegaEq,rho0,P0,zbuffer,eta_buffer,eta_body
  logical :: balanced,isentropic,cooling,drag,deep,sponge
  !
  !Newtonian cooling variables
  real(dp) :: Psurf,Tsurf,Dtrop,dTstra,zstra,t_rad,t_drag,lyTherm
  logical :: coolSimple

  !Shallow water variables
  logical :: forcing
  real(dp) :: Q0

  ! Initial magnetic field (B=B0 for P>Pcrit initially)
  real(dp) :: B0,Pcrit 

  ! Jet properties (for jet test problem)
  real(dp) :: U0,U1
  logical :: forcing_jet

end module atmosphere
