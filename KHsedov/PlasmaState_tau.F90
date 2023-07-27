!!****if* source/physics/utilities/PlasmaState/PlasmaState_tau
!!
!! NAME
!!
!!  PlasmaState_tau
!!
!! SYNOPSIS
!!
!!  call PlasmaState_tau(integer(IN)  :: comp,
!!                          real(IN)  :: temp,
!!                          real(IN)  :: ll,
!!                          real(IN)  :: Z,
!!                          real(IN)  :: n,
!!                          real(OUT) :: tau,
!!                 real(IN), OPTIONAL :: mion)
!!
!! DESCRIPTION
!!
!! Computes collision time (s) for 
!! electron-electron (ee), electron-ion (ei), or ion-ion (ii) collisions.
!! tau here is equivalent to Braginskii
!!
!!
!! ARGUMENTS
!!
!!   comp      :   component (ee, ei, or ii)
!!   temp      :   temperature
!!   ll        :   log of Lambda
!!   Z         :   atomic number
!!   n         :   number density
!!   tau       :   collision time
!!   mion      :   ion mass
!!
!!***

#include "constants.h"

subroutine PlasmaState_tau(comp, temp, ll, z, n, tau, mion)

  implicit none

  integer, intent(in)        :: comp   ! logLambda component (ee, ei, or ii)
  real, intent(in)           :: temp   ! temperature [K]
  real, intent(in)           :: ll     ! the coulomb logarithm [unitless]
  real, intent(in)           :: Z      ! atomic number
  real, intent(in)           :: n      ! number density [cm^-3]
  real, intent(out)          :: tau    ! collision time [s]
  real, intent(in),optional  :: mion   ! ion mass [g]


  !! Stub
  tau =0.0 
  return

end subroutine PlasmaState_tau

