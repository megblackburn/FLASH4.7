!!****if* source/physics/utilities/PlasmaState/PlasmaState_logLambda
!!
!! NAME
!!
!!  PlasmaState_logLambda
!!
!! SYNOPSIS
!!
!!  call PlasmaState_logLambda(integer(IN)  :: comp,
!!                                real(IN)  :: tele,
!!                                real(IN)  :: tion,
!!                                real(IN)  :: nele,
!!                                real(IN)  :: mion,
!!                                real(IN)  :: zbar,
!!                                real(OUT) :: ll
!!                    logical(IN), optional :: ionscrn)
!!
!! DESCRIPTION
!!
!! Computes the log of Lambda (Coulomb logarithm) for 
!! electron-electron (ee), electron-ion (ei), or ion-ion (ii) collisions.
!!
!!
!! ARGUMENTS
!!
!!   comp      :   logLambda component (ee, ei, or ii)
!!   tele      :   electron temperature
!!   tion      :   ion temperature
!!   nele      :   electron number density
!!   mion      :   ion mass
!!   zbar      :   Zbar (average ionization state)
!!   ll        :   log of Lambda
!!   ionscrn   :   optional logical to include ion screening in bmax
!!
!!***

#include "constants.h"

subroutine PlasmaState_logLambda(comp, tele, tion, nele, mion, zbar, ll, ionscrn)

  implicit none

  integer, intent(in) :: comp ! logLambda component (ee, ei, or ii)
  real, intent(in)  :: tele   ! electron temperature [K]
  real, intent(in)  :: tion   ! ion temperature [K]
  real, intent(in)  :: nele   ! electron number density [cm^-3]
  real, intent(in)  :: mion   ! ion mass (g)
  real, intent(in)  :: zbar   ! the average ionization state [unitless]
  real, intent(out) :: ll     ! the coulomb logarithm [unitless]
  logical, optional, intent(in) :: ionscrn  ! include ion screening in debye length

  real :: mu     !! reduced mass
  real :: vtherm !! a thermal velocity
  real :: q2     !! charge of particle 1 * charge of particle 2
  real :: bmax, bmin, bmin_classic, bmin_quantum
  real, parameter :: ll_floor = 1.0

  !! Stub

  ll=0.0
  return

end subroutine PlasmaState_logLambda

