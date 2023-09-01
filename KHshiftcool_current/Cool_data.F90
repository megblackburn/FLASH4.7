module Cool_data

!==============================================================================

  implicit none
#include "Flash.h"

  logical, save :: useCool

  double precision, save :: cl_rho, cl_gamma, cl_smalle
  double precision, save :: cl_Boltzmann

  !! for Townsend
  integer, parameter :: cl_N   = 11
  integer, parameter :: cl_Nm1 = 10

  double precision, save, dimension(cl_Nm1) :: alpha_
  double precision, save, dimension(cl_N)   :: T_, Lambda_, Y_

  double precision, save :: cl_temperature_floor

  !! for heating
  double precision, save, dimension(0:7) :: alphaH_
  double precision, save, dimension(0:8) :: Th_, LambdaH_
  double precision, save, dimension(0:7) :: Yh_ 

  logical, save :: cl_useCoolExactInt, cl_useHeatExactInt
  double precision, save :: cl_temperature_ceiling

  double precision, parameter :: cl_amu = 1.66e-24
  double precision, parameter :: cl_mu = 0.6 
  double precision, save :: cl_kb

!==============================================================================

end module Cool_data

