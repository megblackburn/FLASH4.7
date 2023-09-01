!from RF
subroutine Cool_init(myPe)

  use Simulation_data, ONLY : sim_cloudTemp, sim_coreTemp, sim_tempWind
  use Cool_data
  use Driver_data, ONLY : dr_globalMe

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

#include "constants.h"

  implicit none

  integer,intent(IN) :: MyPE

  !! Townsend
  double precision :: kelvinPerKeV

  !! Heating
  integer, parameter :: Nh = size(Th_) - 1 !! T_ is zero-indexed
  integer, parameter :: Nm1h = Nh - 1
  double precision, dimension(Nm1h) :: cH_
  double precision :: tempWind,densContrast
  integer :: k
  integer :: ii

  double precision :: coreCloudMixTemp, cloudWindMixTemp

!==============================================================================

  call PhysicalConstants_get("Boltzmann", cl_Boltzmann)
  call PhysicalConstants_get("Boltzmann", cl_kb)
  call RuntimeParameters_get("smalle", cl_smalle)
  call RuntimeParameters_get("useCool",useCool)
  call RuntimeParameters_get('sim_cloudTemp',      cl_temperature_floor)
  call RuntimeParameters_get("cl_useCoolExactInt", cl_useCoolExactInt)
  call RuntimeParameters_get("cl_useHeatExactInt", cl_useHeatExactInt)
  call RuntimeParameters_get('sim_densContrast',   densContrast)

  tempWind = cl_temperature_floor * densContrast
  cl_temperature_ceiling = tempWind

  !cl_temperature_floor = cl_temperature_floor / densContrast ! Tcore

  !! 1.602e-12 erg/eV so 1.602e-9 erg/keV
  kelvinPerKeV = 1.602e-9 / cl_Boltzmann 
  cl_gamma = 5.0/3.0

  if (.not. useCool .and. dr_globalMe == MASTER_PE) then
     write(6,*)'WARNING:  You have included the Cool unit but have set '
     write(6,*)'   the runtime parameter useCool to FALSE'
     write(6,*)'   No cooling will occur but Cool_init will continue.'
  end if


!********************** TOWNSEND VARIABLES ****************************
  !! fit from Mike McCourt (sent to me by Max Gronke) to the Sutherland
  !! & Doipita (1993), Z~1 cooling curve (used in Gronke & Oh 2018).

  ! in K:  300,       2e3,     8e3,     1e4,          ~2e4, ~22e4
  T_ = (/ 2.584e-5, 1.723e-4, 6.89e-4, 8.614e-4, 0.0017235, 0.02, &
          0.13, 0.7, 5.0, 100.0, 1e5 /)

  !! that is in keV, need to convert to Kelvin
  T_ = T_ * kelvinPerKeV
  !! Lambda_(N) gets filled in later
  Lambda_ = (/ 1e-3, 1.5e-3, 3e-3, 2.4e-1, 15.438249, &
               66.831473, 2.773501, 1.195229, 1.842056, 6.10541, 0.01 /) 
  !! last one was 0.0
  Lambda_ = Lambda_ * 1e-23

  !alpha_ = (/ 2.0, 1.5, 2.867, 6.0, 0.6, -1.7, -0.5, 0.22, 0.4, 0.4 /)
  !! this guarantees a continuous cooling curve
  do ii = 1, cl_Nm1
    alpha_(ii) = log10(Lambda_(ii+1)/Lambda_(ii)) / log10(T_(ii+1) / T_(ii))
  enddo
  alpha_(cl_Nm1) = alpha_(cl_Nm1-1)

  !! To not have to change things too much I leave the above as
  !! the names for cooling and add here the intervals for heating
  coreCloudMixTemp = sqrt( sim_coreTemp*sim_cloudTemp)
  cloudWindMixTemp = sqrt(sim_cloudTemp*sim_tempWind )

  Th_     = (/  250.0, 300.0, 400.0,   1e3,      4e3,   1e4,      4e4, 8e4, 1e19 /)

  LambdaH_ = (/ 1e-24, 1e-25, 1.03e-26, 4e-27, 2.12e-26, 1e-21, 2.34e-22, 0.001,   0.001 /) 
  !! last 2 were 0.0

  do ii = 0, Nm1h
    alphaH_(ii) = log10(LambdaH_(ii+1)/LambdaH_(ii)) / log10(TH_(ii+1) / TH_(ii))
  enddo
  alphaH_(Nm1h-1) = -10
  alphaH_(Nm1h) = alphaH_(Nm1h-1)

  !alphaH_ = (/     -12.629, -9.61,   -3,   5.06,   8.416,   -1.0477,  -10, -10 /)



  !! Heating
  Yh_(0) = 0.0
  do k=1,Nm1h
    Yh_(k) = Yh_(k-1) - 1.0/(1.0-alphaH_(k-1)) * (LambdaH_(0)/LambdaH_(k-1)) * &
       (Th_(k-1) / Th_(0)) * ( 1.0 - ( Th_(k-1) / Th_(k) )**(alphaH_(k-1)-1.0) )
  end do


!********************** END TOWNSEND VARIABLES ************************

   if (dr_globalMe .EQ. MASTER_PE) print *, "init_cool:  successfully read file."

  if (dr_globalME == MASTER_PE) then
    write(*,*) "T_: ", T_
    write(*,*) "TH_: ", TH_
    write(*,*) "LambdaH_: ", LambdaH_
    write(*,*) "Lambda_: ", Lambda_
    write(*,*) "alpha_: ", alpha_
    write(*,*) "alphaH_: ", alphaH_
    write(*,*) "YH_: ", YH_
  endif

  return
end subroutine Cool_init

