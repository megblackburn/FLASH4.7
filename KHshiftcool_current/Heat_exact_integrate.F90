subroutine Heat_exact_integrate (eint, dt) 

  use Cool_data,        ONLY : cl_kb, cl_mu, cl_amu, cl_rho, cl_gamma, &
                               cl_smalle, alphaH_, Th_, LambdaH_, Yh_, &
                               cl_temperature_floor, cl_temperature_ceiling

  use Cool_interface,   ONLY : Cool_compute_eint
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"

  !! arguments
  double precision, intent(INOUT) :: eint
  double precision, intent(IN)    :: dt

  !! parameters
  integer, parameter :: N = size(Th_) - 1 !! Th_ is zero-indexed
  integer, parameter :: Nm1 = N - 1

  !! arrays
  integer :: temp_loc(1), k ! array index

  !! Physical variables
  double precision :: Temperature
  double precision :: tcool

  !! Townsend variables
  double precision :: Y, Lambda ! Lambda(T)
  double precision :: Y_invarg, Yinv
  double precision :: Lambda_ref, Lambda_k, T_ref, T_k, Y_k, alpha_k
  double precision :: temperature_ceiling


!********************** GET T **************************************
  Temperature = (cl_gamma - 1d0) * eint * cl_mu * cl_amu / (cl_kb * cl_rho)

  if     (Temperature < cl_temperature_floor  ) then
    Temperature = cl_temperature_floor
  elseif (Temperature > cl_temperature_ceiling) then
    Temperature = cl_temperature_ceiling
  endif

  !write(*,*) "HEI init T,e: ", Temperature, eint

!********************** COMPUTE k ********************************
  !! This is the k of Y_k, alpha_k, etc.

  temp_loc = minloc( abs(Th_ - Temperature) )
  k = temp_loc(1) - 1 ! minus 1 to account for zero-indexing
                      ! (minloc assumes 1-indexed)

  if (Temperature < Th_(k)) then
    k = k - 1 ! above looks for closest whereas we want closest below
  endif

  if ((k < 0) .or. (k >= N)) then
    write(*,*) "ERROR! Bad k which is: ", k
    call Driver_abortFlash("line 86 Heat_exact_integrate, aborting...")
  end if

  T_k      = Th_(k)
  alpha_k  = alphaH_(k)
  Lambda_k = LambdaH_(k)


  if (Lambda_k == 0.0) then ! no heating
    return
  endif

!********************** COMPUTE LAMBDA ********************************
  !! Lambda(T) = Lambda_(k) * (T/T_(k))**alpha_(k)    T_k <= T <= T_k+1

  Lambda = Lambda_k * ( Temperature / T_k )**alpha_k

  !write(*,*) "HEI Lambda: ", Lambda
  !write(*,*) "HEI Lambda_k: ", Lambda_k
  !write(*,*) "Temperature: ", Temperature

!********************** COMPUTE TCOOL *********************************  
  tcool = 1.5d0*cl_kb * Temperature * cl_mu * cl_amu / (cl_rho * Lambda)

!********************** COMPUTE Y *************************************
  Y_k = Yh_(k)
  T_ref = Th_(0)
  Lambda_ref = LambdaH_(0)

  temperature_ceiling = cl_temperature_ceiling
  do while (k < Nm1)
    k = k + 1
    if (LambdaH_(k) == 0.0) then
      temperature_ceiling = Th_(k)
    endif
  enddo

  Y = Y_k - 1d0/(1d0-alpha_k) * (Lambda_ref/Lambda_k) * &
      (T_k / T_ref) * ( 1d0 - ( T_k / Temperature )**(alpha_k-1d0) )
!********************** COMPUTE Y_invarg ******************************
  Y_invarg = Y + (Temperature/T_ref)*(Lambda_ref/Lambda)*(dt/tcool)

!********************** COMPUTE Y^{-1} ********************************
  !! T^{n+1} = Y^{-1}

  Y = Y_invarg

  temp_loc = minloc( abs(Yh_ - Y) )
  k = temp_loc(1) - 1 ! minus 1 to account for zero-indexing
                      ! (minloc assumes 1-indexed)
  if (Y < Yh_(k)) then
    k = k - 1
    if (k == N) then
      !! temperature would rise above ceiling, just set it to ceiling
      Temperature = temperature_ceiling
      call Cool_compute_eint(Temperature, eint)
      return
    endif
  endif
  Y_k = Yh_(k)
  T_k = Th_(k)
  alpha_k  = alphaH_(k)
  Lambda_k = LambdaH_(k)

  Yinv = T_k * (1d0 + (1d0 - alpha_k) * (Lambda_k/Lambda_ref) * &
        (T_ref/T_k)*(Y - Y_k))**(1d0/(1d0-alpha_k))

  Temperature = Yinv

  if (Temperature < cl_temperature_floor) then
    Temperature = cl_temperature_floor
  elseif (Temperature > temperature_ceiling) then
    Temperature = temperature_ceiling
  endif

  call Cool_compute_eint(Temperature, eint)

  if     (eint >= 0.0) then
    CONTINUE
  elseif (eint < 0) then
    write(*,*) "ERROR! Negative internal energy: ", eint
    call Driver_abortFlash("line 187 Heat_exact_integrate, aborting...")
  else ! catches NaNs
    write(*,*) "ERROR! Bad internal energy: ", eint
    call Driver_abortFlash("line 190 Heat_exact_integrate, aborting...")
  end if

  !write(*,*) "HEI final T,e: ", Temperature, eint

  return
end subroutine Heat_exact_integrate

