! uses - from RF
subroutine Cool_exact_integrate (e, dt, blockID, ii,jj,kk)

  use Cool_data, ONLY : cl_Boltzmann, cl_rho, cl_gamma, alpha_, T_, &
                        Lambda_, Y_, cl_temperature_floor, cl_N, cl_Nm1

  use Simulation_data, ONLY : sim_amu, sim_wmu, sim_tempWind

  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"

  !! arguments
  real, intent(INOUT) :: e
  real, intent(IN)    :: dt
  integer, intent(IN) :: blockID, ii,jj,kk

  !! arrays
  integer :: k ! array index

  !! Physical variables
  double precision :: T
  double precision :: tcool

  !! Townsend variables
  double precision :: Lambda ! Lambda(T)
  double precision :: Y, Y_invarg, Yinv

  !! convenience variables 
  double precision :: e_dp, cl_rho_dp, KB_dp, gamma_dp, mu_dp, mp_dp
  double precision :: alpha_k, T_k, T_N, Lambda_k, Lambda_N, Y_k

  !! Setup convenience variables
  !!!! double precision
  e_dp         = DBLE(e)
  KB_dp        = DBLE(cl_Boltzmann)
  gamma_dp     = DBLE(cl_gamma)
  cl_rho_dp    = DBLE(cl_rho)
  mu_dp        = DBLE(sim_wmu)
  mp_dp        = DBLE(sim_amu)

  T_N = T_(cl_N)

  Lambda_N = Lambda_(cl_Nm1) * (T_N / T_(cl_Nm1))**alpha_(cl_Nm1)
  Lambda_(cl_N) = Lambda_N

  Y_(cl_N) = 0.0
  do k=cl_Nm1,1,-1
    Y_(k) = Y_(k+1) - 1.0/(1.0-alpha_(k)) * (Lambda_N/Lambda_(k)) * &
         (T_(k) / T_N) * ( 1.0 - ( T_(k) / T_(k+1) )**(alpha_(k)-1.0) )
  end do

  ! Temperature 

  T = (gamma_dp - 1.0) * e_dp * mu_dp * mp_dp / (KB_dp * cl_rho_dp)
  ! Changed 1d0 to 1.0

  ! k 
  !! for Y_k, alpha_k, etc.
  if ( T < cl_temperature_floor ) then
    T = cl_temperature_floor
    e = REAL( KB_dp*T*cl_rho_dp / ((gamma_dp-1d0)*mu_dp*mp_dp) )
    return
  else if (T > 0.6*sim_tempWind) then
    return !! no cooling of hot ambient medium
  else
    if      (T >= T_(1) .AND. T < T_(2)) then
      k = 1
    else if (T >= T_(2) .AND. T < T_(3)) then
      k = 2
    else if (T >= T_(3) .AND. T < T_(4)) then
      k = 3
    else if (T >= T_(4) .AND. T < T_(5)) then
      k = 4
    else if (T >= T_(5) .AND. T < T_(6)) then
      k = 5
    else if (T >= T_(6) .AND. T < T_(7)) then
      k = 6
    else if (T >= T_(7)) then
      k = 7
    else
      write(*,*) "ERROR! Bad T which is: ", T
      write(*,*) "for blockID, ii,jj,kk: ", blockID, ii,jj,kk
      write(*,*) "cl_gamma: ", cl_gamma
      write(*,*) "gamma_dp: ", gamma_dp
      write(*,*) "gamma_dp-1d0: ", gamma_dp - 1d0
      write(*,*) "e: ", e
      write(*,*) "e_dp: ", e_dp
      write(*,*) "k: ", k
      write(*,*) "KB_dp: ", KB_dp
      write(*,*) "cl_rho_dp: ", cl_rho_dp

      call Driver_abortFlash("line 124 Cool_exact_integrate")
    end if
  end if

  Y_k = Y_(k)
  T_k = T_(k)
  alpha_k = alpha_(k)
  Lambda_k = Lambda_(k)

!********************** COMPUTE LAMBDA ********************************

  Lambda = Lambda_k * ( T / T_k )**alpha_k


!********************** COMPUTE Y *************************************
  Y = Y_k + 1d0/(1d0-alpha_k) * (Lambda_N/Lambda_k) * &
      (T_k / T_N) * ( 1d0 - ( T_k / T )**(alpha_k-1d0) )

!********************** COMPUTE TCOOL *********************************  
  tcool = 1.5d0*KB_dp*T * mu_dp * mp_dp / (cl_rho_dp * Lambda)

!********************** COMPUTE Y_invarg ******************************
  Y_invarg = Y + (T/T_N)*(Lambda_N/Lambda)*(DBLE(dt)/tcool)

!********************** COMPUTE Y^{-1} ********************************

  Y = Y_invarg

  if      (Y >= Y_(2) .AND. Y < Y_(1)) then
    k = 1
  else if (Y >= Y_(3) .AND. Y < Y_(2)) then
    k = 2
  else if (Y >= Y_(4) .AND. Y < Y_(3)) then
    k = 3
  else if (Y >= Y_(5) .AND. Y < Y_(4)) then
    k = 4
  else if (Y >= Y_(6) .AND. Y < Y_(5)) then
    k = 5
  else if (Y >= Y_(7) .AND. Y < Y_(6)) then
    k = 6
  else if (Y >= Y_(8) .AND. Y < Y_(7)) then
    k = 7
  else if (Y >= Y_(1)) then
    !! temperature would drop below floor, just set it to the floor.
    T = cl_temperature_floor
    e = REAL( KB_dp*T*cl_rho_dp / ((gamma_dp-1d0)*mu_dp*mp_dp) )
    return
  else
    write(*,*) "ERROR! Y out of range! Y is: ", Y
    write(*,*) "for blockID, ii,jj,kk: ", blockID, ii,jj,kk
    write(*,*) "Y_: ", Y_
    write(*,*) "Temp regime, T_k: ", T_k
    write(*,*) "k based on T: ", k
    write(*,*) "initial temperature T: ", T
    write(*,*) "Crashing now..."
    call Driver_abortFlash("line 180 Cool_exact_integrate")
  end if
  Y_k = Y_(k)
  T_k = T_(k)
  alpha_k = alpha_(k)
  Lambda_k = Lambda_(k)

  Yinv = T_k * (1d0 - (1d0 - alpha_k) * (Lambda_k/Lambda_N) * &
        (T_N/T_k)*(Y - Y_k))**(1d0/(1d0-alpha_k))

  T = Yinv

  if (T < cl_temperature_floor) then
    T = cl_temperature_floor
  endif

  e = REAL( KB_dp*T*cl_rho_dp / ((gamma_dp-1d0)*mu_dp*mp_dp) )

  if (e < 0) then
    write(*,*) "ERROR! Negative e: ", e
    write(*,*) "T: ", T
    write(*,*) "Yinv: ", Yinv
    write(*,*) "T_k: ", T_k
    write(*,*) "alpha_k: ", alpha_k
    write(*,*) "Lambda_k: ", Lambda_k
    write(*,*) "Lambda_N: ", Lambda_N
    write(*,*) "T_N: ", T_N
    write(*,*) "T_k: ", T_k
    write(*,*) "Y: ", Y
    write(*,*) "Y_k: ", Y_k

    call Driver_abortFlash("line 211 Cool_exact_integrate")
  else if (e >= 0.0) then
    CONTINUE
  else
    write(*,*) "ERROR! Bad e: ", e
    write(*,*) "T: ", T
    write(*,*) "Yinv: ", Yinv
    write(*,*) "T_k: ", T_k
    write(*,*) "alpha_k: ", alpha_k
    write(*,*) "Lambda_k: ", Lambda_k
    write(*,*) "Lambda_N: ", Lambda_N
    write(*,*) "T_N: ", T_N
    write(*,*) "T_k: ", T_k
    write(*,*) "Y: ", Y
    write(*,*) "Y_k: ", Y_

    call Driver_abortFlash("line 230 Cool_exact_integrate")
  end if


  return
end subroutine Cool_exact_integrate

