!!****if* source/Simulation/SimulationMain/KHshiftcool_current
!! uses a fourth-order Runge-Kutta method to integrate the Sutherland-Dopita cooling function
!! and compute the new internal energy. This was handled by a third-order method provided 
!! by SVODE in Flash 2
!! from RF


subroutine cool_integrate (e,t,dt)
  implicit none
  real, intent(INOUT) :: e
  real, intent(IN) :: t
  real, intent(IN) :: dt

  return 

end subroutine cool_integrate
