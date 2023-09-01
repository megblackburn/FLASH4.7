subroutine compute_tcool(density, Temperature, Lambda, tcool)

  use Cool_data,   ONLY : cl_Boltzmann
  use Simulation_data,  ONLY : sim_wmu, sim_amu
  implicit none
  double precision, intent(IN) :: density, Temperature, Lambda
  double precision, intent(OUT) :: tcool

  tcool = 1.5d0*cl_Boltzmann * Temperature * sim_wmu * sim_amu &
          / (density * Lambda)

end subroutine compute_tcool
