! computes the specific internal energy in erg/g - from RF

subroutine Cool_compute_eint(Temperature, eint)
  use Cool_data,       ONLY: cl_kb, cl_mu, cl_amu, cl_rho, cl_gamma
  implicit none
  double precision, intent(IN) :: Temperature
  double precision, intent(OUT) :: eint
  eint = cl_kb * Temperature * cl_rho / ((cl_gamma-1.0)*cl_mu*cl_amu)

end subroutine Cool_compute_eint
