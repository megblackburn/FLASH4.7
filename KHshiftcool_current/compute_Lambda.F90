subroutine compute_Lambda(Temperature, kk, Lambda)

  use Cool_data,       ONLY : Lambda_, T_, alpha_

  implicit none
  double precision, intent( IN) :: Temperature
  integer, intent(OUT) :: kk
  double precision, intent(OUT) :: Lambda

  call compute_index(Temperature,kk)

  Lambda = Lambda_(kk) * ( Temperature / T_(kk) )**alpha_(kk)

  return
end subroutine compute_Lambda

