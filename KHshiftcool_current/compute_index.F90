subroutine compute_index(Temperature, kk)

  use Cool_data,    ONLY : T_
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  double precision, intent(IN) :: Temperature
  integer, intent(OUT) :: kk
  integer :: temp_loc(1)
  integer, parameter :: N = size(T_) - 1
  
  temp_loc = minloc(abs(T_-Temperature))
  kk = temp_loc(1)
  if (Temperature < T_(kk)) then 
    kk = kk - 1
  endif 

  if ((kk < 0) .or. (kk >= N)) then 
    write(*,*) "ERROR! k: ", kk
    write(*,*) "Temperature: ", Temperature
    write(*,*) "T_(1)", T_(1)
    call Driver_abortFlash("line 25 compute_index, aborting...")
  endif

  return

end subroutine compute_index


  
