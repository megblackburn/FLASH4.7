! Computes timestep limiter for cooling source term solver - from RF

subroutine Cool_computeDt (blockID, blkLimits, blkLimitsGC, solnData, &
                           dt_check, dt_minloc)

#include "Flash.h"
#include "constants.h"

  implicit none
  integer, intent(IN) :: blockID
  integer, intent(IN), dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, intent(INOUT) :: dt_check
  integer, intent(INOUT) :: dt_minloc(5)
  real, pointer, dimension(:,:,:,:) :: solnData

  return
end subroutine Cool_computeDt
