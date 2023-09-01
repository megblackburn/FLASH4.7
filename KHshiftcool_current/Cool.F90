subroutine Cool (blockCount, blockList, dt, time)

  use Cool_data,       ONLY : useCool, cl_rho, cl_smalle, cl_useCoolExactInt, &
                              cl_useHeatExactInt


  use Cool_interface,  ONLY : Cool_exact_integrate, Heat_exact_integrate
  use Grid_interface,  ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits, &
                              Grid_releaseBlkPtr, Grid_getListOfBlocks
  use Eos_interface,   ONLY : Eos_wrapped

#include "constants.h"
#include "Flash.h"

  implicit none
 
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real, intent(IN) :: dt, time
 
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: blockID
  real :: eid, eidold, ei, ek, dedt
  integer :: i, j, k, lb


!====================================================================================================


  if(.not.useCool) return
 
  do lb=1, blockCount
     blockID=blockList(lb)
     call Grid_getBlkPtr(blockID,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     do k=blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j=blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          
           cl_rho = solnData(DENS_VAR,i,j,k)
           ei = solnData(EINT_VAR,i,j,k)
           
           eid = ei * cl_rho
           eidold = ei * cl_rho
           ek = 0.5 * (solnData(VELX_VAR,i,j,k)**2 + &
                      solnData(VELY_VAR,i,j,k)**2 + &
                      solnData(VELZ_VAR,i,j,k)**2)
           if (cl_useCoolExactInt) then
             call Cool_exact_integrate(eid,dt,blockID,i,j,k)
           endif
           if (cl_useHeatExactInt) then
             call Heat_exact_integrate(eid,dt)
           endif
#ifdef COOL_VAR
           solnData(COOL_VAR,i,j,k) = (eid - eidold)
#endif

           ei = eid/cl_rho
           if (ei < cl_smalle) then
             write(*,*) "WARNING applying internal energy floor"
             ei = max(ei, cl_smalle)
           endif
           solnData(ENER_VAR,i,j,k)= ei + ek
           solnData(EINT_VAR,i,j,k) = ei

#ifdef ENUC_VAR
           solnData(ENUC_VAR,i,j,k) = (eid - eidold) / dt
#endif
           enddo
        enddo 
     enddo  
     
     call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
     call Grid_releaseBlkPtr(blockID,solnData)
  end do
  return


end subroutine Cool

