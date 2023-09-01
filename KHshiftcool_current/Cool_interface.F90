! from RF
Module Cool_interface
#include "constants.h"
#include "Flash.h"
  
  interface 
     subroutine Cool_computeDt (block_no, &
          blkLimits,blkLimitsGC,  &
          solnData,   &   
          dt_check, dt_minloc )
    
       integer, intent(IN) :: block_no
       integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
       real,INTENT(INOUT)    :: dt_check
       integer,INTENT(INOUT)    :: dt_minloc(5)
       real, pointer, dimension(:,:,:,:) :: solnData
     end subroutine Cool_computeDt
  end interface

  interface
     subroutine Cool_compute_eint(Temperature, eint)
       implicit none
       double precision, intent( IN) :: Temperature
       double precision, intent(OUT) :: eint
     end subroutine Cool_compute_eint
  end interface

  interface 
     subroutine Cool_exact_integrate(e, dt, blockID, i,j,k)
       real, intent(INOUT)  :: e
       real, intent(IN)     :: dt
       integer, intent(IN)  :: blockID
       integer, intent(IN)  :: i, j, k
     end subroutine Cool_exact_integrate
  end interface

  interface
    subroutine Heat_exact_integrate(eint, dt) 
      double precision, intent(INOUT) :: eint
      double precision, intent(IN)    :: dt
    end subroutine Heat_exact_integrate
  end interface

  interface 
     subroutine Cool(blockCount,blockList,dt, time)
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount), intent(IN) :: blockList
       real,intent(IN) :: dt, time
     end subroutine Cool
  end interface

  interface
     subroutine Cool_finalize()
     end subroutine Cool_finalize
  end interface

  interface
     subroutine Cool_init()

     end subroutine Cool_init
  end interface

  interface
     subroutine Cool_unitTest( fileUnit, perfect)
         integer, intent(IN) ::  fileUnit
         logical, intent(INOUT) :: perfect
     end subroutine Cool_unitTest
  end interface

end Module Cool_interface

