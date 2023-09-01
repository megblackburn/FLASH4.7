subroutine Driver_shift()
 !! use Driver_data, ONLY : dr_globalNumProcs, dr_globalComm, dr_globalMe
  use Simulation_data, ONLY : sim_kb_mu_mp, sim_shift_start, sim_tempHot, sim_tempCold, &
                            & sim_rhoHot, sim_rhoCold, sim_smallP, sim_zMin, sim_gamma, &
                            & sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr,   &
                             Grid_releaseBlkPtr, Grid_getCellCoords, &
                             Grid_getBlkIndexLimits
 ! use Timers_interface, ONLY :Timers_start, Timers_stop
  use Driver_data, ONLY : dr_dt, dr_simTime, dr_globalComm
  implicit none


#include "Flash_mpi.h"
!!#include "mpif.h"
#include "Flash.h"
#include "constants.h"
 !! integer,intent(IN) ::  blockId
  
  integer :: bb, blkID, blkCount
  integer :: blkList(MAXBLOCKS)
  integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM)
  
  real :: local_cold_mass, global_cold_mass
  real :: local_v2_sum, global_v2_sum
  real :: vy_avg
  !!real :: sizeX, sizeY, sizeZ
  real :: xx,yy,zz
  integer :: i,j,k, ierr

  real :: rho, p, temp
  real :: cell_dimx, cell_dimy, cell_dimz
  real :: cell_vol, global_vol

  real :: front_posn_new
  real :: front_posn_old
  real :: v_shift_t_new
  real :: front_velocity
  !!real :: sim_dt
  !!real :: g
  real :: c_s
  real :: c_s_cap
  real :: e, eint, ek
  real :: time

!!  real :: v3_av
  

  !!double precision :: velx_avg,vely_avg,velz_avg, ekinZone, &
    !!                  v_shiftX, v_shiftY, v_shiftZ
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, pointer, dimension(:,:,:,:) :: solnData
  logical :: gcell = .true.
!===============================================================================================

  local_cold_mass = 0.0
  global_cold_mass = 0.0
  local_v2_sum = 0.0
  global_v2_sum = 0.0
  
  front_posn_old = sim_yMin
  c_s = sqrt(sim_gamma*(sim_tempHot)/sim_kb_mu_mp)
  c_s_cap = 100.0
  call Grid_getListOfBlocks(LEAF, blkList, blkCount)

  do bb = 1, blkCount
    blkID = blkList(bb)
    call Grid_getBlkPtr(blkID, solnData)
    call Grid_getBlkIndexLimits(blkID, blkLimits, blkLimitsGC)

    sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
    allocate(xCoord(sizeX)); xCoord = 0.0
    sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
    allocate(yCoord(sizeY)); yCoord = 0.0
    sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
    allocate(zCoord(sizeZ)); zCoord = 0.0

    if (NDIM == 3) call Grid_getCellCoords&
                        (KAXIS, blkId, CENTER, gcell, zCoord, sizeZ)

    if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blkId, CENTER,gcell, yCoord, sizeY)

    call Grid_getCellCoords(IAXIS, blkId, CENTER, gcell, xCoord, sizeX)

    global_vol = ((sim_xMax-sim_xMin)*(sim_xMax-sim_xMin))+((sim_yMax-sim_yMin)*(sim_xMax-sim_xMin))+((sim_zMax-sim_Zmin)*(sim_xMax-sim_xMin))

  !! cl_kB = boltzmann
    do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
      zz=zCoord(k)
      do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        yy=yCoord(j)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
          xx = xCoord(i)
     !!     call Grid_getBlkPtr(blkID, solnData)
          rho = solnData(DENS_VAR,i,j,k)
          p = solnData(PRES_VAR,i,j,k)
          temp = (p/rho)*sim_kb_mu_mp
  
          if (i>1) then
            if (j>1) then
              if (k>1) then 
                cell_dimx = xx-xCoord(i-1)
                cell_dimy = yy-yCoord(j-1)
                cell_dimz = zz-zCoord(k-1)
                cell_vol = cell_dimx*cell_dimy*cell_dimz
         !! global_vol = sizeX*sizeY*sizeZ
                if (temp <= sim_tempHot) then
                  local_cold_mass = rho*cell_vol
                endif
                local_v2_sum = local_v2_sum +solnData(VELY_VAR,i,j,k)
              endif
            endif  
          endif
        enddo 
      enddo 
    enddo
  
    call MPI_ALLREDUCE(local_cold_mass, global_cold_mass, 1, FLASH_REAL, MPI_SUM, dr_globalComm, ierr)
    call MPI_ALLREDUCE(local_v2_sum, global_v2_sum, 1, FLASH_REAL, MPI_SUM, dr_globalComm, ierr)
  
    vy_avg = global_v2_sum/global_vol
    front_posn_new = sim_yMin + global_cold_mass/(sizeX*sizeY*sim_rhoHot*sim_tempHot/sim_tempCold)
    v_shift_t_new = -2.0 * (front_posn_new-front_posn_old)/dr_dt; !! use Driver_data, ONLY : dr_dt

    if (abs(v_shift_t_new) > (c_s_cap*c_s)) then
      v_shift_t_new = v_shift_t_new / abs(v_shift_t_new);
      v_shift_t_new = v_shift_t_new * c_s * c_s_cap
    endif

    front_velocity = -1.0 * v_shift_t_new;

    if (dr_simTime <= sim_shift_start) then !!use Driver_data, ONLY : dr_dt, dr_simTime
      front_posn_old = front_posn_new
    else
      do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        zz=zCoord(k)
        do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
          yy=yCoord(j)
          do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
            xx = xCoord(i)
          !!dE_kin = ((solnData(VELZ_VAR,i,j,k)+v_shift_t_new)*(solnData(VELZ_VAR,i,j,k)+v_shift_t_new))-(solnData(VELZ_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k))
            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) - v_shift_t_new
            ek  = 0.5*((solnData(VELX_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)) + (solnData(VELY_VAR,i,j,k) &
                     * solnData(VELY_VAR,i,j,k)) + (solnData(VELZ_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)))
            e   = p/(sim_gamma-1.0)
            eint= e/rho
            e   = e/rho + ek
            solnData(ENER_VAR,i,j,k) = max (e, sim_smallP)
            solnData(EINT_VAR,i,j,k) = eint
          enddo 
        enddo
      enddo
    endif
    call Grid_releaseBlkPtr(blkID, solnData)
    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)
  enddo



end subroutine Driver_shift














