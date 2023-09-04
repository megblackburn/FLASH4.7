subroutine Driver_s()
 !! use Driver_data, ONLY : dr_globalNumProcs, dr_globalComm, dr_globalMe
  use Simulation_data, ONLY : sim_kb_mu_mp, sim_shift_start, sim_tempHot, sim_tempCold, &
                            & sim_rhoHot, sim_rhoCold, sim_smallP, sim_zMin, sim_gamma, &
                            & sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr,   &
                             Grid_releaseBlkPtr, Grid_getCellCoords, &
                             Grid_getDeltas, &
                             Grid_getBlkIndexLimits
  use Grid_data, ONLY : gr_eosModeInit
 ! use Timers_interface, ONLY :Timers_start, Timers_stop
  use Driver_data, ONLY : dr_dt, dr_simTime, dr_globalComm, dr_globalMe
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Eos_interface, ONLY : Eos_wrapped
  use Driver_interface, ONLY : Driver_abortFlash


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
  integer :: i,j,k, ierr !! , ierr

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
  real :: Lx, Ly, Lz
  real :: dx, dy, dz
  real, dimension(MDIM) :: del
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
  v_shift_t_new = 0.0
  front_posn_old = sim_yMin
  !front_posn_new = sim_yMin
  c_s = sqrt(sim_gamma*(sim_tempCold)/sim_kb_mu_mp)
  c_s_cap = 100.0
  if (dr_globalMe == MASTER_PE) then
    write(*,*) "Driver_s initiated"
  endif
!!  call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  if (dr_simTime > sim_shift_start) then
    call Grid_getListOfBlocks(LEAF, blkList, blkCount)
    if (dr_globalMe == MASTER_PE) then
      write(*,*) "Driver_s in shift time"
    endif
    do bb = 1, blkCount
      if (dr_globalMe == MASTER_PE) then
        write(*,*) "Driver_s bb loop"
      endif
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
    
      Lx = sim_xMax - sim_xMin
      Ly = sim_yMax - sim_yMin
      Lz = sim_zMax - sim_zMin

      call Grid_getDeltas(blkID, del)
      dx = del(1)
      dy = del(2)
      dz = del(3)

    !global_vol = ((sim_xMax-sim_xMin)*(sim_xMax-sim_xMin))+((sim_yMax-sim_yMin)*(sim_xMax-sim_xMin))+((sim_zMax-sim_Zmin)*(sim_xMax-sim_xMin))
      global_vol = Lx * Ly

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
            temp = (p/rho)/sim_kb_mu_mp
           !! if (dr_globalMe == MASTER_PE) then
              !!write(*,*) "rho", rho
              !!write(*,*) "pressure", p
              !!write(*,*) " temperature", temp
            !!endif  
            if (i>1 .and. j>1) then
      !      if (j>1) then
              !if (k>=1) then 
            !cell_dimx = xx-xCoord(i-1) 
            !cell_dimy = yy-yCoord(j-1)
            
            !if (dr_globalMe == MASTER_PE) then
             ! write(*,*) "cell_dimx", cell_dimx
              !write(*,*) "cell_dimy", cell_dimy
              !write(*,*) "dz", dz
              !write(*,*) "dy", dy
              !write(*,*) "dz", dz
            !endif

                !cell_dimz = zz-zCoord(k-1)
              cell_vol = dx * dy
            !cell_vol = cell_dimx*cell_dimy
              if (temp <= sim_tempHot) then
                local_cold_mass = rho*cell_vol 
              endif
              if (temp > 1.e10 .or. temp < 1.0) then
                call Driver_abortFlash("invalid temp, check units/algebra")
              endif

              local_v2_sum = local_v2_sum +solnData(VELY_VAR,i,j,k)
              call MPI_ALLREDUCE(local_cold_mass, global_cold_mass, 1, FLASH_REAL, MPI_SUM, dr_globalComm, ierr) !!ierr
              call MPI_ALLREDUCE(local_v2_sum, global_v2_sum, 1, FLASH_REAL, MPI_SUM, dr_globalComm, ierr)  !!ierr
              vy_avg = global_v2_sum/global_vol
              front_posn_new = global_cold_mass/(Lx*Ly*sim_rhoHot*sim_tempHot/sim_tempCold)
              v_shift_t_new = (-2.0 * (front_posn_new-front_posn_old)/dr_dt) !! use Driver_data, ONLY : dr_dt

              if (abs(v_shift_t_new) > (c_s_cap*c_s)) then
                v_shift_t_new = v_shift_t_new / abs(v_shift_t_new)
                v_shift_t_new = v_shift_t_new * c_s * c_s_cap
              endif
              front_velocity = -1.0 * v_shift_t_new
!              endif
            else
              front_posn_new = front_posn_old
              v_shift_t_new = (-2.0 * (front_posn_new-front_posn_old)/dr_dt)
           ! endif  
            endif

            !!if (dr_globalMe == MASTER_PE) then
              !!write(*,*) "front_posn_new", front_posn_new
            !!endif
  
            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) - v_shift_t_new
            ek  = 0.5*((solnData(VELX_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)) + (solnData(VELY_VAR,i,j,k) &
                     * solnData(VELY_VAR,i,j,k)) + (solnData(VELZ_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)))
            e   = p/(sim_gamma-1.0)
            eint= e/rho
            e   = e/rho + ek
            solnData(ENER_VAR,i,j,k) = max (e, sim_smallP)
            solnData(EINT_VAR,i,j,k) = eint

            if (dr_simTime >= sim_shift_start) then !!use Driver_data, ONLY : dr_dt, dr_simTime
              front_posn_old = front_posn_new
            endif

          enddo 
        enddo 
      enddo
    
    !if (dr_globalMe == MASTER_PE) then
   !   write(*,*) "v_shift_t_new", v_shift_t_new
    !  write(*,*) "energy", energy
     ! write(*,*) "energy_int", energy_int
      !write(*,*) "vely", vely
 !     write(*,*) "front_posn_old", front_posn_old
  !    write(*,*) "front_posn_new", front_posn_new
   !   write(*,*) "temp", temp
    !  write(*,*) "rho", rho
     ! write(*,*) "pressure", p
      !write(*,*) "global_cold_mass", global_cold_mass
      !write(*,*) "global_vol", global_vol
      !write(*,*) "cell_vol", cell_vol
    !endif

      !!if (dr_globalMe == MASTER_PE) then 
        !!write(*,*) "v_shift_t_new", v_shift_t_new
        !! write(*,*) "front_posn_old", front_posn_old    
   !   write(*,*) "front_posn_new", front_posn_new
     !! endif 
      call Eos_wrapped(gr_eosModeInit, blkLimits, blkList(bb))  
 
      call Grid_releaseBlkPtr(blkID, solnData)
      deallocate(xCoord)
      deallocate(yCoord)
      deallocate(zCoord)
!    call Eos_wrapped(gr_eosModeInit, blkLimits, blkList(bb))  
    enddo
    if (dr_globalMe == MASTER_PE) then
      write(*,*) "Driver_s finished"
    endif
  endif


end subroutine Driver_s














