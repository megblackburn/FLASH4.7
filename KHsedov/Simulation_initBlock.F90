!!****if* source/Simulation/SimulationMain/Sedov/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer(IN) :: blockId)
!!                       
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sedov spherical
!!  explosion problem.
!!
!!  References:  Sedov, L. I., 1959, Similarity and Dimensional Methods
!!                 in Mechanics (New York:  Academic)
!!
!!               Landau, L. D., & Lifshitz, E. M., 1987, Fluid Mechanics,
!!                 2d ed. (Oxford:  Pergamon)
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_expEnergy      Explosion energy (distributed over 2^dimen central zones)
!!  sim_minRhoInit     Density floor for initial condition
!!  sim_rInit          Radial position of inner edge of grid (for 1D )
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

!!REORDER(4): solnData


subroutine Simulation_initBlock(blockId)
  use Driver_data, ONLY: dr_globalMe ! NEW
  use Simulation_data, ONLY: sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin, &
     &  sim_nProfile, sim_rProf, sim_vProf, sim_pProf, sim_pExp, sim_rhoProf, &
     &  sim_tInitial, sim_gamma, sim_expEnergy, sim_pAmbient, sim_rhoAmbient, &
     &  sim_useProfileFromFile, sim_profileInitial, &
     &  sim_tempLo, sim_tempUp, sim_rhoLo, sim_rhoUp, & 
     &  sim_smallX, sim_smallRho, sim_minRhoInit, sim_smallP, sim_rInit, &
     &  sim_smallT, &
     &  sim_nSubZones, sim_xCenter, sim_yCenter, sim_zCenter, sim_inSubzones, sim_inszd, &
     sim_threadBlockList, sim_threadWithinBlock
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData, Grid_subcellGeometry
  use ut_interpolationInterface
  
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) ::  blockId
  
  
  integer,parameter :: op = 2
  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk, kat
  real     ::  drProf
  real,allocatable,dimension(:) :: rProf, vProf, rhoProf, pProf
  real     ::  distInv, xDist, yDist, zDist
  real     ::  sumRho, sumP, sumVX, sumVY, sumVZ
  real     ::  vel, diagonal
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  real     ::  vx, vy, vz, p, rho, e, ek, eint
  real     ::  dist
  real     ::  vSub, rhoSub, pSub, errIgnored
  real     :: temp, R
  
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData
  real :: scale_length ! NEW

!!$  real     :: dvSub(0:sim_nSubZones-1,0:(sim_nSubZones-1)*K2D)
  real,allocatable :: dvSub(:,:)
  real     :: dvc, quotinv

  logical :: gcell = .true.

  if (sim_useProfileFromFile) then
     ! lazy initialization - should already have been done from Simulation_init
     if (sim_tinitial > 0.0) call sim_scaleProfile(sim_tinitial)
  end if

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
! get the coordinate information for the current block from the database
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX)); xCoord = 0.0
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY)); yCoord = 0.0
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ)); zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)

  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)

  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  !
  !     For each cell
  !  

  R=8.314

!!!! triple loop calling zz xx and yy but nothing else endo 
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    zz=zCoord(k)
    do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
      yy=yCoord(j)
      do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
        xx = xCoord(i)

        ! Density, Temperature and x-velocity in different regions

        if (yy>0.5) then
          rho=sim_rhoUp
          vx=-0.5
          temp=sim_tempUp
        else
          rho=sim_rhoLo
          vx=0.5
          temp=sim_tempLo
        end if

         ! Velocities
        vy=0.01*SIN(4.0*PI*xx)
        vz=0.0
        p=R*rho*temp
        ! Energies
        ek  = 0.5*(vx*vx + vy*vy)
        e   = p/(sim_gamma-1.0)
        eint= e/rho
        e   = e/rho + ek
        e   = max (e, p)
        call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
        call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
        call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)
        call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eint)
        call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
        call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
        call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
        call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
        call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)

      enddo
    enddo
  enddo


#ifdef FL_NON_PERMANENT_GUARDCELLS
#endif
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

  return
end subroutine Simulation_initBlock




