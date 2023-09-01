!!****if*source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_computeDtExtMHD
!!
!! NAME
!!
!!  hy_uhd_computeDtExtMHD
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_computeDtExtMHD(real,pointer   :: U(:,:,:,:),
!!                         interger(IN)   :: i,
!!                         interger(IN)   :: j,
!!                         interger(IN)   :: k,
!!                         real(IN)       :: idx,
!!                         real(IN)       :: idy,
!!                         real(IN)       :: idz,
!!                         real(IN)       :: xCenter,
!!                         real(IN)       :: uxgrid,
!!                         real(IN)       :: uygrid,
!!                         real(IN)       :: uzgrid,
!!                         real(IN)       :: sndspd2,
!!                         real(IN)       :: b2,
!!                         real(INOUT)    :: cfx2,
!!                         real(INOUT)    :: cfy2,
!!                         real(INOUT)    :: cfz2,
!!                         real(OUT)      :: dt)
!!  
!!
!! DESCRIPTION
!!
!!  This routine computes the timestep limiter for the Unsplit Hydro solver.
!!  The Courant-Fredrichs-Lewy criterion is used.  The sound
!!  speed is computed and together with the velocities, is used to constrain
!!  the timestep such that no information can propagate more than one zone
!!  per timestep.
!!  Note that this routine only accounts for computing advection time step in hyperbolic
!!  system of equations.
!!
!!  This part of computeDt handles extended MHD terms (e.g., Hall, Biermann, etc.)
!!
!! ARGUMENTS
!!
!!  U             -  the physical, solution data from grid
!!  i,j,k         -  indices of current cell
!!  idx,idy,idz   -  inverse distances in each {*=x, y z} directions
!!  xCenter       -  r-position of current cell (used for CYLINDICAL)
!!  uxgrid        -  velocity of grid expansion in x directions
!!  uygrid        -  velocity of grid expansion in y directions
!!  uzgrid        -  velocity of grid expansion in z directions
!!  sndspd2       -  sound speed^2
!!  b2            -  Alfven speed^2
!!  cfx2          -  starts as fast Alven wave speed^2 in x-direction, becomes generally fastest wave speed^2
!!  cfy2          -  starts as fast Alven wave speed^2 in y-direction, becomes generally fastest wave speed^2
!!  cfz2          -  starts as fast Alven wave speed^2 in z-direction, becomes generally fastest wave speed^2
!!  dt       -  temporary timestep constraint for this cell
!!
!!***

Subroutine hy_uhd_computeDtExtMHD( U,       &
                                   i, j, k, &
                                   idx, idy, idz, &
                                   xCenter, & 
                                   uxgrid, uygrid, uzgrid, &
                                   sndspd2, b2, &
                                   cfx2, cfy2, cfz2) !,  &
                                   !dt)


#include "Flash.h"
#include "constants.h"

#ifdef FLASH_MULTISPECIES
#include "Multispecies.h"
#endif

  use Hydro_data, ONLY: hy_meshMe, hy_useHydro, hy_geometry, &
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       hy_useHall, hy_useBiermann3T, hy_useMagneticResistivity, &
       hy_useCrossMagRes, hy_useAnisoMagRes, &
       hy_useNernst, hy_useCrossField, &
       hy_nernstFlMode, hy_nernstFlCoef, &
       hy_crossFieldFlMode, hy_crossFieldFlCoef, &
#endif
       hy_avogadro, hy_qele, hy_mele, hy_boltz
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY: Eos_getAbarZbar
  use MagneticResistivity_interface, ONLY : MagneticResistivity
  use Thermoelectric_interface, ONLY : Thermoelectric
  use PlasmaState_interface, ONLY : PlasmaState_logLambda, PlasmaState_tau
#ifdef FLASH_MULTISPECIES
  use Multispecies_interface, ONLY: Multispecies_getSumFrac
#else
  use Eos_data, ONLY : eos_singleSpeciesZ
#endif
  implicit none

  !! Arguments type declaration ------------------------------------------
  real, pointer       :: U(:,:,:,:)
  integer, intent(IN) :: i, j, k 
  real,    intent(IN) :: idx, idy, idz, xCenter
  real,    intent(IN) :: uxgrid, uygrid, uzgrid
  real,    intent(IN) :: sndspd2, b2
  real,    intent(INOUT) :: cfx2, cfy2, cfz2
 ! real,    intent(OUT)   :: dt

  !! ----------------------------------------------------------------------

  integer :: tempVar
  real    :: i_inerLen, abar, zbar, dxNe, dyNe, dzNe, ye, nele_right, nele_left, nele, Lne
  real    :: nion, mion, zfullFrac, z2fullFrac, zfull, z2full, zeff, llee, llei
  real    :: whistlerSpeed, HallDriftSpeed, ThermalMagSpeed
  real    :: etaparUnused, etapar, etaperp, etacross, tau, tauCent, tele_right, tele_left, tele
  real    :: dxEtaTau, dyEtaTau, dzEtaTau, etatau_right, etatau_left, etatau, Letatau
  real    :: crossWhistlerSpeed, crossDriftSpeed
  real    :: betapar, betaperp, betacross, betanew, dxTe, dyTe, dzTe, chi, Bmag
  real    :: vEff_X, vEff_Y, vEff_Z   !! effective electron advection velocity
  real    :: vMax_X, vMax_Y, vMax_Z   !! maximum advection velocity
  real    :: vMag
  real    :: jx, jy, jz, jr, vHall_X, vHall_Y, vHall_Z
  real    :: vResist_X, vResist_Y, vResist_Z
  real    :: vNernst_X, vNernst_Y, vNernst_Z
  real    :: vCross_X, vCross_Y, vCross_Z
  real    :: BcrossgradTe_X, BcrossgradTe_Y, BcrossgradTe_Z, crossCoeff
  real    :: Bcrossj_X, Bcrossj_Y, Bcrossj_Z, resistCoeff
  logical :: needCoeffs, useBetaNew
  integer :: gfact

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
#ifdef TELE_VAR
  tempVar = TELE_VAR
#else
  tempVar = TEMP_VAR
#endif

  !! Define gfact (minus sign used for some cross products)
  if (hy_geometry == CARTESIAN) then
     gfact = 1.0
  elseif (hy_geometry == CYLINDRICAL) then
     gfact = -1.0
  end if


  !! Initialize effective electron velocity
  vEff_X = U(VELX_VAR,i,j,k)
  if (NDIM > 1) vEff_Y = U(VELY_VAR,i,j,k)
  if (NDIM > 2) vEff_Z = U(VELZ_VAR,i,j,k)


  !!-----------------------------------------------------------------------
  !!    Compute various quantities needed for extended MHD terms
  !!-----------------------------------------------------------------------
  !! Compute nele (used by all terms)
  call Eos_getAbarZbar(U(:,i,j,k),abar=abar,zbar=zbar, Ye=ye) 
  nele =  ye * hy_avogadro * U(DENS_VAR,i,j,k)         
  nion = nele/zbar

  !! Compute collision time tau and any needed coefficients
  etatau = 0.
  needCoeffs = (hy_useCrossMagRes .or. hy_useNernst .or. hy_useCrossField .or. hy_useAnisoMagRes)
  if (needCoeffs) then
     tele = U(tempVar,i,j,k)
     mion = abar / hy_avogadro
     call PlasmaState_logLambda(ELE_ELE, tele, tele, nele, mion, zbar, llee)
     call PlasmaState_logLambda(ELE_ION, tele, tele, nele, mion, zbar, llei)

     !! "effective" Z as defined by DaviesWen
#ifdef FLASH_MULTISPECIES
     !! These are "full" values (i.e., average atomic number, not ionization state)
     call Multispecies_getSumFrac(Z,zfullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j,k))
     zfull = abar*zfullFrac
     call Multispecies_getSumFrac(Z2,z2fullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j,k))
     z2full = abar*z2fullFrac
#else
     zfull = eos_singleSpeciesZ 
     z2full = zfull**2
#endif
     zeff = (z2full*llei) / (zfull*llee)

#ifdef LLEI_VAR
     U(LLEI_VAR,i,j,k) = llei
#endif
#ifdef LLEE_VAR
     U(LLEE_VAR,i,j,k) = llee
#endif
#ifdef ZEFF_VAR
     U(ZEFF_VAR,i,j,k) = zeff
#endif

     call PlasmaState_tau(ELE_ION, tele, llei, zbar, nion, tau)
     Bmag = sqrt(dot_product(U(MAGX_VAR:MAGZ_VAR,i,j,k),U(MAGX_VAR:MAGZ_VAR,i,j,k)))
     chi = tauCent*hy_qele*Bmag/hy_mele

#ifdef CHIE_VAR
     U(CHIE_VAR,i,j,k) = chi
#endif

     !! resistivity coefficients
     if (hy_useAnisoMagRes .and. hy_useCrossMagRes) then
        call MagneticResistivity(U(:,i,j,k), etapar, etaperp, etacross)
     elseif (hy_useAnisoMagRes .and. (.not. hy_useCrossMagRes)) then
        call MagneticResistivity(U(:,i,j,k), etapar, etaperp)
     elseif ((.not. hy_useAnisoMagRes) .and. hy_useCrossMagRes) then
        call MagneticResistivity(U(:,i,j,k), etaparUnused, resCross=etacross)
        etatau = etacross*tauCent
     end if

     !! thermoelectric coefficients
     if (hy_useNernst .and. hy_useCrossField) then
        call Thermoelectric(U(:,i,j,k), betapar, betaperp, betacross, betanew)
     elseif (hy_useNernst .and. (.not. hy_useCrossField)) then
        call Thermoelectric(U(:,i,j,k), betaCross=betacross)
     elseif ((.not. hy_useNernst) .and. hy_useCrossField) then
        call Thermoelectric(U(:,i,j,k), betapar, betaperp, betaNew=betanew)
     end if

     !! specific to cross-field
     if (hy_useCrossField) then
        if (zeff >= 1.) then
           useBetaNew = .true.
        else
           useBetaNew = .false.
        end if
     end if
  end if


  !! Compute gradients used for Hall and Biermann
  if (hy_useHall .or. hy_useBiermann3T) then
     i_inerLen = sqrt(U(DENS_VAR,i,j,k)) * zbar * hy_qele * hy_avogadro / abar
     i_inerLen = 1.0/i_inerLen 

     !! Compute gradient of electron density for Hall drift waves and thermal magnetic waves
     !! For Hall term, drift speed ~ 1/(qele*nele**2)*grad(nele), and we compute
     !! Lne = 1/nele*grad(ne) to use in a dispersion relation. The other factor of
     !! 1/(qele*nele) is "hidden" inside i_inerLen*sqrt(b2).
     call Eos_getAbarZbar(U(:,i+1,j,k),abar=abar,zbar=zbar, Ye=ye)
     nele_right = ye * hy_avogadro * U(DENS_VAR,i+1,j,k) 

     !! Also compute gradient of eta*tau needed for cross resistivity drift speed
     !! Here cross drift speed ~ -grad(etacross*chi/B) = -grad(etacross*tau)*qele/(mele*clight)
     !! Multiply by (etacross*tau)/(etacross*tau) and use Letatau = 1/(etacross*tau)*grad(etacross*tau)
     !! Now cross drift speed ~ -etacross*tau*qele/(mele*clight)*Letatau, and we again use
     !! Letatau in a dispersion relation. We need to pull out a factor of 1/(qele*nele) in
     !! order to use i_inerLen*sqrt(b2) properly, thus cross drift speed = 
     !! -etacross*tau*qele**2*nele/(mele*clight)*Letatau*i_inerLen*sqrt(b2)
     if (hy_useCrossMagRes) then
        tele_right = U(tempVar,i+1,j,k)
        call PlasmaState_logLambda(ELE_ELE, tele_right, tele_right, &
                                   nele_right, abar/hy_avogadro, zbar, llee)
        call PlasmaState_logLambda(ELE_ION, tele_right, tele_right, &
                                   nele_right, abar/hy_avogadro, zbar, llei)
#ifdef FLASH_MULTISPECIES
        call Multispecies_getSumFrac(Z,zfullFrac,U(SPECIES_BEGIN:SPECIES_END,i+1,j,k))
        zfull = abar*zfullFrac
        call Multispecies_getSumFrac(Z2,z2fullFrac,U(SPECIES_BEGIN:SPECIES_END,i+1,j,k))
        z2full = abar*z2fullFrac
#else
        zfull = eos_singleSpeciesZ
        z2full = zfull**2
#endif
        zeff = (z2full*llei) / (zfull*llee)
        call PlasmaState_tau(ELE_ION, tele_right, llei, zbar, nele_right/zbar, tau)
        call MagneticResistivity(U(:,i+1,j,k), etaparUnused, resCross=etacross)
        etatau_right = etacross*tau
     end if

     call Eos_getAbarZbar(U(:,i-1,j,k),abar=abar,zbar=zbar, Ye=ye)
     nele_left = ye * hy_avogadro * U(DENS_VAR,i-1,j,k) 

     dxNe = 0.5*(nele_right - nele_left)*idx
     dyNe = 0.0
     dzNe = 0.0

     if (hy_useCrossMagRes) then
        tele_left = U(tempVar,i-1,j,k)
        call PlasmaState_logLambda(ELE_ELE, tele_left, tele_left, &
                                   nele_left, abar/hy_avogadro, zbar, llee)
        call PlasmaState_logLambda(ELE_ION, tele_left, tele_left, &
                                   nele_left, abar/hy_avogadro, zbar, llei)
#ifdef FLASH_MULTISPECIES
        call Multispecies_getSumFrac(Z,zfullFrac,U(SPECIES_BEGIN:SPECIES_END,i-1,j,k))
        zfull = abar*zfullFrac
        call Multispecies_getSumFrac(Z2,z2fullFrac,U(SPECIES_BEGIN:SPECIES_END,i-1,j,k))
        z2full = abar*z2fullFrac
#else
        zfull = eos_singleSpeciesZ 
        z2full = zfull**2
#endif
        zeff = (z2full*llei) / (zfull*llee)
        call PlasmaState_tau(ELE_ION, tele_left, llei, zbar, nele_left/zbar, tau)
        call MagneticResistivity(U(:,i-1,j,k), etaparUnused, resCross=etacross)
        etatau_left = etacross*tau

        dxEtaTau = 0.5*(etatau_right - etatau_left)*idx
        dyEtaTau = 0.0
        dzEtaTau = 0.0
     end if

     if (NDIM > 1) then
        call Eos_getAbarZbar(U(:,i,j+1,k),abar=abar,zbar=zbar, Ye=ye)
        nele_right = ye * hy_avogadro * U(DENS_VAR,i,j+1,k) 

        if (hy_useCrossMagRes) then
           tele_right = U(tempVar,i,j+1,k)
           call PlasmaState_logLambda(ELE_ELE, tele_right, tele_right, &
                                      nele_right, abar/hy_avogadro, zbar, llee)
           call PlasmaState_logLambda(ELE_ION, tele_right, tele_right, &
                                      nele_right, abar/hy_avogadro, zbar, llei)
#ifdef FLASH_MULTISPECIES
           call Multispecies_getSumFrac(Z,zfullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j+1,k))
           zfull = abar*zfullFrac
           call Multispecies_getSumFrac(Z2,z2fullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j+1,k))
           z2full = abar*z2fullFrac
#else
           zfull = eos_singleSpeciesZ
           z2full = zfull**2
#endif
           zeff = (z2full*llei) / (zfull*llee)
           call PlasmaState_tau(ELE_ION, tele_right, llei, zbar, nele_right/zbar, tau)
           call MagneticResistivity(U(:,i,j+1,k), etaparUnused, resCross=etacross)
           etatau_right = etacross*tau
        end if

        call Eos_getAbarZbar(U(:,i,j-1,k),abar=abar,zbar=zbar, Ye=ye)
        nele_left = ye * hy_avogadro * U(DENS_VAR,i,j-1,k) 

        dyNe = 0.5*(nele_right - nele_left)*idy
                    
        if (hy_useCrossMagRes) then
           tele_left = U(tempVar,i,j-1,k)
           call PlasmaState_logLambda(ELE_ELE, tele_left, tele_left, &
                                      nele_left, abar/hy_avogadro, zbar, llee)
           call PlasmaState_logLambda(ELE_ION, tele_left, tele_left, &
                                      nele_left, abar/hy_avogadro, zbar, llei)
#ifdef FLASH_MULTISPECIES
           call Multispecies_getSumFrac(Z,zfullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j-1,k))
           zfull = abar*zfullFrac
           call Multispecies_getSumFrac(Z2,z2fullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j-1,k))
           z2full = abar*z2fullFrac
#else
           zfull = eos_singleSpeciesZ 
           z2full = zfull**2
#endif
           zeff = (z2full*llei) / (zfull*llee)
           call PlasmaState_tau(ELE_ION, tele_left, llei, zbar, nele_left/zbar, tau)
           call MagneticResistivity(U(:,i,j-1,k), etaparUnused, resCross=etacross)
           etatau_left = etacross*tau

           dyEtaTau = 0.5*(etatau_right - etatau_left)*idy                    
        end if
     end if

     if (NDIM > 2) then
        call Eos_getAbarZbar(U(:,i,j,k+1),abar=abar,zbar=zbar, Ye=ye)
        nele_right = ye * hy_avogadro * U(DENS_VAR,i,j,k+1) 

        if (hy_useCrossMagRes) then
           tele_right = U(tempVar,i,j,k+1)
           call PlasmaState_logLambda(ELE_ELE, tele_right, tele_right, &
                                      nele_right, abar/hy_avogadro, zbar, llee)
           call PlasmaState_logLambda(ELE_ION, tele_right, tele_right, &
                                      nele_right, abar/hy_avogadro, zbar, llei)
#ifdef FLASH_MULTISPECIES
           call Multispecies_getSumFrac(Z,zfullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j,k+1))
           zfull = abar*zfullFrac
           call Multispecies_getSumFrac(Z2,z2fullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j,k+1))
           z2full = abar*z2fullFrac
#else
           zfull = eos_singleSpeciesZ 
#endif
           zeff = (z2full*llei) / (zfull*llee)
           call PlasmaState_tau(ELE_ION, tele_right, llei, zbar, nele_right/zbar, tau)
           call MagneticResistivity(U(:,i,j,k+1), etaparUnused, resCross=etacross)
           etatau_right = etacross*tau
        end if

        call Eos_getAbarZbar(U(:,i,j,k-1),abar=abar,zbar=zbar, Ye=ye)
        nele_left = ye * hy_avogadro * U(DENS_VAR,i,j,k-1) 

        dzNe = 0.5*(nele_right - nele_left)*idz                   

        if (hy_useCrossMagRes) then
           tele_left = U(tempVar,i,j,k-1)
           call PlasmaState_logLambda(ELE_ELE, tele_left, tele_left, &
                                      nele_left, abar/hy_avogadro, zbar, llee)
           call PlasmaState_logLambda(ELE_ION, tele_left, tele_left, &
                                      nele_left, abar/hy_avogadro, zbar, llei)
#ifdef FLASH_MULTISPECIES
           call Multispecies_getSumFrac(Z,zfullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j,k-1))
           zfull = abar*zfullFrac
           call Multispecies_getSumFrac(Z2,z2fullFrac,U(SPECIES_BEGIN:SPECIES_END,i,j,k-1))
           z2full = abar*z2fullFrac
#else
           zfull = eos_singleSpeciesZ 
           z2full = zfull**2
#endif
           zeff = (z2full*llei) / (zfull*llee)
           call PlasmaState_tau(ELE_ION, tele_left, llei, zbar, nele_left/zbar, tau)
           call MagneticResistivity(U(:,i,j,k-1), etaparUnused, resCross=etacross)
           etatau_left = etacross*tau

           dzEtaTau = 0.5*(etatau_right - etatau_left)*idz                   
        end if
     end if

     Lne = sqrt(dxNe*dxNe + dyNe*dyNe + dzNe*dzNe)/nele !!! in /cm

     if (hy_useCrossMagRes) then
        Letatau = sqrt(dxEtaTau*dxEtaTau + dyEtaTau*dyEtaTau + dzEtaTau*dzEtaTau)/etatau !!! in /cm
     end if 

  end if  !! useHall or useBiermann3T 


  !! Compute electron temperature gradients for thermoelectric terms
  if (hy_useNernst .or. hy_useCrossField) then
     tele_right = U(tempVar,i+1,j,k)
     tele_left = U(tempVar,i-1,j,k)

     dxTe = 0.5*(tele_right - tele_left)*idx
     dyTe = 0.0
     dzTe = 0.0

     if (NDIM > 1) then
        tele_right = U(tempVar,i,j+1,k)
        tele_left = U(tempVar,i,j-1,k)
        dyTe = 0.5*(tele_right - tele_left)*idy
     end if

     if (NDIM > 2) then
        tele_right = U(tempVar,i,j,k+1)
        tele_left = U(tempVar,i,j,k-1)
        dzTe = 0.5*(tele_right - tele_left)*idz
     end if
  end if


  !! Compute currents for Hall and resistive terms
  jx = 0.0
  jy = 0.0
  jz = 0.0
  if (hy_useHall .or. hy_useAnisoMagRes .or. hy_useBiermann3T) then
     jy = -gfact*(U(MAGZ_VAR,i+1,j,k) - U(MAGZ_VAR,i-1,j,k))*0.5*idx
     jz = gfact*(U(MAGY_VAR,i+1,j,k) - U(MAGY_VAR,i-1,j,k))*0.5*idx
     if (hy_geometry == CYLINDRICAL) then
        jy = jy + U(MAGZ_VAR,i,j,k)/xCenter
     end if 
#if NDIM >= 2
     jx = gfact*(U(MAGZ_VAR,i,j+1,k) - U(MAGZ_VAR,i,j-1,k))*0.5*idy
     jz = jz - gfact*(U(MAGX_VAR,i,j+1,k) - U(MAGX_VAR,i,j-1,k))*0.5*idy
#if NDIM == 3
     jx = jx - gfact*(U(MAGY_VAR,i,j,k+1) - U(MAGY_VAR,i,j,k-1))*0.5*idz
     jy = jy + gfact*(U(MAGX_VAR,i,j,k+1) - U(MAGX_VAR,i,j,k-1))*0.5*idz
#endif
#endif
  end if


  !!-----------------------------------------------------------------------
  !!    HALL: Whistler waves + Drift waves + advection velocity
  !!          Cross resistivity affects all speeds 
  !!-----------------------------------------------------------------------
  if (hy_useHall) then
     !! Compute max Whistler wave speed if Hall effect is included
     whistlerSpeed = 2.0*PI*i_inerLen*sqrt(b2)

     !! If including cross magnetic resistivity, alter Whistler speed accordingly
     !! cross whistler ~ etacross*chi/B*(qele*nele/(qele*nele) =
     !! etatau*qele**2*nele/(mele*clight)*i_inerLen*sqrt(b2)
     if (hy_useCrossMagRes) then
        crossWhistlerSpeed = 2.0*PI*i_inerLen*sqrt(b2)*etatau*hy_qele**2*nele/hy_mele
        whistlerSpeed = whistlerSpeed + crossWhistlerSpeed
     end if

     if (NDIM == 1) whistlerSpeed = whistlerSpeed*idx
     if (NDIM == 2) whistlerSpeed = whistlerSpeed*max(idx,idy)
     if (NDIM == 3) whistlerSpeed = whistlerSpeed*max(idx,idy,idz)

     cfx2 = max(cfx2,whistlerSpeed*whistlerSpeed)
     cfy2 = max(cfy2,whistlerSpeed*whistlerSpeed)
     cfz2 = max(cfz2,whistlerSpeed*whistlerSpeed)

     !! Compute Hall drift wave speed if Hall effect is included
     HallDriftSpeed = i_inerLen*Lne*sqrt(b2)

     !! If including cross magnetic resistivity, alter Hall drift speed accordingly
     !! Note that cross drift speed has a minus sign, effectively reduces the Hall drift speed
     !! -etacross*tau*qele**2*nele/(mele*cight)*Letatau*i_inerLen*sqrt(b2)
     if (hy_useCrossMagRes) then
        crossDriftSpeed = -i_inerLen*etatau*hy_qele**2*nele/hy_mele*Letatau*sqrt(b2)
        HallDriftSpeed = HallDriftSpeed + crossDriftSpeed
     end if

     cfx2 = max(cfx2,HallDriftSpeed*HallDriftSpeed)
     cfy2 = max(cfy2,HallDriftSpeed*HallDriftSpeed)
     cfz2 = max(cfz2,HallDriftSpeed*HallDriftSpeed)

  end if  !! useHall

  vHall_X = 0.0
  vHall_Y = 0.0
  vHall_Z = 0.0

  !! vhall = -j/(qele*nele) - etacross*chi/B * j
  vHall_X = -jx*(1.0/(nele*hy_qele) + etatau*hy_qele/hy_mele)
  if (NDIM > 1) vHall_Y = -jy*(1.0/(nele*hy_qele) + etatau*hy_qele/hy_mele)
  if (NDIM > 2) vHall_Z = -jz*(1.0/(nele*hy_qele) + etatau*hy_qele/hy_mele)

  vEff_X = vEff_X + vHall_X
  if (NDIM > 1) vEff_Y = vEff_Y + vHall_Y
  if (NDIM > 2) vEff_Z = vEff_Z + vHall_Z

#ifdef VHLX_VAR
  U(VHLX_VAR,i,j,k) = vHall_X
#endif
#ifdef VHLY_VAR
  U(VHLY_VAR,i,j,k) = vHall_Y
#endif
#ifdef VHLZ_VAR
  U(VHLZ_VAR,i,j,k) = vHall_Z
#endif



  !!-----------------------------------------------------------------------
  !!    BIERMANN: Thermal magnetic waves
  !!-----------------------------------------------------------------------
  if (hy_useBiermann3T) then
     !!! Here we use sndspd2=gamma*Ptot/rho instead of Z*kb*Te/mi which is ok since the former is always larger...
     ThermalMagSpeed = sqrt((U(GAME_VAR,i,j,k) - 1.0)/U(GAMC_VAR,i,j,k))*i_inerLen*Lne*sqrt(sndspd2)

     cfx2 = max(cfx2,ThermalMagSpeed*ThermalMagSpeed)
     cfy2 = max(cfy2,ThermalMagSpeed*ThermalMagSpeed)
     cfz2 = max(cfz2,ThermalMagSpeed*ThermalMagSpeed)
  end if


  !!-----------------------------------------------------------------------
  !!    ANISOTROPIC RESISTIVITY: advection velocity
  !!          This term can be written differently to be diffusive, but in
  !!           this explicit method it is written as an advection term.
  !!-----------------------------------------------------------------------
  if (hy_useAnisoMagRes) then
     !! vresist = (etapar-etaperp)/B^2 * (B x j)
     Bcrossj_X = gfact*(U(MAGY_VAR,i,j,k)*jz - U(MAGZ_VAR,i,j,k)*jy)
     Bcrossj_Y = gfact*(U(MAGZ_VAR,i,j,k)*jx - U(MAGX_VAR,i,j,k)*jz)
     Bcrossj_Z = gfact*(U(MAGX_VAR,i,j,k)*jy - U(MAGY_VAR,i,j,k)*jx)
     if (Bmag == 0.0) then
        resistCoeff = 0.0
     else
        resistCoeff = (etapar-etaperp)*hy_qele*tauCent/(chi*Bmag*hy_mele)
     end if
     vResist_X = resistCoeff*BcrossJ_X
     if (NDIM > 1) vResist_Y = resistCoeff*BcrossJ_Y
     if (NDIM > 2) vResist_Z = resistCoeff*BcrossJ_Z

     vEff_X = vEff_X + vResist_X
     if (NDIM > 1) vEff_Y = vEff_Y + vResist_Y
     if (NDIM > 2) vEff_Z = vEff_Z + vResist_Z
  end if


  !!-----------------------------------------------------------------------
  !!    NERNST: advection velocity
  !!-----------------------------------------------------------------------
  if (hy_useNernst) then
     !! vnernst = -betacross*tau*boltz*gradTe / (mele*c)
     vNernst_X = -hy_boltz*tauCent*betacross/hy_mele*dxTe
     vNernst_Y = 0.0
     vNernst_Z = 0.0
     if (NDIM > 1) vNernst_Y = -hy_boltz*tauCent*betacross/hy_mele*dyTe
     if (NDIM > 2) vNernst_Z = -hy_boltz*tauCent*betacross/hy_mele*dzTe

     vMag = sqrt(vNernst_X**2+vNernst_Y**2+vNernst_Z**2) / & 
         (hy_nernstFlCoef*sqrt(hy_boltz*tele/hy_mele))
     select case(hy_nernstFlMode)
        case(FL_HARMONIC)
           vMag = 1.0 / (1.0 + vMag)
        case(FL_MINMAX)
           vMag = 1.0 / max(1.0, vMag)
        case(FL_LARSEN)
           vMag = 1.0 / sqrt(1.0 + vMag**2)
        case DEFAULT
           call Driver_abortFlash("[hy_uhd_computeDtExtMHD] Invalid Flux limiter type")
     end select
     vNernst_X = vMag * vNernst_X
     vNernst_Y = vMag * vNernst_Y
     vNernst_Z = vMag * vNernst_Z

#ifdef VNRX_VAR
     U(VNRX_VAR,i,j,k) = vNernst_X
#endif
#ifdef VNRY_VAR
     U(VNRY_VAR,i,j,k) = vNernst_Y
#endif
#ifdef VNRZ_VAR
     U(VNRZ_VAR,i,j,k) = vNernst_Z
#endif

     vEff_X = vEff_X + vNernst_X
     if (NDIM > 1) vEff_Y = vEff_Y + vNernst_Y
     if (NDIM > 2) vEff_Z = vEff_Z + vNernst_Z
  end if


  !!-----------------------------------------------------------------------
  !!    Thermoelectric Cross-Field: advection velocity
  !!-----------------------------------------------------------------------
  if (hy_useCrossField) then
     !! vcross = -(betapar-betaperp)*boltz / (qele*B**2) * (B x gradTe)
     BcrossgradTe_X = gfact*(U(MAGY_VAR,i,j,k)*dzTe - U(MAGZ_VAR,i,j,k)*dyTe)
     BcrossgradTe_Y = gfact*(U(MAGZ_VAR,i,j,k)*dxTe - U(MAGX_VAR,i,j,k)*dzTe)
     BcrossgradTe_Z = gfact*(U(MAGX_VAR,i,j,k)*dyTe - U(MAGY_VAR,i,j,k)*dxTe)
     if (Bmag == 0.0) then
        crossCoeff = 0.0
     else
        if (useBetaNew) then
           crossCoeff = -(betanew)*hy_boltz*tauCent/(Bmag*hy_mele)
        else
           crossCoeff = -(betapar-betaperp)*hy_boltz*tauCent/(chi*Bmag*hy_mele)
        endif
     endif
     vCross_X = crossCoeff*BcrossgradTe_X
     if (NDIM > 1) vCross_Y = crossCoeff*BcrossgradTe_Y
     if (NDIM > 2) vCross_Z = crossCoeff*BcrossgradTe_Z

     vMag = sqrt(vCross_X**2+vCross_Y**2+vCross_Z**2) / & 
         (hy_crossFieldFlCoef*sqrt(hy_boltz*tele/hy_mele))
     select case(hy_crossFieldFlMode)
        case(FL_NONE)
           vMag = 1.0
        case(FL_HARMONIC)
           vMag = 1.0 / (1.0 + vMag)
        case(FL_MINMAX)
           vMag = 1.0 / max(1.0, vMag)
        case(FL_LARSEN)
           vMag = 1.0 / sqrt(1.0 + vMag**2)
        case DEFAULT
           call Driver_abortFlash("[hy_uhd_computeDtExtMHD] Invalid Flux limiter type")
     end select
     vCross_X = vMag * vCross_X
     vCross_Y = vMag * vCross_Y
     vCross_Z = vMag * vCross_Z

#ifdef VRGX_VAR
     U(VRGX_VAR,i,j,k) = vCross_X
#endif
#ifdef VRGY_VAR
     U(VRGY_VAR,i,j,k) = vCross_Y
#endif
#ifdef VRGZ_VAR
     U(VRGZ_VAR,i,j,k) = vCross_Z
#endif

     vEff_X = vEff_X + vCross_X
     if (NDIM > 1) vEff_Y = vEff_Y + vCross_Y
     if (NDIM > 2) vEff_Z = vEff_Z + vCross_Z
  end if

#ifdef VEFX_VAR
  U(VEFX_VAR) = vEff_X
#endif
#ifdef VEFY_VAR
  U(VEFY_VAR) = vEff_Y
#endif
#ifdef VEFZ_VAR
  U(VEFZ_VAR) = vEff_Z
#endif

  !! max velocity is max of ion (fluid) velocity and effective electron velocity
  !! (with grid advection included)
  vMax_X = max(abs(U(VELX_VAR,i,j,k) - uxgrid),abs(vEff_X - uxgrid))
  if (NDIM > 1) vMax_Y = max(abs(U(VELY_VAR,i,j,k) - uygrid),abs(vEff_Y - uygrid))
  if (NDIM > 2) vMax_Z = max(abs(U(VELZ_VAR,i,j,k) - uzgrid),abs(vEff_Z - uzgrid))


  !dt = (vMax_X + sqrt(cfx2))*idx
  !if (NDIM > 1) dt = max(dt,(vMax_Y + sqrt(cfy2))*idy)
  !if (NDIM > 2) dt = max(dt,(vMax_Z + sqrt(cfz2))*idz)

#endif
  return                 
   
End Subroutine hy_uhd_computeDtExtMHD

