!!****if* source/Simulation/SimulationMain/Sedov/IO_outputFinal
!!
!! NAME
!!
!!  IO_outputFinal
!!
!!
!! SYNOPSIS
!!
!!  call IO_outputFinal()
!!
!!
!! DESCRIPTION
!!
!!  This routine is called after the code has exited the main timestep
!!  loop.  It outputs the last checkpoint, plotfile, and particle plotfile.
!!  This function should also clean up the IO unit, deallocate memory, etc.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  This implementation is slightly customized for the Sedov setup.
!!  The intent of the customization is to ensure that IO_writeIntegralQuantities
!!  (also customized) is called at least once at the normal end of a run, no
!!  matter why the run ends (nend, tmax, or perhaps a different criterion).
!!
!! SIDE EFFECTS
!!
!!  The state of module level logical variable io_outputInStack.
!!
!! SEE ALSO
!!
!!  IO_writeIntegralQuantities
!!***

subroutine IO_outputFinal()

#include "constants.h"  
  use Driver_interface, ONLY: Driver_getSimTime
  use IO_data, ONLY : io_justCheckpointed, io_outputInStack, &
       io_summaryOutputOnly
  use IO_interface, ONLY : IO_writeIntegralQuantities, IO_writeCheckpoint, IO_writePlotfile, &
    IO_writeParticles

  implicit none

  real :: simTime


  !This setting is used to ensure valid data throughout grid and ancestor blocks
  io_outputInStack = .true.

  call Driver_getSimTime(simTime)
  call IO_writeIntegralQuantities( 0, simTime)

  if (.not.io_summaryOutputOnly) then
     if(.not. io_justCheckpointed) then  
        call IO_writeCheckpoint()
     end if
   
     call IO_writePlotfile(.true.)
     
     call IO_writeParticles(.false.)
  end if

  io_outputInStack = .false.

end subroutine IO_outputFinal
