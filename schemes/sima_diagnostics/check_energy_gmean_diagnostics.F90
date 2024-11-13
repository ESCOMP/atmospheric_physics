! Diagnostic scheme for check_energy_gmean
! This replicates the "debug printout" with the global mean total energy in CAM:
!  nstep, te        0   0.25872620595082335E+10   0.25872620595082335E+10  -0.00000000000000000E+00   0.10126767279027148E+06   0.21940670000000080E+03
! where the numbers correspond to:
! - timestep number
! - global mean input energy using dycore energy formula
! - global mean output energy at end of physics timestep using dycore energy formula
! - global mean heating rate for energy conservation
! - global mean surface pressure
! - global mean model top pressure
! These numbers are very useful for matching models bit-for-bit because they include "everything" in the model.
module check_energy_gmean_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: check_energy_gmean_diagnostics_run  ! main routine

CONTAINS

   !> \section arg_table_check_energy_gmean_diagnostics_run  Argument Table
   !! \htmlinclude check_energy_gmean_diagnostics_run.html
   subroutine check_energy_gmean_diagnostics_run( &
      amIRoot, &
      iulog, &
      nstep, &
      teinp_glob, &
      teout_glob, &
      heat_glob,  &
      psurf_glob, &
      ptopb_glob, &
      errmsg, errflg)

      logical,            intent(in)  :: amIRoot
      integer,            intent(in)  :: iulog         ! log output unit
      integer,            intent(in)  :: nstep         ! current timestep number

      real(kind_phys),    intent(in)  :: teinp_glob    ! global mean energy of input state [J m-2]
      real(kind_phys),    intent(in)  :: teout_glob    ! global mean energy of output state [J m-2]
      real(kind_phys),    intent(in)  :: psurf_glob    ! global mean surface pressure [Pa]
      real(kind_phys),    intent(in)  :: ptopb_glob    ! global mean top boundary pressure [Pa]
      real(kind_phys),    intent(in)  :: heat_glob     ! global mean heating rate [J kg-1 s-1]

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      if (amIRoot) then
         write(iulog,'(1x,a9,1x,i8,5(1x,e25.17))') "nstep, te", nstep, teinp_glob, teout_glob, &
               heat_glob, psurf_glob, ptopb_glob
      endif

   end subroutine check_energy_gmean_diagnostics_run

end module check_energy_gmean_diagnostics
