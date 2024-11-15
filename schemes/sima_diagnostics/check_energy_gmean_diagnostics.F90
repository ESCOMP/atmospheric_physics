! Diagnostic scheme for check_energy_gmean
! This creates a debug printout with:
! - global mean input energy using dycore energy formula
! - global mean output energy at end of physics timestep using dycore energy formula
! - global mean heating rate for energy conservation
! These numbers are very useful for matching models bit-for-bit because they include "everything" in the model.
module check_energy_gmean_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: check_energy_gmean_diagnostics_init
   public :: check_energy_gmean_diagnostics_run  ! main routine

  ! Private module options.
  logical :: print_global_means = .false.

contains

!> \section arg_table_check_energy_gmean_diagnostics_init Argument Table
!! \htmlinclude arg_table_check_energy_gmean_diagnostics_init.html
  subroutine check_energy_gmean_diagnostics_init(print_global_means_in)
    ! Input arguments
    logical,            intent(in)    :: print_global_means_in

    print_global_means = print_global_means_in
  end subroutine check_energy_gmean_diagnostics_init

   !> \section arg_table_check_energy_gmean_diagnostics_run  Argument Table
   !! \htmlinclude check_energy_gmean_diagnostics_run.html
   subroutine check_energy_gmean_diagnostics_run( &
      amIRoot, &
      iulog, &
      teinp_glob, &
      teout_glob, &
      heat_glob,  &
      errmsg, errflg)

      logical,            intent(in)  :: amIRoot       ! are we on the MPI root task?
      integer,            intent(in)  :: iulog         ! log output unit

      real(kind_phys),    intent(in)  :: teinp_glob    ! global mean energy of input state [J m-2]
      real(kind_phys),    intent(in)  :: teout_glob    ! global mean energy of output state [J m-2]
      real(kind_phys),    intent(in)  :: heat_glob     ! global mean heating rate [J kg-1 s-1]

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      if (print_global_means .and. amIRoot) then
         write(iulog,'(1x,a26,1x,e25.17,1x,a15,1x,e25.17)') 'global mean input energy =', teinp_glob, 'output energy =', teout_glob
         write(iulog,'(1x,a38,1x,e25.17)') 'heating rate for energy conservation =', heat_glob
      endif

   end subroutine check_energy_gmean_diagnostics_run

end module check_energy_gmean_diagnostics
