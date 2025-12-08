! Diagnostics for all gravity wave drag parameterizations
! (used to output total diagnostics in the end of all gravity wave drag schemes)
module gravity_wave_drag_common_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: gravity_wave_drag_common_diagnostics_init
   public :: gravity_wave_drag_common_diagnostics_run

contains

   !> \section arg_table_gravity_wave_drag_common_diagnostics_init  Argument Table
   !! \htmlinclude gravity_wave_drag_common_diagnostics_init.html
   subroutine gravity_wave_drag_common_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('EKGW', 'effective_diffusivity_coefficient_at_interfaces_due_to_gravity_wave_drag', 'ilev', 'avg', 'm2 s-1')

   end subroutine gravity_wave_drag_common_diagnostics_init

   !> \section arg_table_gravity_wave_drag_common_diagnostics_run  Argument Table
   !! \htmlinclude gravity_wave_drag_common_diagnostics_run.html
   subroutine gravity_wave_drag_common_diagnostics_run( &
      egwdffi_tot, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input parameters
      real(kind_phys),    intent(in)  :: egwdffi_tot(:,:)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('EKGW', egwdffi_tot)

   end subroutine gravity_wave_drag_common_diagnostics_run

end module gravity_wave_drag_common_diagnostics
