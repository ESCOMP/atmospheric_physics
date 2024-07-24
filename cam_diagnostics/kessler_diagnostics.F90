module kessler_diagnostics

   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: kessler_diagnostics_init ! init routine
   public :: kessler_diagnostics_run  ! main routine

CONTAINS

   !> \section arg_table_kessler_diagnostics_init  Argument Table
   !! \htmlinclude kessler_diagnostics_init.html
   subroutine kessler_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      call history_add_field('PRECT', 'total_precipitation_rate_at_surface', horiz_only, 'avg', 'm s-1')

   end subroutine kessler_diagnostics_init

   !> \section arg_table_kessler_diagnostics_run  Argument Table
   !! \htmlinclude kessler_diagnostics_run.html
   subroutine kessler_diagnostics_run(precl, errmsg, errflg)

      use cam_history, only: history_out_field
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      ! State variables
      real(kind_phys), intent(in) :: precl(:) ! Total precipitation
      ! CCPP error handling variables
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call history_out_field('PRECT', precl)

   end subroutine kessler_diagnostics_run

   !=======================================================================

end module kessler_diagnostics
