! Diagnostics for cloud fraction - convective cloud cover
module convective_cloud_cover_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: convective_cloud_cover_diagnostics_init
   public :: convective_cloud_cover_diagnostics_run

contains

   !> \section arg_table_convective_cloud_cover_diagnostics_init  Argument Table
   !! \htmlinclude convective_cloud_cover_diagnostics_init.html
   subroutine convective_cloud_cover_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('SH_CLD', 'shallow_convective_cloud_area_fraction', 'lev', 'avg', 'fraction')
      call history_add_field('DP_CLD', 'deep_convective_cloud_area_fraction', 'lev', 'avg', 'fraction')
      call history_add_field('CONCLD', 'convective_cloud_area_fraction', 'lev', 'avg', 'fraction')

   end subroutine convective_cloud_cover_diagnostics_init

   !> \section arg_table_convective_cloud_cover_diagnostics_run  Argument Table
   !! \htmlinclude convective_cloud_cover_diagnostics_run.html
   subroutine convective_cloud_cover_diagnostics_run( &
      shallowcu, deepcu, concld, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input parameters
      real(kind_phys),    intent(in)  :: shallowcu(:, :)   ! Shallow convective cloud fraction [fraction]
      real(kind_phys),    intent(in)  :: deepcu(:, :)      ! Deep convective cloud fraction [fraction]
      real(kind_phys),    intent(in)  :: concld(:, :)      ! Convective cloud cover [fraction]


      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('SH_CLD', shallowcu)
      call history_out_field('DP_CLD', deepcu)
      call history_out_field('CONCLD', concld)

   end subroutine convective_cloud_cover_diagnostics_run

end module convective_cloud_cover_diagnostics
