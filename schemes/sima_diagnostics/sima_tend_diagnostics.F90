module sima_tend_diagnostics

   use ccpp_kinds, only:  kind_phys
   use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

   implicit none
   private
   save

   public :: sima_tend_diagnostics_init ! init routine
   public :: sima_tend_diagnostics_run  ! main routine


CONTAINS

   !> \section arg_table_sima_tend_diagnostics_init  Argument Table
   !! \htmlinclude sima_tend_diagnostics_init.html
   subroutine sima_tend_diagnostics_init(errmsg, errflg)
      use cam_history, only: history_add_field
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Add tendency fields
      call history_add_field('TTEND', 'tendency_of_air_temperature_due_to_model_physics', 'lev', 'avg', 'K s-1')
      call history_add_field('UTEND', 'tendency_of_eastward_wind_due_to_model_physics',   'lev', 'avg', 'm s-2')
      call history_add_field('VTEND', 'tendency_of_northward_wind_due_to_model_physics',  'lev', 'avg', 'm s-2')

   end subroutine sima_tend_diagnostics_init

   !> \section arg_table_sima_tend_diagnostics_run  Argument Table
   !! \htmlinclude sima_tend_diagnostics_run.html
   subroutine sima_tend_diagnostics_run(dTdt_total, dudt_total, dvdt_total, errmsg, errflg)
      use cam_history, only: history_out_field
      ! Tendency variables
      real(kind_phys), intent(in) :: dTdt_total(:,:) ! tendency of air temperature due to model physics
      real(kind_phys), intent(in) :: dudt_total(:,:) ! tendency of eastward wind due to model physics
      real(kind_phys), intent(in) :: dvdt_total(:,:) ! tendency of northward wind due to model physics
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Capture tendency fields
      call history_out_field('TTEND', dTdt_total)
      call history_out_field('UTEND', dudt_total)
      call history_out_field('VTEND', dvdt_total)

   end subroutine sima_tend_diagnostics_run
   !=======================================================================
end module sima_tend_diagnostics
