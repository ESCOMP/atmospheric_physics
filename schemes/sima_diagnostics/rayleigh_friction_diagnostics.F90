module rayleigh_friction_diagnostics

   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: rayleigh_friction_diagnostics_init ! init routine
   public :: rayleigh_friction_diagnostics_run  ! main routine

CONTAINS

   !> \section arg_table_rayleigh_friction_diagnostics_init  Argument Table
   !! \htmlinclude rayleigh_friction_diagnostics_init.html
   subroutine rayleigh_friction_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('UTEND_RAYLEIGH', 'Zonal wind tendency due to Rayleigh Friction', 'lev', 'avg', 'm s-2')
      call history_add_field('VTEND_RAYLEIGH', 'Meridional wind tendency due to Rayleigh Friction', 'lev', 'avg', 'm s-2')
      call history_add_field('STEND_RAYLEIGH', 'Dry air enthalpy tendency due to Rayleigh Friction', 'lev', 'avg', 'J kg-1')
      
   end subroutine rayleigh_friction_diagnostics_init

   !> \section arg_table_rayleigh_friction_diagnostics_run  Argument Table
   !! \htmlinclude rayleigh_friction_diagnostics_run.html
   subroutine rayleigh_friction_diagnostics_run(dudt, dvdt, dsdt, errmsg, errflg)

      use cam_history, only: history_out_field
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      ! State variables
      real(kind_phys), intent(in) :: dudt(:,:) !tendency_of_eastward_wind due to RF
      real(kind_phys), intent(in) :: dvdt(:,:) !tendency_of_northward_wind due to RF
      real(kind_phys), intent(in) :: dsdt(:,:) !tendency_of_dry_air_enthalpy_at_constant_pressure due to RF

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('UTEND_RAYLEIGH', dudt)
      call history_out_field('VTEND_RAYLEIGH', dvdt)
      call history_out_field('STEND_RAYLEIGH', dsdt)

   end subroutine rayleigh_friction_diagnostics_run

   !=======================================================================

end module rayleigh_friction_diagnostics
