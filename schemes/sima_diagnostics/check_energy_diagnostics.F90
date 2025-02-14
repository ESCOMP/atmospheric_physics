! Diagnostic scheme for check_energy
! Not all quantities are needed as diagnostics; this module is designed to ease development
module check_energy_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: check_energy_diagnostics_init ! init routine
   public :: check_energy_diagnostics_run  ! main routine

CONTAINS

   !> \section arg_table_check_energy_diagnostics_init  Argument Table
   !! \htmlinclude check_energy_diagnostics_init.html
   subroutine check_energy_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('cp_or_cv_dycore', 'specific_heat_of_air_used_in_dycore', 'lev', 'inst', 'J kg-1 K-1')
      call history_add_field('scaling_dycore', 'ratio_of_specific_heat_of_air_used_in_physics_energy_formula_to_specific_heat_of_air_used_in_dycore_energy_formula', 'lev', 'inst', '1')

      call history_add_field('te_cur_phys', 'vertically_integrated_total_energy_using_physics_energy_formula', horiz_only, 'inst', 'J m-2')
      call history_add_field('tw_cur', 'vertically_integrated_total_water', horiz_only, 'inst', 'kg m-2')

      call history_add_field('tend_te_tnd', 'cumulative_total_energy_boundary_flux_using_physics_energy_formula', horiz_only, 'inst', 'J m-2 s-1')
      call history_add_field('tend_tw_tnd', 'cumulative_total_water_boundary_flux', horiz_only, 'inst', 'kg m-2 s-1')

   end subroutine check_energy_diagnostics_init

   !> \section arg_table_check_energy_diagnostics_run  Argument Table
   !! \htmlinclude check_energy_diagnostics_run.html
   subroutine check_energy_diagnostics_run( &
      cp_or_cv_dycore, scaling_dycore, &
      te_cur_phys, tw_cur, &
      tend_te_tnd, tend_tw_tnd, &
      errmsg, errflg)

      use cam_history, only: history_out_field
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      ! State variables
      real(kind_phys), intent(in) :: cp_or_cv_dycore(:,:)
      real(kind_phys), intent(in) :: scaling_dycore(:,:)
      real(kind_phys), intent(in) :: te_cur_phys(:)
      real(kind_phys), intent(in) :: tw_cur(:)
      real(kind_phys), intent(in) :: tend_te_tnd(:)
      real(kind_phys), intent(in) :: tend_tw_tnd(:)


      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('cp_or_cv_dycore', cp_or_cv_dycore)
      call history_out_field('scaling_dycore', scaling_dycore)
      call history_out_field('te_cur_phys', te_cur_phys)
      call history_out_field('tw_cur', tw_cur)
      call history_out_field('tend_te_tnd', tend_te_tnd)
      call history_out_field('tend_tw_tnd', tend_tw_tnd)

   end subroutine check_energy_diagnostics_run

end module check_energy_diagnostics
