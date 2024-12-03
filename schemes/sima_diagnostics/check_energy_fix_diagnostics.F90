! Diagnostic scheme for check_energy_fix
! This module includes diagnostics that have to be output
! right after the energy fixer has been ran and check_energy_chng has updated energy state
! (before the fluxes for check_energy_chng are zeroed out)
module check_energy_fix_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: check_energy_fix_diagnostics_init
   public :: check_energy_fix_diagnostics_run

contains

   !> \section arg_table_check_energy_fix_diagnostics_init  Argument Table
   !! \htmlinclude check_energy_fix_diagnostics_init.html
   subroutine check_energy_fix_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('teinp', 'vertically_integrated_total_energy_using_dycore_energy_formula_at_start_of_physics_timestep', horiz_only, 'inst', 'J m-2')
      call history_add_field('tefix', 'vertically_integrated_total_energy_using_dycore_energy_formula', horiz_only, 'inst', 'J m-2')
      call history_add_field('teout', 'vertically_integrated_total_energy_at_end_of_physics_timestep', horiz_only, 'inst', 'J m-2')
      call history_add_field('efix', 'net_sensible_heat_flux_through_top_and_bottom_of_atmosphere_column_from_global_total_energy_correction', horiz_only, 'inst', 'J m-2')

   end subroutine check_energy_fix_diagnostics_init

   !> \section arg_table_check_energy_fix_diagnostics_run  Argument Table
   !! \htmlinclude check_energy_diagnostics_run.html
   subroutine check_energy_fix_diagnostics_run( &
      te_ini_dyn, te_cur_dyn, &
      teout, &
      eshflx, &
      errmsg, errflg)

      use cam_history, only: history_out_field
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      ! State variables
      real(kind_phys), intent(in) :: te_ini_dyn(:)
      real(kind_phys), intent(in) :: te_cur_dyn(:)
      real(kind_phys), intent(in) :: teout(:)
      real(kind_phys), intent(in) :: eshflx(:)


      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('teinp', te_ini_dyn)
      call history_out_field('teout', teout)
      call history_out_field('tefix', te_cur_dyn)
      call history_out_field('efix',  eshflx)

   end subroutine check_energy_fix_diagnostics_run

end module check_energy_fix_diagnostics
