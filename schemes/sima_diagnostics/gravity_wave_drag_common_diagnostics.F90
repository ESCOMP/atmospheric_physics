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
      call history_add_field('TTGW', 'tendency_of_air_temperature_due_to_gravity_wave_drag', 'lev', 'avg', 'K s-1')
      call history_add_field('UTGW_TOTAL', 'tendency_of_eastward_wind_due_to_gravity_wave_drag', 'lev', 'avg', 'm s-2')
      call history_add_field('VTGW_TOTAL', 'tendency_of_northward_wind_due_to_gravity_wave_drag', 'lev', 'avg', 'm s-2')

      call history_add_field('QTGW', 'tendency_of_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_gravity_wave_drag', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CLDLIQTGW', 'tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_gravity_wave_drag', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CLDICETGW', 'tendency_of_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_gravity_wave_drag', 'lev', 'avg', 'kg kg-1 s-1')
   end subroutine gravity_wave_drag_common_diagnostics_init

   !> \section arg_table_gravity_wave_drag_common_diagnostics_run  Argument Table
   !! \htmlinclude gravity_wave_drag_common_diagnostics_run.html
   subroutine gravity_wave_drag_common_diagnostics_run( &
      const_props, &
      cpairv, &
      egwdffi_tot, &
      tend_q, & ! all ccpp constituent tendencies
      tend_u, tend_v, tend_s, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      use runtime_obj, only: wv_stdname

      use ccpp_const_utils,     only: ccpp_const_get_idx

      ! framework dependency for const_props
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      ! Input parameters
      type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
      real(kind_phys),    intent(in)  :: cpairv(:,:)
      real(kind_phys),    intent(in)  :: egwdffi_tot(:,:)
      real(kind_phys),    intent(in)  :: tend_q(:,:,:)
      real(kind_phys),    intent(in)  :: tend_u(:,:)
      real(kind_phys),    intent(in)  :: tend_v(:,:)
      real(kind_phys),    intent(in)  :: tend_s(:,:)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables for constituent indices
      integer :: const_wv_idx, const_cldliq_idx, const_cldice_idx

      errmsg = ''
      errflg = 0

      ! Get constituent indices for wv, cldliq, cldice
      call ccpp_const_get_idx(const_props, &
           wv_stdname, &
           const_wv_idx, errmsg, errflg)
      if (errflg /= 0) return

      call ccpp_const_get_idx(const_props, &
           'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', &
           const_cldliq_idx, errmsg, errflg)
      if (errflg /= 0) return

      call ccpp_const_get_idx(const_props, &
           'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', &
           const_cldice_idx, errmsg, errflg)
      if (errflg /= 0) return

      ! History out field calls
      call history_out_field('EKGW', egwdffi_tot)
      call history_out_field('TTGW', tend_s/cpairv)
      call history_out_field('UTGW_TOTAL', tend_u)
      call history_out_field('VTGW_TOTAL', tend_v)

      call history_out_field('QTGW', tend_q(:,:,const_wv_idx))
      if(const_cldliq_idx > 0) then
         call history_out_field('CLDLIQTGW', tend_q(:,:,const_cldliq_idx))
      end if

      if(const_cldice_idx > 0) then
         call history_out_field('CLDICETGW', tend_q(:,:,const_cldice_idx))
      end if

   end subroutine gravity_wave_drag_common_diagnostics_run

end module gravity_wave_drag_common_diagnostics
