! Diagnostics for convective_cloud_water: the CAM conv_water history fields
! (in-cloud and grid-mean convective/total cloud water and ice, and the
! fractional occurrence of cloud with condensate).
module convective_cloud_water_diagnostics
   implicit none
   private

   public :: convective_cloud_water_diagnostics_init
   public :: convective_cloud_water_diagnostics_run

contains

   !> \section arg_table_convective_cloud_water_diagnostics_init  Argument Table
   !! \htmlinclude convective_cloud_water_diagnostics_init.html
   subroutine convective_cloud_water_diagnostics_init(errmsg, errflg)
      use cam_history, only: history_add_field

      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History add field calls (CAM conv_water field names)
      call history_add_field('ICLMRCU',  'in_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_all_convection', 'lev', 'avg', 'kg kg-1')
      call history_add_field('ICIMRCU',  'in_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_all_convection', 'lev', 'avg', 'kg kg-1')
      call history_add_field('ICLMRTOT', 'total_in_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')
      call history_add_field('ICIMRTOT', 'total_in_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')
      call history_add_field('GCLMRSH',  'gridbox_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1')
      call history_add_field('GCLMRDP',  'gridbox_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_deep_convection', 'lev', 'avg', 'kg kg-1')
      call history_add_field('GCIMRSH',  'gridbox_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1')
      call history_add_field('GCIMRDP',  'gridbox_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_deep_convection', 'lev', 'avg', 'kg kg-1')
      call history_add_field('FRESH',    'frequency_of_occurrence_of_shallow_convection_with_condensate', 'lev', 'avg', '1')
      call history_add_field('FREDP',    'frequency_of_occurrence_of_deep_convection_with_condensate', 'lev', 'avg', '1')
      call history_add_field('FRECU',    'frequency_of_occurrence_of_convection_with_condensate', 'lev', 'avg', '1')
      call history_add_field('FRETOT',   'frequency_of_occurrence_of_cloud_with_condensate_above_threshold', 'lev', 'avg', '1')

   end subroutine convective_cloud_water_diagnostics_init

   !> \section arg_table_convective_cloud_water_diagnostics_run  Argument Table
   !! \htmlinclude convective_cloud_water_diagnostics_run.html
   subroutine convective_cloud_water_diagnostics_run( &
      conv_liq, conv_ice, tot_liq, tot_ice, &
      totg_liq_sh, totg_liq_dp, totg_ice_sh, totg_ice_dp, &
      fresh, fredp, frecu, fretot, &
      errmsg, errflg)

      use ccpp_kinds, only: kind_phys
      use cam_history, only: history_out_field

      ! Input parameters
      real(kind_phys),    intent(in)  :: conv_liq(:,:)     ! Convective contributions to in-cloud liquid [kg kg-1]
      real(kind_phys),    intent(in)  :: conv_ice(:,:)     ! Convective contributions to in-cloud ice [kg kg-1]
      real(kind_phys),    intent(in)  :: tot_liq(:,:)      ! Total in-cloud liquid [kg kg-1]
      real(kind_phys),    intent(in)  :: tot_ice(:,:)      ! Total in-cloud ice [kg kg-1]
      real(kind_phys),    intent(in)  :: totg_liq_sh(:,:)  ! Grid-mean liquid from shallow convective cloud [kg kg-1]
      real(kind_phys),    intent(in)  :: totg_liq_dp(:,:)  ! Grid-mean liquid from deep convective cloud [kg kg-1]
      real(kind_phys),    intent(in)  :: totg_ice_sh(:,:)  ! Grid-mean ice from shallow convective cloud [kg kg-1]
      real(kind_phys),    intent(in)  :: totg_ice_dp(:,:)  ! Grid-mean ice from deep convective cloud [kg kg-1]
      real(kind_phys),    intent(in)  :: fresh(:,:)        ! Fractional occurrence of shallow cumulus with condensate [1]
      real(kind_phys),    intent(in)  :: fredp(:,:)        ! Fractional occurrence of deep cumulus with condensate [1]
      real(kind_phys),    intent(in)  :: frecu(:,:)        ! Fractional occurrence of cumulus with condensate [1]
      real(kind_phys),    intent(in)  :: fretot(:,:)       ! Fractional occurrence of cloud with condensate [1]

      ! CCPP error handling variables
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('ICLMRCU',  conv_liq)
      call history_out_field('ICIMRCU',  conv_ice)
      call history_out_field('ICLMRTOT', tot_liq)
      call history_out_field('ICIMRTOT', tot_ice)
      call history_out_field('GCLMRSH',  totg_liq_sh)
      call history_out_field('GCLMRDP',  totg_liq_dp)
      call history_out_field('GCIMRSH',  totg_ice_sh)
      call history_out_field('GCIMRDP',  totg_ice_dp)
      call history_out_field('FRESH',    fresh)
      call history_out_field('FREDP',    fredp)
      call history_out_field('FRECU',    frecu)
      call history_out_field('FRETOT',   fretot)

   end subroutine convective_cloud_water_diagnostics_run

end module convective_cloud_water_diagnostics
