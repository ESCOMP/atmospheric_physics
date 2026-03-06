! Diagnostic scheme for UW shallow convection
module uw_convect_shallow_diagnostics

   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: uw_convect_shallow_diagnostics_init
   public :: uw_convect_shallow_diagnostics_run

contains

   !> \section arg_table_uw_convect_shallow_diagnostics_init  Argument Table
   !! \htmlinclude uw_convect_shallow_diagnostics_init.html
   subroutine uw_convect_shallow_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! --- Fluxes at interfaces ---
      call history_add_field('qtflx_Cu',      'total_water_flux_due_to_shallow_convection_in_mass_units',                       'ilev', 'avg', 'kg m-2 s-1')
      call history_add_field('slflx_Cu',      'liquid_water_static_energy_flux_due_to_shallow_convection',                      'ilev', 'avg', 'W m-2')
      call history_add_field('uflx_Cu',       'zonal_momentum_flux_due_to_shallow_convection',                                 'ilev', 'avg', 'kg m-1 s-2')
      call history_add_field('vflx_Cu',       'meridional_momentum_flux_due_to_shallow_convection',                            'ilev', 'avg', 'kg m-1 s-1')

      ! --- Tendencies at layer centers ---
      call history_add_field('qtten_Cu',      'tendency_of_total_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection',   'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('slten_Cu',      'tendency_of_liquid_water_static_energy_due_to_shallow_convection',               'lev', 'avg', 'J kg-1 s-1')
      call history_add_field('uten_Cu',       'tendency_of_eastward_wind',                                                     'lev', 'avg', 'm s-2')
      call history_add_field('vten_Cu',       'tendency_of_northward_wind',                                                    'lev', 'avg', 'm s-2')
      call history_add_field('qvten_Cu',      'tendency_of_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection',   'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('qlten_Cu',      'tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('qiten_Cu',      'tendency_of_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection',     'lev', 'avg', 'kg kg-1 s-1')

      ! --- Scalar diagnostics (horizontal only) ---
      call history_add_field('cbmf_Cu',         'shallow_convective_cloud_base_mass_flux',                                     horiz_only, 'avg', 'kg m-2 s-1')
      call history_add_field('ufrcinvbase_Cu',  'updraft_area_fraction_at_pbl_top_due_to_shallow_convection',                  horiz_only, 'avg', 'fraction')
      call history_add_field('ufrclcl_Cu',      'updraft_area_fraction_at_lifting_condensation_level_due_to_shallow_convection', horiz_only, 'avg', 'fraction')
      call history_add_field('winvbase_Cu',     'updraft_vertical_velocity_at_pbl_top_due_to_shallow_convection',              horiz_only, 'avg', 'm s-1')
      call history_add_field('wlcl_Cu',         'updraft_vertical_velocity_at_lifting_condensation_level_due_to_shallow_convection', horiz_only, 'avg', 'm s-1')
      call history_add_field('plcl_Cu',         'air_pressure_at_lifting_condensation_level_due_to_shallow_convection',        horiz_only, 'avg', 'Pa')
      call history_add_field('pinv_Cu',         'air_pressure_at_inversion_base_for_shallow_convection',                       horiz_only, 'avg', 'Pa')
      call history_add_field('plfc_Cu',         'air_pressure_at_level_of_free_convection_due_to_shallow_convection',          horiz_only, 'avg', 'Pa')
      call history_add_field('pbup_Cu',         'air_pressure_at_positive_buoyancy_top_due_to_shallow_convection',             horiz_only, 'avg', 'Pa')
      call history_add_field('ppen_Cu',         'air_pressure_at_cumulus_top_due_to_shallow_convection',                       horiz_only, 'avg', 'Pa')
      call history_add_field('qtsrc_Cu',        'source_air_total_water_mixing_ratio_due_to_shallow_convection',               horiz_only, 'avg', 'kg kg-1')
      call history_add_field('thlsrc_Cu',       'source_air_liquid_water_potential_temperature_due_to_shallow_convection',     horiz_only, 'avg', 'K')
      call history_add_field('thvlsrc_Cu',      'source_air_liquid_virtual_potential_temperature_due_to_shallow_convection',   horiz_only, 'avg', 'K')
      call history_add_field('emfkbup_Cu',      'penetrative_entrainment_mass_flux_at_top_level_of_positive_cloud_buoyancy_due_to_shallow_convection', horiz_only, 'avg', 'kg m-2 s-1')
      call history_add_field('cin_Cu',          'convective_inhibition_to_level_of_free_convection_due_to_shallow_convection', horiz_only, 'avg', 'J kg-1')
      call history_add_field('cinlcl_Cu',       'convective_inhibition_to_lifting_condensation_level_due_to_shallow_convection', horiz_only, 'avg', 'J kg-1')
      call history_add_field('cbmflimit_Cu',    'cloud_base_mass_flux_limiter_due_to_shallow_convection',                      horiz_only, 'avg', 'kg m-2 s-1')
      call history_add_field('tkeavg_Cu',       'pbl_averaged_tke_due_to_shallow_convection',                                  horiz_only, 'avg', 'm2 s-2')
      call history_add_field('zinv_Cu',         'inversion_base_height_for_shallow_convection',                                horiz_only, 'avg', 'm')
      call history_add_field('rcwp_Cu',         'cumulus_cloud_water_path_due_to_shallow_convection',                           horiz_only, 'avg', 'kg m-2')
      call history_add_field('rlwp_Cu',         'cumulus_liquid_water_path_due_to_shallow_convection',                          horiz_only, 'avg', 'kg m-2')
      call history_add_field('riwp_Cu',         'cumulus_ice_water_path_due_to_shallow_convection',                             horiz_only, 'avg', 'kg m-2')
      call history_add_field('tophgt_Cu',       'shallow_convective_scale_height',                                             horiz_only, 'avg', 'm')

      ! --- Updraft profiles at interfaces ---
      call history_add_field('wu_Cu',           'updraft_vertical_velocity_due_to_shallow_convection',                         'ilev', 'avg', 'm s-1')
      call history_add_field('ufrc_Cu',         'updraft_area_fraction_due_to_shallow_convection',                             'ilev', 'avg', 'fraction')
      call history_add_field('qtu_Cu',          'updraft_total_water_mixing_ratio_due_to_shallow_convection',                  'ilev', 'avg', 'kg kg-1')
      call history_add_field('thlu_Cu',         'updraft_liquid_water_potential_temperature_due_to_shallow_convection',        'ilev', 'avg', 'K')
      call history_add_field('thvu_Cu',         'updraft_virtual_potential_temperature_due_to_shallow_convection',             'ilev', 'avg', 'K')
      call history_add_field('uu_Cu',           'updraft_eastward_wind_due_to_shallow_convection',                             'ilev', 'avg', 'm s-1')
      call history_add_field('vu_Cu',           'updraft_northward_wind_due_to_shallow_convection',                            'ilev', 'avg', 'm s-1')
      call history_add_field('qtu_emf_Cu',      'penetrative_entrainment_total_water_mixing_ratio_due_to_shallow_convection',  'ilev', 'avg', 'kg kg-1')
      call history_add_field('thlu_emf_Cu',     'penetrative_entrainment_liquid_water_potential_temperature_due_to_shallow_convection', 'ilev', 'avg', 'K')
      call history_add_field('uu_emf_Cu',       'penetrative_entrainment_eastward_wind_due_to_shallow_convection',             'ilev', 'avg', 'm s-1')
      call history_add_field('vu_emf_Cu',       'penetrative_entrainment_northward_wind_due_to_shallow_convection',            'ilev', 'avg', 'm s-1')
      call history_add_field('umf_Cu',          'atmosphere_convective_mass_flux_due_to_shallow_convection',                   'ilev', 'avg', 'kg m-2 s-1')
      call history_add_field('uemf_Cu',         'net_upward_mass_flux_due_to_shallow_convection',                              'ilev', 'avg', 'kg m-2 s-1')

      ! --- In-cumulus cloud properties at layer centers ---
      call history_add_field('qcu_Cu',          'in_cloud_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection',  'lev', 'avg', 'kg kg-1')
      call history_add_field('qlu_Cu',          'in_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1')
      call history_add_field('qiu_Cu',          'in_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection',          'lev', 'avg', 'kg kg-1')
      call history_add_field('cufrc_Cu',        'shallow_convective_cloud_area_fraction_from_shallow_convection',              'lev', 'avg', 'fraction')
      call history_add_field('fer_Cu',          'fractional_lateral_entrainment_rate_due_to_shallow_convection',               'lev', 'avg', 'Pa-1')
      call history_add_field('fdr_Cu',          'fractional_lateral_detrainment_rate_due_to_shallow_convection',               'lev', 'avg', 'Pa-1')

      ! --- Precipitation and flux fields ---
      ! For Hack, PRECSH is output in convect_shallow_diagnostics_after_convective_evaporation
      call history_add_field('PRECSH',          'shallow_convection_precipitation_rate',                                        horiz_only, 'avg', 'm s-1')
      call history_add_field('UWFLXPRC',        'precipitation_flux_at_interface_due_to_uw_shallow_convection',                 'ilev', 'avg', 'kg m-2 s-1')
      call history_add_field('UWFLXSNW',        'snow_flux_at_interface_due_to_uw_shallow_convection',                          'ilev', 'avg', 'kg m-2 s-1')

      ! --- Precipitation microphysics ---
      call history_add_field('dwten_Cu',        'expelled_cloud_liquid_water_tendency_due_to_shallow_convection',               'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('diten_Cu',        'expelled_cloud_ice_tendency_due_to_shallow_convection',                        'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('qrten_Cu',        'tendency_of_rain_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('qsten_Cu',        'tendency_of_snow_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('flxrain_Cu',      'rain_flux_at_interface_due_to_shallow_convection',                             'ilev', 'avg', 'kg m-2 s-1')
      call history_add_field('flxsnow_Cu',      'snow_flux_at_interface_due_to_shallow_convection',                             'ilev', 'avg', 'kg m-2 s-1')
      call history_add_field('ntraprd_Cu',      'tendency_of_precipitation_wrt_moist_air_and_condensed_water_due_to_shallow_convection',        'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('ntsnprd_Cu',      'tendency_of_frozen_precipitation_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1')

      ! --- Buoyancy sorting diagnostics ---
      call history_add_field('excessu_Cu',      'updraft_saturation_excess_due_to_shallow_convection',                         'lev', 'avg', '1')
      call history_add_field('excess0_Cu',      'environmental_saturation_excess_due_to_shallow_convection',                   'lev', 'avg', '1')
      call history_add_field('xc_Cu',           'critical_mixing_fraction_due_to_shallow_convection',                          'lev', 'avg', '1')
      call history_add_field('aquad_Cu',        'buoyancy_sorting_quadratic_coefficient_a_due_to_shallow_convection',          'lev', 'avg', '1')
      call history_add_field('bquad_Cu',        'buoyancy_sorting_quadratic_coefficient_b_due_to_shallow_convection',          'lev', 'avg', '1')
      call history_add_field('cquad_Cu',        'buoyancy_sorting_quadratic_coefficient_c_due_to_shallow_convection',          'lev', 'avg', '1')
      call history_add_field('bogbot_Cu',       'shallow_convective_cloud_buoyancy_at_base_interface',                         'lev', 'avg', '1')
      call history_add_field('bogtop_Cu',       'shallow_convective_cloud_buoyancy_at_top_interface',                          'lev', 'avg', '1')

      ! --- Exit condition flags ---
      call history_add_field('exit_UWCu_Cu',     'exit_flag_for_uwshcu_no_convection',                                         horiz_only, 'avg', '1')
      call history_add_field('exit_conden_Cu',   'exit_flag_for_uwshcu_condensation_failure',                                  horiz_only, 'avg', '1')
      call history_add_field('exit_klclpver_Cu', 'exit_flag_for_uwshcu_lcl_at_model_top',                                      horiz_only, 'avg', '1')
      call history_add_field('exit_klfcpver_Cu', 'exit_flag_for_uwshcu_lfc_at_model_top',                                      horiz_only, 'avg', '1')
      call history_add_field('exit_ufrc_Cu',     'exit_flag_for_uwshcu_updraft_fraction',                                      horiz_only, 'avg', '1')
      call history_add_field('exit_wtw_Cu',      'exit_flag_for_uwshcu_vertical_velocity_squared',                             horiz_only, 'avg', '1')
      call history_add_field('exit_drycore_Cu',  'exit_flag_for_uwshcu_dry_core',                                              horiz_only, 'avg', '1')
      call history_add_field('exit_wu_Cu',       'exit_flag_for_uwshcu_excessive_vertical_velocity',                           horiz_only, 'avg', '1')
      call history_add_field('exit_cufilter_Cu', 'exit_flag_for_uwshcu_cumulus_filter',                                        horiz_only, 'avg', '1')
      call history_add_field('exit_kinv1_Cu',    'exit_flag_for_uwshcu_inversion_base_at_surface',                             horiz_only, 'avg', '1')
      call history_add_field('exit_rei_Cu',      'exit_flag_for_uwshcu_mixing_rate',                                           horiz_only, 'avg', '1')

      ! --- Limiter flags ---
      call history_add_field('limit_shcu_Cu',    'limit_flag_for_uwshcu_forced_shallow_convection',                            horiz_only, 'avg', '1')
      call history_add_field('limit_negcon_Cu',  'limit_flag_for_uwshcu_negative_condensate',                                  horiz_only, 'avg', '1')
      call history_add_field('limit_ufrc_Cu',    'limit_flag_for_uwshcu_updraft_fraction',                                     horiz_only, 'avg', '1')
      call history_add_field('limit_ppen_Cu',    'limit_flag_for_uwshcu_penetration_depth',                                    horiz_only, 'avg', '1')
      call history_add_field('limit_emf_Cu',     'limit_flag_for_uwshcu_entrainment_mass_flux',                                horiz_only, 'avg', '1')
      call history_add_field('limit_cinlcl_Cu',  'limit_flag_for_uwshcu_cin_to_lcl',                                           horiz_only, 'avg', '1')
      call history_add_field('limit_cin_Cu',     'limit_flag_for_uwshcu_cin',                                                  horiz_only, 'avg', '1')
      call history_add_field('limit_cbmf_Cu',    'limit_flag_for_uwshcu_cloud_base_mass_flux',                                 horiz_only, 'avg', '1')
      call history_add_field('limit_rei_Cu',     'limit_flag_for_uwshcu_mixing_rate',                                          horiz_only, 'avg', '1')
      call history_add_field('ind_delcin_Cu',    'indicator_for_uwshcu_explicit_cin_used',                                     horiz_only, 'avg', '1')

   end subroutine uw_convect_shallow_diagnostics_init

   !> \section arg_table_uw_convect_shallow_diagnostics_run  Argument Table
   !! \htmlinclude uw_convect_shallow_diagnostics_run.html
   subroutine uw_convect_shallow_diagnostics_run( &
      ncol, &
      ! Fluxes at interfaces
      cmflq, cmfsl, cmfmc_sh, &
      uflx_diag, vflx_diag, &
      ! Tendencies at layer centers
      qtten_diag, slten_diag, uten, vten, &
      qvten_diag, qlten_diag, qiten_diag, &
      ! Scalar diagnostics
      cbmf, ufrcinvbase_diag, ufrclcl_diag, &
      winvbase_diag, wlcl_diag, &
      plcl_diag, pinv_diag, plfc_diag, pbup_diag, ppen_diag, &
      qtsrc_diag, thlsrc_diag, thvlsrc_diag, &
      emfkbup_diag, cinh_diag, cinlclh_diag, &
      cbmflimit_diag, tkeavg_diag, zinv_diag, &
      rcwp_diag, rlwp_diag, riwp_diag, &
      cush, &
      ! Updraft profiles at interfaces
      wu_diag, ufrc_diag, &
      qtu_diag, thlu_diag, thvu_diag, uu_diag, vu_diag, &
      qtu_emf_diag, thlu_emf_diag, uu_emf_diag, vu_emf_diag, &
      uemf_diag, &
      ! In-cumulus cloud properties at layer centers
      qcu, qlu, qiu, shfrc, &
      fer_out, fdr_out, &
      ! Precipitation and flux fields
      precip_sh, &
      flxrain_diag, flxsnow_diag, &
      ! Precipitation microphysics
      dwten_diag, diten_diag, &
      qrten, qsten, &
      ntraprd_diag, ntsnprd_diag, &
      ! Buoyancy sorting diagnostics
      excessu_arr_diag, excess0_arr_diag, xc_arr_diag, &
      aquad_arr_diag, bquad_arr_diag, cquad_arr_diag, &
      bogbot_arr_diag, bogtop_arr_diag, &
      ! Exit condition flags
      exit_UWCu_diag, exit_conden_diag, &
      exit_klclmkx_diag, exit_klfcmkx_diag, &
      exit_ufrc_diag, exit_wtw_diag, exit_drycore_diag, &
      exit_wu_diag, exit_cufilter_diag, &
      exit_kinv1_diag, exit_rei_diag, &
      ! Limiter flags
      limit_shcu_diag, limit_negcon_diag, limit_ufrc_diag, &
      limit_ppen_diag, limit_emf_diag, limit_cinlcl_diag, &
      limit_cin_diag, limit_cbmf_diag, limit_rei_diag, &
      ind_delcin_diag, &
      ! Latent heat for unit conversion
      latvap, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input parameters
      integer,            intent(in)  :: ncol

      ! Fluxes at interfaces
      real(kind_phys),    intent(in)  :: cmflq(:,:)           ! total water flux due to shallow convection [W m-2]
      real(kind_phys),    intent(in)  :: cmfsl(:,:)           ! liquid water static energy flux due to shallow convection [W m-2]
      real(kind_phys),    intent(in)  :: cmfmc_sh(:,:)        ! convective mass flux due to shallow convection [kg m-2 s-1] (interfaces)
      real(kind_phys),    intent(in)  :: uflx_diag(:,:)       ! zonal momentum flux due to shallow convection [kg m-1 s-2] (interfaces)
      real(kind_phys),    intent(in)  :: vflx_diag(:,:)       ! meridional momentum flux due to shallow convection [kg m-1 s-2] (interfaces)

      ! Tendencies at layer centers
      real(kind_phys),    intent(in)  :: qtten_diag(:,:)      ! tendency of total water mixing ratio [kg kg-1 s-1]
      real(kind_phys),    intent(in)  :: slten_diag(:,:)      ! tendency of liquid water static energy [J kg-1 s-1]
      real(kind_phys),    intent(in)  :: uten(:,:)            ! tendency of eastward wind [m s-2]
      real(kind_phys),    intent(in)  :: vten(:,:)            ! tendency of northward wind [m s-2]
      real(kind_phys),    intent(in)  :: qvten_diag(:,:)      ! tendency of water vapor mixing ratio [kg kg-1 s-1]
      real(kind_phys),    intent(in)  :: qlten_diag(:,:)      ! tendency of cloud liquid water mixing ratio [kg kg-1 s-1]
      real(kind_phys),    intent(in)  :: qiten_diag(:,:)      ! tendency of cloud ice mixing ratio [kg kg-1 s-1]

      ! Scalar diagnostics (horizontal only)
      real(kind_phys),    intent(in)  :: cbmf(:)              ! cloud base mass flux [kg m-2 s-1]
      real(kind_phys),    intent(in)  :: ufrcinvbase_diag(:)  ! updraft area fraction at PBL top [fraction]
      real(kind_phys),    intent(in)  :: ufrclcl_diag(:)      ! updraft area fraction at LCL [fraction]
      real(kind_phys),    intent(in)  :: winvbase_diag(:)     ! updraft vertical velocity at PBL top [m s-1]
      real(kind_phys),    intent(in)  :: wlcl_diag(:)         ! updraft vertical velocity at LCL [m s-1]
      real(kind_phys),    intent(in)  :: plcl_diag(:)         ! LCL pressure [Pa]
      real(kind_phys),    intent(in)  :: pinv_diag(:)         ! PBL top pressure [Pa]
      real(kind_phys),    intent(in)  :: plfc_diag(:)         ! LFC pressure [Pa]
      real(kind_phys),    intent(in)  :: pbup_diag(:)         ! positive buoyancy top pressure [Pa]
      real(kind_phys),    intent(in)  :: ppen_diag(:)         ! cumulus top pressure [Pa]
      real(kind_phys),    intent(in)  :: qtsrc_diag(:)        ! source air total water mixing ratio [kg kg-1]
      real(kind_phys),    intent(in)  :: thlsrc_diag(:)       ! source air liquid water potential temperature [K]
      real(kind_phys),    intent(in)  :: thvlsrc_diag(:)      ! source air liquid virtual potential temperature [K]
      real(kind_phys),    intent(in)  :: emfkbup_diag(:)      ! penetrative entrainment mass flux at kbup [kg m-2 s-1]
      real(kind_phys),    intent(in)  :: cinh_diag(:)         ! CIN to LFC [J kg-1]
      real(kind_phys),    intent(in)  :: cinlclh_diag(:)      ! CIN to LCL [J kg-1]
      real(kind_phys),    intent(in)  :: cbmflimit_diag(:)    ! cloud base mass flux limiter [kg m-2 s-1]
      real(kind_phys),    intent(in)  :: tkeavg_diag(:)       ! PBL-averaged TKE [m2 s-2]
      real(kind_phys),    intent(in)  :: zinv_diag(:)         ! PBL top height [m]
      real(kind_phys),    intent(in)  :: rcwp_diag(:)         ! cumulus cloud water path (LWP+IWP) [kg m-2]
      real(kind_phys),    intent(in)  :: rlwp_diag(:)         ! cumulus liquid water path [kg m-2]
      real(kind_phys),    intent(in)  :: riwp_diag(:)         ! cumulus ice water path [kg m-2]
      real(kind_phys),    intent(in)  :: cush(:)              ! shallow convective scale height [m]

      ! Updraft profiles at interfaces
      real(kind_phys),    intent(in)  :: wu_diag(:,:)         ! updraft vertical velocity [m s-1] (interfaces)
      real(kind_phys),    intent(in)  :: ufrc_diag(:,:)       ! updraft area fraction [fraction] (interfaces)
      real(kind_phys),    intent(in)  :: qtu_diag(:,:)        ! updraft total water mixing ratio [kg kg-1] (interfaces)
      real(kind_phys),    intent(in)  :: thlu_diag(:,:)       ! updraft liquid water potential temperature [K] (interfaces)
      real(kind_phys),    intent(in)  :: thvu_diag(:,:)       ! updraft virtual potential temperature [K] (interfaces)
      real(kind_phys),    intent(in)  :: uu_diag(:,:)         ! updraft eastward wind [m s-1] (interfaces)
      real(kind_phys),    intent(in)  :: vu_diag(:,:)         ! updraft northward wind [m s-1] (interfaces)
      real(kind_phys),    intent(in)  :: qtu_emf_diag(:,:)    ! penetrative entrainment total water [kg kg-1] (interfaces)
      real(kind_phys),    intent(in)  :: thlu_emf_diag(:,:)   ! penetrative entrainment liquid water potential temperature [K] (interfaces)
      real(kind_phys),    intent(in)  :: uu_emf_diag(:,:)     ! penetrative entrainment eastward wind [m s-1] (interfaces)
      real(kind_phys),    intent(in)  :: vu_emf_diag(:,:)     ! penetrative entrainment northward wind [m s-1] (interfaces)
      real(kind_phys),    intent(in)  :: uemf_diag(:,:)       ! net upward mass flux [kg m-2 s-1] (interfaces)

      ! In-cumulus cloud properties at layer centers
      real(kind_phys),    intent(in)  :: qcu(:,:)             ! in-cloud total water [kg kg-1]
      real(kind_phys),    intent(in)  :: qlu(:,:)             ! in-cloud liquid water [kg kg-1]
      real(kind_phys),    intent(in)  :: qiu(:,:)             ! in-cloud ice [kg kg-1]
      real(kind_phys),    intent(in)  :: shfrc(:,:)           ! shallow convective cloud fraction [fraction]
      real(kind_phys),    intent(in)  :: fer_out(:,:)         ! fractional lateral entrainment rate [Pa-1]
      real(kind_phys),    intent(in)  :: fdr_out(:,:)         ! fractional lateral detrainment rate [Pa-1]

      ! Precipitation and flux fields
      real(kind_phys),    intent(in)  :: precip_sh(:)         ! shallow convection precipitation rate [m s-1]
      real(kind_phys),    intent(in)  :: flxrain_diag(:,:)    ! rain flux at interface [kg m-2 s-1] (interfaces)
      real(kind_phys),    intent(in)  :: flxsnow_diag(:,:)    ! snow flux at interface [kg m-2 s-1] (interfaces)

      ! Precipitation microphysics
      real(kind_phys),    intent(in)  :: dwten_diag(:,:)      ! expelled cloud liquid water tendency [kg kg-1 s-1]
      real(kind_phys),    intent(in)  :: diten_diag(:,:)      ! expelled cloud ice tendency [kg kg-1 s-1]
      real(kind_phys),    intent(in)  :: qrten(:,:)           ! tendency of rain [kg kg-1 s-1]
      real(kind_phys),    intent(in)  :: qsten(:,:)           ! tendency of snow [kg kg-1 s-1]
      real(kind_phys),    intent(in)  :: ntraprd_diag(:,:)    ! net rain production tendency [kg kg-1 s-1]
      real(kind_phys),    intent(in)  :: ntsnprd_diag(:,:)    ! net snow production tendency [kg kg-1 s-1]

      ! Buoyancy sorting diagnostics
      real(kind_phys),    intent(in)  :: excessu_arr_diag(:,:)  ! updraft saturation excess [1]
      real(kind_phys),    intent(in)  :: excess0_arr_diag(:,:)  ! environmental saturation excess [1]
      real(kind_phys),    intent(in)  :: xc_arr_diag(:,:)       ! critical mixing fraction [1]
      real(kind_phys),    intent(in)  :: aquad_arr_diag(:,:)    ! buoyancy sorting quadratic coefficient a [1]
      real(kind_phys),    intent(in)  :: bquad_arr_diag(:,:)    ! buoyancy sorting quadratic coefficient b [1]
      real(kind_phys),    intent(in)  :: cquad_arr_diag(:,:)    ! buoyancy sorting quadratic coefficient c [1]
      real(kind_phys),    intent(in)  :: bogbot_arr_diag(:,:)   ! cloud buoyancy at base interface [1]
      real(kind_phys),    intent(in)  :: bogtop_arr_diag(:,:)   ! cloud buoyancy at top interface [1]

      ! Exit condition flags
      real(kind_phys),    intent(in)  :: exit_UWCu_diag(:)
      real(kind_phys),    intent(in)  :: exit_conden_diag(:)
      real(kind_phys),    intent(in)  :: exit_klclmkx_diag(:)
      real(kind_phys),    intent(in)  :: exit_klfcmkx_diag(:)
      real(kind_phys),    intent(in)  :: exit_ufrc_diag(:)
      real(kind_phys),    intent(in)  :: exit_wtw_diag(:)
      real(kind_phys),    intent(in)  :: exit_drycore_diag(:)
      real(kind_phys),    intent(in)  :: exit_wu_diag(:)
      real(kind_phys),    intent(in)  :: exit_cufilter_diag(:)
      real(kind_phys),    intent(in)  :: exit_kinv1_diag(:)
      real(kind_phys),    intent(in)  :: exit_rei_diag(:)

      ! Limiter flags
      real(kind_phys),    intent(in)  :: limit_shcu_diag(:)
      real(kind_phys),    intent(in)  :: limit_negcon_diag(:)
      real(kind_phys),    intent(in)  :: limit_ufrc_diag(:)
      real(kind_phys),    intent(in)  :: limit_ppen_diag(:)
      real(kind_phys),    intent(in)  :: limit_emf_diag(:)
      real(kind_phys),    intent(in)  :: limit_cinlcl_diag(:)
      real(kind_phys),    intent(in)  :: limit_cin_diag(:)
      real(kind_phys),    intent(in)  :: limit_cbmf_diag(:)
      real(kind_phys),    intent(in)  :: limit_rei_diag(:)
      real(kind_phys),    intent(in)  :: ind_delcin_diag(:)

      ! Latent heat for unit conversion
      real(kind_phys),    intent(in)  :: latvap                 ! latent heat of vaporization [J kg-1]

      ! CCPP error handling
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! --- Fluxes at interfaces ---
      ! qtflx_Cu: CAM divides cmflq (W m-2) by latvap to get kg m-2 s-1
      call history_out_field('qtflx_Cu',      cmflq / latvap)
      call history_out_field('slflx_Cu',      cmfsl)
      call history_out_field('uflx_Cu',       uflx_diag)
      call history_out_field('vflx_Cu',       vflx_diag)

      ! --- Tendencies at layer centers ---
      call history_out_field('qtten_Cu',      qtten_diag)
      call history_out_field('slten_Cu',      slten_diag)
      call history_out_field('uten_Cu',       uten)
      call history_out_field('vten_Cu',       vten)
      call history_out_field('qvten_Cu',      qvten_diag)
      call history_out_field('qlten_Cu',      qlten_diag)
      call history_out_field('qiten_Cu',      qiten_diag)

      ! --- Scalar diagnostics ---
      call history_out_field('cbmf_Cu',         cbmf)
      call history_out_field('ufrcinvbase_Cu',  ufrcinvbase_diag)
      call history_out_field('ufrclcl_Cu',      ufrclcl_diag)
      call history_out_field('winvbase_Cu',     winvbase_diag)
      call history_out_field('wlcl_Cu',         wlcl_diag)
      call history_out_field('plcl_Cu',         plcl_diag)
      call history_out_field('pinv_Cu',         pinv_diag)
      call history_out_field('plfc_Cu',         plfc_diag)
      call history_out_field('pbup_Cu',         pbup_diag)
      call history_out_field('ppen_Cu',         ppen_diag)
      call history_out_field('qtsrc_Cu',        qtsrc_diag)
      call history_out_field('thlsrc_Cu',       thlsrc_diag)
      call history_out_field('thvlsrc_Cu',      thvlsrc_diag)
      call history_out_field('emfkbup_Cu',      emfkbup_diag)
      call history_out_field('cin_Cu',          cinh_diag)
      call history_out_field('cinlcl_Cu',       cinlclh_diag)
      call history_out_field('cbmflimit_Cu',    cbmflimit_diag)
      call history_out_field('tkeavg_Cu',       tkeavg_diag)
      call history_out_field('zinv_Cu',         zinv_diag)
      call history_out_field('rcwp_Cu',         rcwp_diag)
      call history_out_field('rlwp_Cu',         rlwp_diag)
      call history_out_field('riwp_Cu',         riwp_diag)
      call history_out_field('tophgt_Cu',       cush)

      ! --- Updraft profiles at interfaces ---
      call history_out_field('wu_Cu',           wu_diag)
      call history_out_field('ufrc_Cu',         ufrc_diag)
      call history_out_field('qtu_Cu',          qtu_diag)
      call history_out_field('thlu_Cu',         thlu_diag)
      call history_out_field('thvu_Cu',         thvu_diag)
      call history_out_field('uu_Cu',           uu_diag)
      call history_out_field('vu_Cu',           vu_diag)
      call history_out_field('qtu_emf_Cu',      qtu_emf_diag)
      call history_out_field('thlu_emf_Cu',     thlu_emf_diag)
      call history_out_field('uu_emf_Cu',       uu_emf_diag)
      call history_out_field('vu_emf_Cu',       vu_emf_diag)
      call history_out_field('umf_Cu',          cmfmc_sh)
      call history_out_field('uemf_Cu',         uemf_diag)

      ! --- In-cumulus cloud properties at layer centers ---
      call history_out_field('qcu_Cu',          qcu)
      call history_out_field('qlu_Cu',          qlu)
      call history_out_field('qiu_Cu',          qiu)
      call history_out_field('cufrc_Cu',        shfrc)
      call history_out_field('fer_Cu',          fer_out)
      call history_out_field('fdr_Cu',          fdr_out)

      ! --- Precipitation and flux fields ---
      ! For Hack, PRECSH is output in convect_shallow_diagnostics_after_convective_evaporation
      call history_out_field('PRECSH',          precip_sh)
      call history_out_field('UWFLXPRC',        flxrain_diag)
      call history_out_field('UWFLXSNW',        flxsnow_diag)

      ! --- Precipitation microphysics ---
      call history_out_field('dwten_Cu',        dwten_diag)
      call history_out_field('diten_Cu',        diten_diag)
      call history_out_field('qrten_Cu',        qrten)
      call history_out_field('qsten_Cu',        qsten)
      call history_out_field('flxrain_Cu',      flxrain_diag)
      call history_out_field('flxsnow_Cu',      flxsnow_diag)
      call history_out_field('ntraprd_Cu',      ntraprd_diag)
      call history_out_field('ntsnprd_Cu',      ntsnprd_diag)

      ! --- Buoyancy sorting diagnostics ---
      call history_out_field('excessu_Cu',      excessu_arr_diag)
      call history_out_field('excess0_Cu',      excess0_arr_diag)
      call history_out_field('xc_Cu',           xc_arr_diag)
      call history_out_field('aquad_Cu',        aquad_arr_diag)
      call history_out_field('bquad_Cu',        bquad_arr_diag)
      call history_out_field('cquad_Cu',        cquad_arr_diag)
      call history_out_field('bogbot_Cu',       bogbot_arr_diag)
      call history_out_field('bogtop_Cu',       bogtop_arr_diag)

      ! --- Exit condition flags ---
      call history_out_field('exit_UWCu_Cu',     exit_UWCu_diag)
      call history_out_field('exit_conden_Cu',   exit_conden_diag)
      call history_out_field('exit_klclpver_Cu', exit_klclmkx_diag)
      call history_out_field('exit_klfcpver_Cu', exit_klfcmkx_diag)
      call history_out_field('exit_ufrc_Cu',     exit_ufrc_diag)
      call history_out_field('exit_wtw_Cu',      exit_wtw_diag)
      call history_out_field('exit_drycore_Cu',  exit_drycore_diag)
      call history_out_field('exit_wu_Cu',       exit_wu_diag)
      call history_out_field('exit_cufilter_Cu', exit_cufilter_diag)
      call history_out_field('exit_kinv1_Cu',    exit_kinv1_diag)
      call history_out_field('exit_rei_Cu',      exit_rei_diag)

      ! --- Limiter flags ---
      call history_out_field('limit_shcu_Cu',    limit_shcu_diag)
      call history_out_field('limit_negcon_Cu',  limit_negcon_diag)
      call history_out_field('limit_ufrc_Cu',    limit_ufrc_diag)
      call history_out_field('limit_ppen_Cu',    limit_ppen_diag)
      call history_out_field('limit_emf_Cu',     limit_emf_diag)
      call history_out_field('limit_cinlcl_Cu',  limit_cinlcl_diag)
      call history_out_field('limit_cin_Cu',     limit_cin_diag)
      call history_out_field('limit_cbmf_Cu',    limit_cbmf_diag)
      call history_out_field('limit_rei_Cu',     limit_rei_diag)
      call history_out_field('ind_delcin_Cu',    ind_delcin_diag)

   end subroutine uw_convect_shallow_diagnostics_run

end module uw_convect_shallow_diagnostics
