! Copyright (C) 2024 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Diagnostics for shallow convection and merged deep + shallow convection
! Haipeng Lin, December 2024
module convect_shallow_diagnostics
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: convect_shallow_diagnostics_init

  ! The shallow convection diagnostics have several phases,
  ! which have to run before the tendency updaters in order to capture
  ! the appropriate tendencies.
  ! 1) after the shallow convective scheme itself
  public :: convect_shallow_diagnostics_after_shallow_scheme_run
  ! 2) after convective evaporation (for schemes using zm_conv_evap)
  public :: convect_shallow_diagnostics_after_convective_evaporation_run
  ! 3) after shallow convective quantities are merged with deep convective quantities
  !    this outputs the final convection diagnostics taking into account both deep and shallow convection
  public :: convect_shallow_diagnostics_after_sum_to_deep_run

contains

!> \section arg_table_convect_shallow_diagnostics_init  Argument Table
!! \htmlinclude convect_shallow_diagnostics_init.html
  subroutine convect_shallow_diagnostics_init( &
    errmsg, errflg)
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    !=======================================================
    ! Convection diagnostics for shallow schemes (common)
    !=======================================================

    call history_add_field('CMFDT', 'tendency_of_air_temperature_at_constant_pressure_due_to_shallow_convection', 'lev', 'avg', 'K s-1') ! T tendency - shallow convection (ftem = ptend%s/cpair)
    call history_add_field('CMFDQ', 'tendency_of_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', '') ! ptend_loc%q(1,1,1)
    call history_add_field('CMFDLIQ', 'tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', '') ! ptend_loc%q(1,1,ixcldliq)
    call history_add_field('CMFDICE', 'tendency_of_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', '') ! ptend_loc%q(1,1,ixcldice)
    call history_add_field('CMFDQR', 'tendency_of_precipitation_wrt_moist_air_and_condensed_water_due_to_shallow_convection_excluding_subcloud_evaporation', 'lev', 'avg', 'kg kg-1 s-1') ! Q tendency - shallow convection rainout

    ! QC and DRP are the same from convect_shallow
    call history_add_field('QC', 'detrainment_of_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1') ! Q tendency - shallow convection LW export
    call history_add_field('DQP', 'detrainment_of_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1') ! Specific humidity tendency due to precipitation

    call history_add_field('ICWMRSH', 'in_cloud_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1') ! Shallow Convection in-cloud water mixing ratio

    call history_add_field('CMFSL', 'liquid_water_static_energy_flux_due_to_shallow_convection_tbd', 'ilev', 'avg', 'W m-2') ! Moist shallow convection liquid water static energy flux
    call history_add_field('CMFLQ', 'total_water_flux_in_energy_unit_due_to_shallow_convection_tbd', 'ilev', 'avg', 'W m-2') ! Moist shallow convection total water flux

    call history_add_field('FREQSH', 'frequency_of_shallow_convection_tbd', horiz_only, 'avg', 'fraction') ! Fractional occurence of shallow convection

    ! UW convection
    ! call history_add_field('CBMF', '', horiz_only, 'avg', 'kg m-2 s-1') ! Cloud base mass flux

    !=======================================================
    ! Convection diagnostics for shallow schemes
    ! with convective evaporation (i.e., Hack)
    ! (convect_shallow_diagnostics_after_convective_evaporation_run)
    !
    ! Standard names are derived from zm_conv_evap
    !=======================================================
    ! These are Hack shallow convection specific, but are just named _due_to_shallow_convection due to
    ! being renamed from generic _convection -> _shallow_convection in scheme set_general_conv_fluxes_to_shallow
    ! FIXME: the history field names are Hack-specific for backwards compatibility.
    call history_add_field('EVAPTCM', 'tendency_of_air_temperature_due_to_convective_evaporation_due_to_shallow_convection', 'lev', 'avg', 'K s-1') ! T tendency - Evaporation/snow prod from Hack convection = tend%s / cpair
    call history_add_field('FZSNTCM', 'tendency_of_air_temperature_due_to_frozen_precipitation_production_due_to_shallow_convection', 'lev', 'avg', 'K s-1') ! T tendency - Rain to snow conversion from Hack convection = tend_snowprd / cpair
    call history_add_field('EVSNTCM', 'tendency_of_air_temperature_due_to_evaporation_and_melting_of_frozen_precipitation_due_to_shallow_convection', 'lev', 'avg', 'K s-1') ! T tendency - Snow to rain prod from Hack convection = tend_snwevmlt / cpair
    call history_add_field('HKNTPRPD', 'tendency_of_precipitation_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1') ! Net precipitation production from HK convection
    call history_add_field('HKNTSNPD', 'tendency_of_frozen_precipitation_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1') ! Net snow production from HK convection
    call history_add_field('HKFLXPRC', 'precipitation_flux_at_interface_due_to_shallow_convection', 'ilev', 'avg', 'kg m-2 s-1') ! Flux of precipitation from HK convection
    call history_add_field('HKFLXSNW', 'frozen_precipitation_flux_at_interface_due_to_shallow_convection', 'ilev', 'avg', 'kg m-2 s-1') ! Flux of snow from HK convection

    ! Fields that need to be saved (duplicated) for history purposes after convective precipitation
    call history_add_field('HKEIHEAT', 'tendency_of_dry_air_enthalpy_at_constant_pressure_due_to_convective_evaporation_due_to_shallow_convection', 'lev', 'avg', 'W kg-1') ! Heating by ice and evaporation in HK convection
    call history_add_field('EVAPQCM', 'tendency_of_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_convective_evaporation_due_to_shallow_convection', 'lev', 'avg', 'K s-1') ! Q tendency - Evaporation from Hack convection (needs to be saved from tend_q output)

    ! For Hack only, the below quantities are modified after zm_conv_evap
    call history_add_field('PRECSH', 'lwe_precipitation_rate_at_surface_due_to_shallow_convection', horiz_only, 'avg', 'm s-1') ! Shallow Convection precipitation rate

    ! call history_add_field('', '', 'lev', 'avg', '')
    ! call history_add_field('', '', horiz_only, 'avg', '')
    ! call history_add_field('', '', 'ilev', 'avg', '')

    ! CMFMC is shallow+deep.
    ! Define a CMFMCSH for shallow only
    call history_add_field('CMFMCSH', 'atmosphere_convective_mass_flux_due_to_shallow_convection', 'ilev', 'avg', 'kg m-2 s-1') ! Moist convection (shallow) mass flux

    !=======================================================
    ! Convection diagnostics for shallow and deep combined
    ! (convect_shallow_diagnostics_after_sum_to_deep_run)
    !=======================================================
    call history_add_field('CLDTOP', 'vertical_index_at_cloud_top_for_all_convection', horiz_only, 'avg', '1') ! Vertical index of cloud top
    call history_add_field('CLDBOT', 'vertical_index_at_cloud_base_for_all_convection', horiz_only, 'avg', '1') ! Vertical index of cloud base
    call history_add_field('PCLDTOP', 'pressure_at_cloud_top_for_all_convection', horiz_only, 'avg', 'Pa') ! Pressure of cloud top
    call history_add_field('PCLDBOT', 'pressure_at_cloud_base_for_all_convection', horiz_only, 'avg', 'Pa') ! Pressure of cloud base

  end subroutine convect_shallow_diagnostics_init

!> \section arg_table_convect_shallow_diagnostics_after_shallow_scheme_run  Argument Table
!! \htmlinclude convect_shallow_diagnostics_after_shallow_scheme_run.html
  subroutine convect_shallow_diagnostics_after_shallow_scheme_run( &
    errmsg, errflg)

    use cam_history, only: history_out_field

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of shallow convection [fraction]

    errmsg = ''
    errflg = 0

    ! Calculate fractional occurrence of shallow convection
    ! FIXME (hplin): this definition looks counter-intuitive but is replicated verbatim from convect_shallow.F90. To check.
    freqsh(:) = 0._kind_phys
    do i = 1, ncol
      if(maxval(cmfmc_sh(i,:)) <= 0._kind_phys) then
        freqsh(i) = 1._kind_phys
      enddo
    enddo

  end subroutine convect_shallow_diagnostics_after_shallow_scheme_run

!> \section arg_table_convect_shallow_diagnostics_after_convective_evaporation_run  Argument Table
!! \htmlinclude convect_shallow_diagnostics_after_convective_evaporation_run.html
  subroutine convect_shallow_diagnostics_after_convective_evaporation_run( &
    errmsg, errflg)

    use cam_history, only: history_out_field

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables

    errmsg = ''
    errflg = 0



  end subroutine convect_shallow_diagnostics_after_convective_evaporation_run

!> \section arg_table_convect_shallow_diagnostics_after_sum_to_deep_run  Argument Table
!! \htmlinclude convect_shallow_diagnostics_after_sum_to_deep_run.html
  subroutine convect_shallow_diagnostics_after_sum_to_deep_run( &
    errmsg, errflg)

    use cam_history, only: history_out_field

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables

    errmsg = ''
    errflg = 0

    ! Calculate fractional occurrence of shallow convection
    ! FIXME (hplin): this definition looks counter-intuitive but is replicated verbatim from convect_shallow.F90. To check.
    freqsh(:) = 0._kind_phys
    do i = 1, ncol
      if(maxval(cmfmc_sh(i,:)) <= 0._kind_phys) then
        freqsh(i) = 1._kind_phys
      enddo
    enddo

  end subroutine convect_shallow_diagnostics_after_sum_to_deep_run

end module convect_shallow_diagnostics
