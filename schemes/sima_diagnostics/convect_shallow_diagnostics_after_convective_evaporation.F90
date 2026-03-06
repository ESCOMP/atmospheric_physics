! Diagnostics specific for shallow convective schemes
! using convective evaporation via ZM (zm_conv_evap)
! (i.e., hack_shallow)
module convect_shallow_diagnostics_after_convective_evaporation
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: convect_shallow_diagnostics_after_convective_evaporation_init

  ! The shallow convection diagnostics have several phases,
  ! which have to run before the tendency updaters in order to capture
  ! the appropriate tendencies.
  ! 2) after convective evaporation
  !    (for schemes using zm_conv_evap, i.e., hack_shallow)
  public :: convect_shallow_diagnostics_after_convective_evaporation_run

contains

!> \section arg_table_convect_shallow_diagnostics_after_convective_evaporation_init  Argument Table
!! \htmlinclude convect_shallow_diagnostics_after_convective_evaporation_init.html
  subroutine convect_shallow_diagnostics_after_convective_evaporation_init( &
    errmsg, errflg)
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

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
    ! and are thus output in convect_shallow_diagnostics_after_convective_evaporation_run.
    call history_add_field('PRECSH', 'lwe_precipitation_rate_at_surface_due_to_shallow_convection', horiz_only, 'avg', 'm s-1') ! Shallow Convection precipitation rate

  end subroutine convect_shallow_diagnostics_after_convective_evaporation_init

!> \section arg_table_convect_shallow_diagnostics_after_convective_evaporation_run  Argument Table
!! \htmlinclude convect_shallow_diagnostics_after_convective_evaporation_run.html
  subroutine convect_shallow_diagnostics_after_convective_evaporation_run( &
    ncol, pver, &
    cpair, &
    tend_s, tend_s_snwprd, tend_s_snwevmlt, &
    tend_q, &
    ntprpd, ntsnprd, flxprec, flxsnow, precc, &
    errmsg, errflg)

    use cam_history, only: history_out_field

    ! Input arguments
    integer,             intent(in)  :: ncol                 ! Number of columns
    integer,             intent(in)  :: pver                 ! Number of model levels
    real(kind_phys),     intent(in)  :: cpair                ! Specific heat of dry air [J kg-1 K-1]
    real(kind_phys),     intent(in)  :: tend_s(:,:)          ! Heating rate from evaporation [J kg-1 s-1]
    real(kind_phys),     intent(in)  :: tend_s_snwprd(:,:)   ! Heating rate from snow production [J kg-1 s-1]
    real(kind_phys),     intent(in)  :: tend_s_snwevmlt(:,:) ! Heating rate from snow evap/melt [J kg-1 s-1]
    real(kind_phys),     intent(in)  :: tend_q(:,:)          ! Water vapor tendency [kg kg-1 s-1]
    real(kind_phys),     intent(in)  :: ntprpd(:,:)          ! Net precipitation production [kg kg-1 s-1]
    real(kind_phys),     intent(in)  :: ntsnprd(:,:)         ! Net snow production [kg kg-1 s-1]
    real(kind_phys),     intent(in)  :: flxprec(:,:)         ! Precipitation flux [kg m-2 s-1]
    real(kind_phys),     intent(in)  :: flxsnow(:,:)         ! Snow flux [kg m-2 s-1]
    real(kind_phys),     intent(in)  :: precc(:)             ! Shallow precipitation rate [m s-1]

    ! Output arguments
    character(len=512),  intent(out) :: errmsg
    integer,             intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Temperature tendencies from energy tendencies
    ! (converted from J kg-1 s-1 to K s-1)
    call history_out_field('EVAPTCM', tend_s(:,:)/cpair)
    call history_out_field('FZSNTCM', tend_s_snwprd(:,:)/cpair)
    call history_out_field('EVSNTCM', tend_s_snwevmlt(:,:)/cpair)

    call history_out_field('HKNTPRPD', ntprpd)
    call history_out_field('HKNTSNPD', ntsnprd)
    call history_out_field('HKFLXPRC', flxprec)
    call history_out_field('HKFLXSNW', flxsnow)
    call history_out_field('HKEIHEAT', tend_s)
    call history_out_field('EVAPQCM',  tend_q)
    call history_out_field('PRECSH',   precc)

  end subroutine convect_shallow_diagnostics_after_convective_evaporation_run

end module convect_shallow_diagnostics_after_convective_evaporation
