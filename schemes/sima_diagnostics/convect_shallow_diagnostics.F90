! Diagnostics for shallow convection and merged deep + shallow convection
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
  !    is now in convect_shallow_diagnostics_after_convective_evaporation
  !    scheme separate from this file.
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
    call history_add_field('CMFDQ', 'tendency_of_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1') ! ptend_loc%q(1,1,1)
    call history_add_field('CMFDLIQ', 'tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1') ! ptend_loc%q(1,1,ixcldliq)
    call history_add_field('CMFDICE', 'tendency_of_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1') ! ptend_loc%q(1,1,ixcldice)
    call history_add_field('CMFDQR', 'tendency_of_precipitation_wrt_moist_air_and_condensed_water_due_to_shallow_convection_excluding_subcloud_evaporation', 'lev', 'avg', 'kg kg-1 s-1') ! Q tendency - shallow convection rainout

    call history_add_field('DQP', 'detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1') ! Specific humidity tendency due to precipitation in shallow convection

    call history_add_field('ICWMRSH', 'in_cloud_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1') ! Shallow Convection in-cloud water mixing ratio

    call history_add_field('CMFSL', 'liquid_water_static_energy_flux_due_to_shallow_convection', 'ilev', 'avg', 'W m-2') ! Moist shallow convection liquid water static energy flux
    call history_add_field('CMFLQ', 'total_water_flux_due_to_shallow_convection', 'ilev', 'avg', 'W m-2') ! Moist shallow convection total water flux

    call history_add_field('FREQSH', 'frequency_of_shallow_convection', horiz_only, 'avg', 'fraction') ! Fractional occurence of shallow convection

    ! CMFMC is shallow+deep.
    ! Define a CMFMCSH for shallow only
    call history_add_field('CMFMCSH', 'atmosphere_convective_mass_flux_due_to_shallow_convection', 'ilev', 'avg', 'kg m-2 s-1') ! Moist convection (shallow) mass flux

    !=======================================================
    ! Convection diagnostics for shallow and deep combined
    ! (convect_shallow_diagnostics_after_sum_to_deep_run)
    !=======================================================
    call history_add_field('CMFMC',  'atmosphere_convective_mass_flux_due_to_all_convection', 'ilev', 'avg', 'kg m-2 s-1') ! Moist convection (shallow+deep) mass flux
    call history_add_field('CLDTOP', 'vertical_index_at_cloud_top_for_all_convection', horiz_only, 'avg', '1') ! Vertical index of cloud top
    call history_add_field('CLDBOT', 'vertical_index_at_cloud_base_for_all_convection', horiz_only, 'avg', '1') ! Vertical index of cloud base
    call history_add_field('PCLDTOP', 'pressure_at_cloud_top_for_all_convection', horiz_only, 'avg', 'Pa') ! Pressure of cloud top
    call history_add_field('PCLDBOT', 'pressure_at_cloud_base_for_all_convection', horiz_only, 'avg', 'Pa') ! Pressure of cloud base

    call history_add_field('ZMDLF', 'detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_all_convection', 'lev', 'avg', 'kg kg-1 s-1')
    call history_add_field('SHDLF', 'detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_shallow_convection', 'lev', 'avg', 'kg kg-1 s-1')

  end subroutine convect_shallow_diagnostics_init

!> \section arg_table_convect_shallow_diagnostics_after_shallow_scheme_run  Argument Table
!! \htmlinclude convect_shallow_diagnostics_after_shallow_scheme_run.html
  subroutine convect_shallow_diagnostics_after_shallow_scheme_run( &
    ncol, pver, pcnst, &
    const_props, &
    cpair, &
    cmfdt, dq, cmfdqr, &
    qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
    errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use cam_history, only: history_out_field

    ! Input arguments
    integer,                          intent(in)  :: ncol, pver, pcnst
    type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
    real(kind_phys),                  intent(in)  :: cpair
    real(kind_phys),                  intent(in)  :: cmfdt (:,:)
    real(kind_phys),                  intent(in)  :: dq    (:,:,:)
    real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
    real(kind_phys),                  intent(in)  :: qc_sh (:,:)
    real(kind_phys),                  intent(in)  :: icwmr (:,:)
    real(kind_phys),                  intent(in)  :: cmfsl (:,:)
    real(kind_phys),                  intent(in)  :: cmflq (:,:)
    real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

    ! Output arguments
    character(len=512),               intent(out) :: errmsg
    integer,                          intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of shallow convection [fraction]
    integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
    character(len=512)              :: const_standard_name
    integer                         :: i, m

    errmsg = ''
    errflg = 0

    ! Find constituent indices for water vapor, cloud liquid, cloud ice
    const_wv_idx     = -1
    const_cldliq_idx = -1
    const_cldice_idx = -1
    const_check_loop: do m = 1, pcnst
      call const_props(m)%standard_name(const_standard_name)
      if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_wv_idx = m
      endif

      if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldliq_idx = m
      endif

      if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldice_idx = m
      endif
    enddo const_check_loop

    ! CMFDT - temperature tendency
    ! Calculate temperature tendency from specific heat flux [K s-1]
    call history_out_field('CMFDT', cmfdt(:,:)/cpair)

    ! Constituent tendencies
    if (const_wv_idx > 0) then
      call history_out_field('CMFDQ',   dq(:,:,const_wv_idx))
    endif
    if (const_cldliq_idx > 0) then
      call history_out_field('CMFDLIQ', dq(:,:,const_cldliq_idx))
    endif
    if (const_cldice_idx > 0) then
      call history_out_field('CMFDICE', dq(:,:,const_cldice_idx))
    endif

    ! CMFDQR
    call history_out_field('CMFDQR', cmfdqr)

    ! Calculate fractional occurrence of shallow convection
    freqsh(:) = 0._kind_phys
    do i = 1, ncol
      if(maxval(cmfmc_sh(i,:)) > 0._kind_phys) then
        freqsh(i) = 1._kind_phys
      endif
    enddo
    call history_out_field('FREQSH', freqsh)

    call history_out_field('DQP',     qc_sh)
    call history_out_field('ICWMRSH', icwmr)
    call history_out_field('CMFSL',   cmfsl)
    call history_out_field('CMFLQ',   cmflq)
    call history_out_field('CMFMCSH', cmfmc_sh)

  end subroutine convect_shallow_diagnostics_after_shallow_scheme_run

!> \section arg_table_convect_shallow_diagnostics_after_sum_to_deep_run  Argument Table
!! \htmlinclude convect_shallow_diagnostics_after_sum_to_deep_run.html
  subroutine convect_shallow_diagnostics_after_sum_to_deep_run( &
    ncol, &
    cmfmc, cnt, cnb, p_cnt, p_cnb, &
    qc_total, qc_sh, &
    errmsg, errflg)

    use cam_history, only: history_out_field

    ! Input arguments
    integer,          intent(in)  :: ncol             ! Number of columns
    real(kind_phys),  intent(in)  :: cmfmc(:,:)       ! Total convective mass flux [kg m-2 s-1]
    integer,          intent(in)  :: cnt(:)           ! Cloud top level index [1]
    integer,          intent(in)  :: cnb(:)           ! Cloud base level index [1]
    real(kind_phys),  intent(in)  :: p_cnt(:)         ! Convective cloud top pressure [Pa]
    real(kind_phys),  intent(in)  :: p_cnb(:)         ! Convective cloud base pressure [Pa]
    real(kind_phys),  intent(in)  :: qc_total(:,:)    !  detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_all_convection [kg kg-1 s-1]
    real(kind_phys),  intent(in)  :: qc_sh(:,:)       !  detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_shallow_convection [kg kg-1 s-1]

    ! Output arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_out_field('CMFMC',   cmfmc)
    call history_out_field('CLDTOP',  real(cnt, kind_phys))
    call history_out_field('CLDBOT',  real(cnb, kind_phys))
    call history_out_field('PCLDTOP', p_cnt)
    call history_out_field('PCLDBOT', p_cnb)

    ! even though history notes this as ZMDLF
    ! (in rk_stratiform_tend)
    ! this appears to be dlf after shallow added
    call history_out_field('ZMDLF',   qc_total)
    call history_out_field('SHDLF',   qc_sh)

  end subroutine convect_shallow_diagnostics_after_sum_to_deep_run

end module convect_shallow_diagnostics
