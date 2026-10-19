module micro_pumas_ccpp_dimensions_pre

  implicit none
  private

  public :: micro_pumas_ccpp_dimensions_pre_init
  public :: micro_pumas_ccpp_dimensions_pre_run

contains

  !> \section arg_table_micro_pumas_ccpp_dimensions_pre_init Argument Table
  !! \htmlinclude micro_pumas_ccpp_dimensions_pre_init.html
  subroutine micro_pumas_ccpp_dimensions_pre_init( ncol, nlev, nlevp1, ndust,       &
                             trop_cloud_top_lev, micro_ncol, micro_nlev,            &
                             micro_nlevp1,                                          &
                             pumas_airT, pumas_airq, pumas_cldliq, pumas_cldice, &
                             pumas_numliq, pumas_numice, pumas_rainliq, &
                             pumas_snowice, pumas_numrain, pumas_numsnow, &
                             pumas_graupice, pumas_numgraup, pumas_relvar, &
                             pumas_accre_enhan, pumas_pmid, pumas_pdel, pumas_pint, &
                             pumas_strat_cldfrc, pumas_strat_liq_cldfrc, &
                             pumas_strat_ice_cldfrc, pumas_qsatfac, pumas_naai, &
                             pumas_npccn, pumas_rndst, pumas_nacon, &
                             pumas_snowice_tend_external, &
                             pumas_numsnow_tend_external, pumas_effi_external, &
                             pumas_frzimm, pumas_frzcnt, pumas_frzdep, &
                             pumas_qcsinksum_rate1ord_out, pumas_airT_tend_out, &
                             pumas_airq_tend_out, pumas_cldliq_tend_out, &
                             pumas_cldice_tend_out, pumas_numliq_tend_out, &
                             pumas_numice_tend_out, pumas_rainliq_tend_out, &
                             pumas_snowice_tend_out, pumas_numrain_tend_out, &
                             pumas_numsnow_tend_out, pumas_graupice_tend_out, &
                             pumas_numgraup_tend_out, pumas_effc_out, &
                             pumas_effc_fn_out, pumas_effi_out, pumas_sadice_out, &
                             pumas_sadsnow_out, pumas_prect_out, pumas_preci_out, &
                             pumas_prec_evap_out, pumas_am_evap_st_out, &
                             pumas_prec_prod_out, pumas_cmeice_out, pumas_deffi_out, &
                             pumas_pgamrad_out, pumas_lamcrad_out, &
                             pumas_snowice_in_prec_out, pumas_scaled_diam_snow_out, &
                             pumas_graupice_in_prec_out, &
                             pumas_numgraup_vol_in_prec_out, &
                             pumas_scaled_diam_graup_out, pumas_lflx_out, &
                             pumas_iflx_out, pumas_gflx_out, pumas_rflx_out, &
                             pumas_sflx_out, pumas_rainliq_in_prec_out, &
                             pumas_reff_rain_out, pumas_reff_snow_out, &
                             pumas_reff_grau_out, pumas_numrain_vol_in_prec_out, &
                             pumas_numsnow_vol_in_prec_out, pumas_refl_out, &
                             pumas_arefl_out, pumas_areflz_out, pumas_frefl_out, &
                             pumas_csrfl_out, pumas_acsrfl_out, pumas_fcsrfl_out, &
                             pumas_refl10cm_out, pumas_reflz10cm_out, &
                             pumas_rercld_out, pumas_ncai_out, pumas_ncal_out, &
                             pumas_rainliq_out, pumas_snowice_out, &
                             pumas_numrain_vol_out, pumas_numsnow_vol_out, &
                             pumas_diam_rain_out, pumas_diam_snow_out, &
                             pumas_graupice_out, pumas_numgraup_vol_out, &
                             pumas_diam_graup_out, pumas_freq_graup_out, &
                             pumas_freq_snow_out, pumas_freq_rain_out, &
                             pumas_frac_ice_out, pumas_frac_cldliq_tend_out, &
                             pumas_rain_evap_out, &
                             errmsg, errcode)

    use pumas_kinds,       only: pumas_r8=>kind_r8

    !Host model dimensions/parameters:
    integer,         intent(in) :: ncol
    integer,         intent(in) :: nlev
    integer,         intent(in) :: nlevp1
    integer,         intent(in) :: ndust               !Number of dust size bins
    integer,         intent(in) :: trop_cloud_top_lev  !Index of the top model level for which
                                                       !cloud physics is applied (1 to nlev)

    !PUMAS dimensions/parameters:
    integer,         intent(out) :: micro_ncol         !Number of horizontal microphysics columns (count)
    integer,         intent(out) :: micro_nlev         !Number of microphysics vertical layers (count)
    integer,         intent(out) :: micro_nlevp1       !Number of microphysics vertical interfaces (count)

    !PUMAS-grid suite variables:
    !These are dimensioned by the microphysics dimensions set above, which are not
    !known when the suite allocates its non-allocatable data (before any init phase),
    !so they are allocated here instead.  This happens only once per run, as an
    !allocatable intent(out) dummy argument in the run phase would be deallocated
    !and reallocated on every call.
    ! microphysics air temperature of new state (K)
    real(pumas_r8), allocatable, intent(out) :: pumas_airT(:, :)
    ! microphysics water vapor mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_airq(:, :)
    ! microphysics cloud liquid mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_cldliq(:, :)
    ! microphysics cloud ice mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_cldice(:, :)
    ! microphysics mass number concentration of cloud liquid wrt moist air and condensed water of new state (kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numliq(:, :)
    ! microphysics mass number concentration of cloud ice wrt moist air and condensed water of new state (kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numice(:, :)
    ! microphysics rain mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_rainliq(:, :)
    ! microphysics snow mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_snowice(:, :)
    ! microphysics mass number concentration of rain wrt moist air and condensed water of new state (kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numrain(:, :)
    ! microphysics mass number concentration of snow wrt moist air and condensed water of new state (kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numsnow(:, :)
    ! microphysics graupel mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_graupice(:, :)
    ! microphysics mass number concentration of graupel wrt moist air and condensed water of new state (kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numgraup(:, :)
    ! microphysics relative variance of cloud water (1)
    real(pumas_r8), allocatable, intent(out) :: pumas_relvar(:, :)
    ! microphysics spatial varying accretion enhancement factor (1)
    real(pumas_r8), allocatable, intent(out) :: pumas_accre_enhan(:, :)
    ! microphysics air pressure (Pa)
    real(pumas_r8), allocatable, intent(out) :: pumas_pmid(:, :)
    ! microphysics air pressure thickness (Pa)
    real(pumas_r8), allocatable, intent(out) :: pumas_pdel(:, :)
    ! microphysics air pressure at interface (Pa)
    real(pumas_r8), allocatable, intent(out) :: pumas_pint(:, :)
    ! microphysics stratiform cloud area fraction (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_strat_cldfrc(:, :)
    ! microphysics stratiform cloud liquid area fraction (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_strat_liq_cldfrc(:, :)
    ! microphysics stratiform cloud ice area fraction (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_strat_ice_cldfrc(:, :)
    ! microphysics subgrid cloud water saturation scaling factor (1)
    real(pumas_r8), allocatable, intent(out) :: pumas_qsatfac(:, :)
    ! microphysics tendency of activated ice nuclei mass number concentration (kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_naai(:, :)
    ! microphysics tendency of activated cloud condensation nuclei mass number concentration (kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_npccn(:, :)
    ! microphysics dust bin radii (m)
    real(pumas_r8), allocatable, intent(out) :: pumas_rndst(:, :, :)
    ! microphysics dust number concentration by size bin (m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_nacon(:, :, :)
    ! microphysics tendency of snow mixing ratio wrt moist air and condensed water from external microphysics (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_snowice_tend_external(:, :)
    ! microphysics tendency of mass number concentration of snow wrt moist air and condensed water from external microphysics (kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numsnow_tend_external(:, :)
    ! microphysics effective radius of stratiform cloud ice particle from external microphysics (m)
    real(pumas_r8), allocatable, intent(out) :: pumas_effi_external(:, :)
    ! microphysics tendency of cloud liquid droplet number concentration due to immersion freezing (cm-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_frzimm(:, :)
    ! microphysics tendency of cloud liquid droplet number concentration due to contact freezing (cm-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_frzcnt(:, :)
    ! microphysics tendency of cloud ice number concentration due to deposition nucleation (cm-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_frzdep(:, :)
    ! microphysics direct conversion rate of stratiform cloud water to precipitation (s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_qcsinksum_rate1ord_out(:, :)
    ! microphysics tendency of dry air enthalpy at constant pressure (J kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_airT_tend_out(:, :)
    ! microphysics tendency of water vapor mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_airq_tend_out(:, :)
    ! microphysics tendency of cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_cldliq_tend_out(:, :)
    ! microphysics tendency of cloud ice mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_cldice_tend_out(:, :)
    ! mass number concentration of cloud liquid water wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numliq_tend_out(:, :)
    ! mass number concentration of cloud ice wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numice_tend_out(:, :)
    ! microphysics tendency of rain mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_rainliq_tend_out(:, :)
    ! microphysics tendency of rain mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_snowice_tend_out(:, :)
    ! microphysics tendency of mass number concentration of rain wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numrain_tend_out(:, :)
    ! microphysics tendency of mass number concentration of snow wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numsnow_tend_out(:, :)
    ! microphysics tendency of graupel water mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_graupice_tend_out(:, :)
    ! microphysics tendency of mass number concentration of graupel wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_numgraup_tend_out(:, :)
    ! microphysics effective radius of stratiform cloud liquid water particle (um)
    real(pumas_r8), allocatable, intent(out) :: pumas_effc_out(:, :)
    ! microphysics effective radius of stratiform cloud liquid water particle assuming cloud water number concentration of 1e8 per kg air (um)
    real(pumas_r8), allocatable, intent(out) :: pumas_effc_fn_out(:, :)
    ! microphysics effective radius of stratiform cloud ice particle (um)
    real(pumas_r8), allocatable, intent(out) :: pumas_effi_out(:, :)
    ! microphysics cloud ice surface area density (cm2 cm-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_sadice_out(:, :)
    ! microphysics snow surface area density (cm2 cm-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_sadsnow_out(:, :)
    ! microphysics total liquid water equivalent large scale precipitation rate at surface (m s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_prect_out(:)
    ! microphysics LWE large scale snowfall rate at surface (m s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_preci_out(:)
    ! microphysics precipitation evaporation rate wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_prec_evap_out(:, :)
    ! microphysics precipitation evaporation area (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_am_evap_st_out(:, :)
    ! microphysics precipitation production rate wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_prec_prod_out(:, :)
    ! microphysics condensation minus evaporation rate of in-cloud ice wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_cmeice_out(:, :)
    ! microphysics effective diameter of stratiform cloud ice particles for radiation (um)
    real(pumas_r8), allocatable, intent(out) :: pumas_deffi_out(:, :)
    ! microphysics cloud particle size distribution shape (gamma) parameter (1)
    real(pumas_r8), allocatable, intent(out) :: pumas_pgamrad_out(:, :)
    ! microphysics cloud particle size distribution slope (lambda) parameter (1)
    real(pumas_r8), allocatable, intent(out) :: pumas_lamcrad_out(:, :)
    ! microphysics snow mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_snowice_in_prec_out(:, :)
    ! microphysics snow scaled diameter (m)
    real(pumas_r8), allocatable, intent(out) :: pumas_scaled_diam_snow_out(:, :)
    ! microphysics graupel water mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_graupice_in_prec_out(:, :)
    ! microphysics graupel number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_numgraup_vol_in_prec_out(:, :)
    ! microphysics graupel diameter (m)
    real(pumas_r8), allocatable, intent(out) :: pumas_scaled_diam_graup_out(:, :)
    ! microphysics cloud liquid sedimentation flux (kg m-2 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_lflx_out(:, :)
    ! microphysics cloud ice sedimentation flux (kg m-2 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_iflx_out(:, :)
    ! microphysics graupel sedimentation flux (kg m-2 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_gflx_out(:, :)
    ! microphysics rain sedimentation flux (kg m-2 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_rflx_out(:, :)
    ! microphysics snow sedimentation flux (kg m-2 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_sflx_out(:, :)
    ! microphysics rain mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_rainliq_in_prec_out(:, :)
    ! microphysics effective radius of stratiform rain particle (um)
    real(pumas_r8), allocatable, intent(out) :: pumas_reff_rain_out(:, :)
    ! microphysics effective radius of stratiform snow particle (um)
    real(pumas_r8), allocatable, intent(out) :: pumas_reff_snow_out(:, :)
    ! microphysics effective radius of stratiform graupel particle (um)
    real(pumas_r8), allocatable, intent(out) :: pumas_reff_grau_out(:, :)
    ! microphysics rain number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_numrain_vol_in_prec_out(:, :)
    ! microphysics snow number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_numsnow_vol_in_prec_out(:, :)
    ! microphysics analytic radar reflectivity at 94 GHz in precipitating fraction of gridcell (dBZ)
    real(pumas_r8), allocatable, intent(out) :: pumas_refl_out(:, :)
    ! microphysics analytic radar reflectivity at 94 GHz (dBZ)
    real(pumas_r8), allocatable, intent(out) :: pumas_arefl_out(:, :)
    ! microphysics analytic radar reflectivity z factor at 94 GHz (mm6 m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_areflz_out(:, :)
    ! microphysics fraction of gridcell with nonzero radar reflectivity (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_frefl_out(:, :)
    ! microphysics analytic radar reflectivity at 94 GHz with CloudSat thresholds in precipitating fraction of gridcell (dBZ)
    real(pumas_r8), allocatable, intent(out) :: pumas_csrfl_out(:, :)
    ! microphysics analytic radar reflectivity at 94 GHz with CloudSat thresholds (dBZ)
    real(pumas_r8), allocatable, intent(out) :: pumas_acsrfl_out(:, :)
    ! microphysics fraction of gridcell with nonzero radar reflectivity with CloudSat thresholds (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_fcsrfl_out(:, :)
    ! microphysics analytic radar reflectivity at 10 cm wavelength (dBZ)
    real(pumas_r8), allocatable, intent(out) :: pumas_refl10cm_out(:, :)
    ! microphysics analytic radar reflectivity z factor at 10 cm wavelength (mm6 m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_reflz10cm_out(:, :)
    ! microphysics effective radius of stratiform cloud liquid plus rain particles (um)
    real(pumas_r8), allocatable, intent(out) :: pumas_rercld_out(:, :)
    ! microphysics available ice nuclei number concentration of new state (m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_ncai_out(:, :)
    ! microphysics available cloud condensation nuclei number concentration of new state (m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_ncal_out(:, :)
    ! microphysics rain mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_rainliq_out(:, :)
    ! microphysics snow mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_snowice_out(:, :)
    ! microphysics rain number concentration of new state (m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_numrain_vol_out(:, :)
    ! microphysics snow number concentration of new state (m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_numsnow_vol_out(:, :)
    ! microphysics average diameter of stratiform rain particle (m)
    real(pumas_r8), allocatable, intent(out) :: pumas_diam_rain_out(:, :)
    ! microphysics average diameter of stratiform snow particle (m)
    real(pumas_r8), allocatable, intent(out) :: pumas_diam_snow_out(:, :)
    ! microphysics graupel mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_graupice_out(:, :)
    ! microphysics graupel number concentration of new state (m-3)
    real(pumas_r8), allocatable, intent(out) :: pumas_numgraup_vol_out(:, :)
    ! microphysics average diameter of stratiform graupel particle (m)
    real(pumas_r8), allocatable, intent(out) :: pumas_diam_graup_out(:, :)
    ! microphysics fraction of gridcell with graupel (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_freq_graup_out(:, :)
    ! microphysics fraction of gridcell with snow (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_freq_snow_out(:, :)
    ! microphysics fraction of gridcell with rain (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_freq_rain_out(:, :)
    ! microphysics fraction of frozen water to total condensed water (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_frac_ice_out(:, :)
    ! microphysics fraction of cloud liquid tendency applied to state (fraction)
    real(pumas_r8), allocatable, intent(out) :: pumas_frac_cldliq_tend_out(:, :)
    ! microphysics rain evaporation rate wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), allocatable, intent(out) :: pumas_rain_evap_out(:, :)

    !CCPP error handling:
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode

    !Initialize error message and error code:
    errmsg  = ''
    errcode = 0

    !Set PUMAS (cloud microphysics) dimensions:
    micro_ncol   = ncol
    micro_nlev   = nlev-trop_cloud_top_lev+1
    micro_nlevp1 = micro_nlev + 1

    !Allocate PUMAS-grid suite variables:
    allocate(pumas_airT(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_airT = 0._pumas_r8
    allocate(pumas_airq(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_airq = 0._pumas_r8
    allocate(pumas_cldliq(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_cldliq = 0._pumas_r8
    allocate(pumas_cldice(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_cldice = 0._pumas_r8
    allocate(pumas_numliq(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numliq = 0._pumas_r8
    allocate(pumas_numice(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numice = 0._pumas_r8
    allocate(pumas_rainliq(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_rainliq = 0._pumas_r8
    allocate(pumas_snowice(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_snowice = 0._pumas_r8
    allocate(pumas_numrain(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numrain = 0._pumas_r8
    allocate(pumas_numsnow(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numsnow = 0._pumas_r8
    allocate(pumas_graupice(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_graupice = 0._pumas_r8
    allocate(pumas_numgraup(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numgraup = 0._pumas_r8
    allocate(pumas_relvar(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_relvar = 0._pumas_r8
    allocate(pumas_accre_enhan(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_accre_enhan = 0._pumas_r8
    allocate(pumas_pmid(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_pmid = 0._pumas_r8
    allocate(pumas_pdel(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_pdel = 0._pumas_r8
    allocate(pumas_pint(micro_ncol, micro_nlevp1), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_pint = 0._pumas_r8
    allocate(pumas_strat_cldfrc(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_strat_cldfrc = 0._pumas_r8
    allocate(pumas_strat_liq_cldfrc(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_strat_liq_cldfrc = 0._pumas_r8
    allocate(pumas_strat_ice_cldfrc(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_strat_ice_cldfrc = 0._pumas_r8
    allocate(pumas_qsatfac(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_qsatfac = 0._pumas_r8
    allocate(pumas_naai(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_naai = 0._pumas_r8
    allocate(pumas_npccn(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_npccn = 0._pumas_r8
    allocate(pumas_rndst(micro_ncol, micro_nlev, ndust), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_rndst = 0._pumas_r8
    allocate(pumas_nacon(micro_ncol, micro_nlev, ndust), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_nacon = 0._pumas_r8
    allocate(pumas_snowice_tend_external(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_snowice_tend_external = 0._pumas_r8
    allocate(pumas_numsnow_tend_external(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numsnow_tend_external = 0._pumas_r8
    allocate(pumas_effi_external(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_effi_external = 0._pumas_r8
    allocate(pumas_frzimm(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_frzimm = 0._pumas_r8
    allocate(pumas_frzcnt(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_frzcnt = 0._pumas_r8
    allocate(pumas_frzdep(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_frzdep = 0._pumas_r8
    allocate(pumas_qcsinksum_rate1ord_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_qcsinksum_rate1ord_out = 0._pumas_r8
    allocate(pumas_airT_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_airT_tend_out = 0._pumas_r8
    allocate(pumas_airq_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_airq_tend_out = 0._pumas_r8
    allocate(pumas_cldliq_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_cldliq_tend_out = 0._pumas_r8
    allocate(pumas_cldice_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_cldice_tend_out = 0._pumas_r8
    allocate(pumas_numliq_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numliq_tend_out = 0._pumas_r8
    allocate(pumas_numice_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numice_tend_out = 0._pumas_r8
    allocate(pumas_rainliq_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_rainliq_tend_out = 0._pumas_r8
    allocate(pumas_snowice_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_snowice_tend_out = 0._pumas_r8
    allocate(pumas_numrain_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numrain_tend_out = 0._pumas_r8
    allocate(pumas_numsnow_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numsnow_tend_out = 0._pumas_r8
    allocate(pumas_graupice_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_graupice_tend_out = 0._pumas_r8
    allocate(pumas_numgraup_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numgraup_tend_out = 0._pumas_r8
    allocate(pumas_effc_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_effc_out = 0._pumas_r8
    allocate(pumas_effc_fn_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_effc_fn_out = 0._pumas_r8
    allocate(pumas_effi_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_effi_out = 0._pumas_r8
    allocate(pumas_sadice_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_sadice_out = 0._pumas_r8
    allocate(pumas_sadsnow_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_sadsnow_out = 0._pumas_r8
    allocate(pumas_prect_out(micro_ncol), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_prect_out = 0._pumas_r8
    allocate(pumas_preci_out(micro_ncol), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_preci_out = 0._pumas_r8
    allocate(pumas_prec_evap_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_prec_evap_out = 0._pumas_r8
    allocate(pumas_am_evap_st_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_am_evap_st_out = 0._pumas_r8
    allocate(pumas_prec_prod_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_prec_prod_out = 0._pumas_r8
    allocate(pumas_cmeice_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_cmeice_out = 0._pumas_r8
    allocate(pumas_deffi_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_deffi_out = 0._pumas_r8
    allocate(pumas_pgamrad_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_pgamrad_out = 0._pumas_r8
    allocate(pumas_lamcrad_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_lamcrad_out = 0._pumas_r8
    allocate(pumas_snowice_in_prec_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_snowice_in_prec_out = 0._pumas_r8
    allocate(pumas_scaled_diam_snow_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_scaled_diam_snow_out = 0._pumas_r8
    allocate(pumas_graupice_in_prec_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_graupice_in_prec_out = 0._pumas_r8
    allocate(pumas_numgraup_vol_in_prec_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numgraup_vol_in_prec_out = 0._pumas_r8
    allocate(pumas_scaled_diam_graup_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_scaled_diam_graup_out = 0._pumas_r8
    allocate(pumas_lflx_out(micro_ncol, micro_nlevp1), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_lflx_out = 0._pumas_r8
    allocate(pumas_iflx_out(micro_ncol, micro_nlevp1), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_iflx_out = 0._pumas_r8
    allocate(pumas_gflx_out(micro_ncol, micro_nlevp1), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_gflx_out = 0._pumas_r8
    allocate(pumas_rflx_out(micro_ncol, micro_nlevp1), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_rflx_out = 0._pumas_r8
    allocate(pumas_sflx_out(micro_ncol, micro_nlevp1), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_sflx_out = 0._pumas_r8
    allocate(pumas_rainliq_in_prec_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_rainliq_in_prec_out = 0._pumas_r8
    allocate(pumas_reff_rain_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_reff_rain_out = 0._pumas_r8
    allocate(pumas_reff_snow_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_reff_snow_out = 0._pumas_r8
    allocate(pumas_reff_grau_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_reff_grau_out = 0._pumas_r8
    allocate(pumas_numrain_vol_in_prec_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numrain_vol_in_prec_out = 0._pumas_r8
    allocate(pumas_numsnow_vol_in_prec_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numsnow_vol_in_prec_out = 0._pumas_r8
    allocate(pumas_refl_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_refl_out = 0._pumas_r8
    allocate(pumas_arefl_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_arefl_out = 0._pumas_r8
    allocate(pumas_areflz_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_areflz_out = 0._pumas_r8
    allocate(pumas_frefl_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_frefl_out = 0._pumas_r8
    allocate(pumas_csrfl_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_csrfl_out = 0._pumas_r8
    allocate(pumas_acsrfl_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_acsrfl_out = 0._pumas_r8
    allocate(pumas_fcsrfl_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_fcsrfl_out = 0._pumas_r8
    allocate(pumas_refl10cm_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_refl10cm_out = 0._pumas_r8
    allocate(pumas_reflz10cm_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_reflz10cm_out = 0._pumas_r8
    allocate(pumas_rercld_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_rercld_out = 0._pumas_r8
    allocate(pumas_ncai_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_ncai_out = 0._pumas_r8
    allocate(pumas_ncal_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_ncal_out = 0._pumas_r8
    allocate(pumas_rainliq_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_rainliq_out = 0._pumas_r8
    allocate(pumas_snowice_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_snowice_out = 0._pumas_r8
    allocate(pumas_numrain_vol_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numrain_vol_out = 0._pumas_r8
    allocate(pumas_numsnow_vol_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numsnow_vol_out = 0._pumas_r8
    allocate(pumas_diam_rain_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_diam_rain_out = 0._pumas_r8
    allocate(pumas_diam_snow_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_diam_snow_out = 0._pumas_r8
    allocate(pumas_graupice_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_graupice_out = 0._pumas_r8
    allocate(pumas_numgraup_vol_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_numgraup_vol_out = 0._pumas_r8
    allocate(pumas_diam_graup_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_diam_graup_out = 0._pumas_r8
    allocate(pumas_freq_graup_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_freq_graup_out = 0._pumas_r8
    allocate(pumas_freq_snow_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_freq_snow_out = 0._pumas_r8
    allocate(pumas_freq_rain_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_freq_rain_out = 0._pumas_r8
    allocate(pumas_frac_ice_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_frac_ice_out = 0._pumas_r8
    allocate(pumas_frac_cldliq_tend_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_frac_cldliq_tend_out = 0._pumas_r8
    allocate(pumas_rain_evap_out(micro_ncol, micro_nlev), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) return
    pumas_rain_evap_out = 0._pumas_r8

  end subroutine micro_pumas_ccpp_dimensions_pre_init

  !> \section arg_table_micro_pumas_ccpp_dimensions_pre_run Argument Table
  !! \htmlinclude micro_pumas_ccpp_dimensions_pre_run.html
  subroutine micro_pumas_ccpp_dimensions_pre_run(ncol,                &
                             trop_cloud_top_lev,                                     &
                             airT_in, pumas_airT, airq_in, pumas_airq,              &
                             cldliq_in, pumas_cldliq,                               &
                             cldice_in, pumas_cldice,                               &
                             micro_mg_num_steps, dtime, pumas_timestep,             &
                             numliq_in, pumas_numliq,                               &
                             numice_in, pumas_numice,                               &
                             rainliq_in, pumas_rainliq,                             &
                             snowice_in, pumas_snowice,                             &
                             numrain_in, pumas_numrain,                             &
                             numsnow_in, pumas_numsnow,                             &
                             graupice_in, pumas_graupice,                           &
                             numgraup_in, pumas_numgraup,                           &
                             relvar_in, pumas_relvar,                               &
                             accre_enhan_in, pumas_accre_enhan,           &
                             pmid_in, pumas_pmid,                                   &
                             pdel_in, pumas_pdel,                                   &
                             pint_in, pumas_pint,                                   &
                             strat_cldfrc_in, pumas_strat_cldfrc,                   &
                             pumas_strat_liq_cldfrc,                                &
                             pumas_strat_ice_cldfrc,                                &
                             qsatfac_in, pumas_qsatfac,                             &
                             naai_in, pumas_naai,                                   &
                             npccn_in, pumas_npccn,                                 &
                             rndst_in, pumas_rndst,                                 &
                             nacon_in, pumas_nacon,                                 &
                             snowice_tend_external_in, pumas_snowice_tend_external, &
                             numsnow_tend_external_in, pumas_numsnow_tend_external, &
                             effi_external_in, pumas_effi_external,                 &
                             frzimm_in, pumas_frzimm,                               &
                             frzcnt_in, pumas_frzcnt,                               &
                             frzdep_in, pumas_frzdep,                               &
                             errmsg, errcode)

    use pumas_kinds,       only: pumas_r8=>kind_r8
    use ccpp_kinds,        only: kind_phys

    !Host model dimensions/parameters:
    integer,         intent(in) :: ncol
    integer,         intent(in) :: trop_cloud_top_lev  !Index of the top model level for which
                                                       !cloud physics is applied (1 to nlev)

    ! Air temperature (K)
    real(kind_phys), intent(in)  :: airT_in(:, :)
    real(pumas_r8), intent(out) :: pumas_airT(:, :)
    ! Water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: airq_in(:, :)
    real(pumas_r8), intent(out) :: pumas_airq(:, :)
    ! Cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: cldliq_in(:, :)
    real(pumas_r8), intent(out) :: pumas_cldliq(:, :)
    ! Cloud ice mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: cldice_in(:, :)
    real(pumas_r8), intent(out) :: pumas_cldice(:, :)
    ! Substepping in microphysics variables
    integer, intent(in)  :: micro_mg_num_steps
    real(kind_phys), intent(in)  :: dtime
    real(pumas_r8), intent(out) :: pumas_timestep
    ! Mass number concentration of cloud liquid water wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numliq_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numliq(:, :)
    ! Mass number concentration of cloud ice wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numice_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numice(:, :)
    ! Rain mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: rainliq_in(:, :)
    real(pumas_r8), intent(out) :: pumas_rainliq(:, :)
    ! Snow mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: snowice_in(:, :)
    real(pumas_r8), intent(out) :: pumas_snowice(:, :)
    ! Mass number concentration of rain wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numrain_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numrain(:, :)
    ! Mass number concentration of snow wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numsnow_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numsnow(:, :)
    ! Graupel mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: graupice_in(:, :)
    real(pumas_r8), intent(out) :: pumas_graupice(:, :)
    ! Mass number concentration of graupel wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numgraup_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numgraup(:, :)
    ! Relative variance of cloud water (1)
    real(kind_phys), intent(in)  :: relvar_in(:, :)
    real(pumas_r8), intent(out) :: pumas_relvar(:, :)
    ! Accretion enhancement factor (1)
    real(kind_phys), intent(in)  :: accre_enhan_in(:, :)
    real(pumas_r8),  intent(out) :: pumas_accre_enhan(:, :)
    ! Air pressure (Pa)
    real(kind_phys), intent(in)  :: pmid_in(:, :)
    real(pumas_r8), intent(out) :: pumas_pmid(:, :)
    ! Air pressure thickness (Pa)
    real(kind_phys), intent(in)  :: pdel_in(:, :)
    real(pumas_r8), intent(out) :: pumas_pdel(:, :)
    ! Air pressure at interfaces (Pa)
    real(kind_phys), intent(in)  :: pint_in(:, :)
    real(pumas_r8), intent(out) :: pumas_pint(:, :)
    ! Stratiform cloud area fraction (fraction)
    real(kind_phys), intent(in)  :: strat_cldfrc_in(:, :)
    real(pumas_r8), intent(out) :: pumas_strat_cldfrc(:, :)
    ! Stratiform cloud liquid area fraction (fraction)
    real(pumas_r8), intent(out) :: pumas_strat_liq_cldfrc(:, :)
    ! Stratiform cloud ice area fraction (fraction)
    real(pumas_r8), intent(out) :: pumas_strat_ice_cldfrc(:, :)
    ! Subgrid cloud water saturation scaling factor (1)
    real(kind_phys), intent(in)  :: qsatfac_in(:, :)
    real(pumas_r8), intent(out) :: pumas_qsatfac(:, :)
    ! Tendency of activated ice nuclei mass number concentration (kg-1 s-1)
    real(kind_phys), intent(in)  :: naai_in(:, :)
    real(pumas_r8), intent(out) :: pumas_naai(:, :)
    ! Tendency of activated cloud condensation nuclei mass number concentration (kg-1 s-1)
    real(kind_phys), intent(in)  :: npccn_in(:, :)
    real(pumas_r8), intent(out) :: pumas_npccn(:, :)
    ! Dust radii by size bin  (m)
    real(kind_phys), intent(in)  :: rndst_in(:, :, :)
    real(pumas_r8), intent(out) :: pumas_rndst(:, :, :)
    ! Dust number concentration by size bin (m-3)
    real(kind_phys), intent(in)  :: nacon_in(:, :, :)
    real(pumas_r8), intent(out) :: pumas_nacon(:, :, :)
    ! Tendency of snow mixing ratio wrt moist air and condensed water from external microphysics (kg kg-1 s-1)
    real(kind_phys), intent(in)  :: snowice_tend_external_in(:, :)
    real(pumas_r8), intent(out) :: pumas_snowice_tend_external(:, :)
    ! Tendency of mass number concentration of snow wrt moist air and condensed water from external microphysics (kg-1 s-1)
    real(kind_phys), intent(in)  :: numsnow_tend_external_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numsnow_tend_external(:, :)
    ! Effective radius of stratiform cloud ice particle from external microphysics (m)
    real(kind_phys), intent(in)  :: effi_external_in(:, :)
    real(pumas_r8), intent(out) :: pumas_effi_external(:, :)
    ! Tendency of cloud liquid droplet number concentration due to immersion freezing (cm-3)
    real(kind_phys), intent(in)  :: frzimm_in(:, :)
    real(pumas_r8), intent(out) :: pumas_frzimm(:, :)
    ! Tendency of cloud liquid droplet number concentration due to contact freezing (cm-3)
    real(kind_phys), intent(in)  :: frzcnt_in(:, :)
    real(pumas_r8), intent(out) :: pumas_frzcnt(:, :)
    ! Tendency of cloud ice number concentration due to deposition nucleation (cm-3)
    real(kind_phys), intent(in)  :: frzdep_in(:, :)
    real(pumas_r8), intent(out) :: pumas_frzdep(:, :)

    character(len=*), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    !Initialize error message and error code:
    errmsg  = ''
    errcode = 0

    pumas_timestep = dtime/micro_mg_num_steps


! Microphysics runs only over levels trop_cloud_top_lev:nlev (the bottom micro_nlev
! layers of the top-down indexed column). Slice each host input onto that range, so
! pumas_xxx(:ncol,:) receives host levels trop_cloud_top_lev:nlev. Mirrors CAM
! micro_pumas_cam.F90, which passes every field as xxx(:ncol,top_lev:).
    pumas_airT(:ncol,:) = real(airT_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_airq(:ncol,:) = real(airq_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_cldliq(:ncol,:) = real(cldliq_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_cldice(:ncol,:) = real(cldice_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_numliq(:ncol,:) = real(numliq_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_numice(:ncol,:) = real(numice_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_rainliq(:ncol,:) = real(rainliq_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_snowice(:ncol,:) = real(snowice_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_numrain(:ncol,:) = real(numrain_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_numsnow(:ncol,:) = real(numsnow_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_graupice(:ncol,:) = real(graupice_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_numgraup(:ncol,:) = real(numgraup_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_relvar(:ncol,:) = real(relvar_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_accre_enhan(:ncol,:) = real(accre_enhan_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_pmid(:ncol,:) = real(pmid_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_pdel(:ncol,:) = real(pdel_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_pint(:ncol,:) = real(pint_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_strat_cldfrc(:ncol,:) = real(strat_cldfrc_in(:,trop_cloud_top_lev:), pumas_r8)
    ! PUMAS uses total stratiform fraction for both liquid and ice stratus.
    ! Mirrors CAM micro_pumas_cam.F90: alst_mic => ast; aist_mic => ast.
    pumas_strat_liq_cldfrc(:ncol,:) = real(strat_cldfrc_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_strat_ice_cldfrc(:ncol,:) = real(strat_cldfrc_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_qsatfac(:ncol,:) = real(qsatfac_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_naai(:ncol,:) = real(naai_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_npccn(:ncol,:) = real(npccn_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_rndst(:ncol,:,:) = real(rndst_in(:,trop_cloud_top_lev:,:), pumas_r8)
    pumas_nacon(:ncol,:,:) = real(nacon_in(:,trop_cloud_top_lev:,:), pumas_r8)
    pumas_snowice_tend_external(:ncol,:) = real(snowice_tend_external_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_numsnow_tend_external(:ncol,:) = real(numsnow_tend_external_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_effi_external(:ncol,:) = real(effi_external_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_frzcnt(:ncol,:) = real(frzcnt_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_frzdep(:ncol,:) = real(frzdep_in(:,trop_cloud_top_lev:), pumas_r8)
    pumas_frzimm(:ncol,:) = real(frzimm_in(:,trop_cloud_top_lev:), pumas_r8)

  end subroutine micro_pumas_ccpp_dimensions_pre_run

end module micro_pumas_ccpp_dimensions_pre
