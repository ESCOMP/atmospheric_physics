module micro_pumas_ccpp_dimensions_post

  implicit none
  private

  public :: micro_pumas_ccpp_dimensions_post_run

! This will be replaced by subcolumn averaging when subcolumns are enabled

contains
  !> \section arg_table_micro_pumas_ccpp_dimensions_post_run Argument Table
  !! \htmlinclude micro_pumas_ccpp_dimensions_post_run.html
  subroutine micro_pumas_ccpp_dimensions_post_run(ncol, &
                   trop_cloud_top_lev, qcsinksum_rate1ord, pumas_qcsinksum_rate1ord, airT_tend,       &
                   pumas_airT_tend, airq_tend, pumas_airq_tend,                                         &
                   cldliq_tend, pumas_cldliq_tend, cldice_tend, pumas_cldice_tend, numliq_tend,         &
                   pumas_numliq_tend, numice_tend, pumas_numice_tend, rainliq_tend, pumas_rainliq_tend, &
                   snowice_tend, pumas_snowice_tend, numrain_tend, pumas_numrain_tend, numsnow_tend,    &
                   pumas_numsnow_tend, graupice_tend, pumas_graupice_tend, numgraup_tend,               &
                   pumas_numgraup_tend, effc, pumas_effc, effc_fn, pumas_effc_fn, effi, pumas_effi,     &
                   sadice, pumas_sadice, sadsnow, pumas_sadsnow, prect, pumas_prect, preci, pumas_preci,&
                   prec_evap, pumas_prec_evap, am_evap_st, pumas_am_evap_st, prec_prod, pumas_prec_prod,&
                   cmeice, pumas_cmeice, deffi, pumas_deffi, pgamrad, pumas_pgamrad, lamcrad,           &
                   pumas_lamcrad, snowice_in_prec, pumas_snowice_in_prec, scaled_diam_snow,             &
                   pumas_scaled_diam_snow, graupice_in_prec, pumas_graupice_in_prec,                    &
                   numgraup_vol_in_prec, pumas_numgraup_vol_in_prec, scaled_diam_graup,                 &
                   pumas_scaled_diam_graup, lflx, pumas_lflx, iflx, pumas_iflx, gflx, pumas_gflx, rflx, &
                   pumas_rflx, sflx, pumas_sflx, rainliq_in_prec, pumas_rainliq_in_prec, reff_rain,     &
                   pumas_reff_rain, reff_snow, pumas_reff_snow, reff_grau, pumas_reff_grau,             &
                   numrain_vol_in_prec, pumas_numrain_vol_in_prec, numsnow_vol_in_prec,                 &
                   pumas_numsnow_vol_in_prec, refl, pumas_refl, arefl, pumas_arefl, areflz, pumas_areflz,&
                   frefl, pumas_frefl, csrfl, pumas_csrfl, acsrfl, pumas_acsrfl, fcsrfl, pumas_fcsrfl,  &
                   refl10cm, pumas_refl10cm, reflz10cm, pumas_reflz10cm, rercld, pumas_rercld, ncai,    &
                   pumas_ncai, ncal, pumas_ncal, rainliq, pumas_rainliq, snowice, pumas_snowice,        &
                   numrain_vol, pumas_numrain_vol, numsnow_vol, pumas_numsnow_vol, diam_rain,           &
                   pumas_diam_rain, diam_snow, pumas_diam_snow, graupice, pumas_graupice, numgraup_vol, &
                   pumas_numgraup_vol, diam_graup, pumas_diam_graup, freq_graup, pumas_freq_graup,      &
                   freq_snow, pumas_freq_snow, freq_rain, pumas_freq_rain, frac_ice, pumas_frac_ice,    &
                   frac_cldliq_tend, pumas_frac_cldliq_tend, rain_evap, pumas_rain_evap,                &
                   errmsg, errcode)

    use pumas_kinds,       only: pumas_r8=>kind_r8
    use ccpp_kinds,        only: kind_phys

    ! horizontal dimension
    integer, intent(in) :: ncol
    ! vertical layer dimension
    ! index of the top model level for which cloud physics is applied (1 to nlev)
    integer, intent(in) :: trop_cloud_top_lev

    !microphysics direct conversion rate of stratiform cloud water to precipitation (s-1)
    real(pumas_r8), intent(in) :: pumas_qcsinksum_rate1ord(:, :)
    !microphysics tendency of dry air enthalpy at constant pressure (J kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_airT_tend(:, :)
    !microphysics tendency of water vapor mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_airq_tend(:, :)
    !microphysics tendency of cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_cldliq_tend(:, :)
    !microphysics tendency of cloud ice mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_cldice_tend(:, :)
    !microphysics tendency of mass number concentration of cloud liquid water wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_numliq_tend(:, :)
    !microphysics tendency of mass number concentration of cloud ice wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_numice_tend(:, :)
    !microphysics tendency of rain mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_rainliq_tend(:, :)
    !microphysics tendency of snow mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_snowice_tend(:, :)
    !microphysics tendency of mass number concentration of rain wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_numrain_tend(:, :)
    !microphysics tendency of mass number concentration of snow wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_numsnow_tend(:, :)
    !microphysics tendency of graupel mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_graupice_tend(:, :)
    !microphysics tendency of mass number concentration of graupel wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_numgraup_tend(:, :)
    !microphysics effective radius of stratiform cloud liquid water particle (um)
    real(pumas_r8), intent(in) :: pumas_effc(:, :)
    !microphysics effective radius of stratiform cloud liquid water particle assuming droplet number concentration of 1e8 kg-1 (um)
    real(pumas_r8), intent(in) :: pumas_effc_fn(:, :)
    !microphysics effective radius of stratiform cloud ice particle (um)
    real(pumas_r8), intent(in) :: pumas_effi(:, :)
    !microphysics cloud ice surface area density (cm2 cm-3)
    real(pumas_r8), intent(in) :: pumas_sadice(:, :)
    !microphysics snow surface area density (cm2 cm-3)
    real(pumas_r8), intent(in) :: pumas_sadsnow(:, :)
    !microphysics LWE large scale precipitation rate at surface (m s-1)
    real(pumas_r8), intent(in) :: pumas_prect(:)
    !microphysics LWE large scale snowfall rate at surface (m s-1)
    real(pumas_r8), intent(in) :: pumas_preci(:)
    !microphysics precipitation evaporation rate wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_prec_evap(:, :)
    !microphysics precipitation evaporation area (fraction)
    real(pumas_r8), intent(in) :: pumas_am_evap_st(:, :)
    !microphysics precipitation production rate wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_prec_prod(:, :)
    !microphysics condensation minus evaporation rate of in-cloud ice wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_cmeice(:, :)
    !microphysics effective diameter of stratiform cloud ice particles for radiation (um)
    real(pumas_r8), intent(in) :: pumas_deffi(:, :)
    !microphysics cloud particle size distribution shape parameter (1)
    real(pumas_r8), intent(in) :: pumas_pgamrad(:, :)
    !microphysics cloud particle size distribution slope parameter (1)
    real(pumas_r8), intent(in) :: pumas_lamcrad(:, :)
    !microphysics snow mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_snowice_in_prec(:, :)
    !microphysics snow scaled diameter (m)
    real(pumas_r8), intent(in) :: pumas_scaled_diam_snow(:, :)
    !microphysics graupel mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_graupice_in_prec(:, :)
    !microphysics graupel number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), intent(in) :: pumas_numgraup_vol_in_prec(:, :)
    !microphysics graupel scaled diameter (m)
    real(pumas_r8), intent(in) :: pumas_scaled_diam_graup(:, :)
    !microphysics cloud liquid sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(in) :: pumas_lflx(:, :)
    !microphysics cloud ice sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(in) :: pumas_iflx(:, :)
    !microphysics graupel sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(in) :: pumas_gflx(:, :)
    !microphysics rain sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(in) :: pumas_rflx(:, :)
    !microphysics snow sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(in) :: pumas_sflx(:, :)
    !microphysics rain mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_rainliq_in_prec(:, :)
    !microphysics effective radius of stratiform rain particle (um)
    real(pumas_r8), intent(in) :: pumas_reff_rain(:, :)
    !microphysics effective radius of stratiform snow particle (um)
    real(pumas_r8), intent(in) :: pumas_reff_snow(:, :)
    !microphysics effective radius of stratiform graupel particle (um)
    real(pumas_r8), intent(in) :: pumas_reff_grau(:, :)
    !microphysics rain number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), intent(in) :: pumas_numrain_vol_in_prec(:, :)
    !microphysics snow number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), intent(in) :: pumas_numsnow_vol_in_prec(:, :)
    !microphysics analytic radar reflectivity at 94 GHz in precipitating fraction of gridcell (dBZ)
    real(pumas_r8), intent(in) :: pumas_refl(:, :)
    !microphysics analytic radar reflectivity at 94 GHz (dBZ)
    real(pumas_r8), intent(in) :: pumas_arefl(:, :)
    !microphysics analytic radar reflectivity z factor at 94 GHz (mm6 m-3)
    real(pumas_r8), intent(in) :: pumas_areflz(:, :)
    !microphysics fraction of gridcell with nonzero radar reflectivity (fraction)
    real(pumas_r8), intent(in) :: pumas_frefl(:, :)
    !microphysics analytic radar reflectivity at 94 GHz with CloudSat thresholds in precipitating fraction of gridcell (dBZ)
    real(pumas_r8), intent(in) :: pumas_csrfl(:, :)
    !microphysics analytic radar reflectivity at 94 GHz with CloudSat thresholds (dBZ)
    real(pumas_r8), intent(in) :: pumas_acsrfl(:, :)
    !microphysics fraction of gridcell with nonzero radar reflectivity with CloudSat thresholds (fraction)
    real(pumas_r8), intent(in) :: pumas_fcsrfl(:, :)
    !microphysics analytic radar reflectivity at 10 cm wavelength (dBZ)
    real(pumas_r8), intent(in) :: pumas_refl10cm(:, :)
    !microphysics analytic radar reflectivity z factor at 10 cm wavelength (mm6 m-3)
    real(pumas_r8), intent(in) :: pumas_reflz10cm(:, :)
    !microphysics effective radius of stratiform cloud liquid plus rain particles (m)
    real(pumas_r8), intent(in) :: pumas_rercld(:, :)
    !microphysics available ice nuclei number concentration of new state (m-3)
    real(pumas_r8), intent(in) :: pumas_ncai(:, :)
    !microphysics available cloud condensation nuclei number concentration of new state (m-3)
    real(pumas_r8), intent(in) :: pumas_ncal(:, :)
    !microphysics rain mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_rainliq(:, :)
    !microphysics snow mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_snowice(:, :)
    !microphysics rain number concentration of new state (m-3)
    real(pumas_r8), intent(in) :: pumas_numrain_vol(:, :)
    !microphysics snow number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), intent(in) :: pumas_numsnow_vol(:, :)
    !microphysics average diameter of stratiform rain particle (m)
    real(pumas_r8), intent(in) :: pumas_diam_rain(:, :)
    !microphysics average diameter of stratiform snow particle (m)
    real(pumas_r8), intent(in) :: pumas_diam_snow(:, :)
    !microphysics graupel mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_graupice(:, :)
    !microphysics graupel number concentration of new state (m-3)
    real(pumas_r8), intent(in) :: pumas_numgraup_vol(:, :)
    !microphysics average diameter of stratiform graupel particle (m)
    real(pumas_r8), intent(in) :: pumas_diam_graup(:, :)
    !microphysics fraction of gridcell with graupel (fraction)
    real(pumas_r8), intent(in) :: pumas_freq_graup(:, :)
    !microphysics fraction of gridcell with snow (fraction)
    real(pumas_r8), intent(in) :: pumas_freq_snow(:, :)
    !microphysics fraction of gridcell with rain (fraction)
    real(pumas_r8), intent(in) :: pumas_freq_rain(:, :)
    !microphysics fraction of frozen water to total condensed water (fraction)
    real(pumas_r8), intent(in) :: pumas_frac_ice(:, :)
    !microphysics fraction of cloud liquid tendency applied to state (fraction)
    real(pumas_r8), intent(in) :: pumas_frac_cldliq_tend(:, :)
    !microphysics rain evaporation rate wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_rain_evap(:, :)

    !direct conversion rate of stratiform cloud water to precipitation (s-1)
    real(kind_phys), intent(out) :: qcsinksum_rate1ord(:, :)
    !tendency of dry air enthalpy at constant pressure (J kg-1 s-1)
    real(kind_phys), intent(out) :: airT_tend(:, :)
    !tendency of water vapor mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: airq_tend(:, :)
    !tendency of cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: cldliq_tend(:, :)
    !tendency of cloud ice mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: cldice_tend(:, :)
    !tendency of mass number concentration of cloud liquid water wrt moist air and condensed water (kg-1 s-1)
    real(kind_phys), intent(out) :: numliq_tend(:, :)
    !tendency of mass number concentration of cloud ice wrt moist air and condensed water (kg-1 s-1)
    real(kind_phys), intent(out) :: numice_tend(:, :)
    !tendency of rain mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: rainliq_tend(:, :)
    !tendency of snow mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: snowice_tend(:, :)
    !tendency of mass number concentration of rain wrt moist air and condensed water (kg-1 s-1)
    real(kind_phys), intent(out) :: numrain_tend(:, :)
    !tendency of mass number concentration of snow wrt moist air and condensed water (kg-1 s-1)
    real(kind_phys), intent(out) :: numsnow_tend(:, :)
    !tendency of graupel mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: graupice_tend(:, :)
    !tendency of mass number concentration of graupel wrt moist air and condensed water (kg-1 s-1)
    real(kind_phys), intent(out) :: numgraup_tend(:, :)
    !effective radius of stratiform cloud liquid water particle (um)
    real(kind_phys), intent(out) :: effc(:, :)
    !effective radius of stratiform cloud liquid water particle assuming droplet number concentration of 1e8 kg-1 (um)
    real(kind_phys), intent(out) :: effc_fn(:, :)
    !effective radius of stratiform cloud ice particle (um)
    real(kind_phys), intent(out) :: effi(:, :)
    !cloud ice surface area density (cm2 cm-3)
    real(kind_phys), intent(out) :: sadice(:, :)
    !snow surface area density (cm2 cm-3)
    real(kind_phys), intent(out) :: sadsnow(:, :)
    !LWE large scale precipitation rate at surface (m s-1)
    real(kind_phys), intent(out) :: prect(:)
    !LWE large scale snowfall rate at surface (m s-1)
    real(kind_phys), intent(out) :: preci(:)
    !precipitation evaporation rate wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: prec_evap(:, :)
    !precipitation evaporation area (fraction)
    real(kind_phys), intent(out) :: am_evap_st(:, :)
    !precipitation production rate wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: prec_prod(:, :)
    !condensation minus evaporation rate of in-cloud ice wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: cmeice(:, :)
    !effective diameter of stratiform cloud ice particles for radiation (um)
    real(kind_phys), intent(out) :: deffi(:, :)
    !cloud particle size distribution shape parameter (1)
    real(kind_phys), intent(out) :: pgamrad(:, :)
    !cloud particle size distribution slope parameter (1)
    real(kind_phys), intent(out) :: lamcrad(:, :)
    !snow mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(kind_phys), intent(out) :: snowice_in_prec(:, :)
    !snow scaled diameter (m)
    real(kind_phys), intent(out) :: scaled_diam_snow(:, :)
    !graupel mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(kind_phys), intent(out) :: graupice_in_prec(:, :)
    !graupel number concentration of new state in precipitating fraction of gridcell (m-3)
    real(kind_phys), intent(out) :: numgraup_vol_in_prec(:, :)
    !graupel scaled diameter (m)
    real(kind_phys), intent(out) :: scaled_diam_graup(:, :)
    !cloud liquid sedimentation flux (kg m-2 s-1)
    real(kind_phys), intent(out) :: lflx(:, :)
    !cloud ice sedimentation flux (kg m-2 s-1)
    real(kind_phys), intent(out) :: iflx(:, :)
    !graupel sedimentation flux (kg m-2 s-1)
    real(kind_phys), intent(out) :: gflx(:, :)
    !rain sedimentation flux (kg m-2 s-1)
    real(kind_phys), intent(out) :: rflx(:, :)
    !snow sedimentation flux (kg m-2 s-1)
    real(kind_phys), intent(out) :: sflx(:, :)
    !rain mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(kind_phys), intent(out) :: rainliq_in_prec(:, :)
    !effective radius of stratiform rain particle (um)
    real(kind_phys), intent(out) :: reff_rain(:, :)
    !effective radius of stratiform snow particle (um)
    real(kind_phys), intent(out) :: reff_snow(:, :)
    !effective radius of stratiform graupel particle (um)
    real(kind_phys), intent(out) :: reff_grau(:, :)
    !rain number concentration of new state in precipitating fraction of gridcell (m-3)
    real(kind_phys), intent(out) :: numrain_vol_in_prec(:, :)
    !snow number concentration of new state in precipitating fraction of gridcell (m-3)
    real(kind_phys), intent(out) :: numsnow_vol_in_prec(:, :)
    !analytic radar reflectivity at 94 GHz in precipitating fraction of gridcell (dBZ)
    real(kind_phys), intent(out) :: refl(:, :)
    !analytic radar reflectivity at 94 GHz (dBZ)
    real(kind_phys), intent(out) :: arefl(:, :)
    !analytic radar reflectivity z factor at 94 GHz (mm6 m-3)
    real(kind_phys), intent(out) :: areflz(:, :)
    !fraction of gridcell with nonzero radar reflectivity (fraction)
    real(kind_phys), intent(out) :: frefl(:, :)
    !analytic radar reflectivity at 94 GHz with CloudSat thresholds in precipitating fraction of gridcell (dBZ)
    real(kind_phys), intent(out) :: csrfl(:, :)
    !analytic radar reflectivity at 94 GHz with CloudSat thresholds (dBZ)
    real(kind_phys), intent(out) :: acsrfl(:, :)
    !fraction of gridcell with nonzero radar reflectivity with CloudSat thresholds (fraction)
    real(kind_phys), intent(out) :: fcsrfl(:, :)
    !analytic radar reflectivity at 10 cm wavelength (dBZ)
    real(kind_phys), intent(out) :: refl10cm(:, :)
    !analytic radar reflectivity z factor at 10 cm wavelength (mm6 m-3)
    real(kind_phys), intent(out) :: reflz10cm(:, :)
    !effective radius of stratiform cloud liquid plus rain particles (m)
    real(kind_phys), intent(out) :: rercld(:, :)
    !available ice nuclei number concentration of new state (m-3)
    real(kind_phys), intent(out) :: ncai(:, :)
    !available cloud condensation nuclei number concentration of new state (m-3)
    real(kind_phys), intent(out) :: ncal(:, :)
    !rain mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(kind_phys), intent(out) :: rainliq(:, :)
    !snow mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(kind_phys), intent(out) :: snowice(:, :)
    !rain number concentration of new state (m-3)
    real(kind_phys), intent(out) :: numrain_vol(:, :)
    !snow number concentration of new state in precipitating fraction of gridcell (m-3)
    real(kind_phys), intent(out) :: numsnow_vol(:, :)
    !average diameter of stratiform rain particle (m)
    real(kind_phys), intent(out) :: diam_rain(:, :)
    !average diameter of stratiform snow particle (m)
    real(kind_phys), intent(out) :: diam_snow(:, :)
    !graupel mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(kind_phys), intent(out) :: graupice(:, :)
    !graupel number concentration of new state (m-3)
    real(kind_phys), intent(out) :: numgraup_vol(:, :)
    !average diameter of stratiform graupel particle (m)
    real(kind_phys), intent(out) :: diam_graup(:, :)
    !fraction of gridcell with graupel (fraction)
    real(kind_phys), intent(out) :: freq_graup(:, :)
    !fraction of gridcell with snow (fraction)
    real(kind_phys), intent(out) :: freq_snow(:, :)
    !fraction of gridcell with rain (fraction)
    real(kind_phys), intent(out) :: freq_rain(:, :)
    !fraction of frozen water to total condensed water (fraction)
    real(kind_phys), intent(out) :: frac_ice(:, :)
    !fraction of cloud liquid tendency applied to state (fraction)
    real(kind_phys), intent(out) :: frac_cldliq_tend(:, :)
    !rain evaporation rate wrt moist air and condensed water (kg kg-1 s-1)
    real(kind_phys), intent(out) :: rain_evap(:, :)

    character(len=*), intent(out) :: errmsg  !PUMAS/CCPP error message (none)
    integer,            intent(out) :: errcode !CCPP error code (1)

    ! Set error variables
    errcode = 0
    errmsg = ''

    ! Copy PUMAS outputs (levels 1:micro_nlev) back onto host levels trop_cloud_top_lev:nlev.
    qcsinksum_rate1ord(:,trop_cloud_top_lev:)   = real(pumas_qcsinksum_rate1ord(:ncol,:), kind_phys)
    airT_tend(:,trop_cloud_top_lev:)            = real(pumas_airT_tend(:ncol,:), kind_phys)
    airq_tend(:,trop_cloud_top_lev:)            = real(pumas_airq_tend(:ncol,:), kind_phys)
    cldliq_tend(:,trop_cloud_top_lev:)          = real(pumas_cldliq_tend(:ncol,:), kind_phys)
    cldice_tend(:,trop_cloud_top_lev:)          = real(pumas_cldice_tend(:ncol,:), kind_phys)
    numliq_tend(:,trop_cloud_top_lev:)          = real(pumas_numliq_tend(:ncol,:), kind_phys)
    numice_tend(:,trop_cloud_top_lev:)          = real(pumas_numice_tend(:ncol,:), kind_phys)
    rainliq_tend(:,trop_cloud_top_lev:)         = real(pumas_rainliq_tend(:ncol,:), kind_phys)
    snowice_tend(:,trop_cloud_top_lev:)         = real(pumas_snowice_tend(:ncol,:), kind_phys)
    numrain_tend(:,trop_cloud_top_lev:)         = real(pumas_numrain_tend(:ncol,:), kind_phys)
    numsnow_tend(:,trop_cloud_top_lev:)         = real(pumas_numsnow_tend(:ncol,:), kind_phys)
    graupice_tend(:,trop_cloud_top_lev:)        = real(pumas_graupice_tend(:ncol,:), kind_phys)
    numgraup_tend(:,trop_cloud_top_lev:)        = real(pumas_numgraup_tend(:ncol,:), kind_phys)
    effc(:,trop_cloud_top_lev:)                 = real(pumas_effc(:ncol,:), kind_phys)
    effc_fn(:,trop_cloud_top_lev:)              = real(pumas_effc_fn(:ncol,:), kind_phys)
    effi(:,trop_cloud_top_lev:)                 = real(pumas_effi(:ncol,:), kind_phys)
    sadice(:,trop_cloud_top_lev:)               = real(pumas_sadice(:ncol,:), kind_phys)
    sadsnow(:,trop_cloud_top_lev:)              = real(pumas_sadsnow(:ncol,:), kind_phys)
    prect(:)                                    = real(pumas_prect(:ncol), kind_phys)
    preci(:)                                    = real(pumas_preci(:ncol), kind_phys)
    prec_evap(:,trop_cloud_top_lev:)            = real(pumas_prec_evap(:ncol,:), kind_phys)
    am_evap_st(:,trop_cloud_top_lev:)           = real(pumas_am_evap_st(:ncol,:), kind_phys)
    prec_prod(:,trop_cloud_top_lev:)            = real(pumas_prec_prod(:ncol,:), kind_phys)
    cmeice(:,trop_cloud_top_lev:)               = real(pumas_cmeice(:ncol,:), kind_phys)
    deffi(:,trop_cloud_top_lev:)                = real(pumas_deffi(:ncol,:), kind_phys)
    pgamrad(:,trop_cloud_top_lev:)              = real(pumas_pgamrad(:ncol,:), kind_phys)
    lamcrad(:,trop_cloud_top_lev:)              = real(pumas_lamcrad(:ncol,:), kind_phys)
    snowice_in_prec(:,trop_cloud_top_lev:)      = real(pumas_snowice_in_prec(:ncol,:), kind_phys)
    scaled_diam_snow(:,trop_cloud_top_lev:)     = real(pumas_scaled_diam_snow(:ncol,:), kind_phys)*1.0e6_kind_phys
    graupice_in_prec(:,trop_cloud_top_lev:)     = real(pumas_graupice_in_prec(:ncol,:), kind_phys)
    numgraup_vol_in_prec(:,trop_cloud_top_lev:) = real(pumas_numgraup_vol_in_prec(:ncol,:), kind_phys)
    scaled_diam_graup(:,trop_cloud_top_lev:)    = real(pumas_scaled_diam_graup(:ncol,:), kind_phys)
    lflx(:,trop_cloud_top_lev:)                 = real(pumas_lflx(:ncol,:), kind_phys)
    iflx(:,trop_cloud_top_lev:)                 = real(pumas_iflx(:ncol,:), kind_phys)
    gflx(:,trop_cloud_top_lev:)                 = real(pumas_gflx(:ncol,:), kind_phys)
    rflx(:,trop_cloud_top_lev:)                 = real(pumas_rflx(:ncol,:), kind_phys)
    sflx(:,trop_cloud_top_lev:)                 = real(pumas_sflx(:ncol,:), kind_phys)
    rainliq_in_prec(:,trop_cloud_top_lev:)      = real(pumas_rainliq_in_prec(:ncol,:), kind_phys)
    reff_rain(:,trop_cloud_top_lev:)            = real(pumas_reff_rain(:ncol,:), kind_phys)
    reff_snow(:,trop_cloud_top_lev:)            = real(pumas_reff_snow(:ncol,:), kind_phys)
    reff_grau(:,trop_cloud_top_lev:)            = real(pumas_reff_grau(:ncol,:), kind_phys)
    numrain_vol_in_prec(:,trop_cloud_top_lev:)  = real(pumas_numrain_vol_in_prec(:ncol,:), kind_phys)
    numsnow_vol_in_prec(:,trop_cloud_top_lev:)  = real(pumas_numsnow_vol_in_prec(:ncol,:), kind_phys)
    refl(:,trop_cloud_top_lev:)                 = real(pumas_refl(:ncol,:), kind_phys)
    arefl(:,trop_cloud_top_lev:)                = real(pumas_arefl(:ncol,:), kind_phys)
    areflz(:,trop_cloud_top_lev:)               = real(pumas_areflz(:ncol,:), kind_phys)
    frefl(:,trop_cloud_top_lev:)                = real(pumas_frefl(:ncol,:), kind_phys)
    csrfl(:,trop_cloud_top_lev:)                = real(pumas_csrfl(:ncol,:), kind_phys)
    acsrfl(:,trop_cloud_top_lev:)               = real(pumas_acsrfl(:ncol,:), kind_phys)
    fcsrfl(:,trop_cloud_top_lev:)               = real(pumas_fcsrfl(:ncol,:), kind_phys)
    refl10cm(:,trop_cloud_top_lev:)             = real(pumas_refl10cm(:ncol,:), kind_phys)
    reflz10cm(:,trop_cloud_top_lev:)            = real(pumas_reflz10cm(:ncol,:), kind_phys)
    rercld(:,trop_cloud_top_lev:)               = real(pumas_rercld(:ncol,:), kind_phys)
    ncai(:,trop_cloud_top_lev:)                 = real(pumas_ncai(:ncol,:), kind_phys)
    ncal(:,trop_cloud_top_lev:)                 = real(pumas_ncal(:ncol,:), kind_phys)
    rainliq(:,trop_cloud_top_lev:)              = real(pumas_rainliq(:ncol,:), kind_phys)
    snowice(:,trop_cloud_top_lev:)              = real(pumas_snowice(:ncol,:), kind_phys)
    numrain_vol(:,trop_cloud_top_lev:)          = real(pumas_numrain_vol(:ncol,:), kind_phys)
    numsnow_vol(:,trop_cloud_top_lev:)          = real(pumas_numsnow_vol(:ncol,:), kind_phys)
    diam_rain(:,trop_cloud_top_lev:)            = real(pumas_diam_rain(:ncol,:), kind_phys)
    diam_snow(:,trop_cloud_top_lev:)            = real(pumas_diam_snow(:ncol,:), kind_phys)
    graupice(:,trop_cloud_top_lev:)             = real(pumas_graupice(:ncol,:), kind_phys)
    numgraup_vol(:,trop_cloud_top_lev:)         = real(pumas_numgraup_vol(:ncol,:), kind_phys)
    diam_graup(:,trop_cloud_top_lev:)           = real(pumas_diam_graup(:ncol,:), kind_phys)
    freq_graup(:,trop_cloud_top_lev:)           = real(pumas_freq_graup(:ncol,:), kind_phys)
    freq_snow(:,trop_cloud_top_lev:)            = real(pumas_freq_snow(:ncol,:), kind_phys)
    freq_rain(:,trop_cloud_top_lev:)            = real(pumas_freq_rain(:ncol,:), kind_phys)
    frac_ice(:,trop_cloud_top_lev:)             = real(pumas_frac_ice(:ncol,:), kind_phys)
    frac_cldliq_tend(:,trop_cloud_top_lev:)     = real(pumas_frac_cldliq_tend(:ncol,:), kind_phys)
    rain_evap(:,trop_cloud_top_lev:)            = real(pumas_rain_evap(:ncol,:), kind_phys)

    ! Zero the levels above the microphysics range (1:trop_cloud_top_lev-1)
    qcsinksum_rate1ord(:,:trop_cloud_top_lev-1)   = 0._kind_phys
    airT_tend(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    airq_tend(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    cldliq_tend(:,:trop_cloud_top_lev-1)          = 0._kind_phys
    cldice_tend(:,:trop_cloud_top_lev-1)          = 0._kind_phys
    numliq_tend(:,:trop_cloud_top_lev-1)          = 0._kind_phys
    numice_tend(:,:trop_cloud_top_lev-1)          = 0._kind_phys
    rainliq_tend(:,:trop_cloud_top_lev-1)         = 0._kind_phys
    snowice_tend(:,:trop_cloud_top_lev-1)         = 0._kind_phys
    numrain_tend(:,:trop_cloud_top_lev-1)         = 0._kind_phys
    numsnow_tend(:,:trop_cloud_top_lev-1)         = 0._kind_phys
    graupice_tend(:,:trop_cloud_top_lev-1)        = 0._kind_phys
    numgraup_tend(:,:trop_cloud_top_lev-1)        = 0._kind_phys
    effc(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    effc_fn(:,:trop_cloud_top_lev-1)              = 0._kind_phys
    effi(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    sadice(:,:trop_cloud_top_lev-1)               = 0._kind_phys
    sadsnow(:,:trop_cloud_top_lev-1)              = 0._kind_phys
    prec_evap(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    am_evap_st(:,:trop_cloud_top_lev-1)           = 0._kind_phys
    prec_prod(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    cmeice(:,:trop_cloud_top_lev-1)               = 0._kind_phys
    deffi(:,:trop_cloud_top_lev-1)                = 0._kind_phys
    pgamrad(:,:trop_cloud_top_lev-1)              = 0._kind_phys
    lamcrad(:,:trop_cloud_top_lev-1)              = 0._kind_phys
    snowice_in_prec(:,:trop_cloud_top_lev-1)      = 0._kind_phys
    scaled_diam_snow(:,:trop_cloud_top_lev-1)     = 0._kind_phys
    graupice_in_prec(:,:trop_cloud_top_lev-1)     = 0._kind_phys
    numgraup_vol_in_prec(:,:trop_cloud_top_lev-1) = 0._kind_phys
    scaled_diam_graup(:,:trop_cloud_top_lev-1)    = 0._kind_phys
    lflx(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    iflx(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    gflx(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    rflx(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    sflx(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    rainliq_in_prec(:,:trop_cloud_top_lev-1)      = 0._kind_phys
    reff_rain(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    reff_snow(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    reff_grau(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    numrain_vol_in_prec(:,:trop_cloud_top_lev-1)  = 0._kind_phys
    numsnow_vol_in_prec(:,:trop_cloud_top_lev-1)  = 0._kind_phys
    refl(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    arefl(:,:trop_cloud_top_lev-1)                = 0._kind_phys
    areflz(:,:trop_cloud_top_lev-1)               = 0._kind_phys
    frefl(:,:trop_cloud_top_lev-1)                = 0._kind_phys
    csrfl(:,:trop_cloud_top_lev-1)                = 0._kind_phys
    acsrfl(:,:trop_cloud_top_lev-1)               = 0._kind_phys
    fcsrfl(:,:trop_cloud_top_lev-1)               = 0._kind_phys
    refl10cm(:,:trop_cloud_top_lev-1)             = 0._kind_phys
    reflz10cm(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    rercld(:,:trop_cloud_top_lev-1)               = 0._kind_phys
    ncai(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    ncal(:,:trop_cloud_top_lev-1)                 = 0._kind_phys
    rainliq(:,:trop_cloud_top_lev-1)              = 0._kind_phys
    snowice(:,:trop_cloud_top_lev-1)              = 0._kind_phys
    numrain_vol(:,:trop_cloud_top_lev-1)          = 0._kind_phys
    numsnow_vol(:,:trop_cloud_top_lev-1)          = 0._kind_phys
    diam_rain(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    diam_snow(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    graupice(:,:trop_cloud_top_lev-1)             = 0._kind_phys
    numgraup_vol(:,:trop_cloud_top_lev-1)         = 0._kind_phys
    diam_graup(:,:trop_cloud_top_lev-1)           = 0._kind_phys
    freq_graup(:,:trop_cloud_top_lev-1)           = 0._kind_phys
    freq_snow(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    freq_rain(:,:trop_cloud_top_lev-1)            = 0._kind_phys
    frac_ice(:,:trop_cloud_top_lev-1)             = 0._kind_phys
    frac_cldliq_tend(:,:trop_cloud_top_lev-1)     = 0._kind_phys
    rain_evap(:,:trop_cloud_top_lev-1)            = 0._kind_phys

  end subroutine micro_pumas_ccpp_dimensions_post_run

end module micro_pumas_ccpp_dimensions_post
