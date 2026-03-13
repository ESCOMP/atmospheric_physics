!Common Community Physics Package (CCPP) wrapper for PUMAS.
module micro_pumas_ccpp

  implicit none
  private
  save

  public :: micro_pumas_ccpp_init
  public :: micro_pumas_ccpp_run

contains

  !> \section arg_table_micro_pumas_ccpp_init Argument Table
  !! \htmlinclude micro_pumas_ccpp_init.html
  subroutine micro_pumas_ccpp_init(micro_ncol, micro_nlev,                               &
                                   Gravit, rair, rh2o, cpair, tmelt, latvap, latice,     &
                                   rhmini, iulog, micro_mg_do_hail, micro_mg_do_graupel, &
                                   microp_uniform, do_cldice, use_hetfrz_classnuc,       &
                                   remove_supersat, micro_mg_evap_sed_off,               &
                                   micro_mg_icenuc_rh_off, micro_mg_icenuc_use_meyers,   &
                                   micro_mg_evap_scl_ifs, micro_mg_evap_rhthrsh_ifs,     &
                                   micro_mg_rainfreeze_ifs, micro_mg_ifs_sed,            &
                                   micro_mg_precip_fall_corr, micro_mg_accre_sees_auto,  &
                                   micro_mg_implicit_fall, micro_mg_nccons,              &
                                   micro_mg_nicons, micro_mg_ngcons, micro_mg_nrcons,    &
                                   micro_mg_nscons, micro_mg_precip_frac_method,         &
                                   micro_mg_warm_rain,                                   &
                                   stochastic_emulated_filename_quantile,                &
                                   stochastic_emulated_filename_input_scale,             &
                                   stochastic_emulated_filename_output_scale,            &
                                   micro_mg_dcs_in,                                      &
                                   micro_mg_berg_eff_factor_in, micro_mg_accre_enhan_fact_in, &
                                   micro_mg_autocon_fact_in, micro_mg_autocon_nd_exp_in,      &
                                   micro_mg_autocon_lwp_exp_in, micro_mg_homog_size_in,       &
                                   micro_mg_vtrmi_factor_in,    micro_mg_vtrms_factor_in,     &
                                   micro_mg_effi_factor_in,     micro_mg_iaccr_factor_in,     &
                                   micro_mg_max_nicons_in, micro_mg_ncnst_in,                 &
                                   micro_mg_ninst_in, micro_mg_ngnst_in, micro_mg_nrnst_in,   &
                                   micro_mg_nsnst_in, micro_proc_rates_out, errmsg, errcode)

  !External dependencies:
  use ccpp_kinds,         only: kind_phys
  use micro_pumas_v1,     only: micro_pumas_init
  use pumas_kinds,        only: pumas_r8=>kind_r8
  use micro_pumas_diags,  only: proc_rates_type

  use pumas_stochastic_collect_tau, only: ncd

  !Subroutine (dummy) arguments:

  !Host model constants:
  integer,         intent(in) :: micro_ncol         !Number of horizontal microphysics columns (count)
  integer,         intent(in) :: micro_nlev         !Number of microphysics vertical layers (count)
  real(kind_phys), intent(in) :: gravit !standard gravitational acceleration                    (m s-2)
  real(kind_phys), intent(in) :: rair   !gas constant for dry air                               (J kg-1 K-1)
  real(kind_phys), intent(in) :: rh2o   !gas constat for water vapor                            (J kg-1 K-1)
  real(kind_phys), intent(in) :: cpair  !specific heat of dry air at constant pressure          (J kg-1 K-1)
  real(kind_phys), intent(in) :: tmelt  !freezing point of water                                (K)
  real(kind_phys), intent(in) :: latvap !latent heat of vaporization of water at 0 degrees C    (J kg-1)
  real(kind_phys), intent(in) :: latice !latent heat of fusion of water at 0 degrees C          (J kg-1)
  real(kind_phys), intent(in) :: rhmini !Minimum RH for ice cloud fraction > 0                  (fraction)

  !Host model variables:
  integer, intent(in) :: iulog          !Log output unit number (1)

  !PUMAS-specific parameters:
  !-------------------------
  logical, intent(in) :: micro_mg_do_hail           !flag for PUMAS to simulate hail                              (flag)
  logical, intent(in) :: micro_mg_do_graupel        !flag for PUMAS to simulate graupel                           (flag)
  logical, intent(in) :: microp_uniform             !flag for PUMAS to perform uniform calc.                      (flag)
  logical, intent(in) :: do_cldice                  !flag for PUMAS to simulate cloud ice                         (flag)
  logical, intent(in) :: use_hetfrz_classnuc        !flag to turn on PUMAS heterogeneous freezing                 (flag)
  logical, intent(in) :: remove_supersat            !flag to remove supersaturation after sedimentation loop      (flag)
  logical, intent(in) :: micro_mg_evap_sed_off      !flag to turn off condensate evap. after sedimentation        (flag)
  logical, intent(in) :: micro_mg_icenuc_rh_off     !flag to turn off RH threshold for ice nucleation             (flag)
  logical, intent(in) :: micro_mg_icenuc_use_meyers !flag to use Meyers 1992 temp. dependent ice nucleation       (flag)
  logical, intent(in) :: micro_mg_evap_scl_ifs      !flag to apply IFS precipitation evap. scaling                (flag)
  logical, intent(in) :: micro_mg_evap_rhthrsh_ifs  !flag to use IFS precipitation evap. RH threshold             (flag)
  logical, intent(in) :: micro_mg_rainfreeze_ifs    !flag to freeze rain at 0 degrees C as is done in IFS         (flag)
  logical, intent(in) :: micro_mg_ifs_sed           !flag to use IFS sedimentation fall speeds                    (flag)
  logical, intent(in) :: micro_mg_precip_fall_corr  !flag to ensure non-zero precip fall speed if precip above    (flag)
  logical, intent(in) :: micro_mg_accre_sees_auto   !flag to add autoconverted liuqid to rain before accretion    (flag)
  logical, intent(in) :: micro_mg_implicit_fall     !flag to use implicit fall speed routine for all hydrometeors (flag)
  logical, intent(in) :: micro_mg_nccons            !flag to have PUMAS hold cloud droplet number constant        (flag)
  logical, intent(in) :: micro_mg_nicons            !flag to have PUMAS hold cloud ice number constant            (flag)
  logical, intent(in) :: micro_mg_ngcons            !flag to have PUMAS hold cloud graupel number constant        (flag)
  logical, intent(in) :: micro_mg_nrcons            !flag to have PUMAS hold cloud rain number constant           (flag)
  logical, intent(in) :: micro_mg_nscons            !flag to have PUMAS hold cloud snow number constant           (flag)

  !type of precipitation fraction method (none):
  character(len=*), intent(in) :: micro_mg_precip_frac_method
  !type of warm rain autoconversion/accr.method to use (none):
  character(len=*), intent(in) :: micro_mg_warm_rain
  !neural net file for warm_rain machine learning (none):
  character(len=*), intent(in) :: stochastic_emulated_filename_quantile
  !neural net input scaling values files for warm_rain machine learning (none):
  character(len=*), intent(in) :: stochastic_emulated_filename_input_scale
  !Neural net output scaling values file for warm_rain machine learning (none):
  character(len=*), intent(in) :: stochastic_emulated_filename_output_scale

  real(kind_phys), intent(in) :: micro_mg_dcs_in              !autoconversion size threshold                      (um)
  real(kind_phys), intent(in) :: micro_mg_berg_eff_factor_in  !efficienty factor for Bergeron process             (1)
  real(kind_phys), intent(in) :: micro_mg_accre_enhan_fact_in !accretion enhancement factor                       (1)
  real(kind_phys), intent(in) :: micro_mg_autocon_fact_in     !autoconverion enhancement prefactor                (1)
  real(kind_phys), intent(in) :: micro_mg_autocon_nd_exp_in   !autconversion cloud liquid exponent factor         (1)
  real(kind_phys), intent(in) :: micro_mg_autocon_lwp_exp_in  !autoconversion LWP exponent factor                 (1)
  real(kind_phys), intent(in) :: micro_mg_homog_size_in       !mean volume radius of homoegenous freezing ice     (m)
  real(kind_phys), intent(in) :: micro_mg_vtrmi_factor_in     !ice fall velocity enhancement factor               (1)
  real(kind_phys), intent(in) :: micro_mg_vtrms_factor_in     !snow fall velocity enhancement factor              (1)
  real(kind_phys), intent(in) :: micro_mg_effi_factor_in      !ice effective radius enhancement factor            (1)
  real(kind_phys), intent(in) :: micro_mg_iaccr_factor_in     !ice accretion factor                               (1)
  real(kind_phys), intent(in) :: micro_mg_max_nicons_in       !max allowed ice number concentration               (m-3)

  !In-cloud droplet number concentration if micro_mg_nccons is True (m-3):
  real(kind_phys), intent(in) :: micro_mg_ncnst_in
  !In-cloud ice number concentration if micro_mg_nicons is True     (m-3):
  real(kind_phys), intent(in) :: micro_mg_ninst_in
  !In-cloud graupel number concentration if micro_mg_ngcons is True (m-3):
  real(kind_phys), intent(in) :: micro_mg_ngnst_in
  !In-cloud rain number concentration when micro_mg_nrcons is True  (m-3):
  real(kind_phys), intent(in) :: micro_mg_nrnst_in
  !In-cloud snow number concentration when micro_mg_nscons is True  (m-3):
  real(kind_phys), intent(in) :: micro_mg_nsnst_in
  !-------------------------

  !Output variables:
  character(len=512), intent(out) :: errmsg  !PUMAS/CCPP error message (none)
  integer,            intent(out) :: errcode !CCPP error code (1)

  !Local variables:
  real(pumas_r8) :: micro_mg_dcs              !autoconversion size threshold                      (um)
  real(pumas_r8) :: micro_mg_berg_eff_factor  !efficienty factor for Bergeron process             (1)
  real(pumas_r8) :: micro_mg_accre_enhan_fact !accretion enhancement factor                       (1)
  real(pumas_r8) :: micro_mg_autocon_fact     !autoconverion enhancement prefactor                (1)
  real(pumas_r8) :: micro_mg_autocon_nd_exp   !autconversion cloud liquid exponent factor         (1)
  real(pumas_r8) :: micro_mg_autocon_lwp_exp  !autoconversion LWP exponent factor                 (1)
  real(pumas_r8) :: micro_mg_homog_size       !mean volume radius of homoegenous freezing ice     (m)
  real(pumas_r8) :: micro_mg_vtrmi_factor     !ice fall velocity enhancement factor               (1)
  real(pumas_r8) :: micro_mg_vtrms_factor     !snow fall velocity enhancement factor              (1)
  real(pumas_r8) :: micro_mg_effi_factor      !ice effective radius enhancement factor            (1)
  real(pumas_r8) :: micro_mg_iaccr_factor     !ice accretion factor                               (1)
  real(pumas_r8) :: micro_mg_max_nicons       !max allowed ice number concentration               (m-3)

  !In-cloud droplet number concentration if micro_mg_nccons is True (m-3):
  real(pumas_r8) :: micro_mg_ncnst
  !In-cloud ice number concentration if micro_mg_nicons is True     (m-3):
  real(pumas_r8) :: micro_mg_ninst
  !In-cloud graupel number concentration if micro_mg_ngcons is True (m-3):
  real(pumas_r8) :: micro_mg_ngnst
  !In-cloud rain number concentration when micro_mg_nrcons is True  (m-3):
  real(pumas_r8) :: micro_mg_nrnst
  !In-cloud snow number concentration when micro_mg_nscons is True  (m-3):
  real(pumas_r8) :: micro_mg_nsnst

  !Local PUMAS error message
  character(len=128) :: pumas_errstring

  !microphysics process rates (none)
  type(proc_rates_type), intent(out) :: micro_proc_rates_out

  !Initialize error message and error code:
  errmsg  = ''
  errcode = 0

  !Convert real-type input fields into appropriate kind:
  micro_mg_dcs              = real(micro_mg_dcs_in, pumas_r8)
  micro_mg_berg_eff_factor  = real(micro_mg_berg_eff_factor_in, pumas_r8)
  micro_mg_accre_enhan_fact = real(micro_mg_accre_enhan_fact_in, pumas_r8)
  micro_mg_autocon_fact     = real(micro_mg_autocon_fact_in, pumas_r8)
  micro_mg_autocon_nd_exp   = real(micro_mg_autocon_nd_exp_in, pumas_r8)
  micro_mg_autocon_lwp_exp  = real(micro_mg_autocon_lwp_exp_in, pumas_r8)
  micro_mg_homog_size       = real(micro_mg_homog_size_in, pumas_r8)
  micro_mg_vtrmi_factor     = real(micro_mg_vtrmi_factor_in, pumas_r8)
  micro_mg_vtrms_factor     = real(micro_mg_vtrms_factor_in, pumas_r8)
  micro_mg_effi_factor      = real(micro_mg_effi_factor_in, pumas_r8)
  micro_mg_iaccr_factor     = real(micro_mg_iaccr_factor_in, pumas_r8)
  micro_mg_max_nicons       = real(micro_mg_max_nicons_in, pumas_r8)
  micro_mg_ncnst            = real(micro_mg_ncnst_in, pumas_r8)
  micro_mg_ninst            = real(micro_mg_ninst_in, pumas_r8)
  micro_mg_ngnst            = real(micro_mg_ngnst_in, pumas_r8)
  micro_mg_nrnst            = real(micro_mg_nrnst_in, pumas_r8)
  micro_mg_nsnst            = real(micro_mg_nsnst_in, pumas_r8)

  !Call PUMAS initialization routine:
  call micro_pumas_init( &
           pumas_r8, gravit, rair, rh2o, cpair, &
           tmelt, latvap, latice, rhmini, &
           micro_mg_dcs,                  &
           micro_mg_do_hail,micro_mg_do_graupel, &
           microp_uniform, do_cldice, use_hetfrz_classnuc, &
           micro_mg_precip_frac_method, micro_mg_berg_eff_factor, &
           micro_mg_accre_enhan_fact , &
           micro_mg_autocon_fact , micro_mg_autocon_nd_exp, micro_mg_autocon_lwp_exp, micro_mg_homog_size, &
           micro_mg_vtrmi_factor, micro_mg_vtrms_factor, micro_mg_effi_factor, &
           micro_mg_iaccr_factor, micro_mg_max_nicons, &
           remove_supersat, micro_mg_warm_rain, &
           micro_mg_evap_sed_off, micro_mg_icenuc_rh_off, micro_mg_icenuc_use_meyers, &
           micro_mg_evap_scl_ifs, micro_mg_evap_rhthrsh_ifs, &
           micro_mg_rainfreeze_ifs,  micro_mg_ifs_sed, micro_mg_precip_fall_corr,&
           micro_mg_accre_sees_auto, micro_mg_implicit_fall, &
           micro_mg_nccons, micro_mg_nicons, micro_mg_ncnst, &
           micro_mg_ninst, micro_mg_ngcons, micro_mg_ngnst, &
           micro_mg_nrcons, micro_mg_nrnst, micro_mg_nscons, micro_mg_nsnst, &
           stochastic_emulated_filename_quantile, stochastic_emulated_filename_input_scale, &
           stochastic_emulated_filename_output_scale, iulog, pumas_errstring)

  !Set error code to non-zero value if PUMAS returns an error message:
  if (trim(pumas_errstring) /= "") then
    errcode = 1
    errmsg = trim(pumas_errstring)
    return
  end if

  ! Allocate the proc_rates DDT
  call micro_proc_rates_out%allocate(micro_ncol, micro_nlev, ncd, micro_mg_warm_rain, pumas_errstring)

  !Set error code to non-zero value if PUMAS returns an error message:
  if (trim(pumas_errstring) /= "") then
    errcode = 1
    errmsg = trim(pumas_errstring)
  end if

  end subroutine micro_pumas_ccpp_init

  !> \section arg_table_micro_pumas_ccpp_run Argument Table
  !! \htmlinclude micro_pumas_ccpp_run.html
  subroutine micro_pumas_ccpp_run(micro_ncol, micro_nlev, micro_nlevp1,             &
                                  micro_dust_nbins, pumas_timestep,              &
                                  pumas_airT, pumas_airq, pumas_cldliq,    &
                                  pumas_cldice,   pumas_numliq,               &
                                  pumas_numice,   pumas_rainliq,              &
                                  pumas_snowice,  pumas_numrain,              &
                                  pumas_numsnow,  pumas_graupice,             &
                                  pumas_numgraup, pumas_relvar,               &
                                  pumas_accre_enhan, pumas_pmid,              &
                                  pumas_pdel, pumas_pint,                     &
                                  pumas_strat_cldfrc, pumas_strat_liq_cldfrc, &
                                  pumas_strat_ice_cldfrc, pumas_qsatfac,      &
                                  pumas_naai, pumas_npccn,                    &
                                  pumas_rndst, pumas_nacon,                   &
                                  pumas_snowice_tend_external,                   &
                                  pumas_numsnow_tend_external,                   &
                                  pumas_effi_external, pumas_frzimm,          &
                                  pumas_frzcnt, pumas_frzdep,                 &
! output vars
                                  pumas_qcsinksum_rate1ord_out,                     &
                                  pumas_airT_tend_out, pumas_airq_tend_out,         &
                                  pumas_cldliq_tend_out, pumas_cldice_tend_out,     &
                                  pumas_numliq_tend_out, pumas_numice_tend_out,     &
                                  pumas_rainliq_tend_out, pumas_snowice_tend_out,   &
                                  pumas_numrain_tend_out, pumas_numsnow_tend_out,   &
                                  pumas_graupice_tend_out, pumas_numgraup_tend_out, &
                                  pumas_effc_out, pumas_effc_fn_out,                &
                                  pumas_effi_out, pumas_sadice_out,                 &
                                  pumas_sadsnow_out, pumas_prect_out,               &
                                  pumas_preci_out, pumas_prec_evap_out,             &
                                  pumas_am_evap_st_out, pumas_prec_prod_out,        &
                                  pumas_cmeice_out, pumas_deffi_out,                &
                                  pumas_pgamrad_out, pumas_lamcrad_out,             &
                                  pumas_snowice_in_prec_out,                        &
                                  pumas_scaled_diam_snow_out,                       &
                                  pumas_graupice_in_prec_out,                       &
                                  pumas_numgraup_vol_in_prec_out,                   &
                                  pumas_scaled_diam_graup_out,                      &
                                  pumas_lflx_out, pumas_iflx_out, pumas_gflx_out,   &
                                  pumas_rflx_out, pumas_sflx_out,                   &
                                  pumas_rainliq_in_prec_out, pumas_reff_rain_out,   &
                                  pumas_reff_snow_out, pumas_reff_grau_out,         &
                                  pumas_numrain_vol_in_prec_out,                    &
                                  pumas_numsnow_vol_in_prec_out,                    &
                                  pumas_refl_out, pumas_arefl_out,                  &
                                  pumas_areflz_out, pumas_frefl_out,                &
                                  pumas_csrfl_out, pumas_acsrfl_out,                &
                                  pumas_fcsrfl_out, pumas_refl10cm_out,             &
                                  pumas_reflz10cm_out, pumas_rercld_out,            &
                                  pumas_ncai_out, pumas_ncal_out,                   &
                                  pumas_rainliq_out, pumas_snowice_out,             &
                                  pumas_numrain_vol_out, pumas_numsnow_vol_out,     &
                                  pumas_diam_rain_out, pumas_diam_snow_out,         &
                                  pumas_graupice_out, pumas_numgraup_vol_out,       &
                                  pumas_diam_graup_out, pumas_freq_graup_out,       &
                                  pumas_freq_snow_out, pumas_freq_rain_out,         &
                                  pumas_frac_ice_out, pumas_frac_cldliq_tend_out,   &
                                  pumas_rain_evap_out, micro_proc_rates_inout,      &
                                  errmsg, errcode)

    !External dependencies:
    use ccpp_kinds,        only: kind_phys
    use micro_pumas_v1,    only: micro_pumas_tend
    use micro_pumas_diags, only: proc_rates_type
    use pumas_kinds,       only: pumas_r8=>kind_r8

    !Subroutine (dummy) input arguments:

    !Host model dimensions/parameters:
    integer,         intent(in) :: micro_ncol         !Number of horizontal microphysics columns (count)
    integer,         intent(in) :: micro_nlev         !Number of microphysics vertical layers (count)
    integer,         intent(in) :: micro_nlevp1       !Number of microphysics vertical interfaces (count)
    integer,         intent(in) :: micro_dust_nbins   !Number of dust size bins

    real(pumas_r8), intent(in) :: pumas_timestep  !Microphysics time step (s)

    !Host model state variables:

    !Microphysics Air temperature (K)
    real(pumas_r8), intent(in) :: pumas_airT(:,:)
    !Microphysics Water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_airq(:,:)
    !Microphysics cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_cldliq(:,:)
    !Microphysics cloud ice mixing ratio wrt moist air and condensed water (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_cldice(:,:)
    !microphysics mass number concentration of cloud liquid water wrt moist air and condensed water (kg-1)
    real(pumas_r8), intent(in) :: pumas_numliq(:,:)
    !microphysics mass number concentration of cloud ice wrt moist air and condensed water (kg-1)
    real(pumas_r8), intent(in) :: pumas_numice(:,:)
    !microphysics rain mixing ratio wrt moist air and condensed water (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_rainliq(:,:)
    !microphysics snow mixing ratio wrt moist air and condensed water (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_snowice(:,:)
    !microphysics mass number concentration of rain wrt moist air and condensed water (kg-1)
    real(pumas_r8), intent(in) :: pumas_numrain(:,:)
    !microphysics mass number concentration of snow wrt moist air and condensed water (kg-1)
    real(pumas_r8), intent(in) :: pumas_numsnow(:,:)
    !microphysics graupel mixing ratio wrt moist air and condensed water (kg kg-1)
    real(pumas_r8), intent(in) :: pumas_graupice(:,:)
    !microphysics mass number concentration of graupel wrt moist air and condensed water (kg-1)
    real(pumas_r8), intent(in) :: pumas_numgraup(:,:)
    !microphysics relative variance of cloud water (1)
    real(pumas_r8), intent(in) :: pumas_relvar(:,:)
    !microphysics accretion enhancement factor (1)
    real(pumas_r8), intent(in) :: pumas_accre_enhan(:,:)
    !microphysics air pressure (Pa)
    real(pumas_r8), intent(in) :: pumas_pmid(:,:)
    !microphysics air pressure thickness (Pa)
    real(pumas_r8), intent(in) :: pumas_pdel(:,:)
    !microphysics air pressure at interfaces (Pa)
    real(pumas_r8), intent(in) :: pumas_pint(:,:)
    !microphysics stratiform cloud area fraction (fraction)
    real(pumas_r8), intent(in) :: pumas_strat_cldfrc(:,:)
    !microphysics stratiform cloud liquid area fraction (fraction)
    real(pumas_r8), intent(in) :: pumas_strat_liq_cldfrc(:,:)
    !microphysics stratiform cloud ice area fraction (fraction)
    real(pumas_r8), intent(in) :: pumas_strat_ice_cldfrc(:,:)
    !microphysics subgrid cloud water saturation scaling factor (1)
    real(pumas_r8), intent(in) :: pumas_qsatfac(:,:)
    !microphysics tendency of activated ice nuclei mass number concentration (kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_naai(:,:)
    !microphysics tendency of activated cloud condensation nuclei mass number concentration (kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_npccn(:,:)
    !microphysics dust radii by size bin  (m)
    real(pumas_r8), intent(in) :: pumas_rndst(:,:,:)
    !microphysics dust number concentration by size bin (m-3)
    real(pumas_r8), intent(in) :: pumas_nacon(:,:,:)
    !microphysics tendency of snow mixing ratio wrt moist air and condensed water from external microphysics (kg kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_snowice_tend_external(:,:)
    !microphysics tendency of mass number concentration of snow wrt moist air and condensed water from external microphysics
    !(kg-1 s-1)
    real(pumas_r8), intent(in) :: pumas_numsnow_tend_external(:,:)
    !microphysics effective radius of stratiform cloud ice particle from external microphysics (m)
    real(pumas_r8), intent(in) :: pumas_effi_external(:,:)
    !microphysics tendency of cloud liquid droplet number concentration due to immersion freezing (cm-3)
    real(pumas_r8), intent(in) :: pumas_frzimm(:,:)
    !microphysics tendency of cloud liquid droplet number concentration due to contact freezing (cm-3)
    real(pumas_r8), intent(in) :: pumas_frzcnt(:,:)
    !microphysics tendency of cloud ice number concentration due to deposition nucleation (cm-3)
    real(pumas_r8), intent(in) :: pumas_frzdep(:,:)

    !Subroutine output arguments:

    !microphysics direct conversion rate of stratiform cloud water to precipitation (s-1)
    real(pumas_r8), intent(out) :: pumas_qcsinksum_rate1ord_out(:,:)
    !microphysics tendency of dry air enthalpy at constant pressure (J kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_airT_tend_out(:,:)
    !microphysics tendency of water vapor mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_airq_tend_out(:,:)
    !microphysics tendency of cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_cldliq_tend_out(:,:)
    !microphysics tendency of cloud ice mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_cldice_tend_out(:,:)
    !microphysics tendency of mass number concentration of cloud liquid water wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_numliq_tend_out(:,:)
    !microphysics tendency of mass number concentration of cloud ice wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_numice_tend_out(:,:)
    !microphysics tendency of rain mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_rainliq_tend_out(:,:)
    !microphysics tendency of snow mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_snowice_tend_out(:,:)
    !microphysics tendency of mass number concentration of rain wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_numrain_tend_out(:,:)
    !microphysics tendency of mass number concentration of snow wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_numsnow_tend_out(:,:)
    !microphysics tendency of graupel mixing ratio wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_graupice_tend_out(:,:)
    !microphysics tendency of mass number concentration of graupel wrt moist air and condensed water (kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_numgraup_tend_out(:,:)
    !microphysics effective radius of stratiform cloud liquid water particle (um)
    real(pumas_r8), intent(out) :: pumas_effc_out(:,:)
    !microphysics effective radius of stratiform cloud liquid water particle assuming droplet number concentration of 1e8 kg-1 (um)
    real(pumas_r8), intent(out) :: pumas_effc_fn_out(:,:)
    !microphysics effective radius of stratiform cloud ice particle (um)
    real(pumas_r8), intent(out) :: pumas_effi_out(:,:)
    !microphysics cloud ice surface area density (cm2 cm-3)
    real(pumas_r8), intent(out) :: pumas_sadice_out(:,:)
    !microphysics snow surface area density (cm2 cm-3)
    real(pumas_r8), intent(out) :: pumas_sadsnow_out(:,:)
    !microphysics LWE large scale precipitation rate at surface (m s-1)
    real(pumas_r8), intent(out) :: pumas_prect_out(:)
    !microphysics LWE large scale snowfall rate at surface (m s-1)
    real(pumas_r8), intent(out) :: pumas_preci_out(:)
    !microphysics precipitation evaporation rate wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_prec_evap_out(:,:)
    !microphysics precipitation evaporation area (fraction)
    real(pumas_r8), intent(out) :: pumas_am_evap_st_out(:,:)
    !microphysics precipitation production rate wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_prec_prod_out(:,:)
    !microphysics condensation minus evaporation rate of in-cloud ice wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_cmeice_out(:,:)
    !microphysics effective diameter of stratiform cloud ice particles for radiation (um)
    real(pumas_r8), intent(out) :: pumas_deffi_out(:,:)
    !microphysics cloud particle size distribution shape parameter (1)
    real(pumas_r8), intent(out) :: pumas_pgamrad_out(:,:)
    !microphysics cloud particle size distribution slope parameter (1)
    real(pumas_r8), intent(out) :: pumas_lamcrad_out(:,:)
    !microphysics snow mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(pumas_r8), intent(out) :: pumas_snowice_in_prec_out(:,:)
    !microphysics snow scaled diameter (m)
    real(pumas_r8), intent(out) :: pumas_scaled_diam_snow_out(:,:)
    !microphysics graupel mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(pumas_r8), intent(out) :: pumas_graupice_in_prec_out(:,:)
    !microphysics graupel number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), intent(out) :: pumas_numgraup_vol_in_prec_out(:,:)
    !microphysics graupel scaled diameter (m)
    real(pumas_r8), intent(out) :: pumas_scaled_diam_graup_out(:,:)
    !microphysics cloud liquid sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(out) :: pumas_lflx_out(:,:)
    !microphysics cloud ice sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(out) :: pumas_iflx_out(:,:)
    !microphysics graupel sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(out) :: pumas_gflx_out(:,:)
    !microphysics rain sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(out) :: pumas_rflx_out(:,:)
    !microphysics snow sedimentation flux (kg m-2 s-1)
    real(pumas_r8), intent(out) :: pumas_sflx_out(:,:)
    !microphysics rain mixing ratio wrt moist air and condensed water of new state in precipitating fraction of gridcell (kg kg-1)
    real(pumas_r8), intent(out) :: pumas_rainliq_in_prec_out(:,:)
    !microphysics effective radius of stratiform rain particle (um)
    real(pumas_r8), intent(out) :: pumas_reff_rain_out(:,:)
    !microphysics effective radius of stratiform snow particle (um)
    real(pumas_r8), intent(out) :: pumas_reff_snow_out(:,:)
    !microphysics effective radius of stratiform graupel particle (um)
    real(pumas_r8), intent(out) :: pumas_reff_grau_out(:,:)
    !microphysics rain number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), intent(out) :: pumas_numrain_vol_in_prec_out(:,:)
    !microphysics snow number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), intent(out) :: pumas_numsnow_vol_in_prec_out(:,:)
    !microphysics analytic radar reflectivity at 94 GHz in precipitating fraction of gridcell (dBZ)
    real(pumas_r8), intent(out) :: pumas_refl_out(:,:)
    !microphysics analytic radar reflectivity at 94 GHz (dBZ)
    real(pumas_r8), intent(out) :: pumas_arefl_out(:,:)
    !microphysics analytic radar reflectivity z factor at 94 GHz (mm6 m-3)
    real(pumas_r8), intent(out) :: pumas_areflz_out(:,:)
    !microphysics fraction of gridcell with nonzero radar reflectivity (fraction)
    real(pumas_r8), intent(out) :: pumas_frefl_out(:,:)
    !microphysics analytic radar reflectivity at 94 GHz with CloudSat thresholds in precipitating fraction of gridcell (dBZ)
    real(pumas_r8), intent(out) :: pumas_csrfl_out(:,:)
    !microphysics analytic radar reflectivity at 94 GHz with CloudSat thresholds (dBZ)
    real(pumas_r8), intent(out) :: pumas_acsrfl_out(:,:)
    !microphysics fraction of gridcell with nonzero radar reflectivity with CloudSat thresholds (fraction)
    real(pumas_r8), intent(out) :: pumas_fcsrfl_out(:,:)
    !microphysics analytic radar reflectivity at 10 cm wavelength (dBZ)
    real(pumas_r8), intent(out) :: pumas_refl10cm_out(:,:)
    !microphysics analytic radar reflectivity z factor at 10 cm wavelength (mm6 m-3)
    real(pumas_r8), intent(out) :: pumas_reflz10cm_out(:,:)
    !microphysics effective radius of stratiform cloud liquid plus rain particles (m)
    real(pumas_r8), intent(out) :: pumas_rercld_out(:,:)
    !microphysics available ice nuclei number concentration of new state (m-3)
    real(pumas_r8), intent(out) :: pumas_ncai_out(:,:)
    !microphysics available cloud condensation nuclei number concentration of new state (m-3)
    real(pumas_r8), intent(out) :: pumas_ncal_out(:,:)
    !microphysics rain mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), intent(out) :: pumas_rainliq_out(:,:)
    !microphysics snow mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), intent(out) :: pumas_snowice_out(:,:)
    !microphysics rain number concentration of new state (m-3)
    real(pumas_r8), intent(out) :: pumas_numrain_vol_out(:,:)
    !microphysics snow number concentration of new state in precipitating fraction of gridcell (m-3)
    real(pumas_r8), intent(out) :: pumas_numsnow_vol_out(:,:)
    !microphysics average diameter of stratiform rain particle (m)
    real(pumas_r8), intent(out) :: pumas_diam_rain_out(:,:)
    !microphysics average diameter of stratiform snow particle (m)
    real(pumas_r8), intent(out) :: pumas_diam_snow_out(:,:)
    !microphysics graupel mixing ratio wrt moist air and condensed water of new state (kg kg-1)
    real(pumas_r8), intent(out) :: pumas_graupice_out(:,:)
    !microphysics graupel number concentration of new state (m-3)
    real(pumas_r8), intent(out) :: pumas_numgraup_vol_out(:,:)
    !microphysics average diameter of stratiform graupel particle (m)
    real(pumas_r8), intent(out) :: pumas_diam_graup_out(:,:)
    !microphysics fraction of gridcell with graupel (fraction)
    real(pumas_r8), intent(out) :: pumas_freq_graup_out(:,:)
    !microphysics fraction of gridcell with snow (fraction)
    real(pumas_r8), intent(out) :: pumas_freq_snow_out(:,:)
    !microphysics fraction of gridcell with rain (fraction)
    real(pumas_r8), intent(out) :: pumas_freq_rain_out(:,:)
    !microphysics fraction of frozen water to total condensed water (fraction)
    real(pumas_r8), intent(out) :: pumas_frac_ice_out(:,:)
    !microphysics fraction of cloud liquid tendency applied to state (fraction)
    real(pumas_r8), intent(out) :: pumas_frac_cldliq_tend_out(:,:)
    !microphysics rain evaporation rate wrt moist air and condensed water (kg kg-1 s-1)
    real(pumas_r8), intent(out) :: pumas_rain_evap_out(:,:)
    !microphysics process rates (none)
    type(proc_rates_type), intent(inout) :: micro_proc_rates_inout

    character(len=512), intent(out) :: errmsg  !PUMAS/CCPP error message (none)
    integer,            intent(out) :: errcode !CCPP error code (1)

    !Local PUMAS error message
    character(len=128) :: pumas_errstring

    !Initialize error message and error code:
    errmsg  = ''
    errcode = 0


    !Call main PUMAS run routine:
    !---------------------------

    call micro_pumas_tend( &
        micro_ncol,             micro_nlev,     pumas_timestep,      &
        pumas_airT,                   pumas_airq,                                &
        pumas_cldliq,                 pumas_cldice,                              &
        pumas_numliq,                 pumas_numice,                              &
        pumas_rainliq,                pumas_snowice,                             &
        pumas_numrain,                pumas_numsnow,                             &
        pumas_graupice,               pumas_numgraup,                            &
        pumas_relvar,                 pumas_accre_enhan,                   &
        pumas_pmid,                   pumas_pdel, pumas_pint,                          &
        pumas_strat_cldfrc,           pumas_strat_liq_cldfrc,                    &
        pumas_strat_ice_cldfrc,       pumas_qsatfac,                             &
        pumas_qcsinksum_rate1ord_out,                                          &
        pumas_naai,                   pumas_npccn,                               &
        pumas_rndst,                  pumas_nacon,                               &
        pumas_airT_tend_out,              pumas_airq_tend_out,                           &
        pumas_cldliq_tend_out,            pumas_cldice_tend_out,                         &
        pumas_numliq_tend_out,            pumas_numice_tend_out,                         &
        pumas_rainliq_tend_out,           pumas_snowice_tend_out,                        &
        pumas_numrain_tend_out,           pumas_numsnow_tend_out,                        &
        pumas_graupice_tend_out,          pumas_numgraup_tend_out,                       &
        pumas_effc_out,                   pumas_effc_fn_out,        pumas_effi_out,                &
        pumas_sadice_out,                 pumas_sadsnow_out,                             &
        pumas_prect_out,                  pumas_preci_out,                               &
        pumas_prec_evap_out,              pumas_am_evap_st_out,                          &
        pumas_prec_prod_out,                                                   &
        pumas_cmeice_out,                 pumas_deffi_out,                               &
        pumas_pgamrad_out,                pumas_lamcrad_out,                             &
        pumas_snowice_in_prec_out,    pumas_scaled_diam_snow_out,                &
        pumas_graupice_in_prec_out,   pumas_numgraup_vol_in_prec_out,            &
        pumas_scaled_diam_graup_out,                                       &
        pumas_lflx_out,                   pumas_iflx_out,                                &
        pumas_gflx_out,                                                        &
        pumas_rflx_out,                   pumas_sflx_out,           pumas_rainliq_in_prec_out, &
        pumas_reff_rain_out,              pumas_reff_snow_out,      pumas_reff_grau_out,           &
        pumas_numrain_vol_in_prec_out,    pumas_numsnow_vol_in_prec_out,         &
        pumas_refl_out,                   pumas_arefl_out,          pumas_areflz_out,              &
        pumas_frefl_out,                  pumas_csrfl_out,          pumas_acsrfl_out,              &
        pumas_fcsrfl_out,   pumas_refl10cm_out,     pumas_reflz10cm_out,      pumas_rercld_out,              &
        pumas_ncai_out,                   pumas_ncal_out,                                &
        pumas_rainliq_out,            pumas_snowice_out,                         &
        pumas_numrain_vol_out,        pumas_numsnow_vol_out,                     &
        pumas_diam_rain_out,          pumas_diam_snow_out,                       &
        pumas_graupice_out,           pumas_numgraup_vol_out, pumas_diam_graup_out,    &
        pumas_freq_graup_out,             pumas_freq_snow_out,        pumas_freq_rain_out,         &
        pumas_frac_ice_out,               pumas_frac_cldliq_tend_out,                    &
        micro_proc_rates_inout, pumas_errstring,                     &
        pumas_snowice_tend_external,  pumas_numsnow_tend_external,               &
        pumas_effi_external,          pumas_rain_evap_out,                     &
        pumas_frzimm,                 pumas_frzcnt,           pumas_frzdep           )

     !---------------------------


    !Set error code to non-zero value if PUMAS returns an error message:
    if (trim(errmsg) /= "") then
      errcode = 1
      errmsg  = trim(pumas_errstring)
    end if

  end subroutine micro_pumas_ccpp_run

end module micro_pumas_ccpp
