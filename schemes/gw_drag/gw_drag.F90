module gw_drag

!--------------------------------------------------------------------------
! CAM and WACCM gravity wave parameterizations were merged by Sean Patrick
! Santos in Summer 2013, and at the same time, gw_drag was split into
! various modules. This is the CAM interface and driver module. The below
! notes are for the old CAM and WACCM versions of gw_drag.
!--------------------------------------------------------------------------
! This file came from wa17 and was modified by Fabrizio: 07-02-2004
! Standard gw_drag with modification (6) of latitude profile of gw spectrum
!--------------------------------------------------------------------------
! Purpose:
!
! Module to compute the forcing due to parameterized gravity waves. Both an
! orographic and an internal source spectrum are considered.
!
! Author: Byron Boville
!
!--------------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys
  use gw_common, only: GWBand, handle_err
  use interpolate_data, only: lininterp

  implicit none

  save
  private

  !
  ! PUBLIC: interfaces
  !
  public :: gw_drag_init                  ! Initialization
  public :: gw_drag_run                   ! interface to actual parameterization

  !
  ! PRIVATE: Rest of the data and interfaces are private to this module
  !

  ! Maximum wave number and width of spectrum bins.
  integer                              :: iulog
  integer                              :: pgwv
  real(kind_phys)                      :: gw_dc
  integer                              :: pgwv_long
  real(kind_phys)                      :: gw_dc_long

  ! Whether or not to enforce an upper boundary condition of tau = 0.
  ! (Like many variables, this is only here to hold the value between
  ! the readnl phase and the init phase of the CAM physics; only gw_common
  ! should actually use it.)
  logical                              :: tau_0_ubc = .false.
  ! Beres (deep convection).
  real(kind_phys)                      :: effgw_beres_dp
  ! Beres (shallow convection).
  real(kind_phys)                      :: effgw_beres_sh
  ! C&M scheme.
  real(kind_phys)                      :: effgw_cm
  ! C&M scheme (inertial waves).
  real(kind_phys)                      :: effgw_cm_igw
  ! Orography.
  real(kind_phys)                      :: effgw_oro
  ! fcrit2 for the mid-scale waves has been made a namelist variable to
  ! facilitate backwards compatibility with the CAM3 version of this
  ! parameterization.  In CAM3, fcrit2=0.5.
  real(kind_phys)                      :: rearth   ! earth radius
  real(kind_phys)                      :: fcrit2   ! critical froude number squared

  ! Ridge scheme.
  logical                              :: use_gw_rdg_beta
  integer                              :: n_rdg_beta
  real(kind_phys)                      :: effgw_rdg_beta
  real(kind_phys)                      :: effgw_rdg_beta_max
  real(kind_phys)                      :: rdg_beta_cd_llb  ! Low-level obstacle drag coefficient Ridge scheme.
  logical                              :: trpd_leewv_rdg_beta

  logical                              :: use_gw_rdg_gamma
  integer                              :: n_rdg_gamma
  real(kind_phys)                      :: effgw_rdg_gamma
  real(kind_phys)                      :: effgw_rdg_gamma_max
  real(kind_phys)                      :: rdg_gamma_cd_llb
  logical                              :: trpd_leewv_rdg_gamma

  character(len=256)                   :: bnd_rdggm ! full pathname for meso-Gamma ridge dataset
  character(len=256)                   :: bnd_topo  ! full pathname for topo dataset
  ! Whether or not to limit tau *before* applying any efficiency factors.
  logical                              :: gw_limit_tau_without_eff = .false.
  logical                              :: gw_lndscl_sgh = .true. ! scale SGH by land frac
  real(kind_phys)                      :: gw_prndl = 0.25_kind_phys
  ! Whether or not to apply tendency max
  real(kind_phys)                      :: gw_qbo_hdepth_scaling = 1._kind_phys ! heating depth scaling factor
  logical                              :: gw_top_taper = .false.
  ! Width of gaussian used to create frontogenesis tau profile [m s-1].
  real(kind_phys)                      :: front_gaussian_width = -huge(1._kind_phys)
  real(kind_phys)                      :: alpha_gw_movmtn
  real(kind_phys)                      :: effgw_rdg_resid
  real(kind_phys)                      :: effgw_movmtn_pbl
  integer                              :: movmtn_source
  real(kind_phys)                      :: movmtn_psteer
  real(kind_phys)                      :: movmtn_plaunch

  real(kind_phys), parameter :: unset_kind_phys = huge(1._kind_phys)

  ! A mid-scale "band" with only stationary waves (l = 0).
  type(GWBand) :: band_oro
  ! Medium scale waves.
  type(GWBand) :: band_mid
  ! Long scale waves for IGWs.
  type(GWBand) :: band_long
  ! Medium scale waves for moving mountain
  type(GWBand) :: band_movmtn

  ! Top level for gravity waves.
  integer :: ktop = 1

  ! Bottom level for frontal waves.
  integer :: kbot_front

  ! Factor for SH orographic waves.
  real(kind_phys) :: gw_oro_south_fac = 1._kind_phys

  ! Frontogenesis function critical threshold.
  real(kind_phys) :: frontgfc = unset_kind_phys

  ! Files to read Beres source spectra from.
  character(len=256) :: gw_drag_file
  character(len=256) :: gw_drag_file_sh
  character(len=256) :: gw_drag_file_mm
  real(kind_phys)   :: gravit          ! gravitational acceleration (m s-2)
  real(kind_phys)   :: rair            ! Dry air gas constant     (J K-1 kg-1)
  real(kind_phys)   :: pi
  real(kind_phys)   :: al0
  real(kind_phys)   :: dlat0
  real(kind_phys), allocatable   :: pref_edge(:), pref_mid(:)
  real(kind_phys)   :: degree2radian
  integer           :: ncid_topo
  logical           :: masterproc

  logical         ::  use_gw_oro
  logical         ::  use_gw_front
  logical         ::  use_gw_rdg_resid
  logical         ::  use_gw_front_igw
  logical         ::  use_gw_convect_dp
  logical         ::  use_gw_convect_sh
  logical         ::  use_simple_phys
  logical         ::  use_gw_movmtn_pbl
  logical         ::  do_molec_diff
  integer         ::  nbot_molec
  ! Horzontal wavelengths [m].
  real(kind_phys) :: wavelength_mid
  real(kind_phys) :: wavelength_long

  ! Background stress source strengths.
  real(kind_phys) :: taubgnd = unset_kind_phys
  real(kind_phys) :: taubgnd_igw = unset_kind_phys

  ! Whether or not to use a polar taper for frontally generated waves.
  logical :: gw_polar_taper = .false.

  ! Whether or not to apply tendency max
  logical :: gw_apply_tndmax = .true.

  ! Water constituent indices for budget
  integer :: ixcldliq = -1
  integer :: ixcldice = -1

  ! Prefixes for history field names
  character(len=1), parameter :: cm_pf = " "
  character(len=1), parameter :: cm_igw_pf = "I"
  character(len=1), parameter :: beres_dp_pf = "B"
  character(len=1), parameter :: beres_sh_pf = "S"

  ! namelist
  logical          :: history_amwg                   ! output the variables used by the AMWG diag package

  real(kind_phys), pointer :: vramp(:) => null()

!==========================================================================
contains
!==========================================================================

!> \section arg_table_gw_drag_init Argument Table
!! \htmlinclude gw_drag_init.html
  subroutine gw_drag_init( &
    iulog_in, &
    ktop_in, &
    masterproc_in, &
    ncol, &
    pver, &
    gravit_in, &
    rair_in, &
    pi_in, &
    fcrit2_in, &
    rearth_in, &
    pref_edge_in, &
    pref_mid_in, &
    pgwv_nl, &
    gw_dc_nl, &
    pgwv_long_nl, &
    gw_dc_long_nl, &
    tau_0_ubc_nl, &
    effgw_beres_dp_nl, &
    effgw_beres_sh_nl, &
    effgw_cm_nl, &
    effgw_cm_igw_nl, &
    effgw_oro_nl, &
    frontgfc_nl, &
    gw_drag_file_nl, &
    gw_drag_file_sh_nl, &
    gw_drag_file_mm_nl, &
    taubgnd_nl, &
    taubgnd_igw_nl, &
    gw_polar_taper_nl, &
    use_gw_rdg_beta_nl, &
    n_rdg_beta_nl, &
    effgw_rdg_beta_nl, &
    effgw_rdg_beta_max_nl, &
    rdg_beta_cd_llb_nl, &
    trpd_leewv_rdg_beta_nl, &
    use_gw_rdg_gamma_nl, &
    n_rdg_gamma_nl, &
    effgw_rdg_gamma_nl, &
    effgw_rdg_gamma_max_nl, &
    rdg_gamma_cd_llb_nl, &
    trpd_leewv_rdg_gamma_nl, &
    bnd_topo_nl, &
    bnd_rdggm_nl, &
    gw_oro_south_fac_nl, &
    gw_limit_tau_without_eff_nl, &
    gw_lndscl_sgh_nl, &
    gw_prndl_nl, &
    gw_apply_tndmax_nl, &
    gw_qbo_hdepth_scaling_nl, &
    gw_top_taper_nl, &
    front_gaussian_width_nl, &
    alpha_gw_movmtn_nl, &
    use_gw_rdg_resid_in, &
    effgw_rdg_resid_in, &
    effgw_movmtn_pbl_in, &
    movmtn_source_in, &
    movmtn_psteer_in, &
    movmtn_plaunch_in, &
    gw_rdg_do_divstream_nl, gw_rdg_C_BetaMax_DS_nl, gw_rdg_C_GammaMax_nl, &
    gw_rdg_Frx0_nl, gw_rdg_Frx1_nl, gw_rdg_C_BetaMax_SM_nl, gw_rdg_Fr_c_nl, &
    gw_rdg_do_smooth_regimes_nl, gw_rdg_do_adjust_tauoro_nl, &
    gw_rdg_do_backward_compat_nl, gw_rdg_orohmin_nl, gw_rdg_orovmin_nl, &
    gw_rdg_orostratmin_nl, gw_rdg_orom2min_nl, gw_rdg_do_vdiff_nl, &
    use_gw_oro_in, &
    use_gw_front_in, &
    use_gw_front_igw_in, &
    use_gw_convect_dp_in, &
    use_gw_convect_sh_in, &
    use_simple_phys_in, &
    use_gw_movmtn_pbl_in, &
    do_molec_diff_in, &
    nbot_molec_in, &
    wavelength_mid_in, &
    wavelength_long_in, &
    errmsg, &
    errflg)

    use ref_pres, only: press_lim_idx

    use gw_common, only: gw_common_init, gw_prof
    use gw_rdg, only: gw_rdg_init
    use gw_front, only: gw_front_init
    use gw_movmtn, only: gw_movmtn_init
    use gw_convect, only: gw_beres_init
    !-----------------------------------------------------------------------
    ! Time independent initialization for multiple gravity wave
    ! parameterization.
    !-----------------------------------------------------------------------

    integer, intent(in)             :: ncol
    integer, intent(in)             :: pver
    real(kind_phys), intent(in)     :: gravit_in          ! gravitational acceleration (m s-2)
    real(kind_phys), intent(in)     :: rair_in            ! Dry air gas constant     (J K-1 kg-1)
    real(kind_phys), intent(in)     :: pi_in
    ! Maximum wave number and width of spectrum bins.
    integer, intent(in)             :: pgwv_nl
    real(kind_phys), intent(in)     :: gw_dc_nl
    integer, intent(in)             :: pgwv_long_nl
    real(kind_phys), intent(in)     :: gw_dc_long_nl
    ! Whether or not to enforce an upper boundary condition of tau = 0.
    ! (Like many variables, this is only here to hold the value between
    ! the readnl phase and the init phase of the CAM physics; only gw_common
    ! should actually use it.)
    logical, intent(in)             :: tau_0_ubc_nl
    real(kind_phys), intent(in)     :: pref_edge_in(:)
    real(kind_phys), intent(in)     :: pref_mid_in(:)
    ! Beres (deep convection).
    real(kind_phys), intent(in)     :: effgw_beres_dp_nl
    ! Beres (shallow convection).
    real(kind_phys), intent(in)     :: effgw_beres_sh_nl
    ! C&M scheme.
    real(kind_phys), intent(in)     :: effgw_cm_nl
    ! C&M scheme (inertial waves).
    real(kind_phys), intent(in)     :: effgw_cm_igw_nl
    ! Orography.
    real(kind_phys), intent(in)     :: effgw_oro_nl
    ! fcrit2 for the mid-scale waves has been made a namelist variable to
    ! facilitate backwards compatibility with the CAM3 version of this
    ! parameterization.  In CAM3, fcrit2=0.5.
    real(kind_phys), intent(in)             :: fcrit2_in   ! critical froude number squared
    real(kind_phys), intent(in)             :: rearth_in   ! earth radius
    ! Frontogenesis function critical threshold.
    real(kind_phys), intent(in)             :: frontgfc_nl
    ! Files to read Beres source spectra from.
    character(len=256), intent(in)             :: gw_drag_file_nl
    character(len=256), intent(in)             :: gw_drag_file_sh_nl
    character(len=256), intent(in)             :: gw_drag_file_mm_nl
    ! Background stress source strengths.
    real(kind_phys), intent(in)             :: taubgnd_nl
    real(kind_phys), intent(in)             :: taubgnd_igw_nl
    ! Whether or not to use a polar taper for frontally generated waves.
    logical, intent(in)             :: gw_polar_taper_nl
    ! Ridge scheme.
    logical, intent(in)              :: use_gw_rdg_beta_nl
    integer, intent(in)             :: n_rdg_beta_nl
    real(kind_phys), intent(in)             :: effgw_rdg_beta_nl
    real(kind_phys), intent(in)             :: effgw_rdg_beta_max_nl
    real(kind_phys), intent(in)             :: rdg_beta_cd_llb_nl  ! Low-level obstacle drag coefficient Ridge scheme.
    logical, intent(in)             :: trpd_leewv_rdg_beta_nl
    logical, intent(in)             :: use_gw_rdg_gamma_nl
    integer, intent(in)             :: n_rdg_gamma_nl
    real(kind_phys), intent(in)             :: effgw_rdg_gamma_nl
    real(kind_phys), intent(in)             :: effgw_rdg_gamma_max_nl
    real(kind_phys), intent(in)             :: rdg_gamma_cd_llb_nl
    logical, intent(in)             :: trpd_leewv_rdg_gamma_nl
    character(len=256), intent(in)    :: bnd_topo_nl ! full pathname for topo file
    character(len=256), intent(in)    :: bnd_rdggm_nl ! full pathname for meso-Gamma ridge dataset
    ! Factor for SH orographic waves.
    real(kind_phys), intent(in)             :: gw_oro_south_fac_nl
    ! Whether or not to limit tau *before* applying any efficiency factors.
    logical, intent(in)             :: gw_limit_tau_without_eff_nl
    logical, intent(in)              :: gw_lndscl_sgh_nl
    real(kind_phys), intent(in)             :: gw_prndl_nl
    ! Whether or not to apply tendency max
    logical, intent(in)             :: gw_apply_tndmax_nl
    real(kind_phys), intent(in)             :: gw_qbo_hdepth_scaling_nl
    logical, intent(in)             :: gw_top_taper_nl
    ! Width of gaussian used to create frontogenesis tau profile [m s-1].
    real(kind_phys), intent(in)             :: front_gaussian_width_nl
    real(kind_phys), intent(in)             :: alpha_gw_movmtn_nl
    logical, intent(in)                     :: use_gw_rdg_resid_in

    real(kind_phys), intent(in)             :: effgw_rdg_resid_in
    real(kind_phys), intent(in)             :: effgw_movmtn_pbl_in
    integer, intent(in)                     :: movmtn_source_in
    real(kind_phys), intent(in)             :: movmtn_psteer_in
    real(kind_phys), intent(in)             :: movmtn_plaunch_in

    logical, intent(in) :: gw_rdg_do_divstream_nl, &
                           gw_rdg_do_smooth_regimes_nl, &
                           gw_rdg_do_adjust_tauoro_nl, &
                           gw_rdg_do_backward_compat_nl, &
                           gw_rdg_do_vdiff_nl
    real(kind_phys), intent(in) :: &
      gw_rdg_C_BetaMax_DS_nl, gw_rdg_C_GammaMax_nl, &
      gw_rdg_Frx0_nl, gw_rdg_Frx1_nl, gw_rdg_C_BetaMax_SM_nl, gw_rdg_Fr_c_nl, &
      gw_rdg_orohmin_nl, gw_rdg_orovmin_nl, gw_rdg_orostratmin_nl, gw_rdg_orom2min_nl

    logical, intent(in)             ::  use_gw_oro_in
    logical, intent(in)             ::  use_gw_front_in
    logical, intent(in)             ::  use_gw_front_igw_in
    logical, intent(in)             ::  use_gw_convect_dp_in
    logical, intent(in)             ::  use_gw_convect_sh_in
    logical, intent(in)             ::  use_simple_phys_in
    logical, intent(in)             ::  use_gw_movmtn_pbl_in
    integer, intent(in)             :: iulog_in
    integer, intent(in)             :: ktop_in
    logical, intent(in)             :: masterproc_in
    integer, intent(in)             :: nbot_molec_in
    logical, intent(in)             :: do_molec_diff_in
    real(kind_phys), intent(in)     :: wavelength_mid_in
    real(kind_phys), intent(in)     :: wavelength_long_in

    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errflg

    ! Local variables
    character(len=*), parameter :: sub = 'gw_drag_init'

    integer :: k, l
    character(len=128) :: errstring

    ! Index for levels at specific pressures.
    integer :: topndx
    integer :: botndx
    integer :: stat

    ! Interpolated Newtonian cooling coefficients.
    real(kind_phys) :: alpha(pver + 1)

    ! Levels of pre-calculated Newtonian cooling (1/day).
    ! The following profile is digitized from:
    ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5

    integer, parameter :: nalph = 71
    real(kind_phys) :: alpha0(nalph) = [ &
                       0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, &
                       0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, &
                       0.1_kind_phys, 0.1_kind_phys, 0.10133333_kind_phys, 0.104_kind_phys, &
                       0.108_kind_phys, 0.112_kind_phys, 0.116_kind_phys, 0.12066667_kind_phys, &
                       0.126_kind_phys, 0.132_kind_phys, 0.138_kind_phys, 0.144_kind_phys, &
                       0.15133333_kind_phys, 0.16_kind_phys, 0.17_kind_phys, 0.18_kind_phys, &
                       0.19_kind_phys, 0.19933333_kind_phys, 0.208_kind_phys, 0.216_kind_phys, &
                       0.224_kind_phys, 0.232_kind_phys, 0.23466667_kind_phys, 0.232_kind_phys, &
                       0.224_kind_phys, 0.216_kind_phys, 0.208_kind_phys, 0.20133333_kind_phys, &
                       0.196_kind_phys, 0.192_kind_phys, 0.188_kind_phys, 0.184_kind_phys, &
                       0.18266667_kind_phys, 0.184_kind_phys, 0.188_kind_phys, 0.192_kind_phys, &
                       0.196_kind_phys, 0.19333333_kind_phys, 0.184_kind_phys, 0.168_kind_phys, &
                       0.152_kind_phys, 0.136_kind_phys, 0.12133333_kind_phys, 0.108_kind_phys, &
                       0.096_kind_phys, 0.084_kind_phys, 0.072_kind_phys, 0.061_kind_phys, &
                       0.051_kind_phys, 0.042_kind_phys, 0.033_kind_phys, 0.024_kind_phys, &
                       0.017666667_kind_phys, 0.014_kind_phys, 0.013_kind_phys, 0.012_kind_phys, &
                       0.011_kind_phys, 0.010333333_kind_phys, 0.01_kind_phys, 0.01_kind_phys, &
                       0.01_kind_phys, 0.01_kind_phys, 0.01_kind_phys &
                       ]

    ! Pressure levels that were used to calculate alpha0 (hPa).
    real(kind_phys) :: palph(nalph) = [ &
                       2.06115E-06_kind_phys, 2.74280E-06_kind_phys, 3.64988E-06_kind_phys, 4.85694E-06_kind_phys, &
                       6.46319E-06_kind_phys, 8.60065E-06_kind_phys, 1.14450E-05_kind_phys, 1.52300E-05_kind_phys, &
                       2.02667E-05_kind_phys, 2.69692E-05_kind_phys, 3.58882E-05_kind_phys, 4.77568E-05_kind_phys, &
                       6.35507E-05_kind_phys, 8.45676E-05_kind_phys, 0.000112535_kind_phys, 0.000149752_kind_phys, &
                       0.000199277_kind_phys, 0.000265180_kind_phys, 0.000352878_kind_phys, 0.000469579_kind_phys, &
                       0.000624875_kind_phys, 0.000831529_kind_phys, 0.00110653_kind_phys, 0.00147247_kind_phys, &
                       0.00195943_kind_phys, 0.00260744_kind_phys, 0.00346975_kind_phys, 0.00461724_kind_phys, &
                       0.00614421_kind_phys, 0.00817618_kind_phys, 0.0108801_kind_phys, 0.0144783_kind_phys, &
                       0.0192665_kind_phys, 0.0256382_kind_phys, 0.0341170_kind_phys, 0.0453999_kind_phys, &
                       0.0604142_kind_phys, 0.0803939_kind_phys, 0.106981_kind_phys, 0.142361_kind_phys, &
                       0.189442_kind_phys, 0.252093_kind_phys, 0.335463_kind_phys, 0.446404_kind_phys, &
                       0.594036_kind_phys, 0.790490_kind_phys, 1.05192_kind_phys, 1.39980_kind_phys, &
                       1.86273_kind_phys, 2.47875_kind_phys, 3.29851_kind_phys, 4.38936_kind_phys, &
                       5.84098_kind_phys, 7.77266_kind_phys, 10.3432_kind_phys, 13.7638_kind_phys, &
                       18.3156_kind_phys, 24.3728_kind_phys, 32.4332_kind_phys, 43.1593_kind_phys, &
                       57.4326_kind_phys, 76.4263_kind_phys, 101.701_kind_phys, 135.335_kind_phys, &
                       180.092_kind_phys, 239.651_kind_phys, 318.907_kind_phys, 424.373_kind_phys, &
                       564.718_kind_phys, 751.477_kind_phys, 1000._kind_phys &
                       ]

    !-----------------------------------------------------------------------
    errmsg = ''
    errflg = 0

    gravit = gravit_in
    rair = rair_in
    pi = pi_in
    allocate (pref_edge(pver + 1))
    pref_edge = pref_edge_in
    allocate (pref_mid(pver))
    pref_mid = pref_mid_in
    masterproc = masterproc_in
    iulog = iulog_in
    pgwv = pgwv_nl
    gw_dc = gw_dc_nl
    pgwv_long = pgwv_long_nl
    gw_dc_long = gw_dc_long_nl
    tau_0_ubc = tau_0_ubc_nl
    effgw_beres_dp = effgw_beres_dp_nl
    effgw_beres_sh = effgw_beres_sh_nl
    effgw_cm = effgw_cm_nl
    effgw_cm_igw = effgw_cm_igw_nl
    effgw_oro = effgw_oro_nl
    fcrit2 = fcrit2_in
    rearth = rearth_in
    frontgfc = frontgfc_nl
    gw_drag_file = trim(gw_drag_file_nl)
    gw_drag_file_sh = trim(gw_drag_file_sh_nl)
    gw_drag_file_mm = trim(gw_drag_file_mm_nl)
    taubgnd = taubgnd_nl
    taubgnd_igw = taubgnd_igw_nl
    gw_polar_taper = gw_polar_taper_nl
    use_gw_rdg_beta = use_gw_rdg_beta_nl
    n_rdg_beta = n_rdg_beta_nl
    effgw_rdg_beta = effgw_rdg_beta_nl
    effgw_rdg_beta_max = effgw_rdg_beta_max_nl
    rdg_beta_cd_llb = rdg_beta_cd_llb_nl
    trpd_leewv_rdg_beta = trpd_leewv_rdg_beta_nl
    use_gw_rdg_gamma = use_gw_rdg_gamma_nl
    n_rdg_gamma = n_rdg_gamma_nl
    effgw_rdg_gamma = effgw_rdg_gamma_nl
    effgw_rdg_gamma_max = effgw_rdg_gamma_max_nl
    rdg_gamma_cd_llb = rdg_gamma_cd_llb_nl
    trpd_leewv_rdg_gamma = trpd_leewv_rdg_gamma_nl
    bnd_topo = trim(bnd_topo_nl)
    bnd_rdggm = trim(bnd_rdggm_nl)
    gw_oro_south_fac = gw_oro_south_fac_nl
    gw_limit_tau_without_eff = gw_limit_tau_without_eff_nl
    gw_lndscl_sgh = gw_lndscl_sgh_nl
    gw_prndl = gw_prndl_nl
    gw_apply_tndmax = gw_apply_tndmax_nl
    gw_qbo_hdepth_scaling = gw_qbo_hdepth_scaling_nl
    gw_top_taper = gw_top_taper_nl
    front_gaussian_width = front_gaussian_width_nl
    alpha_gw_movmtn = alpha_gw_movmtn_nl

    use_gw_rdg_resid = use_gw_rdg_resid_in
    effgw_rdg_resid = effgw_rdg_resid_in
    effgw_movmtn_pbl = effgw_movmtn_pbl_in
    movmtn_source = movmtn_source_in
    movmtn_psteer = movmtn_psteer_in
    movmtn_plaunch = movmtn_plaunch_in

    use_gw_oro = use_gw_oro_in
    use_gw_front = use_gw_front_in
    use_gw_front_igw = use_gw_front_igw_in
    use_gw_convect_dp = use_gw_convect_dp_in
    use_gw_convect_sh = use_gw_convect_sh_in
    use_simple_phys = use_simple_phys_in
    use_gw_movmtn_pbl = use_gw_movmtn_pbl_in
    do_molec_diff = do_molec_diff_in
    nbot_molec = nbot_molec_in
    wavelength_mid = wavelength_mid_in
    wavelength_long = wavelength_long_in
    ktop = ktop_in

    band_oro = GWBand(0, gw_dc, 1.0_kind_phys, wavelength_mid)
    band_mid = GWBand(pgwv, gw_dc, 1.0_kind_phys, wavelength_mid)
    band_long = GWBand(pgwv_long, gw_dc_long, 1.0_kind_phys, wavelength_long)
    band_movmtn = GWBand(0, gw_dc, 1.0_kind_phys, wavelength_mid)

    if (masterproc) then
      write (iulog, *) ' '
      write (iulog, *) "GW_DRAG: band_mid%ngwv = ", band_mid%ngwv
      do l = -band_mid%ngwv, band_mid%ngwv
        write (iulog, '(A,I0,A,F7.2)') &
          "GW_DRAG: band_mid%cref(", l, ") = ", band_mid%cref(l)
      end do
      write (iulog, *) 'GW_DRAG: band_mid%kwv = ', band_mid%kwv
      write (iulog, *) 'GW_DRAG: band_mid%fcrit2 = ', band_mid%fcrit2
      write (iulog, *) ' '
      write (iulog, *) "GW_DRAG: band_long%ngwv = ", band_long%ngwv
      do l = -band_long%ngwv, band_long%ngwv
        write (iulog, '(A,I2,A,F7.2)') &
          "GW_DRAG: band_long%cref(", l, ") = ", band_long%cref(l)
      end do
      write (iulog, *) 'GW_DRAG: band_long%kwv = ', band_long%kwv
      write (iulog, *) 'GW_DRAG: band_long%fcrit2 = ', band_long%fcrit2
      write (iulog, *) ' '
    end if

    ! pre-calculated newtonian damping:
    !     * convert to s-1
    !     * ensure it is not smaller than 1e-6
    !     * convert palph from hpa to pa

    do k = 1, nalph
      alpha0(k) = alpha0(k)/86400._kind_phys
      alpha0(k) = max(alpha0(k), 1.e-6_kind_phys)
      palph(k) = palph(k)*1.e2_kind_phys
    end do

    call lininterp(alpha0, palph, nalph, alpha, pref_edge, pver + 1)
    if (masterproc) then
      write (iulog, *) 'gw_init: newtonian damping (1/day):'
      write (iulog, fmt='(a4,a12,a10)') ' k  ', '  pref_edge      ', &
        '  alpha   '
      do k = 1, pver + 1
        write (iulog, fmt='(i4,1e12.5,1f10.2)') k, pref_edge(k), &
          alpha(k)*86400._kind_phys
      end do
    end if

    if (masterproc) then
      write (iulog, *) 'gw_init: ktop = ', ktop
    end if

    ! Initialize subordinate modules.
    call gw_common_init(pver, &
                        tau_0_ubc, ktop, gravit, rair, alpha, gw_prndl, &
                        gw_qbo_hdepth_scaling, errstring)

    call gw_rdg_init( &
      ncol=ncol, &
      band=band_oro, &
      rearth_c=rearth, &
      effgw_rdg_beta=effgw_rdg_beta, &
      effgw_rdg_gamma=effgw_rdg_gamma, &
      use_gw_rdg_beta_in=use_gw_rdg_beta, &
      use_gw_rdg_gamma_in=use_gw_rdg_gamma, &
      bnd_topo_file_in=bnd_topo, &
      bnd_rdg_file_in=bnd_rdggm, &
      gw_rdg_do_divstream_nl=gw_rdg_do_divstream_nl, &
      gw_rdg_C_BetaMax_DS_nl=gw_rdg_C_BetaMax_DS_nl, &
      gw_rdg_C_GammaMax_nl=gw_rdg_C_GammaMax_nl, &
      gw_rdg_Frx0_nl=gw_rdg_Frx0_nl, &
      gw_rdg_Frx1_nl=gw_rdg_Frx1_nl, &
      gw_rdg_C_BetaMax_SM_nl=gw_rdg_C_BetaMax_SM_nl, &
      gw_rdg_Fr_c_nl=gw_rdg_Fr_c_nl, &
      gw_rdg_do_smooth_regimes_nl=gw_rdg_do_smooth_regimes_nl, &
      gw_rdg_do_adjust_tauoro_nl=gw_rdg_do_adjust_tauoro_nl, &
      gw_rdg_do_backward_compat_nl=gw_rdg_do_backward_compat_nl, &
      gw_rdg_orohmin_nl=gw_rdg_orohmin_nl, &
      gw_rdg_orovmin_nl=gw_rdg_orovmin_nl, &
      gw_rdg_orostratmin_nl=gw_rdg_orostratmin_nl, &
      gw_rdg_orom2min_nl=gw_rdg_orom2min_nl, &
      gw_rdg_do_vdiff_nl=gw_rdg_do_vdiff_nl, &
      masterproc=masterproc, &
      iulog=iulog, &
      errmsg=errmsg, &
      errflg=errflg)
    if (errflg /= 0) return

    call gw_front_init(pver, pref_edge, frontgfc, band_mid, band_long, &
                       taubgnd, taubgnd_igw, &
                       effgw_cm, effgw_cm_igw, use_gw_front, use_gw_front_igw, &
                       front_gaussian_width, masterproc, iulog, errmsg, errflg)
    if (errflg /= 0) return

    call gw_movmtn_init(pver, gw_drag_file_mm, &
                        band_movmtn, &
                        pref_edge, movmtn_psteer, movmtn_plaunch, movmtn_source_in, masterproc, iulog, errmsg, errflg)
    if (errflg /= 0) return

    if (use_gw_convect_dp .or. use_gw_convect_sh) then
      call gw_beres_init(pver, pi, gw_drag_file_sh, gw_drag_file, pref_edge, gw_dc, wavelength_mid, pgwv, &
                         use_gw_convect_dp, use_gw_convect_sh, masterproc, iulog, errmsg, errflg)
    end if
    if (errflg /= 0) return

    if (gw_top_taper) then
      allocate (vramp(pver))
      vramp(:) = 1._kind_phys
      topndx = 1
      botndx = press_lim_idx(0.6E-02_kind_phys, top=.true.)
      if (botndx > 1) then
        do k = botndx, topndx, -1
          vramp(k) = vramp(k + 1)/(pref_edge(k + 1)/pref_edge(k))
        end do
        if (masterproc) then
          write (iulog, '(A)') 'GW taper coef (vramp):'
          do k = 1, pver
            write (iulog, "('k: ',I4,' taper coef,press(Pa): ',F12.8,E12.4)") k, vramp(k), pref_mid(k)
          end do
        end if
      end if
    end if
  end subroutine gw_drag_init

!==========================================================================

!> \section arg_table_gw_drag_run Argument Table
!! \htmlinclude gw_drag_run.html
  subroutine gw_drag_run( &
    ncol, &
    pcnst, &
    pver, &
    cnst_type, &
    dt, &
    cpair, &
    cpairv, &
    pi, &
    frontgf, &
    frontga, &
    degree2radian, &
    al0, &
    dlat0, &
    pint, &
    piln, &
    pdel, &
    pdeldry, &
    zm, &
    zi, &
    lat, &
    landfrac, &
    dse, &
    state_t, &
    state_u, &
    state_v, &
    state_q, &
    vorticity, &
    sgh, &
    kvtt, &
    ttend_dp, &
    ttend_sh, &
    ttend_clubb, &
    thlp2_clubb_gw, &
    wpthlp_clubb_gw, &
    upwp_clubb_gw, &
    vpwp_clubb_gw, &
    s_tend, &
    q_tend, &
    u_tend, &
    v_tend, &
    scheme_name, &
    nbot_molec, &
    egwdffi_tot, &
    flx_heat, &
    errmsg, errflg)
    !-----------------------------------------------------------------------
    ! Interface for multiple gravity wave drag parameterization.
    !-----------------------------------------------------------------------

    use coords_1d, only: Coords1D
    use gw_common, only: gw_prof, gw_drag_prof, calc_taucd
    use gw_common, only: momentum_flux, momentum_fixer, energy_change
    use gw_common, only: energy_fixer, coriolis_speed, adjust_inertial
    use gw_oro, only: gw_oro_src
    use gw_front, only: gw_cm_src, cm_desc, cm_igw_desc
    use gw_convect, only: gw_beres_src, beres_dp_desc, beres_sh_desc
    use gw_movmtn, only: gw_movmtn_run
    use gw_rdg, only: gw_rdg_run

    integer, intent(in)        :: ncol  ! number of atmospheric columns
    integer, intent(in)        :: pcnst ! chunk number
    integer, intent(in)        :: pver  ! number of atmospheric levels
    character*3, intent(in)     :: cnst_type(pcnst) ! wet or dry mixing ratio
    real(kind_phys), intent(in) :: dt          ! physics timestep
    real(kind_phys), intent(in) :: cpair       ! heat capacity of air
    real(kind_phys), intent(in) :: cpairv(:, :) ! location dependent heat capacity of air
    real(kind_phys), intent(in) :: pi          ! pi
    real(kind_phys), intent(in) :: degree2radian
    real(kind_phys), intent(in) :: al0
    real(kind_phys), intent(in) :: dlat0
    real(kind_phys), intent(in) :: pint(:, :)   ! pressure at model interfaces
    real(kind_phys), intent(in) :: piln(:, :)   ! ln pressure at model interfaces
    real(kind_phys), intent(in) :: pdel(:, :)   ! vertical delta-p
    real(kind_phys), intent(in) :: pdeldry(:, :)   ! vertical delta-p
    real(kind_phys), intent(in) :: zm(:, :)
    real(kind_phys), intent(in) :: zi(:, :)
    real(kind_phys), intent(in) :: lat(:)
    real(kind_phys), intent(in) :: landfrac(:)
    real(kind_phys), intent(in) :: dse(:, :)       ! dry static energy
    real(kind_phys), intent(in) :: state_t(:, :)   ! temperature (K)
    real(kind_phys), intent(in) :: state_u(:, :)   ! meridional wind
    real(kind_phys), intent(in) :: state_v(:, :)   ! zonal wind
    real(kind_phys), intent(in) :: state_q(:, :, :) ! constituent array
    real(kind_phys), intent(in) :: vorticity(:, :) ! vorticity
    real(kind_phys), intent(in) :: sgh(:)         !
    real(kind_phys), intent(inout) :: kvtt(:, :)       !
    real(kind_phys), intent(in) :: ttend_dp(:, :)  ! Temperature change due to deep convection.
    real(kind_phys), intent(in) :: ttend_sh(:, :)  ! Temperature change due to shallow convection.
    real(kind_phys), intent(in) :: ttend_clubb(:, :)
    real(kind_phys), intent(in) :: thlp2_clubb_gw(:, :)
    real(kind_phys), intent(in) :: wpthlp_clubb_gw(:, :)
    real(kind_phys), intent(in) :: upwp_clubb_gw(:, :)
    real(kind_phys), intent(in) :: vpwp_clubb_gw(:, :)
    real(kind_phys), intent(inout):: s_tend(:, :)   ! dry air enthalpy tendency
    real(kind_phys), intent(inout):: q_tend(:, :, :)
    real(kind_phys), intent(inout):: u_tend(:, :)
    real(kind_phys), intent(inout):: v_tend(:, :)
    character(len=64), intent(out) :: scheme_name
    integer, intent(in)             :: nbot_molec
    ! Parameterization net tendencies.
    ! sum from the two types of spectral GW
    real(kind_phys), intent(out) :: egwdffi_tot(:, :)
    real(kind_phys), intent(out) :: flx_heat(:)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    !---------------------------Local storage-------------------------------
    character(len=*), parameter :: sub = 'gw_drag_run'

    integer :: stat

    integer :: i, k                   ! loop indices

    type(Coords1D) :: p               ! Pressure coordinates

    real(kind_phys) :: ttgw(ncol, pver) ! temperature tendency
    real(kind_phys) :: utgw(ncol, pver) ! zonal wind tendency
    real(kind_phys) :: vtgw(ncol, pver) ! meridional wind tendency

    real(kind_phys) :: ni(ncol, pver + 1) ! interface Brunt-Vaisalla frequency
    real(kind_phys) :: nm(ncol, pver)   ! midpoint Brunt-Vaisalla frequency
    real(kind_phys) :: rhoi(ncol, pver + 1)     ! interface density
    real(kind_phys), allocatable :: tau(:, :, :)  ! wave Reynolds stress
    real(kind_phys) :: tau0x(ncol)     ! c=0 sfc. stress (zonal)
    real(kind_phys) :: tau0y(ncol)     ! c=0 sfc. stress (meridional)
    real(kind_phys) :: ubi(ncol, pver + 1)! projection of wind at interfaces
    real(kind_phys) :: ubm(ncol, pver)  ! projection of wind at midpoints
    real(kind_phys) :: xv(ncol)        ! unit vector of source wind (x)
    real(kind_phys) :: yv(ncol)        ! unit vector of source wind (y)

    integer :: m                      ! dummy integers
    real(kind_phys) :: qtgw(ncol, pver, pcnst) ! constituents tendencies

    ! Reynolds stress for waves propagating in each cardinal direction.
    real(kind_phys) :: taucd(ncol, pver + 1, 4)

    ! gravity wave wind tendency for each wave
    real(kind_phys), allocatable :: gwut(:, :, :)

    ! Temperature tendencies from diffusion and kinetic energy.
    real(kind_phys) :: dttdf(ncol, pver)
    real(kind_phys) :: dttke(ncol, pver)

    ! Wave phase speeds for each column
    real(kind_phys), allocatable :: phase_speeds(:, :)

    ! Efficiency for a gravity wave source.
    real(kind_phys) :: effgw(ncol)

    ! Coriolis characteristic speed.
    real(kind_phys) :: u_coriolis(ncol)

    ! Adjustment for inertial gravity waves.
    real(kind_phys), allocatable :: ro_adjust(:, :, :)

    ! Frontogenesis
    real(kind_phys), pointer :: frontgf(:, :)
    real(kind_phys), pointer :: frontga(:, :)

    ! gridbox area
    real(kind_phys), pointer :: gbxar(:)

    ! Beta ridges
    ! width of ridges.
    real(kind_phys), pointer :: hwdth(:, :)
    ! length of ridges.
    real(kind_phys), pointer :: clngt(:, :)
    ! Maximum deviations of ridges.
    real(kind_phys), pointer :: mxdis(:, :)
    ! orientation of ridges.
    real(kind_phys), pointer :: angll(:, :)
    ! anisotropy of ridges.
    real(kind_phys), pointer :: anixy(:, :)

    ! Gamma ridges
    ! width of ridges.
    real(kind_phys), pointer :: hwdthg(:, :)
    ! length of ridges.
    real(kind_phys), pointer :: clngtg(:, :)
    ! Maximum deviations of ridges.
    real(kind_phys), pointer :: mxdisg(:, :)
    ! orientation of ridges.
    real(kind_phys), pointer :: angllg(:, :)
    ! anisotropy of ridges.
    real(kind_phys), pointer :: anixyg(:, :)

    ! Indices of gravity wave source and lowest level where wind tendencies
    ! are allowed.
    integer :: src_level(ncol)
    integer :: tend_level(ncol)

    ! Convective source heating depth.
    ! heating depth
    real(kind_phys) :: hdepth(ncol)
    ! maximum heating rate
    real(kind_phys) :: maxq0(ncol)

    ! Scale sgh to account for landfrac.
    real(kind_phys) :: sgh_scaled(ncol)

    ! Parameters for the IGW polar taper.
!!$  real(kind_phys), parameter :: degree2radian = pi/180._kind_phys
!!$  real(kind_phys), parameter :: al0 = 82.5_kind_phys * degree2radian
!!$  real(kind_phys), parameter :: dlat0 = 5.0_kind_phys * degree2radian

    ! effective gw diffusivity at interfaces needed for output
    real(kind_phys) :: egwdffi(ncol, pver + 1)

    ! Momentum fluxes used by fixer.
    real(kind_phys) :: um_flux(ncol), vm_flux(ncol)
    ! Energy change used by fixer.
    real(kind_phys) :: de(ncol)

    ! Which constituents are being affected by diffusion.
    logical  :: lq(pcnst)

    !------------------------------------------------------------------------
    p = Coords1D(pint(:ncol, :))

    ! Profiles of background state variables
    call gw_prof(ncol, p, cpair, state_t, rhoi, nm, ni)

    if (do_molec_diff) then
      !--------------------------------------------------------
      ! Initialize and calculate local molecular diffusivity
      !--------------------------------------------------------

!jt     kvt_in passed in here as kvtt
!!$     call pbuf_get_field(pbuf, kvt_idx, kvt_in)  ! kvt_in(1:pcols,1:pver+1)

!!$     ! Set kvtt from pbuf field; kvtt still needs a factor of 1/cpairv.
!!$     kvtt = kvt_in(:ncol,:)

      ! Use linear extrapolation of cpairv to top interface.
      kvtt(:, 1) = kvtt(:, 1)/ &
                   (1.5_kind_phys*cpairv(:ncol, 1) - &
                    0.5_kind_phys*cpairv(:ncol, 2))

      ! Interpolate cpairv to other interfaces.
      do k = 2, nbot_molec
        kvtt(:, k) = kvtt(:, k)/ &
                     (cpairv(:ncol, k + 1) + cpairv(:ncol, k))*2._kind_phys
      end do

    else

      kvtt = 0._kind_phys

    end if

    if (use_gw_front_igw) then
      u_coriolis = coriolis_speed(band_long, lat(:ncol))
    end if

    ! Totals that accumulate over different sources.
    egwdffi_tot = 0._kind_phys
    flx_heat = 0._kind_phys

    !------------------------------------------------------------------
    ! Convective moving mountain gravity waves (Beres scheme).
    !------------------------------------------------------------------
    if (use_gw_movmtn_pbl) then
      effgw = effgw_movmtn_pbl

      call gw_movmtn_run( &
        ncol                = ncol, &
        band                = band_movmtn, &
        state_t             = state_t(:ncol,:), &
        pcnst               = pcnst, &
        state_u             = state_u(:ncol,:), &
        state_v             = state_v(:ncol,:), &
        p                   = p, &
        ttend_dp            = ttend_dp(:ncol,:), &
        ttend_clubb         = ttend_clubb(:ncol,:), &
        upwp_clubb          = upwp_clubb_gw(:ncol,:), &
        vpwp_clubb          = vpwp_clubb_gw(:ncol,:), &
        vorticity           = vorticity(:ncol,:), &
        zm                  = zm(:ncol,:), &
        alpha_gw_movmtn     = alpha_gw_movmtn, &
        dt                  = dt, &
        vramp               = vramp, &
        pint                = pint(:ncol,:), &
        piln                = piln(:ncol,:), &
        rhoi                = rhoi(:ncol,:), &
        nm                  = nm(:ncol,:), &
        ni                  = ni(:ncol,:), &
        effgw               = effgw(:ncol), &
        kvtt                = kvtt(:ncol,:), &
        state_q             = state_q(:ncol,:,:), &
        dse                 = dse(:ncol,:), &
        gw_apply_tndmax     = gw_apply_tndmax, &
        use_gw_movmtn_pbl   = use_gw_movmtn_pbl, &
        gravit              = gravit, &
        rair                = rair, &
        ! Input/output arguments
        src_level           = src_level(:ncol), &
        tend_level          = tend_level(:ncol), &
        ubm                 = ubm(:ncol,:pver), &
        ubi                 = ubi(:ncol,:pver+1), &
        xv                  = xv(:ncol), &
        yv                  = yv(:ncol), &
        hdepth              = hdepth(:ncol), &
        q_tend              = q_tend(:ncol,:pver,:pcnst), &
        u_tend              = u_tend(:ncol,:pver), &
        v_tend              = v_tend(:ncol,:pver), &
        s_tend              = s_tend(:ncol,:pver), &
        ! Output arguments
        utgw                = utgw(:ncol,:pver), &
        vtgw                = vtgw(:ncol,:pver), &
        ttgw                = ttgw(:ncol,:pver), &
        qtgw                = qtgw(:ncol,:pver,:pcnst), &
        egwdffi             = egwdffi(:ncol,:pver+1), &
        dttdf               = dttdf(:ncol,:pver), &
        dttke               = dttke(:ncol,:pver), &
        flx_heat            = flx_heat(:ncol), &
        errmsg              = errmsg, &
        errflg              = errflg)

      !  add the diffusion coefficients
      do k = 1, pver + 1
        egwdffi_tot(:, k) = egwdffi_tot(:, k) + egwdffi(:, k)
      end do
    end if

    !------------------------------------------------------------------
    ! Convective gravity waves (Beres scheme, deep).
    !------------------------------------------------------------------
    if (use_gw_convect_dp) then

!!$     call gw_beres_run(ncol, band_mid, beres_dp_desc, &
!!$          state_u, state_v, ttend_dp(:ncol,:), zm, src_level, tend_level, tau, &
!!$          ubm, ubi, xv, yv, phase_speeds, hdepth, maxq0, &
!!$          p, dt, frontgf, frontga, &
!!$          state_t, vramp,    &
!!$          piln, rhoi,       nm,   ni,   &
!!$          effgw_dp,  kvtt, state_q,  dse,  utgw,  vtgw, &
!!$          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
!!$          lapply_effgw_in=gw_apply_tndmax, flx_heat_curr)
!!$
!!$     !  add the diffusion coefficients
!!$     do k = 1, pver+1
!!$        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
!!$     end do

      ! Allocate wavenumber fields.
      allocate (tau(ncol, -band_mid%ngwv:band_mid%ngwv, pver + 1), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of tau failed', errmsg)
      allocate (gwut(ncol, pver, -band_mid%ngwv:band_mid%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of gwut failed', errmsg)
      allocate (phase_speeds(ncol, -band_mid%ngwv:band_mid%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of phase_speeds failed', errmsg)

      ! Efficiency of gravity wave momentum transfer.
      ! This is really only to remove the pole points.
      where (pi/2._kind_phys - abs(lat(:ncol)) >= 4*epsilon(1._kind_phys))
        effgw = effgw_beres_dp
      elsewhere
        effgw = 0._kind_phys
      end where

      ! Determine wave sources for Beres deep scheme
      call gw_beres_src( &
        ncol        = ncol, &
        desc        = beres_dp_desc, &
        u           = state_u(:ncol,:), &
        v           = state_v(:ncol,:), &
        netdt       = ttend_dp(:ncol,:), &
        zm          = zm(:ncol,:), &
        ! Output arguments
        src_level   = src_level(:ncol), &
        tend_level  = tend_level(:ncol), &
        tau         = tau(:ncol,-band_mid%ngwv:band_mid%ngwv,:pver+1), &
        ubm         = ubm(:ncol,:pver), &
        ubi         = ubi(:ncol,:pver+1), &
        xv          = xv(:ncol), &
        yv          = yv(:ncol), &
        c           = phase_speeds(:ncol,-band_mid%ngwv:band_mid%ngwv), &
        hdepth      = hdepth(:ncol), &
        maxq0       = maxq0(:ncol))

      ! Solve for the drag profile with Beres source spectrum.
      call gw_drag_prof( &
        ncol                = ncol, &
        band                = band_mid, &
        p                   = p, &
        src_level           = src_level, &
        tend_level          = tend_level, &
        dt                  = dt, &
        t                   = state_t(:ncol,:), &
        vramp               = vramp, &
        piln                = piln(:ncol,:), &
        rhoi                = rhoi(:ncol,:), &
        nm                  = nm(:ncol,:), &
        ni                  = ni(:ncol,:), &
        ubm                 = ubm(:ncol,:), &
        ubi                 = ubi(:ncol,:), &
        xv                  = xv(:ncol), &
        yv                  = yv(:ncol), &
        effgw               = effgw(:ncol), &
        c                   = phase_speeds(:ncol,-band_mid%ngwv:band_mid%ngwv), &
        kvtt                = kvtt(:ncol,:), &
        q                   = state_q(:ncol,:,:), &
        dse                 = dse(:ncol,:), &
        lapply_effgw_in     = gw_apply_tndmax, &
        ! Input/output arguments
        tau                 = tau(:ncol,-band_mid%ngwv:band_mid%ngwv,:pver+1), &
        ! Output arguments
        utgw                = utgw(:ncol,:pver), &
        vtgw                = vtgw(:ncol,:pver), &
        ttgw                = ttgw(:ncol,:pver), &
        qtgw                = qtgw(:ncol,:pver,:pcnst), &
        egwdffi             = egwdffi(:ncol,:pver+1), &
        gwut                = gwut(:ncol,:pver,-band_mid%ngwv:band_mid%ngwv), &
        dttdf               = dttdf(:ncol,:pver), &
        dttke               = dttke(:ncol,:pver))

      ! Project stress into directional components.
      taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

      !  add the diffusion coefficients
      do k = 1, pver + 1
        egwdffi_tot(:, k) = egwdffi_tot(:, k) + egwdffi(:, k)
      end do

      ! Store constituents tendencies
      do m = 1, pcnst
        do k = 1, pver
          q_tend(:ncol, k, m) = q_tend(:ncol, k, m) + qtgw(:, k, m)
        end do
      end do

      ! Find momentum flux, and use it to fix the wind tendencies below
      ! the gravity wave region.
      call momentum_flux(tend_level, taucd, um_flux, vm_flux)
      call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

      ! Add the momentum tendencies to the output tendency arrays.
      do k = 1, pver
        u_tend(:ncol, k) = u_tend(:ncol, k) + utgw(:, k)
        v_tend(:ncol, k) = v_tend(:ncol, k) + vtgw(:, k)
      end do

      ! Find energy change in the current state, and use fixer to apply
      ! the difference in lower levels.
      call energy_change(dt, p, state_u, state_v, u_tend(:ncol, :), &
                         v_tend(:ncol, :), s_tend(:ncol, :) + ttgw, de)
      call energy_fixer(tend_level, p, de - flx_heat(:ncol), ttgw)

      do k = 1, pver
        s_tend(:ncol, k) = s_tend(:ncol, k) + ttgw(:, k)
      end do

      ! Change ttgw to a temperature tendency before outputing it.
!!$     ttgw = ttgw / cpair
!!$     call gw_spec_outflds(beres_dp_pf, ncol, pver, band_mid, phase_speeds, u, v, &
!!$          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
!!$          taucd)

!!$     ! Diagnostic outputs (convert hdepth to km).
!!$     call outfld('NETDT', ttend_dp, pcols, lchnk)
!!$     call outfld('HDEPTH', hdepth/1000._kind_phys, ncol, lchnk)
!!$     call outfld('MAXQ0', maxq0, ncol, lchnk)

      deallocate (tau, gwut, phase_speeds)

    end if

    !------------------------------------------------------------------
    ! Convective gravity waves (Beres scheme, shallow).
    !------------------------------------------------------------------
    if (use_gw_convect_sh) then
!!$
!!$     call gw_beres_run(ncol, band_mid, beres_sh_desc, state_u, state_v, &
!!$     ttend_sh, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
!!$     c, hdepth, maxq0 ,effgw_sh)
!!$
!!$     !  add the diffusion coefficients
!!$     do k = 1, pver+1
!!$        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
!!$     end do

      ! Allocate wavenumber fields.
      allocate (tau(ncol, -band_mid%ngwv:band_mid%ngwv, pver + 1), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of tau failed', errmsg)
      allocate (gwut(ncol, pver, -band_mid%ngwv:band_mid%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of gwut failed', errmsg)
      allocate (phase_speeds(ncol, -band_mid%ngwv:band_mid%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of phase_speeds failed', errmsg)

      ! Efficiency of gravity wave momentum transfer.
      ! This is really only to remove the pole points.
      where (pi/2._kind_phys - abs(lat(:ncol)) >= 4*epsilon(1._kind_phys))
        effgw = effgw_beres_sh
      elsewhere
        effgw = 0._kind_phys
      end where

      ! Determine wave sources for Beres shallow scheme
      call gw_beres_src(ncol, &
                        beres_sh_desc, &
                        state_u, state_v, ttend_sh(:ncol, :), zm, src_level, tend_level, tau, &
                        ubm, ubi, xv, yv, phase_speeds, hdepth, maxq0)

      ! Solve for the drag profile with Beres source spectrum.
      call gw_drag_prof(ncol, band_mid, p, src_level, tend_level, dt, &
                        state_t, vramp, &
                        piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                        effgw, phase_speeds, kvtt, state_q, dse, tau, utgw, vtgw, &
                        ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                        lapply_effgw_in=gw_apply_tndmax)

      ! Project stress into directional components.
      taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

      !  add the diffusion coefficients
      do k = 1, pver + 1
        egwdffi_tot(:, k) = egwdffi_tot(:, k) + egwdffi(:, k)
      end do

      ! Store constituents tendencies
      do m = 1, pcnst
        do k = 1, pver
          q_tend(:ncol, k, m) = q_tend(:ncol, k, m) + qtgw(:, k, m)
        end do
      end do

      ! Add the momentum tendencies to the output tendency arrays.
      ! Don't calculate fixers, since we are too close to the ground to
      ! spread momentum/energy differences across low layers.
      do k = 1, pver
        u_tend(:ncol, k) = u_tend(:ncol, k) + utgw(:, k)
        v_tend(:ncol, k) = v_tend(:ncol, k) + vtgw(:, k)
        s_tend(:ncol, k) = s_tend(:ncol, k) + ttgw(:, k)
      end do

      ! Calculate energy change for output to CAM's energy checker.
      ! This is sort of cheating; we don't have a good a priori idea of the
      ! energy coming from surface stress, so we just integrate what we and
      ! actually have so far and overwrite flx_heat with that.
      call energy_change(dt, p, state_u, state_v, u_tend(:ncol, :), &
                         v_tend(:ncol, :), s_tend(:ncol, :), de)
      flx_heat(:ncol) = de

      ! Change ttgw to a temperature tendency before outputing it.
!!$     ttgw = ttgw / cpair
!!$     call gw_spec_outflds(beres_sh_pf, ncol, pver, band_mid, phase_speeds, u, v, &
!!$          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
!!$          taucd)

!!$     ! Diagnostic outputs (convert SHDEPTH to km).
!!$     call outfld ('SNETDT', ttend_sh, pcols, lchnk)
!!$     call outfld ('SHDEPTH', hdepth/1000._kind_phys, ncol, lchnk)
!!$     call outfld ('SMAXQ0', maxq0, ncol, lchnk)

      deallocate (tau, gwut, phase_speeds)

    end if

    !------------------------------------------------------------------
    ! Frontally generated gravity waves
    !------------------------------------------------------------------
    if (use_gw_front) then
      ! Allocate wavenumber fields.
      allocate (tau(ncol, -band_mid%ngwv:band_mid%ngwv, pver + 1), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of tau failed', errmsg)
      allocate (gwut(ncol, pver, -band_mid%ngwv:band_mid%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of gwut failed', errmsg)
      allocate (phase_speeds(ncol, -band_mid%ngwv:band_mid%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of phase_speeds failed', errmsg)

      ! Efficiency of gravity wave momentum transfer.
      effgw = effgw_cm
      ! Frontogenesis is too high at the poles (at least for the FV
      ! dycore), so introduce a polar taper.
      if (gw_polar_taper) effgw = effgw*cos(lat(:ncol))

      call gw_cm_src(ncol, band_mid, &
                     cm_desc, &
                     state_u, state_v, frontgf(:ncol, :), &
                     src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds)

      ! Solve for the drag profile with C&M source spectrum.
      call gw_drag_prof(ncol, band_mid, p, src_level, tend_level, dt, &
                        state_t, vramp, &
                        piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                        effgw, phase_speeds, kvtt, state_q, dse, tau, utgw, vtgw, &
                        ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                        lapply_effgw_in=gw_apply_tndmax)

      ! Project stress into directional components.
      taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

      ! Add the diffusion coefficients
      do k = 1, pver + 1
        egwdffi_tot(:, k) = egwdffi_tot(:, k) + egwdffi(:, k)
      end do

      ! Add the constituent tendencies
      do m = 1, pcnst
        do k = 1, pver
          q_tend(:ncol, k, m) = q_tend(:ncol, k, m) + qtgw(:, k, m)
        end do
      end do

      ! Find momentum flux, and use it to fix the wind tendencies below
      ! the gravity wave region.
      call momentum_flux(tend_level, taucd, um_flux, vm_flux)
      call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

      ! add the momentum tendencies to the output tendency arrays
      do k = 1, pver
        u_tend(:ncol, k) = u_tend(:ncol, k) + utgw(:, k)
        v_tend(:ncol, k) = v_tend(:ncol, k) + vtgw(:, k)
      end do

      ! Find energy change in the current state, and use fixer to apply
      ! the difference in lower levels.
      call energy_change(dt, p, state_u, state_v, u_tend(:ncol, :), &
                         v_tend(:ncol, :), s_tend(:ncol, :) + ttgw, de)
      call energy_fixer(tend_level, p, de - flx_heat(:ncol), ttgw)

      do k = 1, pver
        s_tend(:ncol, k) = s_tend(:ncol, k) + ttgw(:, k)
      end do

      ! Change ttgw to a temperature tendency before outputing it.
!!$     ttgw = ttgw / cpair
!!$     call gw_spec_outflds(cm_pf, ncol, pver, band_mid, phase_speeds, u, v, &
!!$          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
!!$          taucd)
!!$
!!$     deallocate(tau, gwut, phase_speeds)

    end if

    !------------------------------------------------------------------
    ! Frontally generated inertial gravity waves
    !------------------------------------------------------------------
    if (use_gw_front_igw) then
      ! Allocate wavenumber fields.
      allocate (tau(ncol, -band_long%ngwv:band_long%ngwv, pver + 1), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of tau failed', errmsg)
      allocate (gwut(ncol, pver, -band_long%ngwv:band_long%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of gwut failed', errmsg)
      allocate (phase_speeds(ncol, -band_long%ngwv:band_long%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of phase_speeds failed', errmsg)
      allocate (ro_adjust(ncol, -band_long%ngwv:band_long%ngwv, pver + 1), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of ro_adjust failed', errmsg)

      ! Efficiency of gravity wave momentum transfer.
      effgw = effgw_cm_igw

      ! Frontogenesis is too high at the poles (at least for the FV
      ! dycore), so introduce a polar taper.
      if (gw_polar_taper) then
        where (abs(lat(:ncol)) <= 89._kind_phys*degree2radian)
          effgw = effgw*0.25_kind_phys* &
                  (1._kind_phys + tanh((lat(:ncol) + al0)/dlat0))* &
                  (1._kind_phys - tanh((lat(:ncol) - al0)/dlat0))
        elsewhere
          effgw = 0._kind_phys
        end where
      end if

      call gw_cm_src(ncol, band_long, &
                     cm_igw_desc, &
                     state_u, state_v, frontgf(:ncol, :), &
                     src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds)

      call adjust_inertial(band_long, tend_level, u_coriolis, phase_speeds, ubi, &
                           tau, ro_adjust)

      ! Solve for the drag profile with C&M source spectrum.
      call gw_drag_prof(ncol, band_long, p, src_level, tend_level, dt, &
                        state_t, vramp, &
                        piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                        effgw, phase_speeds, kvtt, state_q, dse, tau, utgw, vtgw, &
                        ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                        ro_adjust=ro_adjust, lapply_effgw_in=gw_apply_tndmax)

      ! Project stress into directional components.
      taucd = calc_taucd(ncol, band_long%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

      !  add the diffusion coefficients
      do k = 1, pver + 1
        egwdffi_tot(:, k) = egwdffi_tot(:, k) + egwdffi(:, k)
      end do

      !Add the constituent tendencies
      do m = 1, pcnst
        do k = 1, pver
          q_tend(:ncol, k, m) = q_tend(:ncol, k, m) + qtgw(:, k, m)
        end do
      end do

      ! Find momentum flux, and use it to fix the wind tendencies below
      ! the gravity wave region.
      call momentum_flux(tend_level, taucd, um_flux, vm_flux)
      call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

      ! add the momentum tendencies to the output tendency arrays
      do k = 1, pver
        u_tend(:ncol, k) = u_tend(:ncol, k) + utgw(:, k)
        v_tend(:ncol, k) = v_tend(:ncol, k) + vtgw(:, k)
      end do

      ! Find energy change in the current state, and use fixer to apply
      ! the difference in lower levels.
      call energy_change(dt, p, state_u, state_v, u_tend(:ncol, :), &
                         v_tend(:ncol, :), s_tend(:ncol, :) + ttgw, de)
      call energy_fixer(tend_level, p, de - flx_heat(:ncol), ttgw)

      do k = 1, pver
        s_tend(:ncol, k) = s_tend(:ncol, k) + ttgw(:, k)
      end do

      ! Change ttgw to a temperature tendency before outputing it.
!!$     ttgw = ttgw / cpair
!!$     call gw_spec_outflds(cm_igw_pf, ncol, pver, band_long, phase_speeds, state_u, state_v, &
!!$          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
!!$          taucd)

      deallocate (tau, gwut, phase_speeds, ro_adjust)

    end if

    !---------------------------------------------------------------------
    ! Orographic stationary gravity waves
    !---------------------------------------------------------------------
    if (use_gw_oro) then

      ! Allocate wavenumber fields.
      allocate (tau(ncol, band_oro%ngwv:band_oro%ngwv, pver + 1), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of tau failed', errmsg)
      allocate (gwut(ncol, pver, band_oro%ngwv:band_oro%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of gwut failed', errmsg)
      allocate (phase_speeds(ncol, band_oro%ngwv:band_oro%ngwv), stat=stat)
      call handle_err(stat, errflg, sub//': Allocate of phase_speeds failed', errmsg)

      if (gw_lndscl_sgh) then
        where (landfrac(:ncol) >= epsilon(1._kind_phys))
          effgw = effgw_oro*landfrac(:ncol)
          sgh_scaled = sgh(:ncol)/sqrt(landfrac(:ncol))
        elsewhere
          effgw = 0._kind_phys
          sgh_scaled = 0._kind_phys
        end where

        ! Determine the orographic wave source
        call gw_oro_src(ncol, band_oro, p, &
                        state_u, state_v, state_t, sgh_scaled, zm, nm, &
                        src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds)
      else
        effgw = effgw_oro

        ! Determine the orographic wave source
        call gw_oro_src(ncol, band_oro, p, &
                        state_u, state_v, state_t, sgh(:ncol), zm, nm, &
                        src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds)
      end if
      do i = 1, ncol
        if (lat(i) < 0._kind_phys) then
          tau(i, :, :) = tau(i, :, :)*gw_oro_south_fac
        end if
      end do

      ! Solve for the drag profile with orographic sources.
      call gw_drag_prof(ncol, band_oro, p, src_level, tend_level, dt, &
                        state_t, vramp, &
                        piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                        effgw, phase_speeds, kvtt, state_q, dse, tau, utgw, vtgw, &
                        ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                        lapply_effgw_in=gw_apply_tndmax)

      ! For orographic waves, don't bother with taucd, since there are no
      ! momentum conservation routines or directional diagnostics.

      !  add the diffusion coefficients
      do k = 1, pver + 1
        egwdffi_tot(:, k) = egwdffi_tot(:, k) + egwdffi(:, k)
      end do

      ! Add the orographic tendencies to the spectrum tendencies.
      ! Don't calculate fixers, since we are too close to the ground to
      ! spread momentum/energy differences across low layers.
      do k = 1, pver
        u_tend(:ncol, k) = u_tend(:ncol, k) + utgw(:, k)
        v_tend(:ncol, k) = v_tend(:ncol, k) + vtgw(:, k)
        s_tend(:ncol, k) = s_tend(:ncol, k) + ttgw(:, k)
        ! Convert to temperature tendency for output.
        ttgw(:, k) = ttgw(:, k)/cpairv(:ncol, k)
      end do

      ! Calculate energy change for output to CAM's energy checker.
      ! This is sort of cheating; we don't have a good a priori idea of the
      ! energy coming from surface stress, so we just integrate what we and
      ! actually have so far and overwrite flx_heat with that.
      call energy_change(dt, p, state_u, state_v, u_tend(:ncol, :), &
                         v_tend(:ncol, :), s_tend(:ncol, :), de)
      flx_heat(:ncol) = de

      do m = 1, pcnst
        do k = 1, pver
          q_tend(:ncol, k, m) = q_tend(:ncol, k, m) + qtgw(:, k, m)
        end do
      end do

!!$     ! Write output fields to history file
!!$     call outfld('TAUAORO', tau(:,0,:),  ncol, lchnk)
!!$     call outfld('UTGWORO', utgw,  ncol, lchnk)
!!$     call outfld('VTGWORO', vtgw,  ncol, lchnk)
!!$     call outfld('TTGWORO', ttgw,  ncol, lchnk)
!!$     call outfld('TTGWSDFORO', dttdf / cpair,  ncol, lchnk)
!!$     call outfld('TTGWSKEORO', dttke / cpair,  ncol, lchnk)
!!$     tau0x = tau(:,0,pver+1) * xv
!!$     tau0y = tau(:,0,pver+1) * yv
!!$     call outfld('TAUGWX', tau0x, ncol, lchnk)
!!$     call outfld('TAUGWY', tau0y, ncol, lchnk)

      deallocate (tau, gwut, phase_speeds)

    end if

    if (use_gw_rdg_beta .or. use_gw_rdg_gamma) then
!!$     ! Save state at top of routine
!!$     ! Useful for unit testing checks
!!$     call outfld('UEGW', state_u ,  ncol, lchnk)
!!$     call outfld('VEGW', state_v ,  ncol, lchnk)
!!$     call outfld('TEGW', state_t ,  ncol, lchnk)
!!$     call outfld('ZEGW', zi , ncol, lchnk)
!!$     call outfld('ZMGW', zm , ncol, lchnk)

      call gw_rdg_run( &
        use_gw_rdg_beta, &
        use_gw_rdg_gamma, &
        vramp, &
        pcnst, pver, ncol, n_rdg_beta, n_rdg_gamma, dt, &
        state_u, state_v, state_t, p, piln, zm, zi, &
        nm, ni, rhoi, kvtt, state_q, dse, &
        effgw_rdg_resid, use_gw_rdg_resid, &
        effgw_rdg_beta, effgw_rdg_beta_max, &
        effgw_rdg_gamma, effgw_rdg_gamma_max, &
        rdg_beta_cd_llb, trpd_leewv_rdg_beta, &
        rdg_gamma_cd_llb, trpd_leewv_rdg_gamma, &
        q_tend, s_tend, u_tend, v_tend, flx_heat, errmsg, errflg)

    end if

!!$  ! Convert the tendencies for the dry constituents to dry air basis.
!!$  do m = 1, pcnst
!!$     if (cnst_type(m).eq.'dry') then
!!$        do k = 1, pver
!!$           do i = 1, ncol
!!$              q_tend(i,k,m) = q_tend(i,k,m)*pdel(i,k)/pdeldry(i,k)
!!$           end do
!!$        end do
!!$     end if
!!$  end do

!!$  ! Write totals to history file.
!!$  call outfld('EKGW', egwdffi_tot , ncol, lchnk)
!!$  call outfld('TTGW', ptend%s/cpairv(:,:),  pcols, lchnk)
!!$
!!$  call outfld('UTGW_TOTAL', ptend%u, pcols, lchnk)
!!$  call outfld('VTGW_TOTAL', ptend%v, pcols, lchnk)
!!$
!!$  call outfld('QTGW', ptend%q(:,:,1), pcols, lchnk)
!!$  call outfld('CLDLIQTGW', ptend%q(:,:,ixcldliq), pcols, lchnk)
!!$  call outfld('CLDICETGW', ptend%q(:,:,ixcldice), pcols, lchnk)

    ! Destroy objects.
    call p%finalize()

  end subroutine gw_drag_run

!==========================================================================

! Add all history fields for a gravity wave spectrum source.
  subroutine gw_spec_addflds(prefix, scheme, band, history_defaults)
    use cam_history, only: addfld, add_default, register_vector_field

    !------------------------------Arguments--------------------------------

    ! One character prefix prepended to output fields.
    character(len=1), intent(in) :: prefix
    ! Gravity wave scheme name prepended to output field descriptions.
    character(len=*), intent(in) :: scheme
    ! Wave speeds.
    type(GWBand), intent(in) :: band
    ! Whether or not to call add_default for fields output by WACCM.
    logical, intent(in) :: history_defaults

    !---------------------------Local storage-------------------------------

    character(len=*), parameter :: sub = 'gw_spec_addflds'
    integer :: l
    ! 7 chars is enough for "-100.00"
    character(len=7)  :: fnum
    ! 10 chars is enough for "BTAUXSn32"
    character(len=10) :: dumc1x, dumc1y
    ! Allow 80 chars for description
    character(len=80) dumc2

    !-----------------------------------------------------------------------

    ! Overall wind tendencies.
    call addfld(trim(prefix)//'UTGWSPEC', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency - gravity wave spectrum')
    call addfld(trim(prefix)//'VTGWSPEC', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' V tendency - gravity wave spectrum')
    call register_vector_field(trim(prefix)//'UTGWSPEC', trim(prefix)//'VTGWSPEC')

    call addfld(trim(prefix)//'TTGWSPEC', (/'lev'/), 'A', 'K s-1', &
                trim(scheme)//' T tendency - gravity wave spectrum')

    ! Wind tendencies broken across five spectral bins.
    call addfld(trim(prefix)//'UTEND1', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency   c < -40')
    call addfld(trim(prefix)//'UTEND2', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency  -40 < c < -15')
    call addfld(trim(prefix)//'UTEND3', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency  -15 < c <  15')
    call addfld(trim(prefix)//'UTEND4', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency   15 < c <  40')
    call addfld(trim(prefix)//'UTEND5', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency   40 < c ')

    ! Reynold's stress toward each cardinal direction, and net zonal stress.
    call addfld(trim(prefix)//'TAUE', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' Eastward Reynolds stress')
    call addfld(trim(prefix)//'TAUW', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' Westward Reynolds stress')
    call addfld(trim(prefix)//'TAUNET', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' E+W Reynolds stress')
    call addfld(trim(prefix)//'TAUN', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' Northward Reynolds stress')
    call addfld(trim(prefix)//'TAUS', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' Southward Reynolds stress')

    ! Momentum flux in each direction.
    call addfld(trim(prefix)//'EMF', (/'lev'/), 'A', 'Pa', &
                trim(scheme)//' Eastward MF')
    call addfld(trim(prefix)//'WMF', (/'lev'/), 'A', 'Pa', &
                trim(scheme)//' Westward MF')
    call addfld(trim(prefix)//'NMF', (/'lev'/), 'A', 'Pa', &
                trim(scheme)//' Northward MF')
    call addfld(trim(prefix)//'SMF', (/'lev'/), 'A', 'Pa', &
                trim(scheme)//' Southward MF')

    ! Temperature tendency terms.
    call addfld(trim(prefix)//'TTGWSDF', (/'lev'/), 'A', 'K s-1', &
                trim(scheme)//' t tendency - diffusion term')
    call addfld(trim(prefix)//'TTGWSKE', (/'lev'/), 'A', 'K s-1', &
                trim(scheme)//' t tendency - kinetic energy conversion term')

    ! Gravity wave source spectra by wave number.
    do l = -band%ngwv, band%ngwv
      ! String containing reference speed.
      write (fnum, fmt='(f7.2)') band%cref(l)

      dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
      dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)
      dumc2 = trim(scheme)//" tau at c= "//trim(fnum)//" m s-1"
      call addfld(trim(dumc1x), (/'lev'/), 'A', 'Pa', dumc2)
      call addfld(trim(dumc1y), (/'lev'/), 'A', 'Pa', dumc2)

    end do

    if (history_defaults) then
      call add_default(trim(prefix)//'UTGWSPEC', 1, ' ')
      call add_default(trim(prefix)//'VTGWSPEC', 1, ' ')
      call add_default(trim(prefix)//'TTGWSPEC', 1, ' ')
      call add_default(trim(prefix)//'TAUE', 1, ' ')
      call add_default(trim(prefix)//'TAUW', 1, ' ')
      call add_default(trim(prefix)//'TAUNET', 1, ' ')
      call add_default(trim(prefix)//'TAUN', 1, ' ')
      call add_default(trim(prefix)//'TAUS', 1, ' ')
    end if

  end subroutine gw_spec_addflds

!==========================================================================

! Outputs for spectral waves.
  subroutine gw_spec_outflds(prefix, ncol, pver, band, phase_speeds, u, v, xv, yv, &
                             gwut, dttdf, dttke, tau, utgw, vtgw, ttgw, taucd)

    use gw_common, only: west, east, south, north

    ! One-character prefix prepended to output fields.
    character(len=1), intent(in) :: prefix
    ! Chunk and number of columns in the chunk.
    integer, intent(in) :: ncol
    integer, intent(in) :: pver
    ! Wave speeds.
    type(GWBand), intent(in) :: band
    ! Wave phase speeds for each column.
    real(kind_phys), intent(in) :: phase_speeds(ncol, -band%ngwv:band%ngwv)
    ! Winds at cell midpoints.
    real(kind_phys), intent(in) :: u(:, :)
    real(kind_phys), intent(in) :: v(:, :)
    ! Unit vector in the direction of wind at source level.
    real(kind_phys), intent(in) :: xv(:)
    real(kind_phys), intent(in) :: yv(:)
    ! Wind tendency for each wave.
    real(kind_phys), intent(in) :: gwut(ncol, pver, -band%ngwv:band%ngwv)
    ! Temperature tendencies from diffusion and kinetic energy.
    real(kind_phys), intent(in) :: dttdf(:, :)
    real(kind_phys), intent(in) :: dttke(:, :)
    ! Wave Reynolds stress.
    real(kind_phys), intent(in) :: tau(ncol, -band%ngwv:band%ngwv, pver)
    ! Zonal and meridional total wind tendencies.
    real(kind_phys), intent(in) :: utgw(:, :)
    real(kind_phys), intent(in) :: vtgw(:, :)
    ! Temperature tendencies.
    real(kind_phys), intent(in) :: ttgw(:, :)
    ! Reynolds stress for waves propagating in each cardinal direction.
    real(kind_phys), intent(in) :: taucd(ncol, pver + 1, 4)

    character(len=*), parameter :: sub = 'gw_spec_outflds'
    ! Indices
    integer :: i, k, l
    integer :: ix(ncol, -band%ngwv:band%ngwv), iy(ncol, -band%ngwv:band%ngwv)
    integer :: iu(ncol), iv(ncol)

    ! Zonal wind tendency, broken up into five bins.
    real(kind_phys) :: utb(ncol, pver, 5)
    ! Definition of the bin boundaries.
    real(kind_phys), parameter :: bounds(4) = (/-40._kind_phys, -15._kind_phys, &
                                                15._kind_phys, 40._kind_phys/)

    ! Momentum flux in the four cardinal directions.
    real(kind_phys) :: mf(ncol, pver, 4)

    ! Wave stress in zonal/meridional direction
    real(kind_phys) :: taux(ncol, -band%ngwv:band%ngwv, pver)
    real(kind_phys) :: tauy(ncol, -band%ngwv:band%ngwv, pver)

    ! Temporaries for output
    real(kind_phys) :: dummyx(ncol, pver)
    real(kind_phys) :: dummyy(ncol, pver)
    ! Variable names
    character(len=10) :: dumc1x, dumc1y

    ! Accumulate wind tendencies binned according to phase speed.

    utb = 0._kind_phys

    ! Find which output bin the phase speed corresponds to.
    ix = find_bin(phase_speeds)

    ! Put the wind tendency in that bin.
    do l = -band%ngwv, band%ngwv
      do k = 1, pver
        do i = 1, ncol
          utb(i, k, ix(i, l)) = utb(i, k, ix(i, l)) + gwut(i, k, l)
        end do
      end do
    end do

    ! Find just the zonal part.
    do l = 1, 5
      do k = 1, pver
        utb(:, k, l) = utb(:, k, l)*xv
      end do
    end do

!!$  call outfld(trim(prefix)//'UTEND1', utb(:,:,1), ncol, lchnk)
!!$  call outfld(trim(prefix)//'UTEND2', utb(:,:,2), ncol, lchnk)
!!$  call outfld(trim(prefix)//'UTEND3', utb(:,:,3), ncol, lchnk)
!!$  call outfld(trim(prefix)//'UTEND4', utb(:,:,4), ncol, lchnk)
!!$  call outfld(trim(prefix)//'UTEND5', utb(:,:,5), ncol, lchnk)
!!$
!!$  ! Output temperature tendencies due to diffusion and from kinetic energy.
!!$  call outfld(trim(prefix)//'TTGWSDF', dttdf / cpair, ncol, lchnk)
!!$  call outfld(trim(prefix)//'TTGWSKE', dttke / cpair, ncol, lchnk)
    ! Output tau broken down into zonal and meridional components.

    taux = 0._kind_phys
    tauy = 0._kind_phys

    ! Project phase_speeds, and convert each component to a wavenumber index.
    ! These are mappings from the wavenumber index of tau to those of taux
    ! and tauy, respectively.
    do l = -band%ngwv, band%ngwv
      ix(:, l) = c_to_l(phase_speeds(:, l)*xv)
      iy(:, l) = c_to_l(phase_speeds(:, l)*yv)
    end do

    ! Find projection of tau.
    do k = 1, pver
      do l = -band%ngwv, band%ngwv
        do i = 1, ncol
          taux(i, ix(i, l), k) = taux(i, ix(i, l), k) &
                                 + abs(tau(i, l, k)*xv(i))
          tauy(i, iy(i, l), k) = tauy(i, iy(i, l), k) &
                                 + abs(tau(i, l, k)*yv(i))
        end do
      end do
    end do

    do l = -band%ngwv, band%ngwv

      dummyx = taux(:, l, :)
      dummyy = tauy(:, l, :)

      dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
      dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)

!jt     call outfld(dumc1x,dummyx,ncol,lchnk)
!jt     call outfld(dumc1y,dummyy,ncol,lchnk)

    end do

    ! Output momentum flux in each cardinal direction.
    mf = 0._kind_phys

    do k = 1, pver

      ! Convert wind speed components to wavenumber indices.
      iu = c_to_l(u(:, k))
      iv = c_to_l(v(:, k))

      ! Sum tau components in each cardinal direction.
      ! Split west/east and north/south based on whether wave speed exceeds
      ! wind speed.
      do l = -band%ngwv, band%ngwv

        where (iu > l)
          mf(:, k, west) = mf(:, k, west) + taux(:, l, k)
        elsewhere
          mf(:, k, east) = mf(:, k, east) + taux(:, l, k)
        end where

        where (iv > l)
          mf(:, k, south) = mf(:, k, south) + tauy(:, l, k)
        elsewhere
          mf(:, k, north) = mf(:, k, north) + tauy(:, l, k)
        end where

      end do

    end do

!!$  call outfld(trim(prefix)//'WMF',mf(:,:,west),ncol,lchnk)
!!$  call outfld(trim(prefix)//'EMF',mf(:,:,east),ncol,lchnk)
!!$  call outfld(trim(prefix)//'SMF',mf(:,:,south),ncol,lchnk)
!!$  call outfld(trim(prefix)//'NMF',mf(:,:,north),ncol,lchnk)
!!$
!!$  ! Simple output fields written to history file.
!!$  ! Total wind tendencies.
!!$  call outfld (trim(prefix)//'UTGWSPEC', utgw , ncol, lchnk)
!!$  call outfld (trim(prefix)//'VTGWSPEC', vtgw , ncol, lchnk)
!!$  call outfld (trim(prefix)//'TTGWSPEC', ttgw , ncol, lchnk)
!!$
!!$  ! Tau in each direction.
!!$  call outfld (trim(prefix)//'TAUE', taucd(:,:,east), ncol, lchnk)
!!$  call outfld (trim(prefix)//'TAUW', taucd(:,:,west), ncol, lchnk)
!!$  call outfld (trim(prefix)//'TAUN', taucd(:,:,north), ncol, lchnk)
!!$  call outfld (trim(prefix)//'TAUS', taucd(:,:,south), ncol, lchnk)
!!$
!!$  call outfld (trim(prefix)//'TAUNET', taucd(:,:,east)+taucd(:,:,west), &
!!$       ncol, lchnk)

  contains

    ! Given a value, finds which bin marked by "bounds" the value falls
    ! into.
    elemental function find_bin(val) result(idx)
      real(kind_phys), intent(in) :: val

      integer :: idx

      ! We just have to count how many bounds are exceeded.
      if (val >= 0._kind_phys) then
        idx = count(val > bounds) + 1
      else
        idx = count(val >= bounds) + 1
      end if

    end function find_bin

    ! Convert a speed to a wavenumber between -ngwv and ngwv.
    elemental function c_to_l(c) result(l)
      real(kind_phys), intent(in) :: c

      integer :: l

      l = min(max(int(c/band%dc), -band%ngwv), band%ngwv)

    end function c_to_l

  end subroutine gw_spec_outflds

!==========================================================================

! Generates names for tau output across the wave spectrum (e.g.
! BTAUXSn01 or TAUYSp05).
! Probably this should use a wavenumber dimension on one field rather
! than creating a ton of numbered fields.
  character(len=9) pure function tau_fld_name(l, prefix, x_not_y)
    ! Wavenumber
    integer, intent(in) :: l
    ! Single-character prefix for output
    character(len=1), intent(in) :: prefix
    ! X or Y?
    logical, intent(in) :: x_not_y

    character(len=2) :: num_str

    tau_fld_name = trim(prefix)

    tau_fld_name = trim(tau_fld_name)//"TAU"

    if (x_not_y) then
      tau_fld_name = trim(tau_fld_name)//"XS"
    else
      tau_fld_name = trim(tau_fld_name)//"YS"
    end if

    if (l < 0) then
      tau_fld_name = trim(tau_fld_name)//"n"
    else
      tau_fld_name = trim(tau_fld_name)//"p"
    end if

    write (num_str, '(I2.2)') abs(l)

    tau_fld_name = trim(tau_fld_name)//num_str

  end function tau_fld_name

end module gw_drag
