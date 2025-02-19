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
  use ccpp_kinds, only:  kind_phys
  use gw_common,      only: GWBand,handle_err
  use gw_convect,     only: BeresSourceDesc
  use gw_movmtn,      only: MovMtnSourceDesc
  use gw_front,       only: CMSourceDesc
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

  character(len=512)                   :: bnd_rdggm ! full pathname for meso-Gamma ridge dataset
  ! Whether or not to limit tau *before* applying any efficiency factors.
  logical                              :: gw_limit_tau_without_eff = .false.
  logical                              :: gw_lndscl_sgh = .true. ! scale SGH by land frac
  real(kind_phys)                      :: gw_prndl = 0.25_kind_phys
  ! Whether or not to apply tendency max
  real(kind_phys)                      :: gw_qbo_hdepth_scaling = 1._kind_phys ! heating depth scaling factor
  logical                              :: gw_top_taper=.false.
  ! Width of gaussian used to create frontogenesis tau profile [m s-1].
  real(kind_phys)                      :: front_gaussian_width = -huge(1._kind_phys)
  real(kind_phys)                      :: alpha_gw_movmtn

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
  real(kind_phys),allocatable   :: pref_edge(:), pref_mid(:)
  real(kind_phys)   :: degree2radian
  integer           :: gw_bot_taper_pres
  integer           :: ncid_topo
  logical           :: masterproc


  logical         ::  use_gw_oro
  logical         ::  use_gw_front
  logical         ::  use_gw_front_igw
  logical         ::  use_gw_convect_dp
  logical         ::  use_gw_convect_sh
  logical         ::  use_simple_phys
  logical         ::  use_gw_movmtn_pbl
  logical         ::  do_molec_diff
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


!!$  ! Beres settings and table.
!!$  type(BeresSourceDesc) :: beres_dp_desc
!!$  type(BeresSourceDesc) :: beres_sh_desc
!!$
!!$  ! Moving mountain settings and table.
!!$  type(MovMtnSourceDesc) :: movmtn_desc
!!$
!!$  ! Frontogenesis wave settings.
!!$  type(CMSourceDesc) :: cm_desc
!!$  type(CMSourceDesc) :: cm_igw_desc
!!$
  real(kind_phys), allocatable, dimension(:), target :: &
     rdg_gbxar

     ! Meso Beta
  real(kind_phys), allocatable, dimension(:,:), target :: &
     rdg_hwdth,  &
     rdg_clngt,  &
     rdg_mxdis,  &
     rdg_anixy,  &
     rdg_angll

  real(kind_phys), allocatable, dimension(:), target :: &
     rdg_gbxarg

  ! Meso Gamma
  real(kind_phys), allocatable, dimension(:,:), target :: &
     rdg_hwdthg, &
     rdg_clngtg, &
     rdg_mxdisg, &
     rdg_anixyg, &
     rdg_angllg

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


  real(kind_phys), allocatable, TARGET :: vramp(:)

  real(kind_phys) :: C_BetaMax_DS, C_GammaMax, &
              Frx0, Frx1, C_BetaMax_SM, Fr_c, &
              orohmin, orovmin, orostratmin, orom2min
  logical ::  do_divstream, do_smooth_regimes, do_adjust_tauoro, &
              do_backward_compat,do_vdiff

  logical                         ::  beres_dp_desc_storm_shift
  integer                         ::  beres_dp_desc_k
  real(kind_phys)                 ::  beres_dp_desc_min_hdepth
  integer                         ::  beres_dp_desc_maxh
  integer                         ::  beres_dp_desc_maxuh
  real(kind_phys),allocatable     ::  beres_dp_desc_hd(:)
  real(kind_phys),allocatable     ::  beres_dp_desc_mfcc(:,:,:)
  logical                         ::  beres_sh_desc_storm_shift
  integer                         ::  beres_sh_desc_k
  real(kind_phys)                 ::  beres_sh_desc_min_hdepth
  integer                         ::  beres_sh_desc_maxh
  integer                         ::  beres_sh_desc_maxuh
  real(kind_phys),allocatable     ::  beres_sh_desc_hd(:)
  real(kind_phys),allocatable     ::  beres_sh_desc_mfcc(:,:,:)
  logical                         ::  movmtn_desc_storm_shift
  integer                         ::  movmtn_desc_k
  real(kind_phys)                 ::  movmtn_desc_min_hdepth
  integer                         ::  movmtn_desc_maxh
  integer                         ::  movmtn_desc_maxuh
  real(kind_phys),allocatable     ::  movmtn_desc_hd(:)
  real(kind_phys),allocatable     ::  movmtn_desc_uh(:)
  real(kind_phys),allocatable     ::  movmtn_desc_mfcc(:,:,:)
  integer                         ::  cm_desc_ksrc
  integer                         ::  cm_desc_kfront
  real(kind_phys)                 ::  cm_desc_frontgfc
  real(kind_phys),allocatable     ::  cm_desc_src_tau(:)
  integer                         ::  cm_igw_desc_ksrc
  integer                         ::  cm_igw_desc_kfront
  real(kind_phys)                 ::  cm_igw_desc_frontgfc
  real(kind_phys),allocatable     ::  cm_igw_desc_src_tau(:)
!==========================================================================
contains
!==========================================================================

!> \section arg_table_gw_drag_init Argument Table
!! \htmlinclude gw_drag_init.html
  subroutine gw_drag_init( &
       pver,  &
       gravit_in,  &
       rair_in,  &
       pi_in,  &
       pgwv_in,  &
       gw_dc_in,  &
       pgwv_long_in,  &
       gw_dc_long_in,  &
       tau_0_ubc_in,  &
       pref_edge_in,  &
       pref_mid_in,  &
       gw_bot_taper_pres,  &
       effgw_beres_dp_in,  &
       effgw_beres_sh_in,  &
       effgw_cm_in,  &
       effgw_cm_igw_in,  &
       effgw_oro_in,  &
       fcrit2_in,  &
       frontgfc_in,  &
       gw_drag_file_in,  &
       gw_drag_file_sh_in,  &
       gw_drag_file_mm_in,  &
       taubgnd_in,  &
       taubgnd_igw_in,  &
       gw_polar_taper_in,  &
       use_gw_rdg_beta_in,  &
       n_rdg_beta_in,  &
       effgw_rdg_beta_in,  &
       effgw_rdg_beta_max_in,  &
       rdg_beta_cd_llb_in,  &
       trpd_leewv_rdg_beta_in,  &
       use_gw_rdg_gamma_in,  &
       n_rdg_gamma_in,  &
       effgw_rdg_gamma_in,  &
       effgw_rdg_gamma_max_in,  &
       rdg_gamma_cd_llb_in,  &
       trpd_leewv_rdg_gamma_in,  &
       bnd_rdggm_in,  &
       gw_oro_south_fac_in,  &
       gw_limit_tau_without_eff_in,  &
       gw_lndscl_sgh_in,  &
       gw_prndl_in,  &
       gw_apply_tndmax_in,  &
       gw_qbo_hdepth_scaling_in,  &
       gw_top_taper_in,  &
       front_gaussian_width_in,  &
       alpha_gw_movmtn_in,  &
       use_gw_front_in,  &
       use_gw_oro_in,  &
       use_gw_front_igw_in,  &
       use_gw_convect_dp_in,  &
       use_gw_convect_sh_in,  &
       use_simple_phys_in,  &
       use_gw_movmtn_pbl_in,  &
       iulog_in,  &
       ktop_in,  &
       masterproc_in, &
       do_molec_diff_in,  &
       wavelength_mid_in,  &
       wavelength_long_in,  &
!!$    gw_rdg_do_divstream,  &
!!$    gw_rdg_C_BetaMax_DS,  &
!!$    gw_rdg_C_GammaMax,  &
!!$    gw_rdg_Frx0,  &
!!$    gw_rdg_Frx1,  &
!!$    gw_rdg_C_BetaMax_SM,  &
!!$    gw_rdg_Fr_c,  &
!!$    gw_rdg_do_smooth_regimes,  &
!!$    gw_rdg_do_adjust_tauoro,  &
!!$    gw_rdg_do_backward_compat,  &
!!$    gw_rdg_orohmin,  &
!!$    gw_rdg_orovmin,  &
!!$    gw_rdg_orostratmin,  &
!!$    gw_rdg_orom2min,  &
!!$    gw_rdg_do_vdiff,  &
!!$       rdg_gbxar_in,  &
!!$       rdg_hwdth_in,  &
!!$       rdg_clngt_in,  &
!!$       rdg_mxdis_in,  &
!!$       rdg_anixy_in,  &
!!$       rdg_angll_in,  &
!!$       rdg_gbxarg_in,  &
!!$       rdg_hwdthg_in,  &
!!$       rdg_clngtg_in,  &
!!$       rdg_mxdisg_in,  &
!!$       rdg_anixyg_in,  &
!!$       rdg_angllg_in,  &
       beres_dp_desc_storm_shift_in, &
       beres_dp_desc_k_in, &
       beres_dp_desc_min_hdepth_in, &
       beres_dp_desc_maxh_in, &
       beres_dp_desc_maxuh_in, &
       beres_dp_desc_hd_in, &
       beres_dp_desc_mfcc_in, &
       beres_sh_desc_storm_shift_in, &
       beres_sh_desc_k_in, &
       beres_sh_desc_min_hdepth_in, &
       beres_sh_desc_maxh_in, &
       beres_sh_desc_maxuh_in, &
       beres_sh_desc_hd_in, &
       beres_sh_desc_mfcc_in, &
       movmtn_desc_storm_shift_in, &
       movmtn_desc_k_in, &
       movmtn_desc_min_hdepth_in, &
       movmtn_desc_maxh_in, &
       movmtn_desc_maxuh_in, &
       movmtn_desc_hd_in, &
       movmtn_desc_uh_in, &
       movmtn_desc_mfcc_in, &
       cm_desc_ksrc_in, &
       cm_desc_kfront_in, &
       cm_desc_frontgfc_in, &
       cm_desc_src_tau_in, &
       cm_igw_desc_ksrc_in, &
       cm_igw_desc_kfront_in, &
       cm_igw_desc_frontgfc_in, &
       cm_igw_desc_src_tau_in, &
       errmsg,  &
       errflg )

    use gw_common,  only: gw_common_init
  !-----------------------------------------------------------------------
  ! Time independent initialization for multiple gravity wave
  ! parameterization.
  !-----------------------------------------------------------------------

  integer, intent(in)             :: pver
  real(kind_phys), intent(in)     :: gravit_in          ! gravitational acceleration (m s-2)
  real(kind_phys), intent(in)     :: rair_in            ! Dry air gas constant     (J K-1 kg-1)
  real(kind_phys), intent(in)     :: pi_in
  ! Maximum wave number and width of spectrum bins.
  integer, intent(in)             :: pgwv_in
  real(kind_phys), intent(in)     :: gw_dc_in
  integer, intent(in)             :: pgwv_long_in
  real(kind_phys), intent(in)     :: gw_dc_long_in
  ! Whether or not to enforce an upper boundary condition of tau = 0.
  ! (Like many variables, this is only here to hold the value between
  ! the readnl phase and the init phase of the CAM physics; only gw_common
  ! should actually use it.)
  logical, intent(in)             :: tau_0_ubc_in
  real(kind_phys), intent(in)     :: pref_edge_in(:)
  real(kind_phys), intent(in)     :: pref_mid_in(:)
  real(kind_phys), intent(in)     :: gw_bot_taper_pres
  ! Beres (deep convection).
  real(kind_phys), intent(in)     :: effgw_beres_dp_in
  ! Beres (shallow convection).
  real(kind_phys), intent(in)     :: effgw_beres_sh_in
  ! C&M scheme.
  real(kind_phys), intent(in)     :: effgw_cm_in
  ! C&M scheme (inertial waves).
  real(kind_phys), intent(in)     :: effgw_cm_igw_in
  ! Orography.
  real(kind_phys), intent(in)     :: effgw_oro_in
  ! fcrit2 for the mid-scale waves has been made a namelist variable to
  ! facilitate backwards compatibility with the CAM3 version of this
  ! parameterization.  In CAM3, fcrit2=0.5.
  real(kind_phys), intent(in)             :: fcrit2_in   ! critical froude number squared
  ! Frontogenesis function critical threshold.
  real(kind_phys), intent(in)             :: frontgfc_in
  ! Files to read Beres source spectra from.
  character(len=256), intent(in)             :: gw_drag_file_in
  character(len=256), intent(in)             :: gw_drag_file_sh_in
  character(len=256), intent(in)             :: gw_drag_file_mm_in
  ! Background stress source strengths.
  real(kind_phys), intent(in)             :: taubgnd_in
  real(kind_phys), intent(in)             :: taubgnd_igw_in
  ! Whether or not to use a polar taper for frontally generated waves.
  logical, intent(in)             :: gw_polar_taper_in
  ! Ridge scheme.
  logical, intent(in)              :: use_gw_rdg_beta_in
  integer , intent(in)             :: n_rdg_beta_in
  real(kind_phys), intent(in)             :: effgw_rdg_beta_in
  real(kind_phys), intent(in)             :: effgw_rdg_beta_max_in
  real(kind_phys), intent(in)             :: rdg_beta_cd_llb_in  ! Low-level obstacle drag coefficient Ridge scheme.
  logical , intent(in)             :: trpd_leewv_rdg_beta_in
  logical , intent(in)             :: use_gw_rdg_gamma_in
  integer , intent(in)             :: n_rdg_gamma_in
  real(kind_phys), intent(in)             :: effgw_rdg_gamma_in
  real(kind_phys), intent(in)             :: effgw_rdg_gamma_max_in
  real(kind_phys), intent(in)             :: rdg_gamma_cd_llb_in
  logical , intent(in)             :: trpd_leewv_rdg_gamma_in
  character(len=256), intent(in)    :: bnd_rdggm_in ! full pathname for meso-Gamma ridge dataset
  ! Factor for SH orographic waves.
  real(kind_phys), intent(in)             :: gw_oro_south_fac_in
  ! Whether or not to limit tau *before* applying any efficiency factors.
  logical, intent(in)             :: gw_limit_tau_without_eff_in
  logical, intent(in)              :: gw_lndscl_sgh_in
  real(kind_phys), intent(in)             :: gw_prndl_in
  ! Whether or not to apply tendency max
  logical, intent(in)             :: gw_apply_tndmax_in
  real(kind_phys), intent(in)             :: gw_qbo_hdepth_scaling_in
  logical, intent(in)             :: gw_top_taper_in
  ! Width of gaussian used to create frontogenesis tau profile [m s-1].
  real(kind_phys), intent(in)             :: front_gaussian_width_in
  real(kind_phys), intent(in)             :: alpha_gw_movmtn_in
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
  logical, intent(in)             :: do_molec_diff_in
  real(kind_phys), intent(in)     :: wavelength_mid_in
  real(kind_phys), intent(in)     :: wavelength_long_in
!!$  logical, intent(in)             :: gw_rdg_do_divstream
!!$  real(kind_phys), intent(in)     ::  gw_rdg_C_BetaMax_DS
!!$  real(kind_phys), intent(in)     ::  gw_rdg_C_GammaMax
!!$  real(kind_phys), intent(in)     ::  gw_rdg_Frx0
!!$  real(kind_phys), intent(in)     ::  gw_rdg_Frx1
!!$  real(kind_phys), intent(in)     ::  gw_rdg_C_BetaMax_SM
!!$  real(kind_phys), intent(in)     ::  gw_rdg_Fr_c
!!$  logical, intent(in)             :: gw_rdg_do_smooth_regimes
!!$  logical, intent(in)             :: gw_rdg_do_adjust_tauoro
!!$  logical, intent(in)             :: gw_rdg_do_backward_compat
!!$  real(kind_phys), intent(in)     ::  gw_rdg_orohmin
!!$  real(kind_phys), intent(in)     ::  gw_rdg_orovmin
!!$  real(kind_phys), intent(in)     ::  gw_rdg_orostratmin
!!$  real(kind_phys), intent(in)     ::  gw_rdg_orom2min
!!$  logical, intent(in)             ::  gw_rdg_do_vdiff
!!$  real(kind_phys), intent(in)     ::  rdg_gbxar_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_hwdth_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_clngt_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_mxdis_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_anixy_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_angll_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_gbxarg_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_hwdthg_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_clngtg_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_mxdisg_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_anixyg_in(:,:)
!!$  real(kind_phys), intent(in)     ::  rdg_angllg_in(:,:)
  logical                         ::  beres_dp_desc_storm_shift_in
  integer                         ::  beres_dp_desc_k_in
  real(kind_phys), intent(in)     ::  beres_dp_desc_min_hdepth_in
  integer                         ::  beres_dp_desc_maxh_in
  integer                         ::  beres_dp_desc_maxuh_in
  real(kind_phys), intent(in)     ::  beres_dp_desc_hd_in(:)
  real(kind_phys), intent(in)     ::  beres_dp_desc_mfcc_in(:,:,:)
  logical                         ::  beres_sh_desc_storm_shift_in
  integer                         ::  beres_sh_desc_k_in
  real(kind_phys), intent(in)     ::  beres_sh_desc_min_hdepth_in
  integer                         ::  beres_sh_desc_maxh_in
  integer                         ::  beres_sh_desc_maxuh_in
  real(kind_phys), intent(in)     ::  beres_sh_desc_hd_in(:)
  real(kind_phys), intent(in)     ::  beres_sh_desc_mfcc_in(:,:,:)
  logical                         ::  movmtn_desc_storm_shift_in
  integer                         ::  movmtn_desc_k_in
  real(kind_phys), intent(in)     ::  movmtn_desc_min_hdepth_in
  integer                         ::  movmtn_desc_maxh_in
  integer                         ::  movmtn_desc_maxuh_in
  real(kind_phys), intent(in)     ::  movmtn_desc_hd_in(:)
  real(kind_phys), intent(in)     ::  movmtn_desc_uh_in(:)
  real(kind_phys), intent(in)     ::  movmtn_desc_mfcc_in(:,:,:)
  integer, intent(in)             ::  cm_desc_ksrc_in
  integer, intent(in)             ::  cm_desc_kfront_in
  real(kind_phys), intent(in)     ::  cm_desc_frontgfc_in
  real(kind_phys), intent(in)     ::  cm_desc_src_tau_in(:)
  integer, intent(in)             ::  cm_igw_desc_ksrc_in
  integer, intent(in)             ::  cm_igw_desc_kfront_in
  real(kind_phys), intent(in)     ::  cm_igw_desc_frontgfc_in
  real(kind_phys), intent(in)     ::  cm_igw_desc_src_tau_in(:)
  character(len=512), intent(out) :: errmsg
  integer, intent(out)            :: errflg

  ! Local variables
  character(len=*), parameter :: sub = 'gw_drag_init'

  integer :: k,l
  character(len=128) :: errstring

  ! Index for levels at specific pressures.
  integer :: topndx
  integer :: botndx
  integer :: stat

  ! Interpolated Newtonian cooling coefficients.
  real(kind_phys) :: alpha(pver+1)

  ! Levels of pre-calculated Newtonian cooling (1/day).
  ! The following profile is digitized from:
  ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5

  integer, parameter :: nalph = 71
  real(kind_phys) :: alpha0(nalph) = [ &
       0.1_kind_phys,         0.1_kind_phys,         0.1_kind_phys,         0.1_kind_phys,         &
       0.1_kind_phys,         0.1_kind_phys,         0.1_kind_phys,         0.1_kind_phys,         &
       0.1_kind_phys,         0.1_kind_phys,         0.10133333_kind_phys,  0.104_kind_phys,       &
       0.108_kind_phys,       0.112_kind_phys,       0.116_kind_phys,       0.12066667_kind_phys,  &
       0.126_kind_phys,       0.132_kind_phys,       0.138_kind_phys,       0.144_kind_phys,       &
       0.15133333_kind_phys,  0.16_kind_phys,        0.17_kind_phys,        0.18_kind_phys,        &
       0.19_kind_phys,        0.19933333_kind_phys,  0.208_kind_phys,       0.216_kind_phys,       &
       0.224_kind_phys,       0.232_kind_phys,       0.23466667_kind_phys,  0.232_kind_phys,       &
       0.224_kind_phys,       0.216_kind_phys,       0.208_kind_phys,       0.20133333_kind_phys,  &
       0.196_kind_phys,       0.192_kind_phys,       0.188_kind_phys,       0.184_kind_phys,       &
       0.18266667_kind_phys,  0.184_kind_phys,       0.188_kind_phys,       0.192_kind_phys,       &
       0.196_kind_phys,       0.19333333_kind_phys,  0.184_kind_phys,       0.168_kind_phys,       &
       0.152_kind_phys,       0.136_kind_phys,       0.12133333_kind_phys,  0.108_kind_phys,       &
       0.096_kind_phys,       0.084_kind_phys,       0.072_kind_phys,       0.061_kind_phys,       &
       0.051_kind_phys,       0.042_kind_phys,       0.033_kind_phys,       0.024_kind_phys,       &
       0.017666667_kind_phys, 0.014_kind_phys,       0.013_kind_phys,       0.012_kind_phys,       &
       0.011_kind_phys,       0.010333333_kind_phys, 0.01_kind_phys,        0.01_kind_phys,        &
       0.01_kind_phys,        0.01_kind_phys,        0.01_kind_phys                         &
       ]

  ! Pressure levels that were used to calculate alpha0 (hPa).
  real(kind_phys) :: palph(nalph) = [ &
       2.06115E-06_kind_phys, 2.74280E-06_kind_phys, 3.64988E-06_kind_phys, 4.85694E-06_kind_phys, &
       6.46319E-06_kind_phys, 8.60065E-06_kind_phys, 1.14450E-05_kind_phys, 1.52300E-05_kind_phys, &
       2.02667E-05_kind_phys, 2.69692E-05_kind_phys, 3.58882E-05_kind_phys, 4.77568E-05_kind_phys, &
       6.35507E-05_kind_phys, 8.45676E-05_kind_phys, 0.000112535_kind_phys, 0.000149752_kind_phys, &
       0.000199277_kind_phys, 0.000265180_kind_phys, 0.000352878_kind_phys, 0.000469579_kind_phys, &
       0.000624875_kind_phys, 0.000831529_kind_phys, 0.00110653_kind_phys,  0.00147247_kind_phys,  &
       0.00195943_kind_phys,  0.00260744_kind_phys,  0.00346975_kind_phys,  0.00461724_kind_phys,  &
       0.00614421_kind_phys,  0.00817618_kind_phys,  0.0108801_kind_phys,   0.0144783_kind_phys,   &
       0.0192665_kind_phys,   0.0256382_kind_phys,   0.0341170_kind_phys,   0.0453999_kind_phys,   &
       0.0604142_kind_phys,   0.0803939_kind_phys,   0.106981_kind_phys,    0.142361_kind_phys,    &
       0.189442_kind_phys,    0.252093_kind_phys,    0.335463_kind_phys,    0.446404_kind_phys,    &
       0.594036_kind_phys,    0.790490_kind_phys,    1.05192_kind_phys,     1.39980_kind_phys,     &
       1.86273_kind_phys,     2.47875_kind_phys,     3.29851_kind_phys,     4.38936_kind_phys,     &
       5.84098_kind_phys,     7.77266_kind_phys,     10.3432_kind_phys,     13.7638_kind_phys,     &
       18.3156_kind_phys,     24.3728_kind_phys,     32.4332_kind_phys,     43.1593_kind_phys,     &
       57.4326_kind_phys,     76.4263_kind_phys,     101.701_kind_phys,     135.335_kind_phys,     &
       180.092_kind_phys,     239.651_kind_phys,     318.907_kind_phys,     424.373_kind_phys,     &
       564.718_kind_phys,     751.477_kind_phys,     1000._kind_phys                        &
       ]

  !-----------------------------------------------------------------------
  errmsg =''
  errflg = 0


  band_oro = GWBand(0, gw_dc_in, fcrit2_in, wavelength_mid_in)
  band_mid = GWBand(pgwv_in, gw_dc_in, 1.0_kind_phys, wavelength_mid_in)
  band_long = GWBand(pgwv_long_in, gw_dc_long_in, 1.0_kind_phys, wavelength_long_in)
  band_movmtn = GWBand(0, gw_dc_in, 1.0_kind_phys, wavelength_mid_in)

  gravit  = gravit_in
  rair    = rair_in
  pi      = pi_in
  allocate(pref_edge(pver))
  pref_edge = pref_edge_in
  allocate(pref_mid(pver))
  pref_mid  = pref_mid_in
  masterproc =  masterproc_in
  iulog  =   iulog_in
  pgwv  =   pgwv_in
  gw_dc =   gw_dc_in
  pgwv_long =   pgwv_long_in
  gw_dc_long =   gw_dc_long_in
  tau_0_ubc =   tau_0_ubc_in
  effgw_beres_dp =   effgw_beres_dp_in
  effgw_beres_sh =   effgw_beres_sh_in
  effgw_cm =   effgw_cm_in
  effgw_cm_igw =   effgw_cm_igw_in
  effgw_oro =   effgw_oro_in
  fcrit2 =   fcrit2_in
  frontgfc =   frontgfc_in
  gw_drag_file    =   trim(gw_drag_file_in)
  gw_drag_file_sh =   trim(gw_drag_file_sh_in)
  gw_drag_file_mm =   trim(gw_drag_file_mm_in)
  taubgnd =   taubgnd_in
  taubgnd_igw =   taubgnd_igw_in
  gw_polar_taper =   gw_polar_taper_in
  use_gw_rdg_beta =   use_gw_rdg_beta_in
  n_rdg_beta =   n_rdg_beta_in
  effgw_rdg_beta =   effgw_rdg_beta_in
  effgw_rdg_beta_max =   effgw_rdg_beta_max_in
  rdg_beta_cd_llb =   rdg_beta_cd_llb_in
  trpd_leewv_rdg_beta =   trpd_leewv_rdg_beta_in
  use_gw_rdg_gamma =   use_gw_rdg_gamma_in
  n_rdg_gamma =   n_rdg_gamma_in
  effgw_rdg_gamma =   effgw_rdg_gamma_in
  effgw_rdg_gamma_max =   effgw_rdg_gamma_max_in
  rdg_gamma_cd_llb =   rdg_gamma_cd_llb_in
  trpd_leewv_rdg_gamma =   trpd_leewv_rdg_gamma_in
  bnd_rdggm =   trim(bnd_rdggm_in)
  gw_oro_south_fac =   gw_oro_south_fac_in
  gw_limit_tau_without_eff =   gw_limit_tau_without_eff_in
  gw_lndscl_sgh =   gw_lndscl_sgh_in
  gw_prndl =   gw_prndl_in
  gw_apply_tndmax =   gw_apply_tndmax_in
  gw_qbo_hdepth_scaling =   gw_qbo_hdepth_scaling_in
  gw_top_taper =   gw_top_taper_in
  front_gaussian_width =   front_gaussian_width_in
  alpha_gw_movmtn =   alpha_gw_movmtn_in
  use_gw_oro =   use_gw_oro_in
  use_gw_front =   use_gw_front_in
  use_gw_front_igw =   use_gw_front_igw_in
  use_gw_convect_dp =   use_gw_convect_dp_in
  use_gw_convect_sh =   use_gw_convect_sh_in
  use_simple_phys =   use_simple_phys_in
  use_gw_movmtn_pbl =   use_gw_movmtn_pbl_in
  do_molec_diff = do_molec_diff_in
  wavelength_mid = wavelength_mid_in
  wavelength_long = wavelength_long_in
  ktop = ktop_in
  if (use_gw_convect_dp) then
  beres_dp_desc_storm_shift = beres_dp_desc_storm_shift_in
  beres_dp_desc_k = beres_dp_desc_k_in
  beres_dp_desc_min_hdepth = beres_dp_desc_min_hdepth_in
  beres_dp_desc_maxh = beres_dp_desc_maxh_in
  beres_dp_desc_maxuh = beres_dp_desc_maxuh_in
  allocate(beres_dp_desc_hd(beres_dp_desc_maxh), stat=stat)
  call handle_err( stat, errflg,sub//': Allocate of beres_dp_desc_hd failed', errmsg )
  beres_dp_desc_hd = beres_dp_desc_hd_in
  allocate(beres_dp_desc_mfcc(beres_dp_desc_maxh,-beres_dp_desc_maxuh:beres_dp_desc_maxuh,&
       -band_mid%ngwv:band_mid%ngwv), stat=stat)
  call handle_err( stat, errflg,sub//': Allocate of beres_dp_desc_mfcc failed', errmsg )
  beres_dp_desc_mfcc = beres_dp_desc_mfcc_in
  end if
  if (use_gw_convect_sh) then
  beres_sh_desc_storm_shift = beres_sh_desc_storm_shift_in
  beres_sh_desc_k = beres_sh_desc_k_in
  beres_sh_desc_min_hdepth = beres_sh_desc_min_hdepth_in
  beres_sh_desc_maxh = beres_sh_desc_maxh_in
  beres_sh_desc_maxuh = beres_sh_desc_maxuh_in
  allocate(beres_sh_desc_hd(beres_sh_desc_maxh))
  call handle_err( stat, errflg,sub//': Allocate of beres_sh_desc_hd failed', errmsg )
  beres_sh_desc_hd = beres_sh_desc_hd_in
  allocate(beres_sh_desc_mfcc(beres_sh_desc_maxh,-beres_sh_desc_maxuh:beres_sh_desc_maxuh,&
       -band_mid%ngwv:band_mid%ngwv), stat=stat)
  call handle_err( stat, errflg,sub//': Allocate of beres_dp_desc_mfcc failed', errmsg )
  beres_sh_desc_mfcc = beres_sh_desc_mfcc_in
  end if
  if (use_gw_movmtn_pbl) then
  movmtn_desc_storm_shift = movmtn_desc_storm_shift_in
  movmtn_desc_k = movmtn_desc_k_in
  movmtn_desc_min_hdepth = movmtn_desc_min_hdepth_in
  movmtn_desc_maxh = movmtn_desc_maxh_in
  movmtn_desc_maxuh = movmtn_desc_maxuh_in
  allocate(movmtn_desc_hd(movmtn_desc_maxh))
  call handle_err( stat, errflg,sub//': Allocate of movmtn_desc_maxh failed', errmsg )
  movmtn_desc_hd = movmtn_desc_hd_in
  allocate(movmtn_desc_uh(movmtn_desc_maxuh))
  call handle_err( stat, errflg,sub//': Allocate of movmtn_desc_maxuh failed', errmsg )
  movmtn_desc_uh = movmtn_desc_uh_in
  allocate(movmtn_desc_mfcc(movmtn_desc_maxh,-movmtn_desc_maxuh:movmtn_desc_maxuh,&
       -band_movmtn%ngwv:band_movmtn%ngwv), stat=stat)
  call handle_err( stat, errflg,sub//': Allocate of movmtn_desc_mfcc failed', errmsg )
  movmtn_desc_mfcc = movmtn_desc_mfcc_in
  end if
  if (use_gw_front) then
  cm_desc_ksrc = cm_desc_ksrc_in
  cm_desc_kfront = cm_desc_kfront_in
  cm_desc_frontgfc = cm_desc_frontgfc_in
  allocate(cm_desc_src_tau(-band_mid%ngwv:band_mid%ngwv))
  call handle_err( stat, errflg,sub//': Allocate of cm_desc_src_tau failed', errmsg )
  cm_desc_src_tau = cm_desc_src_tau_in
  end if
  if (use_gw_front_igw) then
  cm_igw_desc_ksrc = cm_igw_desc_ksrc_in
  cm_igw_desc_kfront = cm_igw_desc_kfront_in
  cm_igw_desc_frontgfc = cm_igw_desc_frontgfc_in
  allocate(cm_igw_desc_src_tau(-band_long%ngwv:band_long%ngwv))
  call handle_err( stat, errflg,sub//': Allocate of cm_igw_desc_src_tau failed', errmsg )
  cm_igw_desc_src_tau = cm_igw_desc_src_tau_in
  end if
!!$  do_divstream = gw_rdg_do_divstream
!!$  C_BetaMax_DS = gw_rdg_C_BetaMax_DS
!!$  C_GammaMax = gw_rdg_C_GammaMax
!!$  Frx0 = gw_rdg_Frx0
!!$  Frx1 = gw_rdg_Frx1
!!$  C_BetaMax_SM = gw_rdg_C_BetaMax_SM
!!$  Fr_c = gw_rdg_Fr_c
!!$  do_smooth_regimes = gw_rdg_do_smooth_regimes
!!$  do_adjust_tauoro = gw_rdg_do_adjust_tauoro
!!$  do_backward_compat = gw_rdg_do_backward_compat
!!$  orohmin = gw_rdg_orohmin
!!$  orovmin = gw_rdg_orovmin
!!$  orostratmin = gw_rdg_orostratmin
!!$  orom2min = gw_rdg_orom2min
!!$  do_vdiff = gw_rdg_do_vdiff



  if (masterproc) then
     write(iulog,*) ' '
     write(iulog,*) "GW_DRAG: band_mid%ngwv = ", band_mid%ngwv
     do l = -band_mid%ngwv, band_mid%ngwv
        write (iulog,'(A,I0,A,F7.2)') &
             "GW_DRAG: band_mid%cref(",l,") = ",band_mid%cref(l)
     enddo
     write(iulog,*) 'GW_DRAG: band_mid%kwv = ', band_mid%kwv
     write(iulog,*) 'GW_DRAG: band_mid%fcrit2 = ', band_mid%fcrit2
     write(iulog,*) ' '
     write(iulog,*) "GW_DRAG: band_long%ngwv = ", band_long%ngwv
     do l = -band_long%ngwv, band_long%ngwv
        write (iulog,'(A,I2,A,F7.2)') &
             "GW_DRAG: band_long%cref(",l,") = ",band_long%cref(l)
     enddo
     write(iulog,*) 'GW_DRAG: band_long%kwv = ', band_long%kwv
     write(iulog,*) 'GW_DRAG: band_long%fcrit2 = ', band_long%fcrit2
     write(iulog,*) ' '
  end if

  ! pre-calculated newtonian damping:
  !     * convert to s-1
  !     * ensure it is not smaller than 1e-6
  !     * convert palph from hpa to pa

  do k=1,nalph
     alpha0(k) = alpha0(k) / 86400._kind_phys
     alpha0(k) = max(alpha0(k), 1.e-6_kind_phys)
     palph(k) = palph(k)*1.e2_kind_phys
  end do

  call lininterp (alpha0  ,palph, nalph , alpha  , pref_edge , pver+1)
  if (masterproc) then
     write (iulog,*) 'gw_init: newtonian damping (1/day):'
     write (iulog,fmt='(a4,a12,a10)') ' k  ','  pref_edge      ', &
          '  alpha   '
     do k = 1, pver+1
        write (iulog,fmt='(i4,1e12.5,1f10.2)') k,pref_edge(k), &
             alpha(k)*86400._kind_phys
     end do
  end if

  if (masterproc) then
     write(iulog,*) 'KTOP        =',ktop
  end if

  ! Initialize subordinate modules.
  call gw_common_init(pver,&
       tau_0_ubc, ktop, gravit, rair, alpha, gw_prndl, &
       gw_qbo_hdepth_scaling, errstring )
!!$  call shr_assert(trim(errstring) == "", "gw_common_init: "//errstring// &
!!$       errorMsg(__FILE__, __LINE__))

  if (gw_top_taper) then
     allocate(vramp(pver))
     call handle_err( stat, errflg,sub//': Allocate of vramp failed', errmsg )
     vramp(:) = 1._kind_phys
     topndx = 1
     botndx = press_lim_idx( gw_bot_taper_pres, pref_mid, top=.true. )
     if (botndx>1) then
        do k=botndx,topndx,-1
           vramp(k) = vramp(k+1)/(pref_edge(k+1)/pref_edge(k))
        end do
        if (masterproc) then
           write(iulog,'(A)') 'GW taper coef (vramp):'
           do k=1,pver
              write(iulog,"('k: ',I4,' taper coef,press(Pa): ',F12.8,E12.4)") k, vramp(k), pref_mid(k)
           enddo
        endif
     endif
  end if
contains
  ! Convert pressure limiters to the appropriate level.
pure function press_lim_idx(p, pref_mid, top) result(k_lim)
  ! Pressure
  real(kind_phys), intent(in) :: p
  real(kind_phys), intent(in) :: pref_mid(:)
  ! Is this a top or bottom limit?
  logical,  intent(in) :: top
  integer :: k_lim, k

  if (top) then
     k_lim = pver+1
     do k = 1, pver
        if (pref_mid(k) > p) then
           k_lim = k
           exit
        end if
     end do
  else
     k_lim = 0
     do k = pver, 1, -1
        if (pref_mid(k) < p) then
           k_lim = k
           exit
        end if
     end do
  end if

end function press_lim_idx

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
     sgh, &
     kvtt, &
     ttend_dp, &
     ttend_sh, &
     ttend_clubb, &
     thlp2_clubb_gw, &
     wpthlp_clubb_gw, &
     upwp_clubb_gw, &
     vpwp_clubb_gw, &
     rdg_gbxar,  &
     rdg_hwdth,  &
     rdg_clngt,  &
     rdg_mxdis,  &
     rdg_anixy,  &
     rdg_angll,  &
     rdg_gbxarg,  &
     rdg_hwdthg,  &
     rdg_clngtg,  &
     rdg_mxdisg,  &
     rdg_anixyg,  &
     rdg_angllg,  &
     s_tend, &
     q_tend, &
     u_tend, &
     v_tend, &
     scheme_name, &
     nbot_molec, &
     flx_heat, &
     errmsg, errflg)
  !-----------------------------------------------------------------------
  ! Interface for multiple gravity wave drag parameterization.
  !-----------------------------------------------------------------------

  use coords_1d,       only: Coords1D
  use gw_common,       only: gw_prof, gw_drag_prof, calc_taucd
  use gw_common,       only: momentum_flux, momentum_fixer, energy_change
  use gw_common,       only: energy_fixer, coriolis_speed, adjust_inertial
  use gw_oro,          only: gw_oro_src
  use gw_front,        only: gw_cm_src
  use gw_convect,      only: gw_beres_src
  use gw_movmtn,       only: gw_movmtn_run

  integer,  intent(in)        :: ncol  ! number of atmospheric columns
  integer,  intent(in)        :: pcnst ! chunk number
  integer,  intent(in)        :: pver  ! number of atmospheric levels
  character*3, intent(in)     :: cnst_type(pcnst) ! wet or dry mixing ratio
  real(kind_phys), intent(in) :: dt          ! physics timestep
  real(kind_phys), intent(in) :: cpair       ! heat capacity of air
  real(kind_phys), intent(in) :: cpairv(:,:) ! location dependent heat capacity of air
  real(kind_phys), intent(in) :: pi          ! pi
  real(kind_phys), intent(in) :: degree2radian
  real(kind_phys), intent(in) :: al0
  real(kind_phys), intent(in) :: dlat0
  real(kind_phys), intent(in) :: pint(:,:)   ! pressure at model interfaces
  real(kind_phys), intent(in) :: piln(:,:)   ! ln pressure at model interfaces
  real(kind_phys), intent(in) :: pdel(:,:)   ! vertical delta-p
  real(kind_phys), intent(in) :: pdeldry(:,:)   ! vertical delta-p
  real(kind_phys), intent(in) :: zm(:,:)
  real(kind_phys), intent(in) :: zi(:,:)
  real(kind_phys), intent(in) :: lat(:)
  real(kind_phys), intent(in) :: landfrac(:)
  real(kind_phys), intent(in) :: dse(:,:)       ! dry static energy
  real(kind_phys), intent(in) :: state_t(:,:)   ! temperature (K)
  real(kind_phys), intent(in) :: state_u(:,:)   ! meridional wind
  real(kind_phys), intent(in) :: state_v(:,:)   ! zonal wind
  real(kind_phys), intent(in) :: state_q(:,:,:) ! constituent array
  real(kind_phys), intent(in) :: sgh(:)         !
  real(kind_phys), intent(in) :: kvtt(:,:)       !
  real(kind_phys), intent(in) :: ttend_dp(:,:)  ! Temperature change due to deep convection.
  real(kind_phys), intent(in) :: ttend_sh(:,:)  ! Temperature change due to shallow convection.
  real(kind_phys), intent(in) :: ttend_clubb(:,:)
  real(kind_phys), intent(in) :: thlp2_clubb_gw(:,:)
  real(kind_phys), intent(in) :: wpthlp_clubb_gw(:,:)
  real(kind_phys), intent(in) :: upwp_clubb_gw(:,:)
  real(kind_phys), intent(in) :: vpwp_clubb_gw(:,:)
  real(kind_phys), intent(in), TARGET :: rdg_gbxar(:)
  real(kind_phys), intent(in), TARGET :: rdg_hwdth(:,:)
  real(kind_phys), intent(in), TARGET :: rdg_clngt(:,:)
  real(kind_phys), intent(inout), TARGET :: rdg_mxdis(:,:)
  real(kind_phys), intent(in), TARGET :: rdg_anixy(:,:)
  real(kind_phys), intent(in), TARGET :: rdg_angll(:,:)
  real(kind_phys), intent(in), TARGET :: rdg_gbxarg(:)
  real(kind_phys), intent(in), TARGET :: rdg_hwdthg(:,:)
  real(kind_phys), intent(in), TARGET :: rdg_clngtg(:,:)
  real(kind_phys), intent(inout), TARGET :: rdg_mxdisg(:,:)
  real(kind_phys), intent(in), TARGET :: rdg_anixyg(:,:)
  real(kind_phys), intent(in), TARGET :: rdg_angllg(:,:)
  real(kind_phys), intent(inout):: s_tend(:,:)   ! dry air enthalpy tendency
  real(kind_phys), intent(inout):: q_tend(:,:,:)
  real(kind_phys), intent(inout):: u_tend(:,:)
  real(kind_phys), intent(inout):: v_tend(:,:)
  character(len=64),  intent(out) :: scheme_name
  integer, intent(in)             :: nbot_molec
  ! Parameterization net tendencies.
  real(kind_phys), intent(inout) :: flx_heat(:)
  character(len=512), intent(out) :: errmsg
  integer,            intent(out) :: errflg

  !---------------------------Local storage-------------------------------
  character(len=*), parameter :: sub = 'gw_drag_run'

  integer :: lchnk                  ! chunk identifier
  integer :: stat

  integer :: i, k                   ! loop indices

  type(Coords1D) :: p               ! Pressure coordinates

  real(kind_phys) :: ttgw(ncol,pver) ! temperature tendency
  real(kind_phys) :: utgw(ncol,pver) ! zonal wind tendency
  real(kind_phys) :: vtgw(ncol,pver) ! meridional wind tendency

  real(kind_phys) :: ni(ncol,pver+1) ! interface Brunt-Vaisalla frequency
  real(kind_phys) :: nm(ncol,pver)   ! midpoint Brunt-Vaisalla frequency
  real(kind_phys) :: rhoi(ncol,pver+1)     ! interface density
  real(kind_phys), allocatable :: tau(:,:,:)  ! wave Reynolds stress
  real(kind_phys) :: tau0x(ncol)     ! c=0 sfc. stress (zonal)
  real(kind_phys) :: tau0y(ncol)     ! c=0 sfc. stress (meridional)
  real(kind_phys) :: ubi(ncol,pver+1)! projection of wind at interfaces
  real(kind_phys) :: ubm(ncol,pver)  ! projection of wind at midpoints
  real(kind_phys) :: xv(ncol)        ! unit vector of source wind (x)
  real(kind_phys) :: yv(ncol)        ! unit vector of source wind (y)

  integer :: m                      ! dummy integers
  real(kind_phys) :: qtgw(ncol,pver,pcnst) ! constituents tendencies

  ! Reynolds stress for waves propagating in each cardinal direction.
  real(kind_phys) :: taucd(ncol,pver+1,4)

  ! gravity wave wind tendency for each wave
  real(kind_phys), allocatable :: gwut(:,:,:)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(kind_phys) :: dttdf(ncol,pver)
  real(kind_phys) :: dttke(ncol,pver)

  ! Wave phase speeds for each column
  real(kind_phys), allocatable :: phase_speeds(:,:)

  ! Efficiency for a gravity wave source.
  real(kind_phys) :: effgw(ncol)

  ! Coriolis characteristic speed.
  real(kind_phys) :: u_coriolis(ncol)

  ! Adjustment for inertial gravity waves.
  real(kind_phys), allocatable :: ro_adjust(:,:,:)

!!$  ! pbuf fields
!!$  ! Molecular diffusivity
!!$  real(kind_phys) :: kvtt(ncol,pver+1)

  ! Frontogenesis
  real(kind_phys), pointer :: frontgf(:,:)
  real(kind_phys), pointer :: frontga(:,:)


  ! gridbox area
  real(kind_phys), pointer :: gbxar(:)

     ! Beta ridges
  ! width of ridges.
  real(kind_phys), pointer :: hwdth(:,:)
  ! length of ridges.
  real(kind_phys), pointer :: clngt(:,:)
  ! Maximum deviations of ridges.
  real(kind_phys), pointer :: mxdis(:,:)
  ! orientation of ridges.
  real(kind_phys), pointer :: angll(:,:)
  ! anisotropy of ridges.
  real(kind_phys), pointer :: anixy(:,:)

     ! Gamma ridges
  ! width of ridges.
  real(kind_phys), pointer :: hwdthg(:,:)
  ! length of ridges.
  real(kind_phys), pointer :: clngtg(:,:)
  ! Maximum deviations of ridges.
  real(kind_phys), pointer :: mxdisg(:,:)
  ! orientation of ridges.
  real(kind_phys), pointer :: angllg(:,:)
  ! anisotropy of ridges.
  real(kind_phys), pointer :: anixyg(:,:)

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
  real(kind_phys) :: egwdffi(ncol,pver+1)
  ! sum from the two types of spectral GW
  real(kind_phys) :: egwdffi_tot(ncol,pver+1)

  ! Momentum fluxes used by fixer.
  real(kind_phys) :: um_flux(ncol), vm_flux(ncol)
  ! Energy change used by fixer.
  real(kind_phys) :: de(ncol)

  ! Which constituents are being affected by diffusion.
  logical  :: lq(pcnst)

  ! Contiguous copies of state arrays.
  real(kind_phys) :: flx_heat_curr(ncol)

  !------------------------------------------------------------------------
  p = Coords1D(pint(:ncol,:))

  ! Profiles of background state variables
  call gw_prof(ncol, p, cpair, state_t, rhoi, nm, ni)

  if (use_gw_front_igw) then
     u_coriolis = coriolis_speed(band_long, lat(:ncol))
  end if

  ! Totals that accumulate over different sources.
  egwdffi_tot = 0._kind_phys
  flx_heat = 0._kind_phys

  if (use_gw_movmtn_pbl) then
     !------------------------------------------------------------------
     !Convective moving mountain gravity waves (Beres scheme).
     !------------------------------------------------------------------
     flx_heat_curr = flx_heat
     effgw=1.0_kind_phys

     call gw_movmtn_run(ncol, &
          band_movmtn , &
          movmtn_desc_storm_shift, &
          movmtn_desc_k, &
          movmtn_desc_min_hdepth, &
          movmtn_desc_maxh, &
          movmtn_desc_maxuh, &
          movmtn_desc_hd, &
          movmtn_desc_uh, &
          movmtn_desc_mfcc, &
          state_t, pcnst, &
          state_u, state_v, p, ttend_dp, ttend_clubb, &
          upwp_clubb_gw, vpwp_clubb_gw, &
          zm, alpha_gw_movmtn, src_level, tend_level, &
          ubm, ubi, xv, yv, hdepth,dt, &
          vramp, pint, piln, rhoi, nm,   ni, &
          effgw,  kvtt, state_q,  dse, utgw,  vtgw, &
          ttgw, qtgw, egwdffi, dttdf, dttke, gw_apply_tndmax, &
          flx_heat_curr, use_gw_movmtn_pbl, &
          rair, gravit,q_tend,u_tend,v_tend,s_tend,errmsg,errflg)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do
  end if

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

     !------------------------------------------------------------------
     ! Convective gravity waves (Beres scheme, deep).
     !------------------------------------------------------------------

     ! Allocate wavenumber fields.
     allocate(tau(ncol,-band_mid%ngwv:band_mid%ngwv,pver+1),stat=stat)
     call handle_err( stat,errflg,sub//': Allocate of tau failed', errmsg )
     allocate(gwut(ncol,pver,-band_mid%ngwv:band_mid%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of gwut failed', errmsg )
     allocate(phase_speeds(ncol,-band_mid%ngwv:band_mid%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of phase_speeds failed', errmsg )

     ! Efficiency of gravity wave momentum transfer.
     ! This is really only to remove the pole points.
     where (pi/2._kind_phys - abs(lat(:ncol)) >= 4*epsilon(1._kind_phys))
        effgw = effgw_beres_dp
     elsewhere
        effgw = 0._kind_phys
     end where

     ! Determine wave sources for Beres deep scheme
     call gw_beres_src(ncol, band_mid, &
          beres_dp_desc_storm_shift, &
          beres_dp_desc_k, &
          beres_dp_desc_min_hdepth, &
          beres_dp_desc_maxh, &
          beres_dp_desc_maxuh, &
          beres_dp_desc_hd, &
          beres_dp_desc_mfcc, &
          state_u, state_v, ttend_dp(:ncol,:), zm, src_level, tend_level, tau, &
          ubm, ubi, xv, yv, phase_speeds, hdepth, maxq0)

     ! Solve for the drag profile with Beres source spectrum.
     call gw_drag_prof(ncol, band_mid, p, src_level, tend_level, dt, &
          state_t, vramp,    &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,   phase_speeds,       kvtt, state_q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          lapply_effgw_in=gw_apply_tndmax)

     ! Project stress into directional components.
     taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     ! Store constituents tendencies
     do m=1, pcnst
        do k = 1, pver
           q_tend(:ncol,k,m) = q_tend(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Find momentum flux, and use it to fix the wind tendencies below
     ! the gravity wave region.
     call momentum_flux(tend_level, taucd, um_flux, vm_flux)
     call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

     ! Add the momentum tendencies to the output tendency arrays.
     do k = 1, pver
        u_tend(:ncol,k) = u_tend(:ncol,k) + utgw(:,k)
        v_tend(:ncol,k) = v_tend(:ncol,k) + vtgw(:,k)
     end do

     ! Find energy change in the current state, and use fixer to apply
     ! the difference in lower levels.
     call energy_change(dt, p, state_u, state_v, u_tend(:ncol,:), &
          v_tend(:ncol,:), s_tend(:ncol,:)+ttgw, de)
     call energy_fixer(tend_level, p, de-flx_heat(:ncol), ttgw)

     do k = 1, pver
        s_tend(:ncol,k) = s_tend(:ncol,k) + ttgw(:,k)
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

     deallocate(tau, gwut, phase_speeds)

  end if

  if (use_gw_convect_sh) then
!!$     !------------------------------------------------------------------
!!$     ! Convective gravity waves (Beres scheme, shallow).
!!$     !------------------------------------------------------------------
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
     allocate(tau(ncol,-band_mid%ngwv:band_mid%ngwv,pver+1),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of tau failed', errmsg )
     allocate(gwut(ncol,pver,-band_mid%ngwv:band_mid%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of gwut failed', errmsg )
     allocate(phase_speeds(ncol,-band_mid%ngwv:band_mid%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of phase_speeds failed', errmsg )

     ! Efficiency of gravity wave momentum transfer.
     ! This is really only to remove the pole points.
     where (pi/2._kind_phys - abs(lat(:ncol)) >= 4*epsilon(1._kind_phys))
        effgw = effgw_beres_sh
     elsewhere
        effgw = 0._kind_phys
     end where

     ! Determine wave sources for Beres shallow scheme
     call gw_beres_src(ncol, band_mid, &
          beres_sh_desc_storm_shift, &
          beres_sh_desc_k, &
          beres_sh_desc_min_hdepth, &
          beres_sh_desc_maxh, &
          beres_sh_desc_maxuh, &
          beres_sh_desc_hd, &
          beres_sh_desc_mfcc, &
          state_u, state_v, ttend_sh(:ncol,:), zm, src_level, tend_level, tau, &
          ubm, ubi, xv, yv, phase_speeds, hdepth, maxq0)

     ! Solve for the drag profile with Beres source spectrum.
     call gw_drag_prof(ncol, band_mid, p, src_level, tend_level,  dt, &
          state_t, vramp,    &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,   phase_speeds,       kvtt, state_q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          lapply_effgw_in=gw_apply_tndmax)

     ! Project stress into directional components.
     taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     ! Store constituents tendencies
     do m=1, pcnst
        do k = 1, pver
           q_tend(:ncol,k,m) = q_tend(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Add the momentum tendencies to the output tendency arrays.
     ! Don't calculate fixers, since we are too close to the ground to
     ! spread momentum/energy differences across low layers.
     do k = 1, pver
        u_tend(:ncol,k) = u_tend(:ncol,k) + utgw(:,k)
        v_tend(:ncol,k) = v_tend(:ncol,k) + vtgw(:,k)
        s_tend(:ncol,k) = s_tend(:ncol,k) + ttgw(:,k)
     end do

     ! Calculate energy change for output to CAM's energy checker.
     ! This is sort of cheating; we don't have a good a priori idea of the
     ! energy coming from surface stress, so we just integrate what we and
     ! actually have so far and overwrite flx_heat with that.
     call energy_change(dt, p, state_u, state_v, u_tend(:ncol,:), &
          v_tend(:ncol,:), s_tend(:ncol,:), de)
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

     deallocate(tau, gwut, phase_speeds)


  end if

  if (use_gw_front) then
     !------------------------------------------------------------------
     ! Frontally generated gravity waves
     !------------------------------------------------------------------
     ! Allocate wavenumber fields.
     allocate(tau(ncol,-band_mid%ngwv:band_mid%ngwv,pver+1),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of tau failed', errmsg )
     allocate(gwut(ncol,pver,-band_mid%ngwv:band_mid%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of gwut failed', errmsg )
     allocate(phase_speeds(ncol,-band_mid%ngwv:band_mid%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of phase_speeds failed', errmsg )

     ! Efficiency of gravity wave momentum transfer.
     effgw = effgw_cm
     ! Frontogenesis is too high at the poles (at least for the FV
     ! dycore), so introduce a polar taper.
     if (gw_polar_taper) effgw = effgw * cos(lat(:ncol))

     call gw_cm_src(ncol, band_mid, &
          cm_desc_ksrc, &
          cm_desc_kfront, &
          cm_desc_frontgfc, &
          cm_desc_src_tau, &
          state_u, state_v, frontgf(:ncol,:), &
          src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds)

     ! Solve for the drag profile with C&M source spectrum.
     call gw_drag_prof(ncol, band_mid, p, src_level, tend_level,  dt, &
          state_t, vramp,   &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,   phase_speeds,       kvtt, state_q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          lapply_effgw_in=gw_apply_tndmax)

     ! Project stress into directional components.
     taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     !Add the constituent tendencies
     do m=1, pcnst
        do k = 1, pver
           q_tend(:ncol,k,m) = q_tend(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Find momentum flux, and use it to fix the wind tendencies below
     ! the gravity wave region.
     call momentum_flux(tend_level, taucd, um_flux, vm_flux)
     call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

     ! add the momentum tendencies to the output tendency arrays
     do k = 1, pver
        u_tend(:ncol,k) = u_tend(:ncol,k) + utgw(:,k)
        v_tend(:ncol,k) = v_tend(:ncol,k) + vtgw(:,k)
     end do

     ! Find energy change in the current state, and use fixer to apply
     ! the difference in lower levels.
     call energy_change(dt, p, state_u, state_v, u_tend(:ncol,:), &
          v_tend(:ncol,:), s_tend(:ncol,:)+ttgw, de)
     call energy_fixer(tend_level, p, de-flx_heat(:ncol), ttgw)

     do k = 1, pver
        s_tend(:ncol,k) = s_tend(:ncol,k) + ttgw(:,k)
     end do

     ! Change ttgw to a temperature tendency before outputing it.
!!$     ttgw = ttgw / cpair
!!$     call gw_spec_outflds(cm_pf, ncol, pver, band_mid, phase_speeds, u, v, &
!!$          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
!!$          taucd)
!!$
!!$     deallocate(tau, gwut, phase_speeds)

  end if

  if (use_gw_front_igw) then
     !------------------------------------------------------------------
     ! Frontally generated inertial gravity waves
     !------------------------------------------------------------------
     ! Allocate wavenumber fields.
     allocate(tau(ncol,-band_long%ngwv:band_long%ngwv,pver+1),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of tau failed', errmsg )
     allocate(gwut(ncol,pver,-band_long%ngwv:band_long%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of gwut failed', errmsg )
     allocate(phase_speeds(ncol,-band_long%ngwv:band_long%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of phase_speeds failed', errmsg )
     allocate(ro_adjust(ncol,-band_long%ngwv:band_long%ngwv,pver+1),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of ro_adjust failed', errmsg )

     ! Efficiency of gravity wave momentum transfer.
     effgw = effgw_cm_igw

     ! Frontogenesis is too high at the poles (at least for the FV
     ! dycore), so introduce a polar taper.
     if (gw_polar_taper) then
        where (abs(lat(:ncol)) <= 89._kind_phys*degree2radian)
           effgw = effgw * 0.25_kind_phys * &
                 (1._kind_phys+tanh((lat(:ncol)+al0)/dlat0)) * &
                 (1._kind_phys-tanh((lat(:ncol)-al0)/dlat0))
        elsewhere
           effgw = 0._kind_phys
        end where
     end if

     call gw_cm_src(ncol, band_long, &
          cm_igw_desc_ksrc, &
          cm_igw_desc_kfront, &
          cm_igw_desc_frontgfc, &
          cm_igw_desc_src_tau, &
          state_u, state_v, frontgf(:ncol,:), &
          src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds)

     call adjust_inertial(band_long, tend_level, u_coriolis, phase_speeds, ubi, &
          tau, ro_adjust)

     ! Solve for the drag profile with C&M source spectrum.
     call gw_drag_prof(ncol, band_long, p, src_level, tend_level,  dt, &
          state_t, vramp,    &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,   phase_speeds,       kvtt, state_q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
          ro_adjust=ro_adjust, lapply_effgw_in=gw_apply_tndmax)

     ! Project stress into directional components.
     taucd = calc_taucd(ncol, band_long%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     !Add the constituent tendencies
     do m=1, pcnst
        do k = 1, pver
           q_tend(:ncol,k,m) = q_tend(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Find momentum flux, and use it to fix the wind tendencies below
     ! the gravity wave region.
     call momentum_flux(tend_level, taucd, um_flux, vm_flux)
     call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

     ! add the momentum tendencies to the output tendency arrays
     do k = 1, pver
        u_tend(:ncol,k) = u_tend(:ncol,k) + utgw(:,k)
        v_tend(:ncol,k) = v_tend(:ncol,k) + vtgw(:,k)
     end do

     ! Find energy change in the current state, and use fixer to apply
     ! the difference in lower levels.
     call energy_change(dt, p, state_u, state_v, u_tend(:ncol,:), &
          v_tend(:ncol,:), s_tend(:ncol,:)+ttgw, de)
     call energy_fixer(tend_level, p, de-flx_heat(:ncol), ttgw)

     do k = 1, pver
        s_tend(:ncol,k) = s_tend(:ncol,k) + ttgw(:,k)
     end do

     ! Change ttgw to a temperature tendency before outputing it.
!!$     ttgw = ttgw / cpair
!!$     call gw_spec_outflds(cm_igw_pf, ncol, pver, band_long, phase_speeds, state_u, state_v, &
!!$          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
!!$          taucd)

     deallocate(tau, gwut, phase_speeds, ro_adjust)

  end if

  if (use_gw_oro) then
     !---------------------------------------------------------------------
     ! Orographic stationary gravity waves
     !---------------------------------------------------------------------

     ! Allocate wavenumber fields.
     allocate(tau(ncol,band_oro%ngwv:band_oro%ngwv,pver+1),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of tau failed', errmsg )
     allocate(gwut(ncol,pver,band_oro%ngwv:band_oro%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of gwut failed', errmsg )
     allocate(phase_speeds(ncol,band_oro%ngwv:band_oro%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of phase_speeds failed', errmsg )

     if (gw_lndscl_sgh) then
        where (landfrac(:ncol) >= epsilon(1._kind_phys))
           effgw = effgw_oro * landfrac(:ncol)
           sgh_scaled = sgh(:ncol) / sqrt(landfrac(:ncol))
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
     endif
     do i = 1, ncol
        if (lat(i) < 0._kind_phys) then
           tau(i,:,:) = tau(i,:,:) * gw_oro_south_fac
        end if
     end do

     ! Solve for the drag profile with orographic sources.
     call gw_drag_prof(ncol, band_oro, p, src_level, tend_level,   dt,   &
          state_t, vramp,   &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw, phase_speeds, kvtt, state_q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          lapply_effgw_in=gw_apply_tndmax)

     ! For orographic waves, don't bother with taucd, since there are no
     ! momentum conservation routines or directional diagnostics.

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     ! Add the orographic tendencies to the spectrum tendencies.
     ! Don't calculate fixers, since we are too close to the ground to
     ! spread momentum/energy differences across low layers.
     do k = 1, pver
        u_tend(:ncol,k) = u_tend(:ncol,k) + utgw(:,k)
        v_tend(:ncol,k) = v_tend(:ncol,k) + vtgw(:,k)
        s_tend(:ncol,k) = s_tend(:ncol,k) + ttgw(:,k)
        ! Convert to temperature tendency for output.
        ttgw(:,k) = ttgw(:,k) / cpairv(:ncol, k)
     end do

     ! Calculate energy change for output to CAM's energy checker.
     ! This is sort of cheating; we don't have a good a priori idea of the
     ! energy coming from surface stress, so we just integrate what we and
     ! actually have so far and overwrite flx_heat with that.
     call energy_change(dt, p, state_u, state_v, u_tend(:ncol,:), &
          v_tend(:ncol,:), s_tend(:ncol,:), de)
     flx_heat(:ncol) = de

     do m = 1, pcnst
        do k = 1, pver
           q_tend(:ncol,k,m) = q_tend(:ncol,k,m) + qtgw(:,k,m)
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

     deallocate(tau, gwut, phase_speeds)

  end if

  if (use_gw_rdg_beta) then
     !---------------------------------------------------------------------
     ! Orographic stationary gravity waves
     !---------------------------------------------------------------------

     ! Efficiency of gravity wave momentum transfer.
     ! Take into account that wave sources are only over land.

     where(rdg_mxdis < 0._kind_phys)
        rdg_mxdis = 0._kind_phys
     end where

!!$     ! Save state at top of routine
!!$     ! Useful for unit testing checks
!!$     call outfld('UEGW', state_u ,  ncol, lchnk)
!!$     call outfld('VEGW', state_v ,  ncol, lchnk)
!!$     call outfld('TEGW', state_t ,  ncol, lchnk)
!!$     call outfld('ZEGW', zi , ncol, lchnk)
!!$     call outfld('ZMGW', zm , ncol, lchnk)

     call gw_rdg_calc(&
        pcnst, pver, 'BETA ', ncol, n_rdg_beta, dt,     &
        state_u, state_v, state_t, p, piln, zm, zi,                 &
        nm, ni, rhoi, kvtt, state_q, dse,               &
        effgw_rdg_beta, effgw_rdg_beta_max,       &
        rdg_hwdth, rdg_clngt, rdg_gbxar, rdg_mxdis, rdg_angll, rdg_anixy, &
        rdg_beta_cd_llb, trpd_leewv_rdg_beta,     &
        q_tend, s_tend, u_tend, v_tend, flx_heat,errmsg, errflg)

  end if

  if (use_gw_rdg_gamma) then
     !---------------------------------------------------------------------
     ! Orographic stationary gravity waves
     !---------------------------------------------------------------------

     ! Efficiency of gravity wave momentum transfer.
     ! Take into account that wave sources are only over land.

     where(rdg_mxdisg < 0._kind_phys)
        rdg_mxdisg = 0._kind_phys
     end where

     call gw_rdg_calc(&
        pcnst, pver, 'GAMMA', ncol, n_rdg_gamma, dt,         &
        state_u, state_v, state_t, p, piln, zm, zi,                      &
        nm, ni, rhoi, kvtt, state_q, dse,                    &
        effgw_rdg_gamma, effgw_rdg_gamma_max,          &
        rdg_hwdthg, rdg_clngtg, rdg_gbxar, rdg_mxdisg, rdg_angllg, rdg_anixyg, &
        rdg_gamma_cd_llb, trpd_leewv_rdg_gamma,        &
        q_tend, s_tend, u_tend, v_tend, flx_heat, errmsg, errflg)

  endif

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

subroutine gw_rdg_calc( &
   pcnst, pver, type, ncol, n_rdg, dt, &
   u, v, t, p, piln, zm, zi, &
   nm, ni, rhoi, kvtt, q, dse, &
   effgw_rdg, effgw_rdg_max, &
   hwdth, clngt, gbxar, &
   mxdis, angll, anixy, &
   rdg_cd_llb, trpd_leewv, &
   q_tend,s_tend,u_tend,v_tend, flx_heat, errmsg, errflg)

   use coords_1d,  only: Coords1D
   use gw_rdg,     only: gw_rdg_src, gw_rdg_belowpeak, gw_rdg_break_trap, gw_rdg_do_vdiff
   use gw_common,  only: gw_drag_prof, energy_change

   integer,          intent(in) :: pcnst        ! number of atmospheric constituents
   integer,          intent(in) :: pver         ! number of atmospheric constituents
   character(len=5), intent(in) :: type         ! BETA or GAMMA
   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: n_rdg
   real(kind_phys),         intent(in) :: dt           ! Time step.

   real(kind_phys),         intent(in) :: u(:,:)    ! Midpoint zonal winds. ( m s-1)
   real(kind_phys),         intent(in) :: v(:,:)    ! Midpoint meridional winds. ( m s-1)
   real(kind_phys),         intent(in) :: t(:,:)    ! Midpoint temperatures. (K)
   type(Coords1D),   intent(in) :: p               ! Pressure coordinates.
   real(kind_phys),         intent(in) :: piln(:,:)  ! Log of interface pressures.
   real(kind_phys),         intent(in) :: zm(:,:)   ! Midpoint altitudes above ground (m).
   real(kind_phys),         intent(in) :: zi(:,:) ! Interface altitudes above ground (m).
   real(kind_phys),         intent(in) :: nm(:,:)   ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real(kind_phys),         intent(in) :: ni(:,:) ! Interface Brunt-Vaisalla frequencies (s-1).
   real(kind_phys),         intent(in) :: rhoi(:,:) ! Interface density (kg m-3).
   real(kind_phys),         intent(in) :: kvtt(:,:) ! Molecular thermal diffusivity.
   real(kind_phys),         intent(in) :: q(:,:,:)        ! Constituent array.
   real(kind_phys),         intent(in) :: dse(:,:)  ! Dry static energy.


   real(kind_phys),         intent(in) :: effgw_rdg       ! Tendency efficiency.
   real(kind_phys),         intent(in) :: effgw_rdg_max
   real(kind_phys),         intent(in) :: hwdth(:,:) ! width of ridges.
   real(kind_phys),         intent(in) :: clngt(:,:) ! length of ridges.
   real(kind_phys),         intent(in) :: gbxar(:)      ! gridbox area

   real(kind_phys),         intent(in) :: mxdis(:,:) ! Height estimate for ridge (m).
   real(kind_phys),         intent(in) :: angll(:,:) ! orientation of ridges.
   real(kind_phys),         intent(in) :: anixy(:,:) ! Anisotropy parameter.

   real(kind_phys),         intent(in) :: rdg_cd_llb      ! Drag coefficient for low-level flow
   logical,          intent(in) :: trpd_leewv

   real(kind_phys), intent(inout):: s_tend(:,:)   ! dry air enthalpy tendency
   real(kind_phys), intent(inout):: q_tend(:,:,:)
   real(kind_phys), intent(inout):: u_tend(:,:)
   real(kind_phys), intent(inout):: v_tend(:,:)
   real(kind_phys),        intent(out) :: flx_heat(:)
   character(len=512), intent(out) :: errmsg
   integer, intent(out)            :: errflg

   !---------------------------Local storage-------------------------------

   character(len=*), parameter :: sub = 'gw_rdg_calc'
   integer :: k, m, nn, stat

   real(kind_phys), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(kind_phys), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(kind_phys), allocatable :: phase_speeds(:,:)

   ! Isotropic source flag [anisotropic orography].
   integer  :: isoflag(ncol)

   ! horiz wavenumber [anisotropic orography].
   real(kind_phys) :: kwvrdg(ncol)

   ! Efficiency for a gravity wave source.
   real(kind_phys) :: effgw(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)
   integer :: bwv_level(ncol)
   integer :: tlb_level(ncol)

   ! Projection of wind at midpoints and interfaces.
   real(kind_phys) :: ubm(ncol,pver)
   real(kind_phys) :: ubi(ncol,pver+1)

   ! Unit vectors of source wind (zonal and meridional components).
   real(kind_phys) :: xv(ncol)
   real(kind_phys) :: yv(ncol)

   ! Averages over source region.
   real(kind_phys) :: ubmsrc(ncol) ! On-ridge wind.
   real(kind_phys) :: usrc(ncol)   ! Zonal wind.
   real(kind_phys) :: vsrc(ncol)   ! Meridional wind.
   real(kind_phys) :: nsrc(ncol)   ! B-V frequency.
   real(kind_phys) :: rsrc(ncol)   ! Density.

   ! normalized wavenumber
   real(kind_phys) :: m2src(ncol)

   ! Top of low-level flow layer.
   real(kind_phys) :: tlb(ncol)

   ! Bottom of linear wave region.
   real(kind_phys) :: bwv(ncol)

   ! Froude numbers for flow/drag regimes
   real(kind_phys) :: Fr1(ncol)
   real(kind_phys) :: Fr2(ncol)
   real(kind_phys) :: Frx(ncol)

   ! Wave Reynolds stresses at source level
   real(kind_phys) :: tauoro(ncol)
   real(kind_phys) :: taudsw(ncol)

   ! Surface streamline displacement height for linear waves.
   real(kind_phys) :: hdspwv(ncol)

   ! Surface streamline displacement height for downslope wind regime.
   real(kind_phys) :: hdspdw(ncol)

   ! Wave breaking level
   real(kind_phys) :: wbr(ncol)

   real(kind_phys) :: utgw(ncol,pver)       ! zonal wind tendency
   real(kind_phys) :: vtgw(ncol,pver)       ! meridional wind tendency
   real(kind_phys) :: ttgw(ncol,pver)       ! temperature tendency
   real(kind_phys) :: qtgw(ncol,pver,pcnst) ! constituents tendencies

   ! Effective gravity wave diffusivity at interfaces.
   real(kind_phys) :: egwdffi(ncol,pver+1)

   ! Temperature tendencies from diffusion and kinetic energy.
   real(kind_phys) :: dttdf(ncol,pver)
   real(kind_phys) :: dttke(ncol,pver)

   ! Wave stress in zonal/meridional direction
   real(kind_phys) :: taurx(ncol,pver+1)
   real(kind_phys) :: taurx0(ncol,pver+1)
   real(kind_phys) :: taury(ncol,pver+1)
   real(kind_phys) :: taury0(ncol,pver+1)
   ! Provisional absolute wave stress from gw_drag_prof
   real(kind_phys) :: tau_diag(ncol,pver+1)

   ! U,V tendency accumulators
   real(kind_phys) :: utrdg(ncol,pver)
   real(kind_phys) :: vtrdg(ncol,pver)
   real(kind_phys) :: ttrdg(ncol,pver)

   ! Energy change used by fixer.
   real(kind_phys) :: de(ncol)

   character(len=1) :: cn
   character(len=9) :: fname(4)
   !----------------------------------------------------------------------------

   errmsg =''
   errflg = 0

   ! Allocate wavenumber fields.
   allocate(tau(ncol,band_oro%ngwv:band_oro%ngwv,pver+1),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of tau failed', errmsg )
   allocate(gwut(ncol,pver,band_oro%ngwv:band_oro%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of gwut failed', errmsg )
   allocate(phase_speeds(ncol,band_oro%ngwv:band_oro%ngwv),stat=stat)
     call handle_err( stat, errflg,sub//': Allocate of phase_speeds failed', errmsg )

   ! initialize accumulated momentum fluxes and tendencies
   taurx = 0._kind_phys
   taury = 0._kind_phys
   ttrdg = 0._kind_phys
   utrdg = 0._kind_phys
   vtrdg = 0._kind_phys
   tau_diag = -9999._kind_phys

   do nn = 1, n_rdg
      kwvrdg  = 0.001_kind_phys / ( hwdth(:,nn) + 0.001_kind_phys ) ! this cant be done every time step !!!
      isoflag = 0
      effgw   = effgw_rdg * ( hwdth(1:ncol,nn)* clngt(1:ncol,nn) ) / gbxar(1:ncol)
      effgw   = min( effgw_rdg_max , effgw )

      call gw_rdg_src(ncol, pver, band_oro, p, &
         u, v, t, mxdis(:,nn), angll(:,nn), anixy(:,nn), kwvrdg, isoflag, zi, nm, &
         src_level, tend_level, bwv_level, tlb_level, tau, ubm, ubi, xv, yv,  &
         ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, phase_speeds)

      call gw_rdg_belowpeak(ncol, pver, band_oro, rdg_cd_llb, &
         t, mxdis(:,nn), anixy(:,nn), kwvrdg, &
         zi, nm, ni, rhoi, &
         src_level, tau, &
         ubmsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, &
         tauoro, taudsw, hdspwv, hdspdw)


      call gw_rdg_break_trap(ncol, pver, band_oro, &
         zi, nm, ni, ubm, ubi, rhoi, kwvrdg , bwv, tlb, wbr, &
         src_level, tlb_level, hdspwv, hdspdw,  mxdis(:,nn), &
         tauoro, taudsw, tau, &
         ldo_trapped_waves=trpd_leewv)


      call gw_drag_prof(ncol, band_oro, p, src_level, tend_level, dt, &
         t, vramp,    &
         piln, rhoi, nm, ni, ubm, ubi, xv, yv,   &
         effgw, phase_speeds, kvtt, q, dse, tau, utgw, vtgw, &
         ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, &
         kwvrdg=kwvrdg, &
         satfac_in = 1._kind_phys, lapply_vdiff=gw_rdg_do_vdiff , tau_diag=tau_diag )

      ! Add the tendencies from each ridge to the totals.
      do k = 1, pver
         ! diagnostics
         utrdg(:,k) = utrdg(:,k) + utgw(:,k)
         vtrdg(:,k) = vtrdg(:,k) + vtgw(:,k)
         ttrdg(:,k) = ttrdg(:,k) + ttgw(:,k)
         ! physics tendencies
         u_tend(:ncol,k) = u_tend(:ncol,k) + utgw(:,k)
         v_tend(:ncol,k) = v_tend(:ncol,k) + vtgw(:,k)
         s_tend(:ncol,k) = s_tend(:ncol,k) + ttgw(:,k)
      end do

      do m = 1, pcnst
         do k = 1, pver
            q_tend(:ncol,k,m) = q_tend(:ncol,k,m) + qtgw(:,k,m)
         end do
      end do

      do k = 1, pver+1
         taurx0(:,k) = tau(:,0,k)*xv
         taury0(:,k) = tau(:,0,k)*yv
         taurx(:,k)  = taurx(:,k) + taurx0(:,k)
         taury(:,k)  = taury(:,k) + taury0(:,k)
      end do

      if (nn == 1) then
!!$         call outfld('BWV_HT1', bwv,     ncol, lchnk)
!!$         call outfld('TLB_HT1', tlb,     ncol, lchnk)
!!$         call outfld('WBR_HT1', wbr,     ncol, lchnk)
!!$         call outfld('TAUDSW1', taudsw,  ncol, lchnk)
!!$         call outfld('TAUORO1', tauoro,  ncol, lchnk)
!!$         call outfld('UBMSRC1', ubmsrc,  ncol, lchnk)
!!$         call outfld('USRC1',   usrc,    ncol, lchnk)
!!$         call outfld('VSRC1',   vsrc,    ncol, lchnk)
!!$         call outfld('NSRC1'  , nsrc,    ncol, lchnk)
!!$         ! Froude numbers
!!$         call outfld('Fr1_DIAG' , Fr1,    ncol, lchnk)
!!$         call outfld('Fr2_DIAG' , Fr2,    ncol, lchnk)
!!$         call outfld('Frx_DIAG' , Frx,    ncol, lchnk)
!!$         ! Ridge quantities - don't change.  Written for convenience
!!$         call outfld('MXDIS1' , mxdis(:,nn) ,  ncol, lchnk)
!!$         call outfld('ANGLL1' , angll(:,nn) ,  ncol, lchnk)
!!$         call outfld('ANIXY1' , anixy(:,nn) ,  ncol, lchnk)
!!$         call outfld('HWDTH1' , hwdth(:,nn) ,  ncol, lchnk)
!!$         call outfld('CLNGT1' , clngt(:,nn) ,  ncol, lchnk)
!!$         call outfld('GBXAR1' , gbxar ,        ncol, lchnk)
!!$         call outfld('TAUM1_DIAG' , tau_diag ,  ncol, lchnk)
!!$         call outfld('TAU1RDG'//trim(type)//'M', tau(:,0,:),  ncol, lchnk)
!!$         call outfld('UBM1'//trim(type),         ubm,         ncol, lchnk)
!!$         call outfld('UBT1RDG'//trim(type),      gwut,        ncol, lchnk)
      end if

      if (nn <= 6) then
!!$         write(cn, '(i1)') nn
!!$         call outfld('TAU'//cn//'RDG'//trim(type)//'X', taurx0,  ncol, lchnk)
!!$         call outfld('TAU'//cn//'RDG'//trim(type)//'Y', taury0,  ncol, lchnk)
!!$         call outfld('UT'//cn//'RDG'//trim(type),       utgw,    ncol, lchnk)
!!$         call outfld('VT'//cn//'RDG'//trim(type),       vtgw,    ncol, lchnk)
      end if

   end do ! end of loop over multiple ridges

   ! Calculate energy change for output to CAM's energy checker.
   call energy_change(dt, p, u, v, u_tend(:ncol,:), &
          v_tend(:ncol,:), s_tend(:ncol,:), de)
   flx_heat(:ncol) = de

!!$   call outfld('TAUARDG'//trim(type)//'X', taurx,  ncol, lchnk)
!!$   call outfld('TAUARDG'//trim(type)//'Y', taury,  ncol, lchnk)

   if (trim(type) == 'BETA') then
      fname(1) = 'TAUGWX'
      fname(2) = 'TAUGWY'
      fname(3) = 'UTGWORO'
      fname(4) = 'VTGWORO'
   else if (trim(type) == 'GAMMA') then
      fname(1) = 'TAURDGGMX'
      fname(2) = 'TAURDGGMY'
      fname(3) = 'UTRDGGM'
      fname(4) = 'VTRDGGM'
   else
      call handle_err( 1, errflg,sub//'gw_rdg_calc: FATAL: type must be either BETA or GAMMA'&
                  //' type= '//type, errmsg )
   end if

!!$   call outfld(fname(1), taurx(:,pver+1), ncol, lchnk)
!!$   call outfld(fname(2), taury(:,pver+1), ncol, lchnk)
!!$   call outfld(fname(3), utrdg,  ncol, lchnk)
!!$   call outfld(fname(4), vtrdg,  ncol, lchnk)
!!$   call outfld('TTGWORO', ttrdg / cpair,  ncol, lchnk)

   deallocate(tau, gwut, phase_speeds)

end subroutine gw_rdg_calc

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
  call addfld (trim(prefix)//'UTGWSPEC',(/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency - gravity wave spectrum')
  call addfld (trim(prefix)//'VTGWSPEC',(/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' V tendency - gravity wave spectrum')
  call register_vector_field(trim(prefix)//'UTGWSPEC',trim(prefix)//'VTGWSPEC')

  call addfld (trim(prefix)//'TTGWSPEC',(/ 'lev' /), 'A','K s-1', &
       trim(scheme)//' T tendency - gravity wave spectrum')

  ! Wind tendencies broken across five spectral bins.
  call addfld (trim(prefix)//'UTEND1',  (/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency   c < -40')
  call addfld (trim(prefix)//'UTEND2',  (/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency  -40 < c < -15')
  call addfld (trim(prefix)//'UTEND3',  (/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency  -15 < c <  15')
  call addfld (trim(prefix)//'UTEND4',  (/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency   15 < c <  40')
  call addfld (trim(prefix)//'UTEND5',  (/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency   40 < c ')

  ! Reynold's stress toward each cardinal direction, and net zonal stress.
  call addfld (trim(prefix)//'TAUE' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Eastward Reynolds stress')
  call addfld (trim(prefix)//'TAUW' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Westward Reynolds stress')
  call addfld (trim(prefix)//'TAUNET' , (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' E+W Reynolds stress')
  call addfld (trim(prefix)//'TAUN' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Northward Reynolds stress')
  call addfld (trim(prefix)//'TAUS' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Southward Reynolds stress')

  ! Momentum flux in each direction.
  call addfld (trim(prefix)//'EMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Eastward MF')
  call addfld (trim(prefix)//'WMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Westward MF')
  call addfld (trim(prefix)//'NMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Northward MF')
  call addfld (trim(prefix)//'SMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Southward MF')

  ! Temperature tendency terms.
  call addfld (trim(prefix)//'TTGWSDF' , (/ 'lev' /), 'A','K s-1', &
       trim(scheme)//' t tendency - diffusion term')
  call addfld (trim(prefix)//'TTGWSKE' , (/ 'lev' /), 'A','K s-1', &
       trim(scheme)//' t tendency - kinetic energy conversion term')

  ! Gravity wave source spectra by wave number.
  do l=-band%ngwv,band%ngwv
     ! String containing reference speed.
     write (fnum,fmt='(f7.2)') band%cref(l)

     dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
     dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)
     dumc2 = trim(scheme)//" tau at c= "//trim(fnum)//" m s-1"
     call addfld (trim(dumc1x),(/ 'lev' /), 'A','Pa',dumc2)
     call addfld (trim(dumc1y),(/ 'lev' /), 'A','Pa',dumc2)

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
  real(kind_phys), intent(in) :: phase_speeds(ncol,-band%ngwv:band%ngwv)
  ! Winds at cell midpoints.
  real(kind_phys), intent(in) :: u(:,:)
  real(kind_phys), intent(in) :: v(:,:)
  ! Unit vector in the direction of wind at source level.
  real(kind_phys), intent(in) :: xv(:)
  real(kind_phys), intent(in) :: yv(:)
  ! Wind tendency for each wave.
  real(kind_phys), intent(in) :: gwut(ncol,pver,-band%ngwv:band%ngwv)
  ! Temperature tendencies from diffusion and kinetic energy.
  real(kind_phys), intent(in) :: dttdf(:,:)
  real(kind_phys), intent(in) :: dttke(:,:)
  ! Wave Reynolds stress.
  real(kind_phys), intent(in) :: tau(ncol,-band%ngwv:band%ngwv,pver)
  ! Zonal and meridional total wind tendencies.
  real(kind_phys), intent(in) :: utgw(:,:)
  real(kind_phys), intent(in) :: vtgw(:,:)
  ! Temperature tendencies.
  real(kind_phys), intent(in) :: ttgw(:,:)
  ! Reynolds stress for waves propagating in each cardinal direction.
  real(kind_phys), intent(in) :: taucd(ncol,pver+1,4)

  character(len=*), parameter :: sub = 'gw_spec_outflds'
  ! Indices
  integer :: i, k, l
  integer :: ix(ncol, -band%ngwv:band%ngwv), iy(ncol, -band%ngwv:band%ngwv)
  integer :: iu(ncol), iv(ncol)

  ! Zonal wind tendency, broken up into five bins.
  real(kind_phys) :: utb(ncol, pver, 5)
  ! Definition of the bin boundaries.
  real(kind_phys), parameter :: bounds(4) = (/ -40._kind_phys, -15._kind_phys, &
       15._kind_phys, 40._kind_phys /)

  ! Momentum flux in the four cardinal directions.
  real(kind_phys) :: mf(ncol, pver, 4)

  ! Wave stress in zonal/meridional direction
  real(kind_phys) :: taux(ncol,-band%ngwv:band%ngwv,pver)
  real(kind_phys) :: tauy(ncol,-band%ngwv:band%ngwv,pver)

  ! Temporaries for output
  real(kind_phys) :: dummyx(ncol,pver)
  real(kind_phys) :: dummyy(ncol,pver)
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
           utb(i,k,ix(i,l)) = utb(i,k,ix(i,l)) + gwut(i,k,l)
        end do
     end do
  end do

  ! Find just the zonal part.
  do l = 1, 5
     do k = 1, pver
        utb(:, k, l) = utb(:, k, l) * xv
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
  do l=-band%ngwv,band%ngwv
     ix(:,l) = c_to_l(phase_speeds(:,l)*xv)
     iy(:,l) = c_to_l(phase_speeds(:,l)*yv)
  end do

  ! Find projection of tau.
  do k = 1, pver
     do l = -band%ngwv,band%ngwv
        do i = 1, ncol
           taux(i,ix(i,l),k) = taux(i,ix(i,l),k) &
                + abs(tau(i,l,k)*xv(i))
           tauy(i,iy(i,l),k) = tauy(i,iy(i,l),k) &
                + abs(tau(i,l,k)*yv(i))
        end do
     end do
  end do

  do l=-band%ngwv,band%ngwv

     dummyx = taux(:,l,:)
     dummyy = tauy(:,l,:)

     dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
     dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)

!!$     call outfld(dumc1x,dummyx,ncol,lchnk)
!!$     call outfld(dumc1y,dummyy,ncol,lchnk)

  enddo


  ! Output momentum flux in each cardinal direction.
  mf = 0._kind_phys

  do k = 1, pver

     ! Convert wind speed components to wavenumber indices.
     iu = c_to_l(u(:,k))
     iv = c_to_l(v(:,k))

     ! Sum tau components in each cardinal direction.
     ! Split west/east and north/south based on whether wave speed exceeds
     ! wind speed.
     do l = -band%ngwv, band%ngwv

        where (iu > l)
           mf(:,k,west) = mf(:,k,west) + taux(:,l,k)
        elsewhere
           mf(:,k,east) = mf(:,k,east) + taux(:,l,k)
        end where

        where (iv > l)
           mf(:,k,south) = mf(:,k,south) + tauy(:,l,k)
        elsewhere
           mf(:,k,north) = mf(:,k,north) + tauy(:,l,k)
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

    l = min( max(int(c/band%dc),-band%ngwv), band%ngwv )

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

  write(num_str,'(I2.2)') abs(l)

  tau_fld_name = trim(tau_fld_name)//num_str

end function tau_fld_name

end module gw_drag
