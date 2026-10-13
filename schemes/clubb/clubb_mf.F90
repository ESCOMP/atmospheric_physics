module clubb_mf

! =============================================================================== !
! Mass-flux module for use with CLUBB                                             !
! Together (CLUBB+MF) they comprise a eddy-diffusivity mass-flux approach (EDMF)  !
! =============================================================================== !
  use ccpp_kinds,    only: kind_phys

  implicit none
  private
  save

  public :: clubb_mf_init, &
            integrate_mf, &
            do_clubb_mf, &
            do_clubb_mf_diag, &
            clubb_mf_nup, &
            do_clubb_mf_rad, &
            clubb_mf_Lopt, &
            clubb_mf_ddalph, &
            clubb_mf_up_ndt, &
            clubb_mf_cp_ndt, &
            do_clubb_mf_cmt, &
            do_clubb_mf_addtke, &
            clubb_mf_cldfrac_fac

  !
  ! Lopt 0 = fixed L0
  !      1 = tke_clubb L0
  !      2 = wpthlp_clubb L0
  !      3 = test plume L0
  !      6 = ztopm1
  !      7 = rel.hum. at 500 hPa
  !      8 = column int. rel.hum.
  integer  :: clubb_mf_Lopt    = 6
  real(kind_phys) :: clubb_mf_a0      = 0.15_kind_phys
  real(kind_phys) :: clubb_mf_b0      = 1.0_kind_phys
  real(kind_phys) :: clubb_mf_L0      = 200._kind_phys
  real(kind_phys) :: clubb_mf_ent0    = 0.2_kind_phys
  real(kind_phys) :: clubb_mf_alphturb= 3.0_kind_phys
  real(kind_phys) :: clubb_mf_max_L0  = 1.e3_kind_phys
  real(kind_phys) :: clubb_mf_fdd     = 0.5_kind_phys
  real(kind_phys) :: clubb_mf_ddalph  = 12.0_kind_phys
  real(kind_phys) :: clubb_mf_ddbeta  = 1.0_kind_phys
  real(kind_phys) :: clubb_mf_pwfac   = 1.0_kind_phys
  real(kind_phys) :: clubb_mf_ddexp   = 3.0_kind_phys
  real(kind_phys) :: clubb_mf_pwmin   = 1.5_kind_phys
  real(kind_phys) :: clubb_mf_pwmax   = 3.0_kind_phys
  real(kind_phys) :: clubb_mf_cldfrac_fac = 10._kind_phys
  integer  :: clubb_mf_up_ndt  = 1
  integer  :: clubb_mf_cp_ndt  = 1
  integer  :: clubb_mf_kseed = 1
  integer, protected :: clubb_mf_nup     = 0
  ! do_clubb_mf is NOT a namelist variable: it is derived in
  ! clubb_mf_readnl from deep_scheme='CLUBB_MF' (phys_ctl_nl),
  ! making the deep-scheme choice and the MF plume ensemble a single switch
  logical, protected :: do_clubb_mf = .false.
  !
  logical, protected :: do_clubb_mf_diag = .false.
  logical, protected :: do_clubb_mf_rad = .true.
  logical, protected :: do_clubb_mf_addtke = .true.
  logical, protected :: do_clubb_mf_coldpool = .true.
  logical, protected :: do_clubb_mf_ustar = .false.
  logical, protected :: do_clubb_mf_mixd = .false.
  logical, protected :: do_clubb_mf_precip = .true.
  logical, protected :: do_clubb_mf_rhtke = .false.
  logical, protected :: do_clubb_mf_cmt = .false.
  logical, protected :: do_clubb_mf_aloft = .true.
  logical, protected :: do_clubb_mf_coldpool_init = .false.
  logical, protected :: do_clubb_mf_coldpool_perplume = .true.
  logical, protected :: do_clubb_mf_lscale_perplume = .true.
  logical, protected :: do_clubb_mf_aspd = .false.
  logical, protected :: do_clubb_mf_pblcull = .false.
  ! shut off the MF plume ensemble where the lower-tropospheric
  ! inversion is strong (thl700 - thl1000 >= 20 K), mirroring the CLUBB-core
  ! expldiff criterion (advance_clubb_core_module).  Strong-inversion
  ! columns are the marine-Sc regime, where CLUBB should own the boundary
  ! layer; unlike do_clubb_mf_rhtke this is a direct on/off trigger rather
  ! than an indirect entrainment enhancement.
  logical, protected :: do_clubb_mf_invswitch = .true.
  logical :: tht_tweaks = .true.
  integer :: mf_num_cin = 5

  contains

  subroutine clubb_mf_init( &
       clubb_mf_Lopt_in, clubb_mf_a0_in, clubb_mf_b0_in, clubb_mf_L0_in,            &
       clubb_mf_ent0_in, clubb_mf_alphturb_in, clubb_mf_nup_in, clubb_mf_max_L0_in, &
       clubb_mf_fdd_in, clubb_mf_ddalph_in, clubb_mf_ddbeta_in, clubb_mf_pwfac_in,  &
       clubb_mf_ddexp_in, clubb_mf_up_ndt_in, clubb_mf_cp_ndt_in, clubb_mf_kseed_in,&
       clubb_mf_pwmin_in, clubb_mf_pwmax_in, clubb_mf_cldfrac_fac_in,               &
       do_clubb_mf_in, do_clubb_mf_diag_in, do_clubb_mf_precip_in, do_clubb_mf_rad_in, &
       do_clubb_mf_coldpool_in, do_clubb_mf_coldpool_init_in,                       &
       do_clubb_mf_coldpool_perplume_in, do_clubb_mf_lscale_perplume_in,            &
       do_clubb_mf_ustar_in, do_clubb_mf_mixd_in, do_clubb_mf_rhtke_in,             &
       do_clubb_mf_cmt_in, do_clubb_mf_addtke_in, do_clubb_mf_aloft_in,             &
       do_clubb_mf_aspd_in, do_clubb_mf_pblcull_in, do_clubb_mf_invswitch_in,       &
       errmsg, errflg)

    ! Real-valued namelist parameters
    real(kind_phys), intent(in) :: clubb_mf_a0_in, clubb_mf_b0_in, clubb_mf_L0_in
    real(kind_phys), intent(in) :: clubb_mf_ent0_in, clubb_mf_alphturb_in, clubb_mf_max_L0_in
    real(kind_phys), intent(in) :: clubb_mf_fdd_in, clubb_mf_ddalph_in, clubb_mf_ddbeta_in
    real(kind_phys), intent(in) :: clubb_mf_pwfac_in, clubb_mf_ddexp_in
    real(kind_phys), intent(in) :: clubb_mf_pwmin_in, clubb_mf_pwmax_in, clubb_mf_cldfrac_fac_in

    ! Integer-valued namelist parameters
    integer, intent(in) :: clubb_mf_Lopt_in, clubb_mf_nup_in, clubb_mf_up_ndt_in
    integer, intent(in) :: clubb_mf_cp_ndt_in, clubb_mf_kseed_in

    ! Logical switches (do_clubb_mf_in is derived on the CAM side from deep_scheme)
    logical, intent(in) :: do_clubb_mf_in, do_clubb_mf_diag_in, do_clubb_mf_precip_in
    logical, intent(in) :: do_clubb_mf_rad_in, do_clubb_mf_coldpool_in, do_clubb_mf_coldpool_init_in
    logical, intent(in) :: do_clubb_mf_coldpool_perplume_in, do_clubb_mf_lscale_perplume_in
    logical, intent(in) :: do_clubb_mf_ustar_in, do_clubb_mf_mixd_in, do_clubb_mf_rhtke_in
    logical, intent(in) :: do_clubb_mf_cmt_in, do_clubb_mf_addtke_in, do_clubb_mf_aloft_in
    logical, intent(in) :: do_clubb_mf_aspd_in, do_clubb_mf_pblcull_in, do_clubb_mf_invswitch_in

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! ---- Validation (formerly endrun calls in clubb_mf_readnl) ----
    if (clubb_mf_Lopt_in < 0 .or. clubb_mf_Lopt_in == 4 .or. &
        clubb_mf_Lopt_in == 5 .or. clubb_mf_Lopt_in > 8) then
       errmsg = 'clubb_mf_init: clubb_mf_Lopt value must be one of 0,1,2,3,6,7,8'
       errflg = 1
       return
    end if

    if ((.not. do_clubb_mf_in) .and. do_clubb_mf_diag_in) then
       errmsg = 'clubb_mf_init: cannot turn on do_clubb_mf_diag without also turning on do_clubb_mf'
       errflg = 1
       return
    end if

    ! ---- Set module variables ----
    clubb_mf_Lopt     = clubb_mf_Lopt_in
    clubb_mf_a0       = clubb_mf_a0_in
    clubb_mf_b0       = clubb_mf_b0_in
    clubb_mf_L0       = clubb_mf_L0_in
    clubb_mf_ent0     = clubb_mf_ent0_in
    clubb_mf_alphturb = clubb_mf_alphturb_in
    clubb_mf_nup      = clubb_mf_nup_in
    clubb_mf_max_L0   = clubb_mf_max_L0_in
    clubb_mf_fdd      = clubb_mf_fdd_in
    clubb_mf_ddalph   = clubb_mf_ddalph_in
    clubb_mf_ddbeta   = clubb_mf_ddbeta_in
    clubb_mf_pwfac    = clubb_mf_pwfac_in
    clubb_mf_ddexp    = clubb_mf_ddexp_in
    clubb_mf_up_ndt   = clubb_mf_up_ndt_in
    clubb_mf_cp_ndt   = clubb_mf_cp_ndt_in
    clubb_mf_kseed    = clubb_mf_kseed_in
    clubb_mf_pwmin    = clubb_mf_pwmin_in
    clubb_mf_pwmax    = clubb_mf_pwmax_in
    clubb_mf_cldfrac_fac = clubb_mf_cldfrac_fac_in

    do_clubb_mf                   = do_clubb_mf_in
    do_clubb_mf_diag              = do_clubb_mf_diag_in
    do_clubb_mf_precip            = do_clubb_mf_precip_in
    do_clubb_mf_rad               = do_clubb_mf_rad_in
    do_clubb_mf_coldpool          = do_clubb_mf_coldpool_in
    do_clubb_mf_coldpool_init     = do_clubb_mf_coldpool_init_in
    do_clubb_mf_coldpool_perplume = do_clubb_mf_coldpool_perplume_in
    do_clubb_mf_lscale_perplume   = do_clubb_mf_lscale_perplume_in
    do_clubb_mf_ustar             = do_clubb_mf_ustar_in
    do_clubb_mf_mixd              = do_clubb_mf_mixd_in
    do_clubb_mf_rhtke             = do_clubb_mf_rhtke_in
    do_clubb_mf_cmt               = do_clubb_mf_cmt_in
    do_clubb_mf_addtke            = do_clubb_mf_addtke_in
    do_clubb_mf_aloft             = do_clubb_mf_aloft_in
    do_clubb_mf_aspd              = do_clubb_mf_aspd_in
    do_clubb_mf_pblcull           = do_clubb_mf_pblcull_in
    do_clubb_mf_invswitch         = do_clubb_mf_invswitch_in

  end subroutine clubb_mf_init

  subroutine integrate_mf( nzm,     nzt,     dtime,                                 & ! input
                           cpair,   gravit,  latvap,  latice,    zvir,              & ! input
                           h2otrip,          pi,                                    & ! input
                           rho_zm,           zm,      p_zm,      iexner_zm,         & ! input
                           rho_zt,  dzt,     zt,      p_zt,      iexner_zt,         & ! input
                           u,       v,       thl,     qt,        thv,               & ! input
                                    w,       th,      qv,        qc,                & ! input
                                             thl_zm,  qt_zm,     thv_zm,            & ! input
                                             th_zm,   qv_zm,     qc_zm,             & ! input
                           ustar,   ths,     wthl_sfc,wqt_sfc,   pblh,              & ! input
                           tke,     tpert,   lts,                                   & ! input
                           wpthlp_env,       wpthvp_env,         wpqtp_env,         & ! input
                           ztopm1,           ddcp,               cbm1,              & ! inout
                           mcape,                                                   & ! output
                           upa,     dna,                                            & ! output
                           upw,     dnw,                                            & ! output
                           upmf,                                                    & ! output
                           upqt,    dnqt,                                           & ! output
                           upthl,   dnthl,                                          & ! output
                           upthv,   dnthv,                                          & ! output
                           upth,    dnth,                                           & ! output
                           upqc,    dnqc,                                           & ! output
                           upbuoy,                                                  & ! output
                           upent,                                                   & ! output
                           updet,                                                   & ! output
                           dry_a,   moist_a,                                        & ! output
                           dry_w,   moist_w,                                        & ! output
                           dry_qt,  moist_qt,                                       & ! output
                           dry_thl, moist_thl,                                      & ! output
                           dry_u,   moist_u,                                        & ! output
                           dry_v,   moist_v,                                        & ! output
                                    moist_qc,                                       & ! output
                           ae,                                                      & ! output
                           ac,      aup,     adn,                                   & ! output
                           aw,      awup,    awdn,                                  & ! output
                           aww,     awwup,   awwdn,                                 & ! output
                           awthlup, awqtup,  awuup, awvup,                          & ! output
                           awthldn, awqtdn,  awudn, awvdn,                          & ! output
                           awthl,   awqt,                                           & ! output
                           awu,     awv,                                            & ! output
                           thlflxup,qtflxup, uflxup, vflxup,                        & ! output
                           thlflxdn,qtflxdn, uflxdn, vflxdn,                        & ! output
                           thlflx,  qtflx,   uflx,   vflx,                          & ! output - variables needed for solver
                           thvflx,                                                  & ! output
                           sqtup,   sthlup,                                         & ! output
                           sqtdn,   sthldn,                                         & ! output
                           sqt,     sthl,                                           & ! output - variables needed for solver
                           sac,     sev,                                            & ! output - plume autoconversion / evaporation
                           precc,                                                   & ! output
                           ztop,    dynamic_L0,                                     & ! output
                           mfup,    entup,   detup,                                 & ! output
                           mfdn,    entdn,   detdn,                                 & ! output
                           kctop )                                                    ! output

  ! ================================================================================= !
  ! Mass-flux algorithm                                                               !
  !                                                                                   !
  ! Provides rtm and thl fluxes due to mass flux ensemble,                            !
  ! which are fed into the mixed explicit/implicit clubb solver as explicit terms     !
  !                                                                                   !
  ! Mass flux variables are computed on edges (i.e. momentum grid):                   !
  ! upa,upw,upqt,...                                                                  !
  ! dry_a,moist_a,dry_w,moist_w, ...                                                  !
  !                                                                                   !
  ! In CLUBB (unlike CAM) nlevs of momentum grid = nlevs of thermodynamic grid,       !
  ! due to a subsurface thermodynamic layer. To avoid confusion, below the variables  !
  ! are grouped by the grid they are on.                                              !
  !                                                                                   !
  ! *note that state on the lowest thermo level is equal to state on the lowest       !
  ! momentum level due to state_zt(1) = state_zt(2), and lowest momentum level        !
  ! is a weighted combination of the lowest two thermodynamic levels.                 !
  !                                                                                   !
  ! ---------------------------------Authors----------------------------------------  !
  ! Marcin Kurowski, JPL                                                              !
  ! Modified heavily by Mikael Witte, UCLA/JPL for implementation in CESM2/E3SM       !
  ! Additional modifications by Adam Herrington, NCAR                                 !
  ! ================================================================================= !

     use wv_saturation,      only : qsat

     integer,  intent(in)                 :: nzm, nzt
     real(kind_phys), intent(in)                 :: dtime                   ! host-model physics timestep [s]

     real(kind_phys), intent(in)          :: cpair                   ! specific heat of dry air at constant pressure [J kg-1 K-1]
     real(kind_phys), intent(in)          :: gravit                  ! gravitational acceleration [m s-2]
     real(kind_phys), intent(in)          :: latvap                  ! latent heat of vaporization [J kg-1]
     real(kind_phys), intent(in)          :: latice                  ! latent heat of fusion of water [J kg-1]
     real(kind_phys), intent(in)          :: zvir                    ! rh2o/rair - 1 [1]
     real(kind_phys), intent(in)          :: h2otrip                 ! triple point temperature of water [K]
     real(kind_phys), intent(in)          :: pi                      ! pi [1]

     real(kind_phys), dimension(nzt), intent(in) :: u,      v,            & ! thermodynamic grid
                                             thl,    thv,          & ! thermodynamic grid
                                             th,     qv,           & ! thermodynamic grid
                                             qt,     qc,           & ! thermodynamic grid
                                             p_zt,   iexner_zt,    & ! thermodynamic grid
                                             dzt,    rho_zt,       & ! thermodynamic grid
                                             zt

     real(kind_phys), dimension(nzm), intent(in) :: thl_zm, thv_zm,       & ! momentum grid
                                             w,                    &
                                             th_zm,  qv_zm,        &
                                             qt_zm,  qc_zm,        & ! momentum grid
                                             p_zm,   iexner_zm,    & ! momentum grid
                                                     rho_zm,       & ! momentum grid
                                             zm,                   & ! momentum grid
                                             tke,    wpthlp_env,   & ! momentum grid
                                             wpthvp_env, wpqtp_env

     real(kind_phys), intent(in)                :: wthl_sfc,wqt_sfc
     real(kind_phys), intent(in)                :: pblh,tpert
     ! lower-tropospheric stability thl(700 hPa) - thl(1000 hPa) (K),
     ! computed in clubb_intr (kept there verbatim for bit-for-bit
     ! reproducibility of the trigger; see the note at the computation site)
     real(kind_phys), intent(in)                :: lts
     ! the rhinv predictor is computed locally from the column inputs
     ! (previously preprocessed in clubb_intr)
     real(kind_phys)                            :: rhinv, rhlev
     real(kind_phys), dimension(nzt)            :: rhprof
     real(kind_phys)                            :: tmpt_rh, es_rh, qs_rh, tmpmq, tmpmqs
     integer                             :: k_rh
     real(kind_phys), intent(in)                :: ths,ustar
     real(kind_phys), intent(inout)             :: cbm1

     real(kind_phys),dimension(clubb_mf_nup), intent(inout)  :: ztopm1,ddcp

     real(kind_phys),dimension(nzm,clubb_mf_nup), intent(out) :: upa,     & ! momentum grid
                                                          upw,     & ! momentum grid
                                                          upmf,    & ! momentum grid
                                                          upqt,    & ! momentum grid
                                                          upthl,   & ! momentum grid
                                                          upthv,   & ! momentum grid
                                                          upth,    & ! momentum grid
                                                          upqc,    & ! momentum grid
                                                          upbuoy,  & ! momentum grid
                                                          upent,   & ! momentum grid
                                                          updet

     real(kind_phys),dimension(nzm,clubb_mf_nup), intent(out) :: dna,     & ! momentum grid
                                                          dnw,     & ! momentum grid
                                                          dnqt,    & ! momentum grid
                                                          dnthl,   & ! momentum grid
                                                          dnthv,   & ! momentum grid
                                                          dnth,    & ! momentum grid
                                                          dnqc

     real(kind_phys),dimension(nzm), intent(out) :: dry_a,   moist_a,     & ! momentum grid
                                             dry_w,   moist_w,     & ! momentum grid
                                             dry_qt,  moist_qt,    & ! momentum grid
                                             dry_thl, moist_thl,   & ! momentum grid
                                             dry_u,   moist_u,     & ! momentum grid
                                             dry_v,   moist_v,     & ! momentum grid
                                                      moist_qc       ! momentum grid

     real(kind_phys),dimension(nzm), intent(out) :: ae,                                &
                                             ac,      aup,     adn,              &
                                             aw,      awup,    awdn,             &
                                             aww,     awwup,  awwdn,             &
                                             awthlup, awqtup, awuup, awvup,      & ! momentum grid
                                             awthldn, awqtdn, awudn, awvdn,      & ! momentum grid
                                             awthl,   awqt,                      & ! momentum grid
                                             awu,     awv,                       & ! momentum grid
                                             thlflxup,qtflxup, uflxup, vflxup,   & ! momentum grid
                                             thlflxdn,qtflxdn, uflxdn, vflxdn,   & ! momentum grid
                                             thlflx,  qtflx,   uflx,   vflx,     & ! momentum grid
                                             thvflx,  precc

     real(kind_phys),dimension(nzt), intent(out) :: sqtup,   sthlup,                   & ! thermodynamic grid
                                             sqtdn,   sthldn,                    & ! thermodynamic grid
                                             sqt,     sthl,                      & ! thermodynamic grid
                                             sac,     sev                          ! thermodynamic grid

     ! ensemble plume mass flux [kg/m2/s], fractional entrainment/detrainment [1/m]
     ! for updrafts/downdrafts (momentum grid), and ensemble plume-top index
     ! (counted in momentum interfaces from the surface, orientation independent)
     ! needed to populate deep pbuf variables for use by convtran2 and convproc_aer
     real(kind_phys),dimension(nzm), intent(out) :: mfup, entup, detup, &
                                             mfdn, entdn, detdn
     real(kind_phys), intent(out)                :: kctop

     real(kind_phys),dimension(clubb_mf_nup), intent(out) :: ztop, dynamic_L0, mcape

     ! =============================================================================== !
     ! INTERNAL VARIABLES
     !

     ! =============================================================================== !
     ! GRID ORIENTATION GENERALIZATION VARIABLES
     ! ------------------------------------------------------------------------------- !
     ! To support both top-down (CAM) and bottom-up (CLUBB) grid orientations without
     ! duplicating code, these variables abstract the vertical loop bounds and slices.
     !
     ! ksfcm / ksfct : Index of the surface for momentum (m) and thermodynamic (t) grids.
     ! ktopm / ktopt : Index of the model top for momentum (m) and thermodynamic (t) grids.
     ! kdir          : Directional step (+1 for moving up, -1 for moving down).
     !
     ! STAGGERED GRID INDEXING
     ! Because momentum (zm) and thermodynamic (zt) grids are staggered, the relative
     ! index of the cell center (zt) to the interface (zm) flips depending on whether
     ! memory is loaded top-down or bottom-up. These variables dynamically map them:
     !
     ! kt    : The active thermodynamic cell center associated with the current step.
     !         Upward Sweep:   kt = k - (1-kdir)/2
     !         Downward Sweep: kt = k - (1+kdir)/2
     !
     ! kn    : The NEXT momentum interface in the direction of the current sweep.
     !         Upward Sweep:   kn = k + kdir
     !         Downward Sweep: kn = k - kdir
     !
     ! kt_up : The thermodynamic cell center physically ABOVE momentum interface k.
     !         kt_up = k - (1-kdir)/2
     !
     ! kt_dn : The thermodynamic cell center physically BELOW momentum interface k.
     !         kt_dn = k - (1+kdir)/2
     ! =============================================================================== !
     integer :: ksfcm, ktopm, ksfct, ktopt, kdir, kt, kn, kt_up, kt_dn

     ! sums over all plumes
     real(kind_phys), dimension(nzm)              :: moist_th,   dry_th,      &
                                              thl_env,    qt_env,      &
                                              thv_env,                 &
                                              thvflxup,   thvflxdn,    &
                                              awthvup,    awthvdn
     ! updraft properties
     real(kind_phys), dimension(nzm,clubb_mf_nup) :: upqv,     upqs,          & ! momentum grid
                                              upql,     upqi,          & ! momentum grid
                                              upu,      upv,           & ! momentum grid
                                              uplmix                     ! momentum grid
     ! downdraft properties
     real(kind_phys), dimension(nzm,clubb_mf_nup) ::           dnqs,          & ! momentum grid
                                              dnql,     dnqi,          & ! momentum grid
                                              dnu,      dnv,           & ! momentum grid
                                              dnlmix                     ! momentum grid
     ! microphyiscs terms
     real(kind_phys), dimension(nzt,clubb_mf_nup) :: supqt,    supthl,        & ! thermodynamic grid
                                              sdnqt,    sdnthl,        & ! thermodynamic grid
                                              upauto,   upevap           ! thermodynamic grid
     ! ensemble updraft rain evaporation (area-weighted sum over plumes)
     real(kind_phys), dimension(nzt)              :: sevup                      ! thermodynamic grid
     ! precipitation rates
     real(kind_phys), dimension(nzm,clubb_mf_nup) :: uprr,     dnrr             ! momentum grid
     !
     ! grid-mean (area-weighted) rain reservoirs. uprr/dnrr above are
     ! per-unit-plume-area and drive the local evaporation physics; uprg/dnrg carry
     ! the same budgets in grid-mean units (weighted by the plume areas used in the
     ! sqtup/sqtdn accumulations) so that the area-weighted evaporation can never
     ! exceed the area-weighted rain supply. This guarantees the column integral of
     ! sqt (and hence surface precc / PRECC) is non-negative.
     real(kind_phys), dimension(nzm,clubb_mf_nup) :: uprg,     dnrg             ! momentum grid
     real(kind_phys)                              :: drytot,   lamk,   sdnmax
     real(kind_phys), parameter                   :: mf_dry_frac_max = 0.9_kind_phys   ! max fractional total water removable per timestep
     !
     ! entrainment profiles
     real(kind_phys), dimension(nzt,clubb_mf_nup) :: entf,     mix              ! thermodynamic grid
     integer,  dimension(nzt,clubb_mf_nup) :: enti                       ! thermodynamic grid
     !
     ! other variables
     integer                              :: k,i,kstart,ddtopm,kcb,kpbl,kmid,nbot
     integer,  dimension(clubb_mf_nup)    :: ddbotm,kcbarr
     real(kind_phys), dimension(clubb_mf_nup)    :: zcb,cpfac
     real(kind_phys)                             :: zcb_unset,                &
                                             wthv_sfc, wthv,   wqt,    &
                                             ddint,   iddcp,   &
                                             wstar,  qstar,   thvstar, &
                                             sigmaw, sigmaqt, sigmathv,&
                                             convh,  wmin,    wmax,    &
                                             wlv,    wtv,     wp,      &
                                             B,                        & ! thermodynamic grid
                                             entexp, entexpu, entw,    & ! thermodynamic grid
                                             Mn,                       & ! momentum grid
                                             eturb,  det,     lmixt,   & ! thermodynamic grid
                                             qtovqs, sevap,   taum1,   & ! thermodynamic grid
                                             sqtint, sthlint, alphint, &
                                             qtmp,   betathl, betaqt,  & ! thermodynamic grid
                                             thln,   thvn,    thn,     & ! momentum grid
                                             qtn,    qsn,              & ! momentum grid
                                             qcn,    qln,     qin,     & ! momentum grid
                                             un,     vn,      wn2,     & ! momentum grid
                                             wn,                       & ! momentum grid
                                             lmixn,   srfarea,         & ! momentum grid
                                             srfwqtu, srfwthvu,        &
                                             facqtu,  facthvu,         &
                                             zsub,    wcb,    rh_L0,   &
                                             dzext

     ! limit convective area
     logical                                :: limarea = .false.
     real(kind_phys),parameter                     :: amax = 0.6_kind_phys

     ! buoyancy sorting variables
     logical                                :: bsort = .false.
     real(kind_phys),parameter                     :: rle = 0.1_kind_phys
     integer                                :: niter_xc = 1
     integer                                :: iter_xc
     real(kind_phys)                               :: es, entn, detn, mfn, ee2, ud2
     real(kind_phys)                               :: thln0,   qtn0,    wn0


     ! aloft trigger flag
     logical                                :: aloft = .false.

     ! strong-inversion (marine-Sc) plume inhibition (do_clubb_mf_invswitch)
     logical                                :: mf_inhibit
     ! rng seed
     real(kind_phys), dimension(4)                 :: u_seed

     ! alpha relates star quantities to stddev Suselj etal 2019 DOI: 10.1175/JAS-D-18-0239.1
     real(kind_phys),parameter                     :: alphw   = 0.572_kind_phys,        &
                                               alphqt  = 2.890_kind_phys,        &
                                               alphthv = 2.890_kind_phys
     ! w' covariance Suselj etal 2019 2019 DOI: 10.1175/JAS-D-18-0239.1
     real(kind_phys),parameter                     :: cwqt  = 0.32_kind_phys,           &
                                               cwthv = 0.58_kind_phys
     ! virtual mass coefficients for w-eqn Suselj etal 2019 DOI: 10.1175/JAS-D-18-0239.1
     real(kind_phys),parameter                     :: wa = 1.0_kind_phys,               &
                                               wb = 1.5_kind_phys
     ! min values to avoid singularities
     real(kind_phys),parameter                     :: wstarmin = 1.e-3_kind_phys,       &
                                               pblhmin  = 100._kind_phys
     ! evaporation efficiency Suselj etal 2019 DOI: 10.1175/JAS-D-18-0239.1
     real(kind_phys),parameter                     :: ke = 2.5e-4_kind_phys
     !
     ! height where downdrafts feel the surface
     real(kind_phys),parameter                     :: z00dn = 1.e3_kind_phys, &
                                               tinynum = 1.e-7_kind_phys
     ! to fix entrainmnet rate
     logical                                :: fixent = .false.
     !
     ! fixed entrainment rate
     real(kind_phys),parameter                     :: fixent_ent = 2.e-4_kind_phys
     !
     ! Lower limit on entrainment length scale
     real(kind_phys),parameter                     :: min_L0 = 0.5_kind_phys
     !
     ! limiter for tke enahnced fractional entrainment
     real(kind_phys),parameter                     :: max_eturb = 10._kind_phys
     !
     ! floor on the ensemble mass flux for the deep-hookup ratio
     ! diagnostics (entup/detup/entdn/detdn): near plume tops the
     ! mass-flux-weighted entrainment divides by a vanishing mfup while
     ! upent ~ 1/upw diverges, producing immense (finite) ratios
     real(kind_phys),parameter                     :: mf_tiny = 1.e-12_kind_phys
     !
     ! use implicit method for plume updraft velocity
     logical, parameter                     :: do_implicit = .false.
     !
     ! to scale surface fluxes
     logical, parameter                     :: scalesrf = .false.
     !
     ! minimum downdraft speed
     real(kind_phys),parameter                     :: mindnw = 1.E-2_kind_phys
     !
     ! limiter on cold pool effects
     real(kind_phys),parameter                     :: max_cpfac = 5._kind_phys
     !
     ! max limiter on cold pool init effects
     real(kind_phys),parameter                     :: max_cpinit = 0.5_kind_phys
     !
     ! minimum depth of plume expressed as number of levels
     integer,parameter                      :: kspan_min = 4
     !
     ! threshold on surface buoyancy flux below which elevated ("aloft")
     ! convection triggering is used instead of surface-based
     real(kind_phys), parameter :: aloft_wthv_thresh = 0.01_kind_phys

     ! reference pressure (Pa) used to locate the mid-troposphere level
     ! for the aloft-triggering search, and as
     ! the target level for RH interpolation in the rhinv predictor
     real(kind_phys), parameter :: p_midtrop = 50000._kind_phys
     !
     ! floor on the RH-interpolation denominator (qsat)
     real(kind_phys), parameter :: qsat_floor = 1.e-30_kind_phys
     !
     ! cap on interpolated relative humidity before inverting it into rhinv,
     ! to avoid a singularity as RH -> 1
     real(kind_phys), parameter :: rh_cap = 0.990_kind_phys
     !
     ! lower-tropospheric-stability threshold (K) above which the marine-Sc
     ! inversion is judged strong enough to inhibit the MF plume ensemble
     ! entirely (mirrors the CLUBB-core expldiff criterion)
     real(kind_phys), parameter :: lts_inhibit_thresh = 20._kind_phys
     !
     ! vertical extension depth (m) used to stretch an aloft-triggered plume
     ! back down toward the surface when its launch level is within this
     ! depth of the surface
     real(kind_phys), parameter :: aloft_ext_depth = 1000._kind_phys
     !
     ! coefficient and exponent in the RH-based dynamic entrainment-length
     ! diagnostic: rh_L0 = rh_L0_coef * rhinv**rh_L0_exp
     real(kind_phys), parameter :: rh_L0_coef = 50._kind_phys
     real(kind_phys), parameter :: rh_L0_exp  = 3._kind_phys
     !
     ! rh_L0 threshold above which do_clubb_mf_rhtke forces eturb to 1
     ! (TKE-entrainment enhancement switched off)
     real(kind_phys), parameter :: rh_L0_thresh = 733.34_kind_phys
     !
     ! momentum entrains at a reduced rate relative to the thermodynamic
     ! entrainment rate entn -- classic mass-flux assumption
     real(kind_phys), parameter :: ent_mom_reduction = 3._kind_phys
     !
     ! sentinel "not yet reached" value for zcb (cloud-base height);
     ! must be far above any physically possible height
     real(kind_phys), parameter :: zcb_unset_val = 9999999._kind_phys

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! DETERMINE GRID ORIENTATION
     if (zm(1) < zm(nzm)) then
        ! Bottom-up
        ksfcm = 1
        ktopm = nzm
        kdir  = 1
     else
        ! Top-down
        ksfcm = nzm
        ktopm = 1
        kdir  = -1
     end if

     if (zt(1) < zt(nzt)) then
        ksfct = 1
        ktopt = nzt
     else
        ksfct = nzt
        ktopt = 1
     end if

     ! INITIALIZE OUTPUT VARIABLES
     dry_a     = 0._kind_phys
     moist_a   = 0._kind_phys
     dry_w     = 0._kind_phys
     moist_w   = 0._kind_phys
     dry_qt    = 0._kind_phys
     moist_qt  = 0._kind_phys
     dry_thl   = 0._kind_phys
     moist_thl = 0._kind_phys
     dry_u     = 0._kind_phys
     moist_u   = 0._kind_phys
     dry_v     = 0._kind_phys
     moist_v   = 0._kind_phys
     moist_qc  = 0._kind_phys

     ! this is the environmental area - by default 1.
     ae        = 1._kind_phys
     ac        = 0._kind_phys
     aup       = 0._kind_phys
     adn       = 0._kind_phys
     aw        = 0._kind_phys
     awup      = 0._kind_phys
     awdn      = 0._kind_phys
     aww       = 0._kind_phys
     awwup     = 0._kind_phys
     awwdn     = 0._kind_phys
     awuup     = 0._kind_phys
     awvup     = 0._kind_phys
     awudn     = 0._kind_phys
     awvdn     = 0._kind_phys
     awthvup   = 0._kind_phys
     awthvdn   = 0._kind_phys
     awthlup   = 0._kind_phys
     awqtup    = 0._kind_phys
     awthldn   = 0._kind_phys
     awqtdn    = 0._kind_phys
     awthl     = 0._kind_phys
     awqt      = 0._kind_phys
     awu       = 0._kind_phys
     awv       = 0._kind_phys
     thvflxup  = 0._kind_phys
     thlflxup  = 0._kind_phys
     qtflxup   = 0._kind_phys
     thvflxdn  = 0._kind_phys
     thlflxdn  = 0._kind_phys
     qtflxdn   = 0._kind_phys
     thvflx    = 0._kind_phys
     thlflx    = 0._kind_phys
     uflxup    = 0._kind_phys
     vflxup    = 0._kind_phys
     uflxdn    = 0._kind_phys
     vflxdn    = 0._kind_phys
     uflx      = 0._kind_phys
     vflx      = 0._kind_phys
     qtflx     = 0._kind_phys
     sqtup     = 0._kind_phys
     sthlup    = 0._kind_phys
     sqtdn     = 0._kind_phys
     sthldn    = 0._kind_phys
     sqt       = 0._kind_phys
     sthl      = 0._kind_phys
     precc     = 0._kind_phys

     sac       = 0._kind_phys
     sev       = 0._kind_phys
     sevup     = 0._kind_phys

     mfup      = 0._kind_phys
     entup     = 0._kind_phys
     detup     = 0._kind_phys
     mfdn      = 0._kind_phys
     entdn     = 0._kind_phys
     detdn     = 0._kind_phys
     kctop     = 1._kind_phys

     mix       = 0._kind_phys
     entf      = 0._kind_phys
     enti      = 0
     det       = 0._kind_phys

     ! START MAIN COMPUTATION
     upw   = 0._kind_phys
     upth  = 0._kind_phys
     upthl = 0._kind_phys
     upthv = 0._kind_phys
     upqt  = 0._kind_phys
     upa   = 0._kind_phys
     upmf  = 0._kind_phys
     upu   = 0._kind_phys
     upv   = 0._kind_phys
     upqc  = 0._kind_phys
     upth  = 0._kind_phys
     upql  = 0._kind_phys
     upqi  = 0._kind_phys
     upqv  = 0._kind_phys
     upqs  = 0._kind_phys
     upbuoy= 0._kind_phys
     uplmix= 0._kind_phys
     uprr  = 0._kind_phys
     uprg  = 0._kind_phys
     supqt = 0._kind_phys
     supthl= 0._kind_phys
     upent = 0._kind_phys
     updet = 0._kind_phys
     upauto= 0._kind_phys
     upevap= 0._kind_phys

     dnw   = 0._kind_phys
     dna   = 0._kind_phys
     dnu   = 0._kind_phys
     dnv   = 0._kind_phys
     dnqt  = 0._kind_phys
     dnthl = 0._kind_phys
     dnthv = 0._kind_phys
     dnrr  = 0._kind_phys
     dnrg  = 0._kind_phys
     dnth  = 0._kind_phys
     dnqc  = 0._kind_phys
     dnql  = 0._kind_phys
     dnqi  = 0._kind_phys
     dnqs  = 0._kind_phys
     dnlmix= 0._kind_phys
     sdnqt = 0._kind_phys
     sdnthl= 0._kind_phys

     dynamic_L0 = 0._kind_phys
     ztop = 0._kind_phys
     mcape = 0._kind_phys
     ddbotm = 0

     if (bsort) then
       niter_xc = 3
       limarea = .true.
     end if

     ! with the ASPD limiter, area growth is only bounded per layer, so the
     ! accumulated column area must be capped by the existing limarea
     ! machinery (total convective area rescaled to amax)
     if (do_clubb_mf_aspd) limarea = .true.

     ! unique identifier
     zcb_unset = zcb_unset_val !zcb_unset_val = 9999999._kind_phys
     zcb       = zcb_unset

     ! surface buoyancy flux
     !wthv = wthl+zvir*ths*wqt
     wthv_sfc = wthl_sfc+zvir*ths*wqt_sfc

     ! PBL-top index
     kpbl = ksfcm
     do while (zm(kpbl) < pblh .and. kpbl /= ktopm)
       kpbl = kpbl + kdir
     end do

     if (do_clubb_mf_aloft .and. wthv_sfc < aloft_wthv_thresh) then ! aloft_wthv_thresh = 0.01_kind_phys
       aloft = .true.

       kmid = ksfcm
       ! Use a pressure-based criterion to locate the mid-level within the troposphere
       do while (p_zm(kmid) > p_midtrop) !p_midtrop = 50000._kind_phys
         kmid = kmid + kdir
       end do

       ! Search absolute bounds by converting relative slice index
       kstart = maxloc(wpthvp_env(kpbl : kmid : kdir), DIM=1)
       kstart = kpbl + (kstart - 1) * kdir

       wthv = wpthvp_env(kstart)
       wqt = wpqtp_env(kstart)

       if (kstart == ktopm) then
         wthv = 0._kind_phys
         wqt  = 0._kind_phys
       end if

     else
       aloft = .false.
       kstart = ksfcm
       wthv = wthv_sfc
       wqt  = wqt_sfc
     end if

     ! ------------------------------------------------------------------- !
     ! local predictors, computed from the column inputs (previously        !
     ! preprocessed in clubb_intr):                                         !
     !                                                                      !
     ! (1) rhinv, consumed by the do_clubb_mf_rhtke entrainment factor and  !
     !     by get_Lscale under clubb_mf_Lopt 7/8 -- computed only when one  !
     !     of those consumers is active.  Lopt 7 (and the rhtke factor)     !
     !     use RH interpolated to 500 hPa; Lopt 8 uses the column-          !
     !     integrated RH (rho_zt*dzt = dp/g, so the mass weighting matches  !
     !     the previous state%pdel/g integral exactly).                     !
     ! (2) the lower-tropospheric stability thl(700 hPa) - thl(1000 hPa)    !
     !     for the marine-Sc inhibition is computed in clubb_intr and       !
     !     passed in as lts -- exactly the CLUBB-core expldiff criterion    !
     !     (advance_clubb_core_module).                                     !
     ! ------------------------------------------------------------------- !
     rhinv = 0._kind_phys
     if (do_clubb_mf_rhtke .or. clubb_mf_Lopt==7 .or. clubb_mf_Lopt==8) then
       if (clubb_mf_Lopt == 8) then
         tmpmq  = 0._kind_phys
         tmpmqs = 0._kind_phys
         do k_rh = 1, nzt
           ! T from the full potential temperature (matches the state%t the
           ! previous clubb_intr preprocessing used)
           tmpt_rh = th(k_rh)/iexner_zt(k_rh)
           call qsat(tmpt_rh, p_zt(k_rh), es_rh, qs_rh)
           tmpmq  = tmpmq  + rho_zt(k_rh)*dzt(k_rh)*qv(k_rh)
           tmpmqs = tmpmqs + rho_zt(k_rh)*dzt(k_rh)*qs_rh
         end do
         rhlev = tmpmq/max(tmpmqs, qsat_floor) !qsat_floor = 1.e-30_kind_phys
       else
         do k_rh = 1, nzt
           tmpt_rh = th(k_rh)/iexner_zt(k_rh)
           call qsat(tmpt_rh, p_zt(k_rh), es_rh, qs_rh)
           rhprof(k_rh) = qv(k_rh)/max(qs_rh, qsat_floor) !qsat_floor = 1.e-30_kind_phys
         end do
         rhlev = mf_pinterp(nzt, p_zt, rhprof, p_midtrop)  !p_midtrop = 50000._kind_phys
       end if
       if (rhlev >= 1._kind_phys) rhlev = rh_cap ! rh_cap = 0.990_kind_phys
       if (rhlev > 0._kind_phys) rhinv = 1._kind_phys/((1._kind_phys/rhlev) - 1._kind_phys)
     end if

     ! marine-Sc inhibition: mirror the CLUBB-core expldiff
     ! criterion (advance_clubb_core_module: expldiff is applied only where
     ! thlm700 - thlm1000 < 20 K).  Where the lower-tropospheric inversion is
     ! at least that strong -- the marine stratocumulus regime -- the MF plume
     ! ensemble is shut off entirely and CLUBB owns the boundary layer.  All
     ! plume outputs were zeroed above, so the inhibited column follows the
     ! same code path as a negatively buoyant (no-convection) one.
     mf_inhibit = do_clubb_mf_invswitch .and. (lts >= lts_inhibit_thresh) ! lts_inhibit_thresh = 20._kind_phys

     ! Early return if no positive buoyancy and no moisture flux
     if (.not. (wthv > 0._kind_phys .and. wqt > 0._kind_phys .and. (.not. mf_inhibit) ) ) then

       ddcp(:) = 0._kind_phys
       ztopm1(:) = zm(ksfcm)
       cbm1 = zm(ksfcm)
       return

    end if

    ! we have positive surface buoyancy and moisture flux lets do mass-flux

       if (do_clubb_mf_mixd) then
         convh = max(cbm1,pblhmin)
       else
         convh = max(pblh,pblhmin)
       end if

       ! --------------------------------------------------------- !
       ! Initialize using Deardorff convective velocity scale      !
       ! --------------------------------------------------------- !
       wstar = max( wstarmin, (gravit/thv(kstart - (1-kdir)/2)*wthv*convh)**(1._kind_phys/3._kind_phys) )

       ! --------------------------------------------------------- !
       ! Compute cold pool feedback parameter                      !
       ! --------------------------------------------------------- !

       cpfac(:) = 1._kind_phys
       if (do_clubb_mf_coldpool) then
         do i=1,clubb_mf_nup
           cpfac(i) = min( (max(ddcp(i)/wstar,1._kind_phys))**clubb_mf_ddbeta, max_cpfac )
         end do
       end if

       ! --------------------------------------------------------- !
       ! Construct tri-variate PDF at the surface from wstar       !
       ! and initialize plume thv, qt, w                           !
       ! --------------------------------------------------------- !

       if (do_clubb_mf_ustar) then
         qstar   = wqt / max(wstarmin,ustar)
         thvstar = wthv / max(wstarmin,ustar)
       else
         qstar   = wqt / wstar
         thvstar = wthv / wstar
       end if

       do i=1,clubb_mf_nup

         if (do_clubb_mf_coldpool_init) then
           sigmaw   = alphw * wstar * (1._kind_phys + max_cpinit*cpfac(i)/max_cpfac)
           sigmaqt  = alphqt * abs(qstar) * (1._kind_phys + max_cpinit*cpfac(i)/max_cpfac)
           sigmathv = alphthv * abs(thvstar) * (1._kind_phys + max_cpinit*cpfac(i)/max_cpfac)
         else
           sigmaw   = alphw * wstar
           sigmaqt  = alphqt * abs(qstar)
           sigmathv = alphthv * abs(thvstar)
         end if

         wmin = sigmaw * clubb_mf_pwmin
         wmax = sigmaw * clubb_mf_pwmax

         wlv = wmin + (wmax-wmin) / (real(clubb_mf_nup,kind_phys)) * (real(i-1, kind_phys))
         wtv = wmin + (wmax-wmin) / (real(clubb_mf_nup,kind_phys)) * real(i,kind_phys)

         upw(kstart,i) = 0.5_kind_phys * (wlv+wtv)
         upa(kstart,i) = 0.5_kind_phys * erf( wtv/(sqrt(2._kind_phys)*sigmaw) ) &
                    - 0.5_kind_phys * erf( wlv/(sqrt(2._kind_phys)*sigmaw) )

         upu(kstart,i) = u(kstart - (1-kdir)/2)
         upv(kstart,i) = v(kstart - (1-kdir)/2)

         upqt(kstart,i)  = cwqt * upw(kstart,i) * sigmaqt/sigmaw
         upthv(kstart,i) = cwthv * upw(kstart,i) * sigmathv/sigmaw
       enddo

       facqtu=1._kind_phys
       facthvu=1._kind_phys
       if (scalesrf) then
         ! scale surface fluxes
         ! (req'd for conservation if not running with zero-flux b.c.'s)
         srfwqtu = 0._kind_phys
         srfwthvu = 0._kind_phys
         srfarea = 0._kind_phys
         do i=1,clubb_mf_nup
             srfwqtu=srfwqtu+upqt(kstart,i)*upw(kstart,i)*upa(kstart,i)
             srfwthvu=srfwthvu+upthv(kstart,i)*upw(kstart,i)*upa(kstart,i)
             srfarea = srfarea+upa(kstart,i)
         end do
         facqtu=srfarea*wqt/srfwqtu
         facthvu=srfarea*wthv/srfwthvu
       end if

       do i=1,clubb_mf_nup

         kt = kstart - (1-kdir)/2
         kn = kstart + kdir

         betaqt = (qt(kt+2*kdir)-qt(kt))/(0.5_kind_phys*(dzt(kt+2*kdir)+2._kind_phys*dzt(kt+kdir)+dzt(kt)))
         betathl = (thv(kt+2*kdir)-thv(kt))/(0.5_kind_phys*(dzt(kt+2*kdir)+2._kind_phys*dzt(kt+kdir)+dzt(kt)))

         if (.not.aloft) then
           upqt(kstart,i)= qt(kt)-betaqt*dzt(kt)+facqtu*upqt(kstart,i)
           upthv(kstart,i)= thv(kt)-betathl*dzt(kt)+facthvu*upthv(kstart,i)
         else
           upqt(kstart,i)= qt(kt)+upqt(kstart,i)
           upthv(kstart,i)= thv(kt)+upthv(kstart,i)
           if (w(kstart) > 0._kind_phys) upw(kstart,i)= w(kstart)+upw(kstart,i)
         end if

         upthl(kstart,i) = upthv(kstart,i) / (1._kind_phys+zvir*upqt(kstart,i))
         upth(kstart,i)  = upthl(kstart,i)
         upmf(kstart,i) = rho_zm(kstart)*upa(kstart,i)*upw(kstart,i)

         ! get cloud, lowest momentum level
         call condensation_mf(cpair, zvir, h2otrip, latvap, latice, &
                              upqt(kstart,i), upthl(kstart,i), p_zm(kstart), iexner_zm(kstart), &
                              thvn, qcn, thn, qln, qin, qsn, lmixn)
         upthv(kstart,i) = thvn
         upqc(kstart,i)  = qcn
         upql(kstart,i)  = qln
         upqi(kstart,i)  = qin
         upqs(kstart,i)  = qsn
         upth(kstart,i)  = thn
         if (qcn > 0._kind_phys) zcb(i) = zm(kstart)
       end do

       ! if aloft extend the mass flux plume below kstart nbot levels
       if (aloft) then
         zsub = zm(kstart)
         dzext = aloft_ext_depth ! aloft_ext_depth = 1000._kind_phys
         !if greater than dzext above the surface
         if (zsub > dzext) then
             ! find nbot levs below kstart
             nbot = 0
             do k = kstart-kdir, ksfcm, -kdir
               if ((zm(kstart)-zm(k)) < dzext) then
                 nbot = nbot + 1
               end if
             end do
         else
           !else set dxext to height above the suface
           dzext = zm(kstart) - zm(ksfcm)
           nbot = abs(kstart-ksfcm)
         end if

         zsub = zm(kstart)
         do i=1,clubb_mf_nup
           wcb  = upw(kstart,i)
           do k = kstart-kdir, kstart-nbot*kdir, -kdir
             upw(k,i) = wcb - (wcb/(dzext**clubb_mf_ddexp))*(zsub - zm(k))**clubb_mf_ddexp
             upa(k,i) = upa(kstart,i)
             upmf(k,i) = rho_zm(k)*upa(k,i)*upw(k,i)

             upu(k,i) = upu(kstart,i)
             upv(k,i) = upv(kstart,i)
             upqt(k,i) = upqt(kstart,i)
             upthv(k,i) = upthv(kstart,i)
             upthl(k,i) = upthl(kstart,i)
             upth(k,i)  = upth(kstart,i)
             upqc(k,i)  = upqc(kstart,i)
             upql(k,i)  = upql(kstart,i)
             upqi(k,i)  = upqi(kstart,i)
             upqs(k,i)  = upqs(kstart,i)
           end do
         end do
       end if

       do i=1,clubb_mf_nup
         ! --------------------------------------------------------- !
         ! Calculate ztop and dynamic_L based on value of namelist   !
         ! --------------------------------------------------------- !
         call get_Lscale (cpair, gravit, zvir, h2otrip, latvap, latice, &
                          nzt, nzm, zm, tke, wpthlp_env, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, th, &
                          wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, wa, wb,  &
                          qv, p_zt, zt, tpert, pblh, convh, rhinv, ztopm1(i), dynamic_L0(i), ztop(i), mcape(i))

         ! cold pool feedback on the entrainmnet length scale
         dynamic_L0(i) = dynamic_L0(i) * cpfac(i)

         ! limit max/min
         dynamic_L0(i) = max(min_L0,dynamic_L0(i))
         dynamic_L0(i) = min(clubb_mf_max_L0,dynamic_L0(i))

         ! --------------------------------------------------------- !
         ! Stochastic entrainmnet calculation                        !
         ! From Suselj et al 2019 DOI: 10.1175/JAS-D-18-0239.1       !
         ! after Romps and Kuang 2010 DOI: 10.1175/2009JAS3307.1     !
         ! --------------------------------------------------------- !
         do k = ksfct, ktopt-kdir, kdir
           ! get entrainment coefficient, dz/L0
           entf(k,i) = dzt(k) / dynamic_L0(i)
         end do
       end do

       ! get poisson, P(dz/L0)
       ! Grab the u-wind from the exact physical near-surface layers to preserve the PRNG seed
       u_seed(1) = u(ksfct)
       u_seed(2) = u(ksfct + (clubb_mf_kseed + 0)* kdir)
       u_seed(3) = u(ksfct + (clubb_mf_kseed + 1) * kdir)
       u_seed(4) = u(ksfct + (clubb_mf_kseed + 2) * kdir)

       ! get poisson, P(dz/L0)
       ! Ideally we would not fill the enti array above plume tops,
       ! but plume tops have not yet been determined.
       call poisson( nzt, clubb_mf_nup, ksfct, ktopt, kdir, entf, enti, u_seed, pi )

       ! --------------------------------------------------------- !
       ! Main upward sweep to compute updraft properties           !
       !                                                           !
       ! --------------------------------------------------------- !
       do i=1,clubb_mf_nup
         do k = kstart, ktopm-kdir, kdir

           kt = k - (1-kdir)/2
           kn = k + kdir

           ! get microphysics, autoconversion
           if (do_clubb_mf_precip .and. upqc(k,i) > 0._kind_phys) then
             call precip_mf(upqs(k,i),upqt(k,i),upw(k,i),dzt(kt),abs(zm(kn)-zcb(i)),supqt(kt,i))
             supthl(kt,i) = -1._kind_phys*lmixn*supqt(kt,i)*iexner_zt(kt)/cpair
           else
             supqt(kt,i)  = 0._kind_phys
             supthl(kt,i) = 0._kind_phys
           end if

           ! compute mixing rate
           if (fixent) then
             mix(kt,i) = fixent_ent
           else
             ! get entrainment, ent=ent0/dz*P(dz/L0)
             mix(kt,i) = real( enti(kt,i),kind_phys)*clubb_mf_ent0/dzt(kt)
           end if

           !ee2 and the plume's post-entrainment state are mutually dependent
           !iterate for ee2 a few times using successive substitution
           do iter_xc = 1, niter_xc

             if (iter_xc==1) then
               qtn  = upqt(k,i)
               thln = upthl(k,i)
               wn   = upw(k,i)
             else
               qtn  = 0.5_kind_phys*(qtn + qtn0)
               thln = 0.5_kind_phys*(thln + thln0)
               wn = 0.5_kind_phys*(wn + wn0)
             end if

             qtn0  = qtn
             thln0 = thln
             wn0   = wn

             if (bsort) then
               !calculate mixing fractions (ee2, ud2) given a plumes thermodynamic state
               call buoyancy_sort_mixfrac( gravit, cpair, zvir, h2otrip, latvap, latice, &
                                           qtn, thln, wn,                      &
                                           thl_zm(kn), qt_zm(kn), thv_zm(kn),  &
                                           iexner_zm(kn), p_zm(kn), ztopm1(i), &
                                           rle, wa,                            &
                                           ee2, ud2 )
               detn  = mix(kt,i) * ud2
             else !simple closure, mixing fractions are constants.
               ee2 = 1._kind_phys
               ud2 = 1._kind_phys
             end if

             ! entrainment rate
             entn = mix(kt,i) * ee2

             ! --------------------------------------------------------- !
             ! TKE enhanced entrainment                                  !
             ! switches off when dynamic_L0 > max_L0                     !
             ! --------------------------------------------------------- !
             eturb = (1._kind_phys + clubb_mf_alphturb*sqrt(tke(k))/upw(k,i))
             if (do_clubb_mf_rhtke) then
               rh_L0 = rh_L0_coef * (rhinv**rh_L0_exp) ! rh_L0_coef = 50._kind_phys, rh_L0_exp  = 3._kind_phys
               if (rh_L0 >= rh_L0_thresh) eturb = 1._kind_phys ! rh_L0_thresh = 733.34_kind_phys
             else
               if (dynamic_L0(i) >= clubb_mf_max_L0) eturb = 1._kind_phys
             end if
             entn = entn * eturb

             ! integrate updraft
             entexp  = exp(-entn*eturb*dzt(kt))
             entexpu = exp(-entn*dzt(kt)/ent_mom_reduction) ! ent_mom_reduction = 3._kind_phys

             qtn  = qt(kt) *(1._kind_phys-entexp ) + upqt (k,i)*entexp + supqt(kt,i)
             thln = thl(kt)*(1._kind_phys-entexp ) + upthl(k,i)*entexp + supthl(kt,i)
             un   = u(kt)  *(1._kind_phys-entexpu) + upu  (k,i)*entexpu
             vn   = v(kt)  *(1._kind_phys-entexpu) + upv  (k,i)*entexpu

             ! convert source terms to a tendency (convert from S*dz/w to S)
             supqt(kt,i) = supqt(kt,i)*upw(k,i)/dzt(kt)
             upauto(kt,i) = supqt(kt,i)
             supthl(kt,i) = supthl(kt,i)*upw(k,i)/dzt(kt)

             ! get cloud, momentum levels
             call condensation_mf(cpair, zvir, h2otrip, latvap, latice, &
                                  qtn, thln, p_zm(kn), iexner_zm(kn), &
                                  thvn, qcn, thn, qln, qin, qsn, lmixn)
             if (zcb(i)==zcb_unset .and. qcn > 0._kind_phys) zcb(i) = zm(kn)

             ! get buoyancy
             B=gravit*(0.5_kind_phys*(thvn + upthv(k,i))/thv(kt)-1._kind_phys)

             if (do_implicit) then
               wp = clubb_mf_alphturb*wb*entn*sqrt(0.5_kind_phys*(tke(kn)+tke(k)))*dzt(kt)
               wn = (-wp + sqrt(wp**2._kind_phys + (1._kind_phys + 2._kind_phys*wb*entn*dzt(kt))* &
                     (upw(k,i)**2._kind_phys + 2._kind_phys*wa*B*dzt(kt))) )/(1._kind_phys + 2._kind_phys*wb*entn*dzt(kt))
             else
               ! get wn2
               wp = wb*entn*eturb
               if (wp==0._kind_phys) then
                 wn2 = upw(k,i)**2._kind_phys+2._kind_phys*wa*B*dzt(kt)
               else
                 entw = exp(-2._kind_phys*wp*dzt(kt))
                 wn2 = entw*upw(k,i)**2._kind_phys+(1._kind_phys-entw)*wa*B/wp
               end if
               wn = sqrt(max(wn2, 0._kind_phys))
             end if

           end do !iter_xc

           if (wn>0._kind_phys) then

             upthv(kn,i) = thvn
             upthl(kn,i) = thln
             upqt(kn,i)  = qtn
             upqc(kn,i)  = qcn
             upqs(kn,i)  = qsn
             upu(kn,i)   = un
             upv(kn,i)   = vn
             upql(kn,i)  = qln
             upqi(kn,i)  = qin
             upqv(kn,i)  = qtn - qcn
             uplmix(kn,i)= lmixn
             upth(kn,i)  = thn

             if (bsort) then
               mfn = upmf(k,i)*exp( dzt(kt)*( entn - detn ))
               upa(kn,i) = mfn/(wn*rho_zm(kn))
             else
               upa(kn,i) = upa(k,i)
               mfn = rho_zm(kn)*upa(kn,i)*wn
               detn = entn - (mfn - rho_zm(k)*upa(k,i)*upw(k,i)) &
                             /(rho_zm(k)*upa(k,i)*upw(k,i)*dzt(kt))
             end if

             upbuoy(kn,i)= B
             upw(kn,i)   = wn
             upmf(kn,i)  = mfn
             upent(kn,i) = entn
             updet(kn,i) = detn

           else
             ! cull plumes terminating at/below the PBL top (do_clubb_mf_pblcull),
             if ( abs(k-kstart+kdir)<kspan_min .or. &
                  (do_clubb_mf_pblcull .and. (k - kpbl)*kdir <= 0) ) then
               supqt(:,i) = 0._kind_phys
               upauto(:,i)= 0._kind_phys
               supthl(:,i)= 0._kind_phys
               upa(:,i)   = 0._kind_phys
               upbuoy(:,i)= 0._kind_phys
               upw(:,i)   = 0._kind_phys
               upmf(:,i)  = 0._kind_phys
               upent(:,i) = 0._kind_phys
               updet(:,i) = 0._kind_phys
               upthv(:,i) = 0._kind_phys
               upthl(:,i) = 0._kind_phys
               upqt(:,i)  = 0._kind_phys
               upqc(:,i)  = 0._kind_phys
               upqs(:,i)  = 0._kind_phys
               upu(:,i)   = 0._kind_phys
               upv(:,i)   = 0._kind_phys
               upql(:,i)  = 0._kind_phys
               upqi(:,i)  = 0._kind_phys
               upqv(:,i)  = 0._kind_phys
               uplmix(:,i)= 0._kind_phys
               upth(:,i)  = 0._kind_phys
             end if
             ! exit updraft integration
             exit
           end if
         enddo
       enddo

       ! ---------------------------------------------------------- !
       ! In-scheme production limiter                               !
       ! Cap the area-weighted plume drying at each level so the    !
       ! grid-mean microphysical sink cannot remove more than       !
       ! mf_dry_frac_max of the grid-cell total water in a timestep.!
       ! Because the cap is applied BEFORE the rain/evaporation     !
       ! sweeps and the downdraft section, all downstream rain      !
       ! accounting is consistent with the limited production. A    !
       ! common per-level factor lamk preserves the relative        !
       ! weighting of the plumes in the ensemble. The limiter       !
       ! establishes the invariant -sqt*dtime <= mf_dry_frac_max*qt !
       ! (asserted in clubb_intr.F90); the margin below 1.0 leaves  !
       ! room for other same-step moisture sinks and puts the       !
       ! assertion safely beyond roundoff.                          !
       ! ---------------------------------------------------------- !
       if (do_clubb_mf_precip) then
         do k = ksfcm+kdir, ktopm, kdir
           kt = k - (1+kdir)/2
           kn = k - kdir
           drytot = 0._kind_phys
           do i=1,clubb_mf_nup
             drytot = drytot + upa(kn,i)*max(-1._kind_phys*supqt(kt,i), 0._kind_phys)
           end do
           if (drytot*dtime > mf_dry_frac_max*max(qt(kt),0._kind_phys)) then
             lamk = mf_dry_frac_max*max(qt(kt),0._kind_phys)/(drytot*dtime)
             do i=1,clubb_mf_nup
               if (supqt(kt,i) < 0._kind_phys) then
                 supqt(kt,i)  = lamk*supqt(kt,i)
                 supthl(kt,i) = lamk*supthl(kt,i)
                 upauto(kt,i) = lamk*upauto(kt,i)
               end if
             end do
           end if
         end do
       end if

       ! --------------------------------------------------------- !
       ! downward sweep for rain evaporation, snow melting         !
       ! --------------------------------------------------------- !
       if (do_clubb_mf_precip) then
         do i=1,clubb_mf_nup
           do k = ktopm, ksfcm+kdir, -kdir
             kt = k - (1+kdir)/2
             kn = k - kdir

             ! get rain evaporation
             if ((upqs(k,i) + upqs(kn,i)) <=0._kind_phys) then
               qtovqs = 0._kind_phys
             else
               qtovqs = (upqt(k,i) + upqt(kn,i))/(upqs(k,i) + upqs(kn,i))
             end if
             qtovqs = min(1._kind_phys,qtovqs)
             sevap = ke*(1._kind_phys - qtovqs)*sqrt(max(uprr(k,i),0._kind_phys))

             ! limit evaporation to available precip
             sevap = min(sevap,( uprr(k,i)/(rho_zt(kt)*dzt(kt)) - supqt(kt,i)*(1._kind_phys-clubb_mf_fdd) ))

             ! further limit the area-weighted evaporation to the
             ! area-weighted (grid-mean) rain supply. Weight is upa
             ! at the interface below the cell, matching the sqtup accumulation
             if (upa(kn,i) > 0._kind_phys) then
               sevap = min(sevap, uprg(k,i)/(rho_zt(kt)*dzt(kt)*upa(kn,i)) &
                                  - supqt(kt,i)*(1._kind_phys-clubb_mf_fdd) )
             end if
             sevap = max(sevap, 0._kind_phys)

             ! get rain rate
             uprr(kn,i) = uprr(k,i) &
                         - rho_zt(kt)*dzt(kt)*( supqt(kt,i)*(1._kind_phys-clubb_mf_fdd) + sevap )

             ! grid-mean rain reservoir (max() only guards roundoff;
             ! the cap above keeps this non-negative by construction)
             uprg(kn,i) = max( uprg(k,i) &
                         - rho_zt(kt)*dzt(kt)*upa(kn,i)*( supqt(kt,i)*(1._kind_phys-clubb_mf_fdd) + sevap ), 0._kind_phys )

             ! update source terms
             lmixt = 0.5_kind_phys*(uplmix(k,i)+uplmix(kn,i))
             supqt(kt,i) = supqt(kt,i) + sevap
             supthl(kt,i) = supthl(kt,i) - lmixt*sevap*iexner_zt(kt)/cpair
             upevap(kt,i) = sevap
           end do
         end do
       end if

       ! --------------------------------------------------------- !
       ! begin computing downdrafts                                !
       ! --------------------------------------------------------- !
       if (do_clubb_mf_precip .and. clubb_mf_fdd > 0._kind_phys) then

         do i=1,clubb_mf_nup

           ! find cloud base
           do k = ksfcm, ktopm, kdir
             if (upqc(k,i) > 0._kind_phys) then
               ddbotm(i) = k
               exit
             end if
           end do

           ! find cloud top
           ddtopm = 0
           do k = ksfcm, ktopm, kdir
             if (uprr(k,i) > 0._kind_phys) ddtopm = k
           end do

           if (ddtopm /= 0) then
             ! initilaize downdrafts

             ! Kay initializes using negative of the updraft velocity
             ! this causes anomalouly large downdrafts at the initializaiton level
             ! I am intializing with zero velocity as that is more physically defensible
             dnw(ddtopm,i)   = -1._kind_phys*mindnw
             dna(ddtopm,i)   = upa(ddtopm,i)
             dnu(ddtopm,i)   = 0.5_kind_phys*(u(ddtopm - (1-kdir)/2)+u(ddtopm - (1+kdir)/2))
             dnv(ddtopm,i)   = 0.5_kind_phys*(v(ddtopm - (1-kdir)/2)+v(ddtopm - (1+kdir)/2))
             dnqt(ddtopm,i)  = qt_zm(ddtopm)

             ! no cloud in downdrafts, set to cloud free thl
             dnthl(ddtopm,i) = thl_zm(ddtopm)
             dnthv(ddtopm,i) = thv_zm(ddtopm) ! includes condensate loading (!)

             ! get rain generated in the updraft, appropriate it to the downdraft
             ! (we are using an upward sweep notation to be consistent with how uprr(ddtop) was computed)
             dnrr(ddtopm,i)  = -1._kind_phys*dzt(ddtopm - (1-kdir)/2)*rho_zt(ddtopm - (1-kdir)/2)*upauto(ddtopm - (1-kdir)/2,i)*clubb_mf_fdd

             ! grid-mean downdraft rain reservoir seed. area weight is
             ! upa at the interface below the generating cell, which is ddtopm
             ! itself (matching the sqtup accumulation weighting of the
             ! production that this rain came from).
             dnrg(ddtopm,i)  = -1._kind_phys*dzt(ddtopm - (1-kdir)/2)*rho_zt(ddtopm - (1-kdir)/2) &
                               *upa(ddtopm,i)*upauto(ddtopm - (1-kdir)/2,i)*clubb_mf_fdd

             if (fixent) then
               entn = fixent_ent
             else
               ! use deterministic mean entrainment
               entn = clubb_mf_ent0/dynamic_L0(i)
             end if

             ! downdraft qsat
             call qsat(dnthl(ddtopm,i)/iexner_zm(ddtopm),p_zm(ddtopm),es,dnqs(ddtopm,i))

             do k = ddtopm, ksfcm+kdir, -kdir

               kt = k - (1+kdir)/2
               kn = k - kdir

               ! assume fixed area with height
               dna(kn,i) = dna(k,i)

               ! get rain evaporation in integrated form
               taum1 = ke*sqrt(dnrr(k,i))/dnqs(k,i)
               alphint = exp(dzt(kt)*taum1/dnw(k,i))
               sqtint = max( (dnqs(k,i) - dnqt(k,i))*(1._kind_phys - alphint) ,0._kind_phys)

               ! limit to available rain
               sqtint = min( sqtint, -1._kind_phys*dnrr(k,i) / (rho_zt(kt)*dzt(kt)*dnw(k,i)) )
               sthlint = -1._kind_phys*latvap*sqtint*iexner_zt(kt)/cpair

               ! get rain evaporation in tendency form
               sdnqt(kt,i) = max( (dnqs(k,i) - dnqt(k,i))*taum1, 0._kind_phys )

               ! cap the tendency-form downdraft evaporation by the
               ! area-weighted rain actually available to the downdraft.
               ! weights: dna at the interface above the cell
               ! (matching the sqtdn accumulation), upa at the interface below
               ! (matching sqtup).
               if (dna(k,i) > 0._kind_phys) then
                 sdnmax = ( dnrg(k,i)/(rho_zt(kt)*dzt(kt)) &
                            - upa(kn,i)*upauto(kt,i)*clubb_mf_fdd ) / dna(k,i)
                 sdnqt(kt,i) = min( sdnqt(kt,i), max(sdnmax, 0._kind_phys) )
               end if
               sdnthl(kt,i) = -1._kind_phys*latvap*sdnqt(kt,i)*iexner_zt(kt)/cpair

               ! compute rain rate (rain above - evaporation + appropriate updraft rain)
               dnrr(kn,i) = max( dnrr(k,i) &
                                - rho_zt(kt)*dzt(kt)*(sdnqt(kt,i) + upauto(kt,i)*clubb_mf_fdd) , 0._kind_phys )

               ! compute the TKE entrainment enhancement locally for the downdraft.
               eturb = 1._kind_phys + clubb_mf_alphturb*sqrt(tke(k))/abs(dnw(k,i))
               if (do_clubb_mf_rhtke) then
                 rh_L0 = rh_L0_coef*(rhinv**rh_L0_exp)  ! rh_L0_coef = 50._kind_phys, rh_L0_exp  = 3._kind_phys
                 if (rh_L0 >= rh_L0_thresh) eturb = 1._kind_phys ! rh_L0_thresh = 733.34_kind_phys
               else
                 if (dynamic_L0(i) >= clubb_mf_max_L0) eturb = 1._kind_phys
               end if
               eturb = min(eturb, max_eturb)

               entexp  = exp(-1._kind_phys*entn*eturb*dzt(kt))
               entexpu = exp(-1._kind_phys*entn*dzt(kt)/ent_mom_reduction) ! ent_mom_reduction = 3._kind_phys

               ! integrate downward
               dnu(kn,i)   = u(kt)  *(1._kind_phys-entexpu) + dnu  (k,i)*entexpu
               dnv(kn,i)   = v(kt)  *(1._kind_phys-entexpu) + dnv  (k,i)*entexpu
               dnqt(kn,i)  = qt(kt) *(1._kind_phys-entexp ) + dnqt (k,i)*entexp + sqtint
               dnthl(kn,i) = thl(kt)*(1._kind_phys-entexp ) + dnthl(k,i)*entexp + sthlint

               ! get qsat
               call qsat(dnthl(kn,i)/iexner_zm(kn),p_zm(kn),es,dnqs(kn,i))

               ! no supersaturation in downdrafts
               if (dnqt(kn,i) > dnqs(kn,i)) then
                 ! set qt to saturation vapor pressure
                 dnqt(kn,i) = dnqs(kn,i)

                 ! find evaporation that gives saturation vapor pressure
                 sqtint = dnqt(kn,i) - (qt(kt) *(1._kind_phys-entexp ) + dnqt (k,i)*entexp)

                 ! limit to available rain
                 sqtint = min( sqtint, -1._kind_phys*dnrr(k,i) / (rho_zt(kt)*dzt(kt)*dnw(k,i)) )
                 sthlint = -1._kind_phys*latvap*sqtint*iexner_zt(kt)/cpair

                 ! find new evap tendency
                 if ((alphint - 1._kind_phys) /= 0._kind_phys) then
                   qtmp = dnqs(k,i) + sqtint/(alphint - 1._kind_phys)
                   sdnqt(kt,i) = max( (dnqs(k,i) - qtmp)*taum1, 0._kind_phys )
                 else
                   sdnqt(kt,i) = 0._kind_phys
                 end if

                 ! apply the same area-weighted available-rain cap to
                 ! the recomputed (saturation-adjusted) evaporation tendency
                 if (dna(k,i) > 0._kind_phys) then
                   sdnmax = ( dnrg(k,i)/(rho_zt(kt)*dzt(kt)) &
                              - upa(kn,i)*upauto(kt,i)*clubb_mf_fdd ) / dna(k,i)
                   sdnqt(kt,i) = min( sdnqt(kt,i), max(sdnmax, 0._kind_phys) )
                 end if
                 sdnthl(kt,i) = -1._kind_phys*latvap*sdnqt(kt,i)*iexner_zt(kt)/cpair

                 ! re-compute thl with new evaporation rate
                 dnthl(kn,i) = thl(kt)*(1._kind_phys-entexp ) + dnthl(k,i)*entexp + sthlint

                 ! adjust rain
                 dnrr(kn,i) = max( dnrr(k,i) &
                                  - rho_zt(kt)*dzt(kt)*(sdnqt(kt,i) + upauto(kt,i)*clubb_mf_fdd) , 0._kind_phys )
               end if

               ! grid-mean downdraft rain reservoir, using the final
               ! (capped) evaporation tendency. max() only guards roundoff; the
               ! sdnqt cap keeps this non-negative by construction.
               dnrg(kn,i) = max( dnrg(k,i) &
                    - rho_zt(kt)*dzt(kt)*( dna(k,i)*sdnqt(kt,i) &
                                           + upa(kn,i)*upauto(kt,i)*clubb_mf_fdd ), 0._kind_phys )

               ! get virtual temperature
               dnthv(kn,i) = dnthl(kn,i)*(1._kind_phys+zvir*dnqt(kn,i))

               if ((kn - ddbotm(i))*kdir > 0) then
                 ! get virtual temperature
                 dnthv(kn,i) = dnthl(kn,i)*(1._kind_phys+zvir*dnqt(kn,i))

                 ! get buoyancy
                 ! (midpoint k is surrounded by interface k and k-1,
                 ! and therefore we can't compute B at the midpoint properly)
                 B = gravit*(dnthv(kn,i)/thv(kt)-1._kind_phys)

                 ! get wn2
                 wp = wb*entn*eturb  &
                      + clubb_mf_pwfac/( 2._kind_phys*zm(kn + kdir)+tinynum ) * max( 1._kind_phys - exp( zm(kn + kdir)/z00dn-1._kind_phys), 0._kind_phys )
                 if (wp==0._kind_phys) then
                   wn2 = dnw(k,i)**2._kind_phys-2._kind_phys*wa*B*dzt(kt)
                 else
                   entw = exp(-2._kind_phys*wp*dzt(kt))
                   wn2 = entw*dnw(k,i)**2._kind_phys-(1._kind_phys-entw)*wa*B/wp
                 end if
                 wn2 = max(wn2,mindnw**2._kind_phys)
                 dnw(kn,i) = -1._kind_phys*sqrt(wn2)
               else
                 ! analytical taper in sub-cloud layer to
                 ! avoid unphysical flux divergence near surface
                 zsub = zm(ddbotm(i)+kdir)
                 wcb  = dnw(ddbotm(i)+kdir,i)
                 dnw(kn,i) = wcb - (wcb/(zsub**clubb_mf_ddexp))*(zsub - zm(kn))**clubb_mf_ddexp
                 if (dnw(kn,i) >= -1._kind_phys*mindnw) then
                   ! sub-cloud profile has decayed to the velocity floor: the downdraft
                   ! ends here. Levels below keep their zero initialization.
                   dnw(kn,i) = 0._kind_phys
                   dna(kn,i) = 0._kind_phys
                   exit
                 end if
               end if

             end do!k

           end if

         end do!i

       end if
       ! end computing downdrafts

       ! --------------------------------------------------------- !
       ! Arakawa-Schubert positive detrainment limiter             !
       ! --------------------------------------------------------- !
       !   Arakawa-Schubert positive detrainment limiter, repaired:
       ! - the scheme's mass-continuity convention is the detn diagnosis in the
       !   ascent loop: 1/M dM/dz = upent - updet, where upent ALREADY contains
       !   the TKE eturb enhancement.  The old code multiplied by the scalar
       !   eturb again, which both double counted the enhancement and used a
       !   stale value (whatever the last ascent/downdraft iteration left in
       !   it) rather than the (k,i) being limited.
       ! - the entrainment used for the pure-entrainment regrowth is capped at
       !   mix*max_eturb, restoring the documented purpose of max_eturb
       !   (eturb = 1 + alphturb*sqrt(tke)/w diverges as w -> 0 and an uncapped
       !   exp(upent*dz) can produce runaway plume areas).
       ! - guard against M(k)=0 (plume base and aloft-launch levels) and start
       !   at kstart so the aloft downward extension is untouched.
       ! - upmf/updet are kept consistent with the adjusted areas (downstream
       !   ensemble sums and the deep-hookup mass fluxes use upa/upw directly
       !   and are computed after this block).
       if (do_clubb_mf_aspd) then
         do k = kstart, ktopm-kdir, kdir
           kt = k - (1-kdir)/2
           kn = k + kdir
           do i=1,clubb_mf_nup
             if (upw(kn,i)>0._kind_phys .and. upa(k,i)*upw(k,i)>0._kind_phys) then
               ! capped effective entrainment for mass continuity
               entn = min(upent(kn,i), mix(kt,i)*max_eturb)
               ! diagnose detrainment
               Mn = rho_zm(k)*upa(k,i)*upw(k,i)
               det = entn - (rho_zm(kn)*upa(kn,i)*upw(kn,i) - Mn) &
                             /(Mn*dzt(kt))
               if (det < 0._kind_phys) then
                 ! eliminate negative detrainment: grow M by pure entrainment
                 ! and absorb the change into the plume area.
                 ! cap the per-layer growth exponent: where dynamic_L0 is small the capped
                 ! entrainment mix*max_eturb can still reach O(1) m^-1 and exp(entn*dz)
                 ! overflows, producing runaway plume areas (found as area >> 1 in the
                 ! edmf_S_AE diagnostics of the first ASPD test runs)
                 Mn = Mn*exp(min(entn*dzt(kt), 0.5_kind_phys))
                 upa(kn,i)   = Mn/(rho_zm(kn)*upw(kn,i))
                 upmf(kn,i)  = Mn
                 updet(kn,i) = 0._kind_phys
               end if
             end if
           end do
         end do
         ! downdraft counterpart: enforce non-negative detrainment along the
         ! descent.  The downdraft mass entrainment rate mirrors its construction
         ! (per-plume constant entn; fixed dna with dnw from the w equation), so
         ! where |Md| grows downward faster than entrainment allows, cap the
         ! growth and absorb it into dna.
         ! Interaction with the sub-cloud exponential dnw boundary condition
         ! (which keeps flux gradients across the lowest levels gentle and
         ! dz-insensitive): the decaying |Md| tail is positive detrainment, so
         ! the limiter never modifies it; sweeping downward with the current
         ! (already-limited) upper-level dna means any in-cloud dna reduction
         ! cascades smoothly across cloud base (bounded by exp(entn*dz) per
         ! layer) instead of leaving a discontinuity there; and since the cap
         ! only ever REDUCES dna, near-surface flux gradients can only become
         ! gentler.  Do NOT restrict this loop to the in-cloud region: that
         ! would strand reduced in-cloud dna against the unmodified sub-cloud
         ! constant and reintroduce a cloud-base flux-gradient discontinuity.
         if (do_clubb_mf_precip .and. clubb_mf_fdd > 0._kind_phys) then
           do i=1,clubb_mf_nup
             if (fixent) then
               entn = fixent_ent
             else if (dynamic_L0(i) > 0._kind_phys) then
               entn = clubb_mf_ent0/dynamic_L0(i)
             else
               cycle
             end if
             do k = ktopm, ksfcm+kdir, -kdir
               kt = k - (1+kdir)/2
               kn = k - kdir
               if (dna(kn,i)*abs(dnw(kn,i))>0._kind_phys .and. &
                   dna(k,i)*abs(dnw(k,i))>0._kind_phys) then
                 ! mass flux magnitude above (k) and below (kn)
                 Mn = rho_zm(k)*dna(k,i)*abs(dnw(k,i))
                 det = entn - (rho_zm(kn)*dna(kn,i)*abs(dnw(kn,i)) - Mn) &
                               /(Mn*dzt(kt))
                 if (det < 0._kind_phys) then
                   ! same per-layer growth cap as the updraft limiter
                   Mn = Mn*exp(min(entn*dzt(kt), 0.5_kind_phys))
                   dna(kn,i) = Mn/(rho_zm(kn)*abs(dnw(kn,i)))
                 end if
               end if
             end do
           end do
         end if
       end if

       ! --------------------------------------------------------- !
       ! integrate for total convective area                       !
       ! --------------------------------------------------------- !
       do k=ksfcm, ktopm-kdir, kdir
         do i=1,clubb_mf_nup
           aup(k) = aup(k) + upa(k,i)
           adn(k) = adn(k) + dna(k,i)
         end do
         ac(k) = aup(k) + adn(k)
         if (limarea .and. ac(k) > amax) then
           upa(k,:) = upa(k,:)*amax/ac(k)
           ac(k) = amax
         end if
         ae(k) = ae(k) - ac(k)
       end do

       ! --------------------------------------------------------- !
       ! updraft properties for output                             !
       ! --------------------------------------------------------- !
       do k=ksfcm, ktopm, kdir

         ! first sum over all i-updrafts
         do i=1,clubb_mf_nup
           if (upqc(k,i)>0._kind_phys) then
             moist_a(k)   = moist_a(k)   + upa(k,i)
             moist_w(k)   = moist_w(k)   + upa(k,i)*upw(k,i)
             moist_qt(k)  = moist_qt(k)  + upa(k,i)*upqt(k,i)
             moist_thl(k) = moist_thl(k) + upa(k,i)*upthl(k,i)
             moist_u(k)   = moist_u(k)   + upa(k,i)*upu(k,i)
             moist_v(k)   = moist_v(k)   + upa(k,i)*upv(k,i)
             moist_qc(k)  = moist_qc(k)  + upa(k,i)*upqc(k,i)
           else
             dry_a(k)     = dry_a(k)     + upa(k,i)
             dry_w(k)     = dry_w(k)     + upa(k,i)*upw(k,i)
             dry_qt(k)    = dry_qt(k)    + upa(k,i)*upqt(k,i)
             dry_thl(k)   = dry_thl(k)   + upa(k,i)*upthl(k,i)
             dry_u(k)     = dry_u(k)     + upa(k,i)*upu(k,i)
             dry_v(k)     = dry_v(k)     + upa(k,i)*upv(k,i)
           endif
         enddo

         if ( dry_a(k) > 0._kind_phys ) then
           dry_w(k)   = dry_w(k)   / dry_a(k)
           dry_qt(k)  = dry_qt(k)  / dry_a(k)
           dry_thl(k) = dry_thl(k) / dry_a(k)
           dry_u(k)   = dry_u(k)   / dry_a(k)
           dry_v(k)   = dry_v(k)   / dry_a(k)
         else
           dry_w(k)   = 0._kind_phys
           dry_qt(k)  = 0._kind_phys
           dry_thl(k) = 0._kind_phys
           dry_u(k)   = 0._kind_phys
           dry_v(k)   = 0._kind_phys
         endif

         if ( moist_a(k) > 0._kind_phys ) then
           moist_w(k)   = moist_w(k)   / moist_a(k)
           moist_qt(k)  = moist_qt(k)  / moist_a(k)
           moist_thl(k) = moist_thl(k) / moist_a(k)
           moist_u(k)   = moist_u(k)   / moist_a(k)
           moist_v(k)   = moist_v(k)   / moist_a(k)
           moist_qc(k)  = moist_qc(k)  / moist_a(k)
         else
           moist_w(k)   = 0._kind_phys
           moist_qt(k)  = 0._kind_phys
           moist_thl(k) = 0._kind_phys
           moist_u(k)   = 0._kind_phys
           moist_v(k)   = 0._kind_phys
           moist_qc(k)  = 0._kind_phys
         endif

       enddo
       ! --------------------------------------------------------- !
       ! get ensemble mean                                         !
       ! --------------------------------------------------------- !

       ! 1. Momentum Grid Accumulations (Interfaces)
       ! Iterates over all interfaces
       do k=ksfcm, ktopm, kdir
         do i=1,clubb_mf_nup
           awup(k) = awup(k) + upa(k,i)*upw(k,i)
           awdn(k) = awdn(k) + dna(k,i)*dnw(k,i)

           awwup(k) = awwup(k) + upa(k,i)*upw(k,i)*upw(k,i)
           awwdn(k) = awwdn(k) + dna(k,i)*dnw(k,i)*dnw(k,i)

           awuup(k) = awuup(k) + upa(k,i)*upw(k,i)*upu(k,i)
           awudn(k) = awudn(k) + dna(k,i)*dnw(k,i)*dnu(k,i)
           awvup(k) = awvup(k) + upa(k,i)*upw(k,i)*upv(k,i)
           awvdn(k) = awvdn(k) + dna(k,i)*dnw(k,i)*dnv(k,i)

           awthvdn(k)= awthvdn(k)+ dna(k,i)*dnw(k,i)*dnthv(k,i)
           awthldn(k)= awthldn(k)+ dna(k,i)*dnw(k,i)*dnthl(k,i)
           awqtdn(k) = awqtdn(k) + dna(k,i)*dnw(k,i)*dnqt(k,i)

           awthvup(k)= awthvup(k)+ upa(k,i)*upw(k,i)*upthv(k,i)
           awthlup(k)= awthlup(k)+ upa(k,i)*upw(k,i)*upthl(k,i)
           awqtup(k) = awqtup(k) + upa(k,i)*upw(k,i)*upqt(k,i)
         enddo

         aw (k) = awup(k)+ awdn(k)
         aww(k) = awwup(k)+ awwdn(k)
         awu(k) = awuup(k)+ awudn(k)
         awv(k) = awvup(k)+ awvdn(k)
       enddo


       ! 2. Thermodynamic Grid Accumulations (Cells)
       ! loop counter k used to derive
       ! kt_dn (the cell) and
       ! kn (the interface immediately below the cell, toward ksfcm)
       ! Loop counter k range:      Top-down: nzm-1 down to 1.  Bottom-up: 2 up to nzm.
       ! Actual cell range (kt_dn): Top-down: nzm-1 down to 1.  Bottom-up: 1 up to nzm-1.
       do k = ksfcm+kdir, ktopm, kdir
         kn = k - kdir
         kt_dn = k - (1+kdir)/2

         do i=1,clubb_mf_nup
           sqtup(kt_dn)  = sqtup(kt_dn)  + upa(kn,i)*supqt(kt_dn,i)
           sthlup(kt_dn) = sthlup(kt_dn) + upa(kn,i)*supthl(kt_dn,i)

           sqtdn(kt_dn)  = sqtdn(kt_dn)  + dna(k,i)*sdnqt(kt_dn,i)
           sthldn(kt_dn) = sthldn(kt_dn) + dna(k,i)*sdnthl(kt_dn,i)

           sac(kt_dn)   = sac(kt_dn)   + upa(kn,i)*upauto(kt_dn,i)
           sevup(kt_dn) = sevup(kt_dn) + upa(kn,i)*upevap(kt_dn,i)
         end do

         sqt(kt_dn)  = sqtup(kt_dn)  + sqtdn(kt_dn)
         sthl(kt_dn) = sthlup(kt_dn) + sthldn(kt_dn)

         sev(kt_dn) = sevup(kt_dn) + sqtdn(kt_dn)
       enddo

       ! --------------------------------------------------------- !
       ! ensemble mass flux, entrainment and detrainment for the   !
       ! deep-convection hookup.                                   !
       ! Total updraft (downdraft) mass flux =                     !
       !   rho * sum_i a_i*w_i = rho*awup (awdn).                  !
       ! Fractional entrainment is the mass-flux-weighted plume    !
       ! entrainment; detrainment is then derived from discrete    !
       ! mass continuity so that (ent-det) is consistent with      !
       ! d(mf)/dz, which keeps the downstream ZM_MU/EU/DU arrays   !
       ! internally consistent.                                    !
       ! --------------------------------------------------------- !
       do k = ksfcm, ktopm, kdir
         mfup(k) = rho_zm(k)*awup(k)
         mfdn(k) = rho_zm(k)*awdn(k)
         do i=1,clubb_mf_nup
           entup(k) = entup(k) + rho_zm(k)*upa(k,i)*upw(k,i)*upent(k,i)
         enddo
         if (mfup(k) > mf_tiny) then
           entup(k) = entup(k)/mfup(k)
         else
           entup(k) = 0._kind_phys
         end if
       enddo
       do k = ksfcm+kdir, ktopm, kdir
         kn = k - kdir
         kt_dn = k - (1+kdir)/2
         ! updraft: det = ent - d(ln mf)/dz, clipped >= 0
         if (mfup(kn) > mf_tiny .and. mfup(k) > mf_tiny) then
           detup(k) = entup(k) - (mfup(k) - mfup(kn))/(mfup(kn)*dzt(kt_dn))
           if (detup(k) < 0._kind_phys) detup(k) = 0._kind_phys
         end if
       enddo
       do k = ktopm, ksfcm+kdir, -kdir
         kt = k - (1+kdir)/2
         kn = k - kdir
         ! downdraft grows downward: ent from d(|mf|)/dz descending, det clipped
         if (abs(mfdn(k)) > mf_tiny .and. abs(mfdn(kn)) > mf_tiny) then
           entdn(kn) = (abs(mfdn(kn)) - abs(mfdn(k)))/(abs(mfdn(k))*dzt(kt))
           if (entdn(kn) < 0._kind_phys) then
             detdn(kn) = -entdn(kn)
             entdn(kn) = 0._kind_phys
           end if
         end if
       enddo
       ! ensemble plume-top index, counted in momentum interfaces from the
       ! surface (1 = surface interface); orientation independent
       do k = ksfcm, ktopm-kdir, kdir
         if (ac(k) > 0._kind_phys) kctop = real(abs(k-ksfcm)+1, kind_phys)
       enddo
       ! --------------------------------------------------------- !
       ! ztopm1 calculation                                        !
       ! --------------------------------------------------------- !
       do i=1,clubb_mf_nup
         do k=kstart, ktopm, kdir
           ! return if no convection at first level above surface
           if (k == (ksfcm+kdir) .and. ac(k) == 0._kind_phys .and. .not.aloft) then
             sqt(k) = 0._kind_phys
             sthl(k) = 0._kind_phys
             ztopm1(:) = zm(ksfcm)
             ddcp(:) = 0._kind_phys
             cbm1 = zm(ksfcm)
             return
           end if
           ! height of the plume ensemble
           if (do_clubb_mf_lscale_perplume) then
             if ((upa(k,i)+dna(k,i)) > 0._kind_phys) ztopm1(i) = zm(k)
           else
             if (ac(k) > 0._kind_phys) ztopm1(:) = zm(k)
           end if
         end do
       end do

       !subtract init level from ztop for aloft plumes
       if (aloft) then
         ztopm1 = ztopm1 - zm(kstart-nbot*kdir)
         do i=1,clubb_mf_nup
           if (ztopm1(i) < zm(ksfcm)) ztopm1(i) = zm(ksfcm)
         end do
       end if

       ! --------------------------------------------------------- !
       ! cloud base / mixing depth calculation                     !
       ! --------------------------------------------------------- !
       cbm1 = 0._kind_phys
       do i=1,clubb_mf_nup
         kcbarr(i) = 0
         do k=kstart, ktopm, kdir
           if (upqc(k,i) > 0._kind_phys) then
             kcbarr(i) = k
             exit
           end if
         end do

         ! find height of dry plumes
         if (kcbarr(i) == 0) then
           do k=kstart, ktopm, kdir
             if (upw(k,i) <= 0._kind_phys) then
               kcbarr(i) = k
               exit
             end if
           end do
         end if

         ! Edge case where kdbarr is still 0, plume persists through the entire column;
         if (kcbarr(i) == 0) then
            kcbarr(i) = ktopm  ! use the model top as the cap
         end if

         cbm1 = cbm1 + zm(kcbarr(i))

       end do
       cbm1 = cbm1/REAL(clubb_mf_nup,kind_phys)

       ! --------------------------------------------------------- !
       ! bulk downdraft velocity for coldpool parameterization     !
       ! --------------------------------------------------------- !

       ! reset ddcp
       ddcp(:) = 0._kind_phys
       if (do_clubb_mf_coldpool .and. clubb_mf_fdd > 0._kind_phys) then
         ! use single level for cold pool param. note this differs from the sub-cloud mean
         ! in the original implementation (Suselj et al. 2019; DOI: 10.1175/JAS-D-18-0239.1)
         ! because we've implemented a different boundary condition for the downdrafts
         ! to ensure they smoothly decay towards the surface using an analytical function.
         do i=1,clubb_mf_nup
           if (ddbotm(i) == 0) then
             continue
           else
             if (do_clubb_mf_coldpool_perplume) then
               ddcp(i) = -1._kind_phys*dnw(ddbotm(i)+kdir,i)
             else
               ddcp(:) = -1._kind_phys*dna(ddbotm(i)+kdir,i)*dnw(ddbotm(i)+kdir,i) + ddcp(:)
             end if
           end if
         end do
       end if

       ! --------------------------------------------------------- !
       ! downward sweep to get ensemble mean precip                !
       ! --------------------------------------------------------- !
       do k = ktopm, ksfcm+kdir, -kdir
         kt_dn = k - (1+kdir)/2
         precc(k-kdir) = precc(k) - rho_zt(kt_dn)*dzt(kt_dn)*sqt(kt_dn)
       end do

       ! this clamp only removes roundoff (upstream limiters conserve mass)
       precc(ksfcm) = max( precc(ksfcm), 0._kind_phys )

       ! --------------------------------------------------------- !
       ! get turbulent fluxes                                      !
       ! --------------------------------------------------------- !
       kstart = ksfcm + kdir
       if (scalesrf) then
         kstart = ksfcm
       end if

       do k=kstart, ktopm-kdir, kdir

         ! Secure boundary cells to prevent array out-of-bounds on zero-flux boundaries
         kt_up = max(1, min(nzt, k - (1-kdir)/2))
         kt_dn = max(1, min(nzt, k - (1+kdir)/2))

         thvflxup(k)= awthvup(k) - awup(k)*thv(kt_up)
         thlflxup(k)= awthlup(k) - awup(k)*thl(kt_up)
         qtflxup (k)= awqtup (k) - awup(k)*qt(kt_up)

         uflxup  (k)= awuup(k) - awup(k)*u(kt_up)
         vflxup  (k)= awvup(k) - awup(k)*v(kt_up)

         thvflxdn(k)= awthvdn(k) - awdn(k)*thv(kt_dn)
         thlflxdn(k)= awthldn(k) - awdn(k)*thl(kt_dn)
         qtflxdn (k)= awqtdn (k) - awdn(k)*qt(kt_dn)

         uflxdn  (k)= awudn(k) - awdn(k)*u(kt_dn)
         vflxdn  (k)= awvdn(k) - awdn(k)*v(kt_dn)

         thvflx(k)  = thvflxup(k) + thvflxdn(k)
         thlflx(k)  = thlflxup(k) + thlflxdn(k)
         qtflx (k)  = qtflxup (k) + qtflxdn (k)

         uflx(k)    = uflxup(k) + uflxdn(k)
         vflx(k)    = vflxup(k) + vflxdn(k)
       enddo

  end subroutine integrate_mf

  subroutine buoyancy_sort_mixfrac( gravit, cpair, zvir, h2otrip, latvap, latice, &
                                    qtn, thln, wn,                     &
                                    thl_zm_kn, qt_zm_kn, thv_zm_kn,    &
                                    iexner_zm_kn, p_zm_kn, ztopm1_i,   &
                                    rle, wa,                           &
                                    ee2, ud2 )
  ! =============================================================================== !
  ! Buoyancy-sorting neutral mixing fraction, after Bretherton et al 2004.          !
  ! DOI: 10.1175/1520-0493(2004)132<0864:ANPFSC>2.0.CO;2                            !
  ! Iterates the trial (qtn,thln,wn) state internally and returns only the final    !
  ! entrainment/detrainment area-fraction terms (ee2, ud2).                         !
  ! =============================================================================== !

    use wv_saturation, only: qsat

    real(kind_phys), intent(in)  :: gravit, cpair, zvir, h2otrip, latvap, latice
    real(kind_phys), intent(in)  :: qtn, thln, wn   ! this iteration's trial plume state -- caller's job to pick it
    real(kind_phys), intent(in)  :: thl_zm_kn, qt_zm_kn, thv_zm_kn
    real(kind_phys), intent(in)  :: iexner_zm_kn, p_zm_kn
    real(kind_phys), intent(in)  :: ztopm1_i
    real(kind_phys), intent(in)  :: rle, wa

    real(kind_phys), intent(out) :: ee2, ud2

    ! local variables
    integer  :: kk, status
    real(kind_phys) :: qtn0, thln0, wn0, xc
    real(kind_phys) :: tlm, es, qsm, excessm, tln, qsn, excessn
    real(kind_phys) :: thvn, qcn, thn, qln, qin, lmixn
    real(kind_phys) :: cridis, xsat, thlxsat, qtxsat, thvxsat
    real(kind_phys) :: thv_x0, thv_x1, aquad, bquad, cquad, xs1, xs2, x_cu, x_en


    ! --------------------------------------------------------- !
    ! Compute excess water to derive neutral mixing fraction    !
    ! after Bretherton et al 2004                               !
    ! DOI: 10.1175/1520-0493(2004)132<0864:ANPFSC>2.0.CO;2      !
    ! --------------------------------------------------------- !

    ! qexcess of the envrionment
    tlm = thl_zm_kn/iexner_zm_kn
    call qsat(tlm,p_zm_kn,es,qsm)
    excessm = qt_zm_kn - qsm

    ! qexcess in plume
    tln  = thln/iexner_zm_kn
    call qsat(tln,p_zm_kn,es,qsn)
    excessn = qtn - qsn

    call condensation_mf(cpair, zvir, h2otrip, latvap, latice, &
                         qtn, thln, p_zm_kn, iexner_zm_kn, &
                         thvn, qcn, thn, qln, qin, qsn, lmixn)

    ! critical stopping distance
    cridis = rle*ztopm1_i

    ! ----------------------------------------------------------------- !
    ! Case 1 : When both cumulus and env. are unsaturated or saturated. !
    ! ----------------------------------------------------------------- !
    if (excessm*excessn > 0._kind_phys) then
      xc = min(1._kind_phys,max(0._kind_phys,1._kind_phys-2._kind_phys*wa*gravit*cridis/wn**2._kind_phys*(1._kind_phys-thvn/thv_zm_kn)))
    else
      ! -------------------------------------------------- !
      ! Case 2 : When either cumulus or env. is saturated. !
      ! -------------------------------------------------- !
      xsat    = excessn / ( excessn - excessm );
      thlxsat = thln + xsat * ( thl_zm_kn - thln );
      qtxsat  = qtn  + xsat * ( qt_zm_kn - qtn );
      call condensation_mf(cpair, zvir, h2otrip, latvap, latice, &
                           qtxsat, thlxsat, p_zm_kn, iexner_zm_kn, &
                           thvxsat, qcn, thn, qln, qin, qsn, lmixn)
      ! -------------------------------------------------- !
      ! kk=1 : Cumulus Segment, kk=2 : Environment Segment !
      ! -------------------------------------------------- !
      do kk = 1, 2
        if( kk == 1 ) then
          thv_x0 = thvn
          thv_x1 = ( 1._kind_phys - 1._kind_phys/xsat ) * thvn + ( 1._kind_phys/xsat ) * thvxsat
        else
          thv_x1 = thv_zm_kn
          thv_x0 = ( xsat / ( xsat - 1._kind_phys ) ) * thv_zm_kn + ( 1._kind_phys/( 1._kind_phys - xsat ) ) * thvxsat
        endif
        aquad =  wn**2
        bquad =  2._kind_phys*wa*gravit*cridis*(thv_x1 - thv_x0)/thv_zm_kn - 2._kind_phys*wn**2
        cquad =  2._kind_phys*wa*gravit*cridis*(thv_x0 - thv_zm_kn)/thv_zm_kn  + wn**2
        if( kk == 1 ) then
          if( ( bquad**2-4._kind_phys*aquad*cquad ) >= 0._kind_phys ) then
            call roots(aquad,bquad,cquad,xs1,xs2,status)
            x_cu = min(1._kind_phys,max(0._kind_phys,min(xsat,min(xs1,xs2))))
          else
            x_cu = xsat
          endif
        else
          if( ( bquad**2-4._kind_phys*aquad*cquad) >= 0._kind_phys ) then
            call roots(aquad,bquad,cquad,xs1,xs2,status)
            x_en = min(1._kind_phys,max(0._kind_phys,max(xsat,min(xs1,xs2))))
          else
            x_en = 1._kind_phys
          endif
        endif
      enddo
      if( x_cu == xsat ) then
        xc = max(x_cu, x_en)
      else
        xc = x_cu
      endif
    endif

    ee2 = xc**2
    ud2 = 1._kind_phys - 2._kind_phys*xc + xc**2

  end subroutine buoyancy_sort_mixfrac

  ! clamped linear-in-pressure interpolation of a column field
  ! to a target pressure level, mirroring the CLUBB-core pvertinterp
  ! (advance_helper_module): outside the sounding the boundary value is used,
  ! so the criterion degrades gracefully over high terrain exactly as the
  ! CLUBB-core expldiff switch does.  Orientation-agnostic (brackets on the
  ! pressure array itself), so it works for either vertical index convention.
  function mf_pinterp(nz, p, fld, pout) result(val)
    integer,  intent(in) :: nz
    real(kind_phys), intent(in) :: p(nz)     ! level pressures (Pa)
    real(kind_phys), intent(in) :: fld(nz)   ! field on the same levels
    real(kind_phys), intent(in) :: pout      ! target pressure (Pa)
    real(kind_phys)             :: val

    integer  :: k, kbot(1), ktop(1)
    real(kind_phys) :: wgt

    kbot = maxloc(p)
    ktop = minloc(p)
    if (pout >= p(kbot(1))) then
      val = fld(kbot(1))
    else if (pout <= p(ktop(1))) then
      val = fld(ktop(1))
    else
      val = fld(kbot(1))
      do k = 1, nz-1
        if ((pout - p(k))*(pout - p(k+1)) <= 0._kind_phys .and. p(k) /= p(k+1)) then
          wgt = (pout - p(k+1))/(p(k) - p(k+1))
          val = wgt*fld(k) + (1._kind_phys - wgt)*fld(k+1)
          exit
        end if
      end do
    end if

  end function mf_pinterp

  subroutine get_Lscale(cpair, gravit, zvir, h2otrip, latvap, latice, &
                        nzt, nzm, zm, tke, wpthlp_env, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, th, &
                        wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, wa, wb,  &
                        qv, p_zt, zt, tpert, pblh, convh, rhinv, ztopm1, dynamic_L0, ztop, mcape )
  ! --------------------------------------------------------- !
  ! Calculate ztop and dynamic_L based on value of namelist   !
  ! --------------------------------------------------------- !
     integer,  intent(in)                 :: nzt, nzm
     real(kind_phys), dimension(nzt), intent(in) :: thl,    thv,          &
                                             th,                   &
                                             qt,     qv,           &
                                             p_zt,   iexner_zt,    &
                                             dzt,    zt

     real(kind_phys), dimension(nzm), intent(in) :: p_zm,   iexner_zm,    &
                                             zm,     tke,    wpthlp_env

     real(kind_phys), intent(in) ::                wmax,   wmin,       tpert, &
                                            sigmaw, sigmaqt, sigmathv, &
                                            cwqt,   cwthv,  zcb_unset, &
                                            wa,     wb,        ztopm1, &
                                            pblh,   convh,     rhinv

     real(kind_phys), intent(in) :: cpair, gravit, zvir, h2otrip, latvap, latice

     real(kind_phys), intent(out) ::               dynamic_L0, ztop, mcape

     ! local variables
     ! =============================================================================== !
     ! GRID ORIENTATION GENERALIZATION VARIABLES
     ! ------------------------------------------------------------------------------- !
     ! To support both top-down (CAM) and bottom-up (CLUBB) grid orientations without
     ! duplicating code, these variables abstract the vertical loop bounds and slices.
     !
     ! ksfcm / ksfct : Index of the surface for momentum (m) and thermodynamic (t) grids.
     ! ktopm / ktopt : Index of the model top for momentum (m) and thermodynamic (t) grids.
     ! kdir          : Directional step (+1 for moving up, -1 for moving down).
     !
     ! STAGGERED GRID INDEXING
     ! Because momentum (zm) and thermodynamic (zt) grids are staggered, the relative
     ! index of the cell center (zt) to the interface (zm) flips depending on whether
     ! memory is loaded top-down or bottom-up. These variables dynamically map them:
     !
     ! kt    : The active thermodynamic cell center associated with the current step.
     !         Upward Sweep:   kt = k - (1-kdir)/2
     !         Downward Sweep: kt = k - (1+kdir)/2
     !
     ! kn    : The NEXT momentum interface in the direction of the current sweep.
     !         Upward Sweep:   kn = k + kdir
     !         Downward Sweep: kn = k - kdir
     !
     ! kt_up : The thermodynamic cell center physically ABOVE momentum interface k.
     !         kt_up = k - (1-kdir)/2
     !
     ! kt_dn : The thermodynamic cell center physically BELOW momentum interface k.
     !         kt_dn = k - (1+kdir)/2
     ! =============================================================================== !
     integer :: ksfc, ktop, kdir


     real(kind_phys), dimension(nzt)               :: t_zt
     real(kind_phys), dimension(nzt)               :: tp,       qstp
     real(kind_phys), dimension(nzt,1)             :: dmpdz
     real(kind_phys), dimension(1)                 :: tl,                     &
                                               cape,     cin
     integer,  dimension(1)                 :: lcl,      lel
     real(kind_phys)                               :: landfrac
     integer                                :: kpbl,     msg,          &
                                               lon,      mx,           &
                                               k

     ! upper search bound (m) for locating the TKE-gradient or heat-flux-
     ! gradient based L0 diagnostic
     real(kind_phys), parameter :: Lscale_search_top = 20000._kind_phys
     !
     ! TKE vertical-gradient threshold (clubb_mf_Lopt==1) marking the
     ! diagnosed top of the turbulent layer
     real(kind_phys), parameter :: tke_grad_thresh = 1.e-5_kind_phys
     !
     ! heat-flux vertical-gradient threshold (clubb_mf_Lopt==2), same role
     ! as tke_grad_thresh but for the heat-flux-based diagnostic
     real(kind_phys), parameter :: hflux_grad_thresh = 1.e-4_kind_phys
     !
     ! intialize local variables
     cape      = 0._kind_phys
     mcape     = 0._kind_phys
     ztop      = 0._kind_phys
     dmpdz     = 0._kind_phys

     if (zt(1) < zt(nzt)) then
        ksfc = 1
        ktop = nzt
        kdir = 1
     else
        ksfc = nzt
        ktop = 1
        kdir = -1
     end if

     if (clubb_mf_Lopt == 0) then
       !Constant L0
       dynamic_L0 = clubb_mf_L0
       ztop = clubb_mf_L0
     else if (clubb_mf_Lopt == 1) then
       !TKE
       do k = ktop-2*kdir, ksfc+kdir, -kdir
         if (zm(k) < Lscale_search_top .and. tke(k) - tke(k+kdir) > tke_grad_thresh) then ! Lscale_search_top = 20000._kind_phys,tke_grad_thresh = 1.e-5_kind_phys
           ztop = zm(k)
           exit
         endif
       enddo
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)

     else if (clubb_mf_Lopt == 2) then
       !Heat flux
       do k = ktop-2*kdir, ksfc+kdir, -kdir
         if (zm(k) < Lscale_search_top .and. abs(abs(wpthlp_env(k))-abs(wpthlp_env(k-kdir))) > hflux_grad_thresh) then !hflux_grad_thresh = 1.e-4_kind_phys
           ztop = zm(k)
           exit
         endif
       enddo
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)

     else if (clubb_mf_Lopt == 3) then
       !Test plume
       call oneplume( cpair, gravit, zvir, h2otrip, latvap, latice, &
                      nzm, nzt, zm, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, &
                      wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, &
                      wa, wb, tke, do_clubb_mf_precip, ztop )

       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)
     else if (clubb_mf_Lopt == 6) then
       ! grab ztop from max height of ensemble in prior time-step(s)
       ztop = ztopm1
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)
     else if (clubb_mf_Lopt == 7 .or. clubb_mf_Lopt == 8) then
       ztop = rhinv
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)
     end if

  end subroutine get_Lscale

  subroutine condensation_mf( cpair, zvir, h2otrip, latvap, latice, &
                              qt, thl, p, iex, thv, qc, th, ql, qi, qs, lmix )
  ! =============================================================================== !
  ! zero or one condensation for edmf: calculates thv and qc                        !
  ! =============================================================================== !
     use wv_saturation,      only : qsat

     real(kind_phys),intent(in):: cpair,zvir,h2otrip,latvap,latice
     real(kind_phys),intent(in) :: qt,thl,p,iex
     real(kind_phys),intent(out):: thv,qc,th,ql,qi,qs,lmix

     !local variables
     integer  :: niter,i
     real(kind_phys) :: diff,t,qstmp,qcold,es,wf
     logical  :: noice = .true.

     !
     ! max fixed-point iterations for the qc/T condensation solve
     integer,  parameter :: condensation_max_iter = 50
     !
     ! convergence tolerance on successive qc iterates
     real(kind_phys), parameter :: condensation_tol = 2.e-5_kind_phys

     ! max number of iterations
     niter = condensation_max_iter  ! condensation_max_iter = 50
     ! minimum difference
     diff = condensation_tol ! condensation_tol = 2.e-5_kind_phys

     qc=0._kind_phys
     t=thl/iex

     !by definition:
     ! T   = Th*Exner, Exner=(p/p0)^(R/cp)   (1)
     ! Thl = Th - L/cp*ql/Exner              (2)
     !so:
     ! Th  = Thl + L/cp*ql/Exner             (3)
     ! T   = Th*Exner=(Thl+L/cp*ql/Exner)*Exner    (4)
     !     = Thl*Exner + L/cp*ql
     do i=1,niter

       if (noice) then
         wf = 1._kind_phys
       else
         wf = get_watf(h2otrip,t)
       end if
       t = thl/iex+get_alhl(latvap,latice,wf)/cpair*qc   !as in (4)

       ! qsat, p is in pascal
       call qsat(t,p,es,qstmp)
       qcold = qc
       qc = max(0.5_kind_phys*qc+0.5_kind_phys*(qt-qstmp),0._kind_phys)
       if (abs(qc-qcold)<diff) exit
     enddo

     if (noice) then
       wf = 1._kind_phys
     else
       wf = get_watf(h2otrip,t)
     end if
     t = thl/iex+get_alhl(latvap,latice,wf)/cpair*qc

     call qsat(t,p,es,qs)
     qc = max(qt-qs,0._kind_phys)
     thv = (thl+get_alhl(latvap,latice,wf)/cpair*iex*qc)*(1._kind_phys+zvir*(qt-qc)-qc)
     lmix = get_alhl(latvap,latice,wf)
     th = t*iex
     qi = qc*(1._kind_phys-wf)
     ql = qc*wf

     contains

     function get_watf(h2otrip,t)
       real(kind_phys), intent(in) :: h2otrip
       real(kind_phys)            :: t,get_watf,tc
       real(kind_phys), parameter :: &
                              tmax=-10._kind_phys, &
                              tmin=-40._kind_phys

       tc=t-h2otrip

       if (tc>tmax) then
         get_watf=1._kind_phys
       else if (tc<tmin) then
         get_watf=0._kind_phys
       else
         get_watf=(tc-tmin)/(tmax-tmin);
       end if

     end function get_watf


     function get_alhl(latvap,latice,wf)
     !latent heat of the mixture based on water fraction
       real(kind_phys), intent(in) :: latvap, latice
       real(kind_phys) :: get_alhl,wf

       get_alhl = wf*latvap+(1._kind_phys-wf)*(latvap+latice)

     end function get_alhl

  end subroutine condensation_mf


  subroutine precip_mf(qs,qt,w,dz,dzcld,Supqt)
  !**********************************************************************
  ! Precipitation microphysics
  ! By Adam Herrington, algorithm from Kay Suselj DOI: 10.1175/JAS-D-18-0239.1
  !**********************************************************************

       real(kind_phys),intent(in)  :: qs,qt,w,dz,dzcld
       real(kind_phys),intent(out) :: Supqt
       ! ! local vars
       real(kind_phys)            :: tauwgt, tau,       & ! time-scale vars
                              qstar                ! excess cloud liquid

       real(kind_phys),parameter  :: tau0  = 15._kind_phys,    & ! base time-scale
                              zmin  = 300._kind_phys,   & ! small cloud thick
                              zmax  = 3000._kind_phys,  & ! large cloud thick
                              qcmin = 0.00125_kind_phys   ! supersat threshold

       qstar = qs+qcmin

       if (qt > qstar) then
         ! get precip efficiency
         tauwgt = (dzcld-zmin)/(zmax-zmin)
         tauwgt = min(max(tauwgt,0._kind_phys),1._kind_phys)
         tau    = tauwgt/tau0

         ! get source for updraft
         Supqt = (qstar-qt)*(1._kind_phys - exp(-1._kind_phys*tau*dz/w))
       else
         Supqt = 0._kind_phys
       end if

  end subroutine precip_mf

  subroutine poisson(nz, nup, ksfc, ktop, kdir, lambda, poi, state, pi)
  !**********************************************************************
  ! Set a unique (but reproducible) seed for the kiss RNG
  ! Call Poisson deviate
  ! By Adam Herrington
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       integer,                     intent(in)  :: nz, nup, ksfc, ktop, kdir
       real(kind_phys), dimension(4),      intent(in)  :: state
       real(kind_phys), dimension(nz,nup), intent(in)  :: lambda
       real(kind_phys), intent(in) :: pi
       integer,  dimension(nz,nup), intent(out) :: poi
       integer,  dimension(1,4)                 :: tmpseed
       integer                                  :: i,j
       type(ShrKissRandGen)                     :: kiss_gen

       ! Compute seed
       tmpseed(1,1) = int((state(1) - int(state(1))) * 1000000000._kind_phys)
       tmpseed(1,2) = int((state(2) - int(state(2))) * 1000000000._kind_phys)
       tmpseed(1,3) = int((state(3) - int(state(3))) * 1000000000._kind_phys)
       tmpseed(1,4) = int((state(4) - int(state(4))) * 1000000000._kind_phys)

       ! Set seed
       kiss_gen = ShrKissRandGen(tmpseed)

       ! Loop from the SURFACE to the TOP to preserve the exact PRNG sequence
       do i = ksfc, ktop, kdir
         do j = 1, nup
           call hybridRNG(kiss_gen, lambda(i,j), poi(i,j), pi)
         enddo
       enddo

  end subroutine poisson


  subroutine hybridRNG(kiss_gen,lambda,kout,pi)
  !**********************************************************************
  ! Interface for the two poisson rng subroutines
  ! chooses the appropriate subroutine based on the value of lambda
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       type(ShrKissRandGen), intent(inout) :: kiss_gen
       real(kind_phys),             intent(in)    :: lambda
       integer,              intent(out)   :: kout
       real(kind_phys), intent(in) :: pi

       if (lambda < 10._kind_phys) then
          call knuth(kiss_gen,lambda,kout)
       else
          call hormann(kiss_gen,lambda,kout,pi)
       end if

  end subroutine hybridRNG


  subroutine knuth(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Discrete random poisson from Knuth
  ! The Art of Computer Programming, v2, 137-138
  ! By Adam Herrington
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       type(ShrKissRandGen), intent(inout) :: kiss_gen
       real(kind_phys),             intent(in)    :: lambda
       integer,              intent(out)   :: kout

       ! Local variables
       real(kind_phys), dimension(1,1) :: tmpuni
       real(kind_phys)                 :: puni, explam
       integer                  :: k

       k = 0
       explam = exp(-1._kind_phys*lambda)
       puni = 1._kind_phys
       do while (puni > explam)
         k = k + 1
         call kiss_gen%random(tmpuni)
         puni = puni*tmpuni(1,1)
       end do
       kout = k - 1

  end subroutine knuth


  subroutine hormann(kiss_gen,lambda,kout,pi)
  !**********************************************************************
  ! Discrete random poisson
  ! Implements Poisson Transformed Rejection with Squeeze (PTRS)
  ! from W. Hormann Insurance: Mathematics and Economics 12, 39-45 (1993)
  ! By Jake Reschke
  ! DOI: 10.1016/0167-6687(93)90997-4
  !**********************************************************************
  use shr_RandNum_mod, only: ShrKissRandGen

      type(ShrKissRandGen), intent(inout) :: kiss_gen
      real(kind_phys),             intent(in)    :: lambda
      integer,              intent(out)   :: kout
      real(kind_phys), intent(in) :: pi

      ! Local variables
      real(kind_phys), dimension(1,1) :: U,V
      real(kind_phys)                 :: a,b,vr,alphinv,us,loggam
      integer                  :: k,i

      b = 0.931_kind_phys + 2.53_kind_phys*sqrt(lambda)
      a = -0.059_kind_phys + 0.02483_kind_phys*b
      vr = 0.9277_kind_phys - 3.6224_kind_phys/(b - 2._kind_phys)
      alphinv = 1.1239_kind_phys + 1.1328_kind_phys/(b - 3.4_kind_phys)

      do
         call kiss_gen%random(U)
         call kiss_gen%random(V)
         U(1,1) = U(1,1) - 0.5_kind_phys
         us = 0.5_kind_phys - abs(U(1,1))
         k = floor( (2._kind_phys*a/us + b)*U(1,1) + lambda + 0.43_kind_phys )
         if (us >= 0.07_kind_phys .and.  V(1,1) <= vr) then
            kout = k
            exit
         end if
         if (k <= 0 .or. (us < 0.013_kind_phys .and. V(1,1) > us)) then
            cycle
         end if
         ! compute log(k!). If k >=10 use stirling's approximation
         if (k < 10) then
            loggam = 0._kind_phys
            do i = 1, k
               loggam = loggam + log(1._kind_phys*i)
            end do
         else
            loggam = log(sqrt(2._kind_phys*pi)) + (k + 0.5_kind_phys)*log(1._kind_phys*k) - k + (1._kind_phys/12._kind_phys - 1._kind_phys/(360._kind_phys*k*k))/k
         end if
         if (log( V(1,1)*alphinv/(a/(us*us) + b) ) <= -1._kind_phys*lambda + k*log(lambda) - loggam) then
            kout = k
            exit
         end if
      end do

  end subroutine hormann

  subroutine roots(a,b,c,r1,r2,status)
  ! --------------------------------------------------------- !
  ! Subroutine to solve the second order polynomial equation. !
  ! after uwshcu.F90                                          !
  ! --------------------------------------------------------- !
    real(kind_phys), intent(in)  :: a
    real(kind_phys), intent(in)  :: b
    real(kind_phys), intent(in)  :: c
    real(kind_phys), intent(out) :: r1
    real(kind_phys), intent(out) :: r2
    integer , intent(out) :: status
    real(kind_phys)              :: q

    status = 0

    if( a == 0._kind_phys ) then                              ! Form b*x + c = 0
        if( b == 0._kind_phys ) then                          ! Failure: c = 0
            status = 1
        else                                           ! b*x + c = 0
            r1 = -c/b
        endif
        r2 = r1
    else
        if( b == 0._kind_phys ) then                          ! Form a*x**2 + c = 0
            if( a*c > 0._kind_phys ) then                     ! Failure: x**2 = -c/a < 0
                status = 2
            else                                       ! x**2 = -c/a
                r1 = sqrt(-c/a)
            endif
            r2 = -r1
       else                                            ! Form a*x**2 + b*x + c = 0
            if( (b**2 - 4._kind_phys*a*c) < 0._kind_phys ) then      ! Failure, no real roots
                 status = 3
            else
                 q  = -0.5_kind_phys*(b + sign(1.0_kind_phys,b)*sqrt(b**2 - 4._kind_phys*a*c))
                 r1 =  q/a
                 r2 =  c/q
            endif
       endif
    endif

    return

  end subroutine roots


  subroutine oneplume( cpair, gravit, zvir, h2otrip, latvap, latice, &
                       nzm, nzt, zm, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, &
                       wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, &
                       wa, wb, tke, do_precip, plumeheight )
  !**********************************************************************
  ! Calculate a single plume with fixed entrainment
  ! to be used for a dynamic mixing length calculation
  ! By Rachel Storer
  !**********************************************************************

    real(kind_phys), intent(in)   :: cpair, gravit, zvir, h2otrip, latvap, latice
    integer,  intent(in)                 :: nzm, nzt
    real(kind_phys), dimension(nzt), intent(in) :: dzt, iexner_zt,  &
                                            qt, thv, thl
    real(kind_phys), dimension(nzm), intent(in) :: zm, p_zm, iexner_zm, tke

    real(kind_phys), intent(in)                :: wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, &
                                           cwthv, zcb_unset, wa, wb
    logical, intent(in)                 :: do_precip

    real(kind_phys), intent(inout)             :: plumeheight

    !local variables

    ! =============================================================================== !
    ! GRID ORIENTATION GENERALIZATION VARIABLES
    ! ------------------------------------------------------------------------------- !
    ! To support both top-down (CAM) and bottom-up (CLUBB) grid orientations without
    ! duplicating code, these variables abstract the vertical loop bounds and slices.
    !
    ! ksfcm / ksfct : Index of the surface for momentum (m) and thermodynamic (t) grids.
    ! ktopm / ktopt : Index of the model top for momentum (m) and thermodynamic (t) grids.
    ! kdir          : Directional step (+1 for moving up, -1 for moving down).
    !
    ! STAGGERED GRID INDEXING
    ! Because momentum (zm) and thermodynamic (zt) grids are staggered, the relative
    ! index of the cell center (zt) to the interface (zm) flips depending on whether
    ! memory is loaded top-down or bottom-up. These variables dynamically map them:
    !
    ! kt    : The active thermodynamic cell center associated with the current step.
    !         Upward Sweep:   kt = k - (1-kdir)/2
    !         Downward Sweep: kt = k - (1+kdir)/2
    !
    ! kn    : The NEXT momentum interface in the direction of the current sweep.
    !         Upward Sweep:   kn = k + kdir
    !         Downward Sweep: kn = k - kdir
    !
    ! kt_up : The thermodynamic cell center physically ABOVE momentum interface k.
    !         kt_up = k - (1-kdir)/2
    !
    ! kt_dn : The thermodynamic cell center physically BELOW momentum interface k.
    !         kt_dn = k - (1+kdir)/2
    ! =============================================================================== !
    integer :: ksfcm, ktopm, kdir, kt, kn
    integer                     :: k
    real(kind_phys)                    :: thvn, qtn, thln, qcn, thn, qln, qin, qsn, lmixn, zcb, B, wn2, pentexp, pturb, pentw, wp
    real(kind_phys), dimension(nzm)     :: upw, upa, upqt, upthv, upthl, upth, upqs, &
                                   upqc, upql, upqi
    real(kind_phys), dimension(nzt)     :: supqt, supthl
    ! ! fractional entrainment rate
    real(kind_phys), parameter         :: pent = 1.E-3_kind_phys
    ! ! use tke enhanced entrainment
    logical                     :: do_tptke = .false.

    if (zm(1) < zm(nzm)) then
       ksfcm = 1
       ktopm = nzm
       kdir = 1
    else
       ksfcm = nzm
       ktopm = 1
       kdir = -1
    end if

    zcb = zcb_unset

    upw(ksfcm) = 0.5_kind_phys * wmax
    upa(ksfcm) = 0.5_kind_phys * erf( wmax/(sqrt(2._kind_phys)*sigmaw) )

    upqt(ksfcm)  = cwqt * upw(ksfcm) * sigmaqt/sigmaw
    upthv(ksfcm) = cwthv * upw(ksfcm) * sigmathv/sigmaw

    upqt(ksfcm) = qt(ksfcm)+upqt(ksfcm)
    upthv(ksfcm) = thv(ksfcm)+upthv(ksfcm)
    upthl(ksfcm) = upthv(ksfcm) / (1._kind_phys+zvir*upqt(ksfcm))
    upth(ksfcm)  = upthl(ksfcm)

    ! get cloud, lowest momentum level
    call condensation_mf(cpair, zvir, h2otrip, latvap, latice, &
                         upqt(ksfcm), upthl(ksfcm), p_zm(ksfcm), iexner_zm(ksfcm), &
                         thvn, qcn, thn, qln, qin, qsn, lmixn)
    upthv(ksfcm) = thvn
    upqc(ksfcm)  = qcn
    upql(ksfcm)  = qln
    upqi(ksfcm)  = qin
    upqs(ksfcm)  = qsn
    upth(ksfcm)  = thn
    if (qcn > 0._kind_phys) zcb = zm(ksfcm)

    do k = ksfcm, ktopm-kdir, kdir
      kt = k - (1-kdir)/2
      kn = k + kdir

      ! get microphysics, autoconversion
      if (do_precip .and. upqc(k) > 0._kind_phys) then
        call precip_mf(upqs(k),upqt(k),upw(k),dzt(kt),abs(zm(kn)-zcb),supqt(kt))
        supthl(kt) = -1._kind_phys*lmixn*supqt(kt)*iexner_zt(kt)/cpair
      else
        supqt(kt)  = 0._kind_phys
        supthl(kt) = 0._kind_phys
      end if
      ! integrate updraft
      if (do_tptke) then
        pturb = (1._kind_phys + clubb_mf_alphturb*sqrt(tke(k))/upw(k))
      else
        pturb = 1._kind_phys
      end if
      pentexp  = exp(-pent*pturb*dzt(kt))
      qtn  = qt(kt) *(1._kind_phys-pentexp ) + upqt (k)*pentexp + supqt(kt)
      thln = thl(kt)*(1._kind_phys-pentexp ) + upthl(k)*pentexp + supthl(kt)

      ! convert source terms to a tendency
      supqt(kt) = supqt(kt)*upw(k)/dzt(kt)
      supthl(kt) = supthl(kt)*upw(k)/dzt(kt)

      ! get cloud, momentum levels
      call condensation_mf(cpair, zvir, h2otrip, latvap, latice, &
                           qtn, thln, p_zm(kn), iexner_zm(kn), &
                           thvn, qcn, thn, qln, qin, qsn, lmixn)
      if (zcb == zcb_unset .and. qcn > 0._kind_phys) zcb = zm(kn)
      ! get buoyancy
      B=gravit*(0.5_kind_phys*(thvn + upthv(k))/thv(kt)-1._kind_phys)

      ! get wn^2
      wp = wb*pent*pturb
      if (wp==0._kind_phys) then
         wn2 = upw(k)**2._kind_phys+2._kind_phys*wa*B*dzt(kt)
      else
         pentw = exp(-2._kind_phys*wp*dzt(kt))
         wn2 = pentw*upw(k)**2._kind_phys+(1._kind_phys-pentw)*wa*B/wp
      end if

      if (wn2>0._kind_phys) then
        upw(kn)   = sqrt(wn2)
        upthv(kn) = thvn
        upthl(kn) = thln
        upqt(kn)  = qtn
        upqc(kn)  = qcn
        upqs(kn)  = qsn
        upa(kn)   = upa(k)
        upql(kn)  = qln
        upqi(kn)  = qin
        upth(kn)  = thn
      else
        plumeheight = zm(k)
        exit
      end if
    enddo

  end subroutine oneplume

end module clubb_mf
