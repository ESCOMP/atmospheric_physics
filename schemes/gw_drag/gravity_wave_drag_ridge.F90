! This module handles gravity waves from orographic sources, and was
! extracted from gw_drag in May 2013,
! and subsequently CCPPized in August 2025.
module gravity_wave_drag_ridge
  use ccpp_kinds,    only: kind_phys
  use gw_common,     only: unset_kind_phys, gwband
  use coords_1d,     only: Coords1D

  implicit none
  private
  save

  ! Public CCPP-compliant interfaces.
  public :: gravity_wave_drag_ridge_init        ! CCPP I/O read Ridge data into state and initialize.

  public :: gravity_wave_drag_ridge_beta_run    ! Meso-Beta.
  public :: gravity_wave_drag_ridge_gamma_run   ! Meso-Gamma.

  ! For CAM compatibility only.
  public :: gw_rdg_init                         ! Init routine (public only for current CAM compat.)

  ! use separate dividing streamlines for downslope wind and flow splitting regimes ("DS" configuration)?
  ! or use single dividing streamline as in Scinocca and McFarlane 2000 ("SM" configuration).
  logical         :: do_divstream

  !===========================================
  ! Parameters for DS2017 (do_divstream=.T.)
  !===========================================
  ! Amplification factor - 1.0 for
  ! high-drag/windstorm regime
  real(kind_phys) :: C_BetaMax_DS

  ! Max Ratio Fr2:Fr1 - 1.0
  real(kind_phys) :: C_GammaMax

  ! Normalized limits for Fr2(Frx) function
  real(kind_phys) :: Frx0
  real(kind_phys) :: Frx1

  !===========================================
  ! Parameters for SM2000
  ! Scinocca, J.F. and McFarlane, N.A. (2000),
  ! The parametrization of drag induced by stratified flow over anisotropic orography.
  ! Q.J.R. Meteorol. Soc., 126: 2353-2393. https://doi.org/10.1002/qj.49712656802
  !===========================================
  ! Amplification factor - 1.0 for
  ! high-drag/windstorm regime
  real(kind_phys) :: C_BetaMax_SM

  ! NOTE: Critical inverse Froude number Fr_c is
  ! 1./(SQRT(2.)~0.707 in SM2000 (should be <= 1), pp. 2388.
  real(kind_phys) :: Fr_c

  logical :: gw_rdg_do_vdiff = .true.

  ! Limiters (min/max values)
  ! min surface displacement height for orographic waves (m)
  real(kind_phys) :: orohmin
  ! min wind speed for orographic waves
  real(kind_phys) :: orovmin
  ! min stratification allowing wave behavior
  real(kind_phys) :: orostratmin
  ! min stratification allowing wave behavior
  real(kind_phys) :: orom2min

  logical            :: do_smooth_regimes
  logical            :: do_adjust_tauoro
  logical            :: do_backward_compat

  ! GWBand for orographic wave drag.
  type(GWBand)       :: band_oro

  ! anisotropic ridge fields
  integer, parameter :: prdg = 16

  logical            :: use_gw_rdg_gamma
  logical            :: use_gw_rdg_beta

contains
  subroutine gw_rdg_init(&
                         gw_delta_c, &
                         effgw_rdg_beta, &
                         effgw_rdg_gamma, &
                         use_gw_rdg_beta_in, &
                         use_gw_rdg_gamma_in, &
                         gw_rdg_do_divstream_nl, gw_rdg_C_BetaMax_DS_nl, gw_rdg_C_GammaMax_nl, &
                         gw_rdg_Frx0_nl, gw_rdg_Frx1_nl, gw_rdg_C_BetaMax_SM_nl, gw_rdg_Fr_c_nl, &
                         gw_rdg_do_smooth_regimes_nl, gw_rdg_do_adjust_tauoro_nl, &
                         gw_rdg_do_backward_compat_nl, gw_rdg_orohmin_nl, gw_rdg_orovmin_nl, &
                         gw_rdg_orostratmin_nl, gw_rdg_orom2min_nl, gw_rdg_do_vdiff_nl, &
                         errmsg, errflg)

    use gw_common, only: wavelength_mid

    ! Gravity wave band parameters
    real(kind_phys), intent(in)      :: gw_delta_c                    ! Width of speed bins (delta c) for gravity wave spectrum [m s-1]

    ! Ridge efficiency parameters
    real(kind_phys), intent(in)      :: effgw_rdg_beta                ! Beta ridge efficiency factor [1]
    real(kind_phys), intent(in)      :: effgw_rdg_gamma               ! Gamma ridge efficiency factor [1]

    ! Ridge scheme control flags
    logical, intent(in)              :: use_gw_rdg_beta_in            ! Enable beta ridge scheme [flag]
    logical, intent(in)              :: use_gw_rdg_gamma_in           ! Enable gamma ridge scheme [flag]

    ! Dividing streamline (DS2017) parameters
    logical, intent(in)              :: gw_rdg_do_divstream_nl        ! Enable dividing streamline parameterization [flag]
    real(kind_phys), intent(in)      :: gw_rdg_C_BetaMax_DS_nl        ! Enhancement factor for downslope wind stress in DS [1]
    real(kind_phys), intent(in)      :: gw_rdg_C_GammaMax_nl          ! Enhancement factor for depth of downslope wind regime in DS configuration [1]
    real(kind_phys), intent(in)      :: gw_rdg_Frx0_nl                ! Lower inverse Froude number limits on linear ramp terminating downslope wind regime for high mountains in DS configuration [1]
    real(kind_phys), intent(in)      :: gw_rdg_Frx1_nl                ! Upper inverse Froude number limits on linear ramp terminating downslope wind regime for high mountains in DS configuration [1]

    ! Scinocca & McFarlane (SM2000) parameters
    real(kind_phys), intent(in)      :: gw_rdg_C_BetaMax_SM_nl        ! Enhancement factor for downslope wind stress in SM configuration. [1]
    real(kind_phys), intent(in)      :: gw_rdg_Fr_c_nl                ! Critical inverse Froude number [1]

    ! Ridge scheme behavior flags
    logical, intent(in)              :: gw_rdg_do_smooth_regimes_nl   ! Enable smooth regime transitions [flag]
    logical, intent(in)              :: gw_rdg_do_adjust_tauoro_nl    ! Enable orographic stress adjustment [flag]
    logical, intent(in)              :: gw_rdg_do_backward_compat_nl  ! Adjust for bit-for-bit answers with the ("N5") configuration [flag]

    ! Ridge scheme physical limits
    real(kind_phys), intent(in)      :: gw_rdg_orohmin_nl             ! Minimum surface displacement height for orographic waves [m]
    real(kind_phys), intent(in)      :: gw_rdg_orovmin_nl             ! Minimum wind speed for orographic waves [m s-1]
    real(kind_phys), intent(in)      :: gw_rdg_orostratmin_nl         ! Minimum stratification allowing wave behavior [s-1]
    real(kind_phys), intent(in)      :: gw_rdg_orom2min_nl            ! Minimum normalized vertical wavenumber squared [1]
    logical, intent(in)              :: gw_rdg_do_vdiff_nl            ! Ridge scheme contribute to vdiff tendencies? [flag]

    character(len=512), intent(out)  :: errmsg
    integer, intent(out)             :: errflg

    character(len=*), parameter :: sub = 'gw_rdg_init'

    ! Initialize gravity wave band based on wavelength
    band_oro = GWBand(0, gw_delta_c, 1.0_kind_phys, wavelength_mid)

    ! Set the local variables from namelist read
    do_divstream = gw_rdg_do_divstream_nl
    C_BetaMax_DS = gw_rdg_C_BetaMax_DS_nl
    C_GammaMax = gw_rdg_C_GammaMax_nl
    Frx0 = gw_rdg_Frx0_nl
    Frx1 = gw_rdg_Frx1_nl
    C_BetaMax_SM = gw_rdg_C_BetaMax_SM_nl
    Fr_c = gw_rdg_Fr_c_nl
    do_smooth_regimes = gw_rdg_do_smooth_regimes_nl
    do_adjust_tauoro = gw_rdg_do_adjust_tauoro_nl
    do_backward_compat = gw_rdg_do_backward_compat_nl
    orohmin = gw_rdg_orohmin_nl
    orovmin = gw_rdg_orovmin_nl
    orostratmin = gw_rdg_orostratmin_nl
    orom2min = gw_rdg_orom2min_nl
    gw_rdg_do_vdiff = gw_rdg_do_vdiff_nl
    use_gw_rdg_beta = use_gw_rdg_beta_in
    use_gw_rdg_gamma = use_gw_rdg_gamma_in

    if (use_gw_rdg_beta) then
      if (effgw_rdg_beta == unset_kind_phys) then
        errmsg = sub//": ERROR: Anisotropic OGW enabled, but effgw_rdg_beta was not set."
        errflg = 1
        return
      end if
    end if

    if (use_gw_rdg_gamma) then
      if (effgw_rdg_gamma == unset_kind_phys) then
        errmsg = sub//": ERROR: Anisotropic OGW enabled, but effgw_rdg_gamma was not set."
        errflg = 1
        return
      end if
    end if
  end subroutine gw_rdg_init

  subroutine gravity_wave_drag_ridge_init( &
    ncol, &
    amIRoot, iulog, &
    rearth, &
    use_gw_rdg_beta_in, &
    bnd_topo_file, &
    use_gw_rdg_gamma_in, &
    bnd_rdggm_file, &
    gbxar, isovar, isowgt, &
    hwdth, clngt, mxdis, anixy, angll, &
    gbxarg, &
    hwdthg, clngtg, mxdisg, anixyg, angllg, &
    gw_delta_c, &
    effgw_rdg_beta, &
    effgw_rdg_gamma, &
    gw_rdg_do_divstream_nl, gw_rdg_C_BetaMax_DS_nl, gw_rdg_C_GammaMax_nl, &
    gw_rdg_Frx0_nl, gw_rdg_Frx1_nl, gw_rdg_C_BetaMax_SM_nl, gw_rdg_Fr_c_nl, &
    gw_rdg_do_smooth_regimes_nl, gw_rdg_do_adjust_tauoro_nl, &
    gw_rdg_do_backward_compat_nl, gw_rdg_orohmin_nl, gw_rdg_orovmin_nl, &
    gw_rdg_orostratmin_nl, gw_rdg_orom2min_nl, gw_rdg_do_vdiff_nl, &
    errmsg, errflg)

    use ccpp_io_reader, only: abstract_netcdf_reader_t, create_netcdf_reader_t

    ! Input arguments
    integer, intent(in)              :: ncol
    !integer, intent(in)              :: prdg                          ! second dimension of the 2-D outs below.
    logical, intent(in)              :: amIRoot
    integer, intent(in)              :: iulog
    real(kind_phys), intent(in)      :: rearth                        ! Earth radius [m]

    logical, intent(in)              :: use_gw_rdg_beta_in            ! Enable Meso-beta ridges [flag]
    character(len=256), intent(in)   :: bnd_topo_file                 ! Filepath of topo file
    logical, intent(in)              :: use_gw_rdg_gamma_in           ! Enable Meso-gamma ridges [flag]
    character(len=256), intent(in)   :: bnd_rdggm_file                ! Filepath of ridge (gamma) file

    ! Gravity wave band parameters
    real(kind_phys), intent(in)      :: gw_delta_c                    ! Width of speed bins (delta c) for gravity wave spectrum [m s-1]

    ! Ridge efficiency parameters
    real(kind_phys), intent(in)      :: effgw_rdg_beta                ! Beta ridge efficiency factor [1]
    real(kind_phys), intent(in)      :: effgw_rdg_gamma               ! Gamma ridge efficiency factor [1]

    ! Dividing streamline (DS2017) parameters
    logical, intent(in)              :: gw_rdg_do_divstream_nl        ! Enable dividing streamline parameterization [flag]
    real(kind_phys), intent(in)      :: gw_rdg_C_BetaMax_DS_nl        ! Enhancement factor for downslope wind stress in DS [1]
    real(kind_phys), intent(in)      :: gw_rdg_C_GammaMax_nl          ! Enhancement factor for depth of downslope wind regime in DS configuration [1]
    real(kind_phys), intent(in)      :: gw_rdg_Frx0_nl                ! Lower inverse Froude number limits on linear ramp terminating downslope wind regime for high mountains in DS configuration [1]
    real(kind_phys), intent(in)      :: gw_rdg_Frx1_nl                ! Upper inverse Froude number limits on linear ramp terminating downslope wind regime for high mountains in DS configuration [1]

    ! Scinocca & McFarlane (SM2000) parameters
    real(kind_phys), intent(in)      :: gw_rdg_C_BetaMax_SM_nl        ! Enhancement factor for downslope wind stress in SM configuration. [1]
    real(kind_phys), intent(in)      :: gw_rdg_Fr_c_nl                ! Critical inverse Froude number [1]

    ! Ridge scheme behavior flags
    logical, intent(in)              :: gw_rdg_do_smooth_regimes_nl   ! Enable smooth regime transitions [flag]
    logical, intent(in)              :: gw_rdg_do_adjust_tauoro_nl    ! Enable orographic stress adjustment [flag]
    logical, intent(in)              :: gw_rdg_do_backward_compat_nl  ! Adjust for bit-for-bit answers with the ("N5") configuration [flag]

    ! Ridge scheme physical limits
    real(kind_phys), intent(in)      :: gw_rdg_orohmin_nl             ! Minimum surface displacement height for orographic waves [m]
    real(kind_phys), intent(in)      :: gw_rdg_orovmin_nl             ! Minimum wind speed for orographic waves [m s-1]
    real(kind_phys), intent(in)      :: gw_rdg_orostratmin_nl         ! Minimum stratification allowing wave behavior [s-1]
    real(kind_phys), intent(in)      :: gw_rdg_orom2min_nl            ! Minimum normalized vertical wavenumber squared [1]
    logical, intent(in)              :: gw_rdg_do_vdiff_nl            ! Ridge scheme contribute to vdiff tendencies? [flag]

    ! Output arguments
    real(kind_phys),    intent(out), pointer :: gbxar (:)
    real(kind_phys),    intent(out), pointer :: isovar(:)
    real(kind_phys),    intent(out), pointer :: isowgt(:)
    real(kind_phys),    intent(out), pointer :: hwdth (:,:)
    real(kind_phys),    intent(out), pointer :: clngt (:,:)
    real(kind_phys),    intent(out), pointer :: mxdis (:,:)
    real(kind_phys),    intent(out), pointer :: anixy (:,:)
    real(kind_phys),    intent(out), pointer :: angll (:,:)

    real(kind_phys),    intent(out), pointer :: gbxarg (:)
    real(kind_phys),    intent(out), pointer :: hwdthg (:,:)
    real(kind_phys),    intent(out), pointer :: clngtg (:,:)
    real(kind_phys),    intent(out), pointer :: mxdisg (:,:)
    real(kind_phys),    intent(out), pointer :: anixyg (:,:)
    real(kind_phys),    intent(out), pointer :: angllg (:,:)

    character(len=512), intent(out)  :: errmsg
    integer, intent(out)             :: errflg

    ! Local variables
    logical :: has_gbxar_from_topo

    ! Temporaries for providing to the I/O reader; data will be copied to the
    ! respective pointer variables.
    real(kind_phys), allocatable :: alloc1D(:)    ! 1-D temporary for I/O reader
    real(kind_phys), allocatable :: alloc2D(:,:)  ! 2-D temporary for I/O reader

    class(abstract_netcdf_reader_t), allocatable :: reader

    errmsg = ''
    errflg = 0

    has_gbxar_from_topo = .false.

    reader = create_netcdf_reader_t()

    if(use_gw_rdg_beta) then
      call reader%open_file(bnd_topo_file, errmsg, errflg)
      if (errflg /= 0) return

      call reader%get_var('GBXAR', alloc1D, errmsg, errflg)
      if (errflg /= 0) return
      gbxar(:) = alloc1D(:)*(rearth/1000._kind_phys)*(rearth/1000._kind_phys) ! transform to km^2
      deallocate(alloc1D, stat=errflg)
      has_gbxar_from_topo = .true.

      call reader%get_var('ISOVAR', alloc1D, errmsg, errflg)
      if (errflg /= 0) then
        ! topo files do not currently contain ISOVAR
        isovar(:) = 0._kind_phys
        errflg = 0
        errmsg = ''
      else
        isovar(:) = alloc1D(:)
        deallocate(alloc1D, stat=errflg)
      endif

      call reader%get_var('ISOWGT', alloc1D, errmsg, errflg)
      if (errflg /= 0) then
        ! topo files do not currently contain ISOWGT
        isowgt(:) = 0._kind_phys
        errflg = 0
        errmsg = ''
      else
        isowgt(:) = alloc1D(:)
        deallocate(alloc1D, stat=errflg)
      endif

      call reader%get_var('HWDTH', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      hwdth(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      call reader%get_var('CLNGT', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      clngt(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      call reader%get_var('MXDIS', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      mxdis(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      call reader%get_var('ANIXY', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      anixy(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      call reader%get_var('ANGLL', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      angll(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      ! close topo file only if it was opened here
      if (len_trim(bnd_topo_file) > 0) then
        call reader%close_file(errmsg, errflg)
        if (errflg /= 0) then
          return
        end if

        if (amIRoot) then
          write(iulog, *) "gravity_wave_drag_ridge_init: Read in ridge beta source file (topo)."
        end if
      end if
    end if

    if(use_gw_rdg_gamma) then
      call reader%open_file(bnd_rdggm_file, errmsg, errflg)
      if (errflg /= 0) return

      ! Note (hplin 8/6/25): the behavior, based on the original gw_drag, for GBXAR, ISOVAR, ISOWGT
      ! are very strange.
      !
      ! If GBXAR is available from topo file, then GBXAR from topo is used;
      ! otherwise, GBXAR is read from rdggm.
      !
      ! ISOVAR and ISOWGT are not populated at all from the rdggm file; the only files
      ! available as of writing were fv files that did not contain these variables.
      ! It appears that the gw_drag.F90 code simply uses dangling pointers for those,
      ! which we have replaced with actual isovar/isowgt from the "topo file" (actually still zeroes)
      ! as they are also not in the topo file, which I believe is closest to the original intent
      ! of the code.
      if(.not. has_gbxar_from_topo) then
        call reader%get_var('GBXAR', alloc1D, errmsg, errflg)
        if (errflg /= 0) return
        gbxarg(:) = alloc1D(:)*(rearth/1000._kind_phys)*(rearth/1000._kind_phys) ! transform to km^2
        deallocate(alloc1D, stat=errflg)
      endif

      ! call reader%get_var('ISOVAR', alloc1D, errmsg, errflg)
      ! if (errflg /= 0) return
      ! isovarg(:) = alloc1D(:)
      ! deallocate(alloc1D, stat=errflg)

      ! call reader%get_var('ISOWGT', alloc1D, errmsg, errflg)
      ! if (errflg /= 0) return
      ! isowgtg(:) = alloc1D(:)
      ! deallocate(alloc1D, stat=errflg)

      call reader%get_var('HWDTH', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      hwdthg(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      call reader%get_var('CLNGT', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      clngtg(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      call reader%get_var('MXDIS', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      mxdisg(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      call reader%get_var('ANIXY', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      anixyg(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      call reader%get_var('ANGLL', alloc2D, errmsg, errflg)
      if (errflg /= 0) return
      angllg(:,:) = alloc2D(:,:)
      deallocate(alloc2D, stat=errflg)

      ! close topo file only if it was opened here
      if (len_trim(bnd_rdggm_file) > 0) then
        call reader%close_file(errmsg, errflg)
        if (errflg /= 0) then
          return
        end if

        if (amIRoot) then
          write(iulog, *) "gravity_wave_drag_gamma_ridge_init: Read in ridge gamma source file."
        end if
      end if
    end if

    ! Call underlying initialization subroutine to populate namelist variables
    !REMOVECAM: once CAM is retired, gw_rdg_init can be collapsed into this subroutine.
    call gw_rdg_init( &
      gw_delta_c                   = gw_delta_c, &
      effgw_rdg_beta               = effgw_rdg_beta, &
      effgw_rdg_gamma              = effgw_rdg_gamma, &
      use_gw_rdg_beta_in           = use_gw_rdg_beta_in, &
      use_gw_rdg_gamma_in          = use_gw_rdg_gamma_in, &
      gw_rdg_do_divstream_nl       = gw_rdg_do_divstream_nl, &
      gw_rdg_C_BetaMax_DS_nl       = gw_rdg_C_BetaMax_DS_nl, &
      gw_rdg_C_GammaMax_nl         = gw_rdg_C_GammaMax_nl, &
      gw_rdg_Frx0_nl               = gw_rdg_Frx0_nl, &
      gw_rdg_Frx1_nl               = gw_rdg_Frx1_nl, &
      gw_rdg_C_BetaMax_SM_nl       = gw_rdg_C_BetaMax_SM_nl, &
      gw_rdg_Fr_c_nl               = gw_rdg_Fr_c_nl, &
      gw_rdg_do_smooth_regimes_nl  = gw_rdg_do_smooth_regimes_nl, &
      gw_rdg_do_adjust_tauoro_nl   = gw_rdg_do_adjust_tauoro_nl, &
      gw_rdg_do_backward_compat_nl = gw_rdg_do_backward_compat_nl, &
      gw_rdg_orohmin_nl            = gw_rdg_orohmin_nl, &
      gw_rdg_orovmin_nl            = gw_rdg_orovmin_nl, &
      gw_rdg_orostratmin_nl        = gw_rdg_orostratmin_nl, &
      gw_rdg_orom2min_nl           = gw_rdg_orom2min_nl, &
      gw_rdg_do_vdiff_nl           = gw_rdg_do_vdiff_nl, &
      errmsg = errmsg, &
      errflg = errflg)
    if (errflg /= 0) return

  end subroutine gravity_wave_drag_ridge_init

  subroutine gravity_wave_drag_ridge_beta_run( &
    ncol, pver, pcnst, dt, pi, cpair, rair, &
    vramp, p, &
    n_rdg_beta, &
    u, v, t, q, dse, &
    piln, zm, zi, &
    nm, ni, rhoi, kvt_gw, &
    use_gw_rdg_resid, effgw_rdg_resid, &
    effgw_rdg_beta, effgw_rdg_beta_max, &
    rdg_beta_cd_llb, trpd_leewv_rdg_beta, &
    gbxar, isovar, isowgt, &
    hwdth, clngt, mxdis, anixy, angll, &
    taurx, taury, &
    tauardgx, tauardgy, &
    utrdg, vtrdg, ttrdg, &
    q_tend, s_tend, u_tend, v_tend, flx_heat, errmsg, errflg)

    use coords_1d, only: Coords1D

    ! Input arguments
    integer,             intent(in)    :: pver
    integer,             intent(in)    :: ncol
    integer,             intent(in)    :: pcnst
    real(kind_phys),     intent(in)    :: dt                           ! Physics timestep [s]
    real(kind_phys),     intent(in)    :: pi                           ! pi_constant [1]
    real(kind_phys),     intent(in)    :: cpair                        ! specific_heat_of_dry_air_at_constant_pressure [J kg-1 K-1]
    real(kind_phys),     intent(in)    :: rair                         ! gas_constant_of_dry_air [J kg-1 K-1]

    real(kind_phys),     intent(in)    :: vramp(:)                     ! Vertical tapering function [1]
    type(coords1d),      intent(in)    :: p                            ! Pressure coordinates, Coords1D

    integer,             intent(in)    :: n_rdg_beta                   ! Number of meso-Beta ridges (per gridbox) to invoke [count]

    real(kind_phys),     intent(in)    :: u(:, :)                      ! Midpoint zonal winds [m s-1]
    real(kind_phys),     intent(in)    :: v(:, :)                      ! Midpoint meridional winds [m s-1]
    real(kind_phys),     intent(in)    :: t(:, :)                      ! Midpoint temperatures [K]
    real(kind_phys),     intent(in)    :: q(:, :, :)                   ! Constituent array [kg kg-1]
    real(kind_phys),     intent(in)    :: dse(:, :)                    ! Dry static energy [J kg-1]

    real(kind_phys),     intent(in)    :: piln(:, :)                   ! Log of interface pressures [ln(Pa)]
    real(kind_phys),     intent(in)    :: zm(:, :)                     ! Midpoint altitudes above ground [m]
    real(kind_phys),     intent(in)    :: zi(:, :)                     ! Interface altitudes above ground [m]

    real(kind_phys),     intent(in)    :: nm(:, :)                     ! Midpoint Brunt-Vaisalla frequencies [s-1]
    real(kind_phys),     intent(in)    :: ni(:, :)                     ! Interface Brunt-Vaisalla frequencies [s-1]
    real(kind_phys),     intent(in)    :: rhoi(:, :)                   ! Interface density [kg m-3]
    real(kind_phys),     intent(in)    :: kvt_gw(:, :)                 ! Molecular thermal diffusivity [m2 s-1]

    ! Ridge scheme parameters
    logical,             intent(in)    :: use_gw_rdg_resid             ! Enable residual (non-ridge) orography GW. [flag]
    real(kind_phys),     intent(in)    :: effgw_rdg_resid              ! Efficiency scaling factor associated with residual non-ridge topo. [1]
    real(kind_phys),     intent(in)    :: effgw_rdg_beta               ! Efficiency scaling factor associated with anisotropic OGW. [1]
    real(kind_phys),     intent(in)    :: effgw_rdg_beta_max           ! Max efficiency associated with anisotropic OGW. [1]
    real(kind_phys),     intent(in)    :: rdg_beta_cd_llb              ! Drag coefficient for obstacles in low-level flow. [1]
    logical,             intent(in)    :: trpd_leewv_rdg_beta          ! Allow trapping for meso-Beta Ridges? [flag]

    ! Beta ridge input data
    real(kind_phys),     intent(in), pointer :: gbxar (:)              ! Grid box area [km2]
    real(kind_phys),     intent(in), pointer :: isovar(:)              ! Isotropic variance [m]
    real(kind_phys),     intent(in), pointer :: isowgt(:)              ! Isotropic weight [1]
    real(kind_phys),     intent(in), pointer :: hwdth (:,:)            ! Ridge half-width [km]
    real(kind_phys),     intent(in), pointer :: clngt (:,:)            ! Ridge length [km]
    real(kind_phys),     intent(in), pointer :: mxdis (:,:)            ! Ridge/obstacle height [m]
    real(kind_phys),     intent(in), pointer :: anixy (:,:)            ! Ridge anisotropy [1]
    real(kind_phys),     intent(in), pointer :: angll (:,:)            ! Ridge clockwise angle wrt north-south [degrees]

    ! Output tendencies
    real(kind_phys),     intent(inout) :: q_tend(:, :, :)              ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys),     intent(inout) :: s_tend(:, :)                 ! Dry static energy tendency [J kg-1 s-1]
    real(kind_phys),     intent(inout) :: u_tend(:, :)                 ! Zonal wind tendency [m s-2]
    real(kind_phys),     intent(inout) :: v_tend(:, :)                 ! Meridional wind tendency [m s-2]
    real(kind_phys),     intent(out)   :: flx_heat(:)                  ! Surface heat flux for check energy [W m-2]

    ! Diagnostic outputs
    real(kind_phys),     intent(out)   :: taurx(:,:)                   ! wave stress in zonal direction, interfaces [N m-2]
    real(kind_phys),     intent(out)   :: taury(:,:)                   ! wave stress in meridional direction, interfaces [N m-2]
    real(kind_phys),     intent(out)   :: tauardgx(:,:)                ! ridge based momentum flux profile, zonal, interfaces [N m-2]
    real(kind_phys),     intent(out)   :: tauardgy(:,:)                ! ridge based momentum flux profile, meridional, interfaces [N m-2]
    real(kind_phys),     intent(out)   :: utrdg(:,:)                   ! U tendency accummulator [m s-1]
    real(kind_phys),     intent(out)   :: vtrdg(:,:)                   ! V tendency accummulator [m s-1]
    real(kind_phys),     intent(out)   :: ttrdg(:,:)                   ! T tendency from orographic gravity wave drag [K s-1]

    ! Error handling
    character(len=512),  intent(out)   :: errmsg
    integer,             intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    call gw_rdg_calc(band_oro          = band_oro, &
                     vramp             = vramp(:), &
                     type              = 'BETA ', &
                     pver              = pver, &
                     ncol              = ncol, &
                     pcnst             = pcnst, &
                     prdg              = prdg, &
                     n_rdg             = n_rdg_beta, &
                     dt                = dt, &
                     pi                = pi, &
                     cpair             = cpair, &
                     rair              = rair, &
                     u                 = u(:, :), &
                     v                 = v(:, :), &
                     t                 = t(:, :), &
                     p                 = p, &
                     piln              = piln(:, :), &
                     zm                = zm(:, :), &
                     zi                = zi(:, :), &
                     nm                = nm(:, :), &
                     ni                = ni(:, :), &
                     rhoi              = rhoi(:, :), &
                     kvtt              = kvt_gw(:, :), &
                     q                 = q(:, :, :), &
                     dse               = dse(:, :), &
                     effgw_rdg         = effgw_rdg_beta, &
                     effgw_rdg_max     = effgw_rdg_beta_max, &
                     effgw_rdg_resid   = effgw_rdg_resid, &
                     luse_gw_rdg_resid = use_gw_rdg_resid, &
                     hwdth             = hwdth(:, :), &
                     clngt             = clngt(:, :), &
                     gbxar             = gbxar(:), &
                     mxdis             = mxdis(:, :), &
                     angll             = angll(:, :), &
                     anixy             = anixy(:, :), &
                     isovar            = isovar(:), &
                     isowgt            = isowgt(:), &
                     rdg_cd_llb        = rdg_beta_cd_llb, &
                     trpd_leewv        = trpd_leewv_rdg_beta, &
                     u_tend            = u_tend(:, :), &
                     v_tend            = v_tend(:, :), &
                     s_tend            = s_tend(:, :), &
                     q_tend            = q_tend(:, :, :), &
                     flx_heat          = flx_heat(:), &
                     taurx             = taurx(:,:), &
                     taury             = taury(:,:), &
                     tauardgx          = tauardgx(:,:), &
                     tauardgy          = tauardgy(:,:), &
                     utrdg             = utrdg(:,:), &
                     vtrdg             = vtrdg(:,:), &
                     ttrdg             = ttrdg(:,:), &
                     errmsg            = errmsg, &
                     errflg            = errflg)

  end subroutine gravity_wave_drag_ridge_beta_run

  subroutine gravity_wave_drag_ridge_gamma_run( &
    ncol, pver, pcnst, dt, pi, cpair, rair, &
    vramp, p, &
    n_rdg_gamma, &
    u, v, t, q, dse, &
    piln, zm, zi, &
    nm, ni, rhoi, kvt_gw, &
    use_gw_rdg_resid, effgw_rdg_resid, &
    effgw_rdg_gamma, effgw_rdg_gamma_max, &
    rdg_gamma_cd_llb, trpd_leewv_rdg_gamma, &
    gbxarg, &
    hwdthg, clngtg, mxdisg, anixyg, angllg, &
    q_tend, s_tend, u_tend, v_tend, flx_heat, &
    taurx, taury, &
    tauardgx, tauardgy, &
    utrdg, vtrdg, ttrdg, &
    errmsg, errflg)

    use coords_1d, only: Coords1D

    ! Input arguments
    integer,             intent(in)    :: pver                         ! Number of vertical layers [count]
    integer,             intent(in)    :: ncol                         ! Number of columns [count]
    integer,             intent(in)    :: pcnst                        ! Number of constituents [count]
    real(kind_phys),     intent(in)    :: dt                           ! Physics timestep [s]
    real(kind_phys),     intent(in)    :: pi                           ! pi_constant [1]
    real(kind_phys),     intent(in)    :: cpair                        ! specific_heat_of_dry_air_at_constant_pressure [J kg-1 K-1]
    real(kind_phys),     intent(in)    :: rair                         ! gas_constant_of_dry_air [J kg-1 K-1]

    real(kind_phys),     intent(in)    :: vramp(:)                     ! Vertical tapering function [1]
    type(coords1d),      intent(in)    :: p                            ! Pressure coordinates, Coords1D

    integer,             intent(in)    :: n_rdg_gamma                  ! Number of meso-gamma ridges (per gridbox) to invoke.

    real(kind_phys),     intent(in)    :: u(:, :)                      ! Midpoint zonal winds [m s-1]
    real(kind_phys),     intent(in)    :: v(:, :)                      ! Midpoint meridional winds [m s-1]
    real(kind_phys),     intent(in)    :: t(:, :)                      ! Midpoint temperatures [K]
    real(kind_phys),     intent(in)    :: q(:, :, :)                   ! Constituent array [kg kg-1]
    real(kind_phys),     intent(in)    :: dse(:, :)                    ! Dry static energy [J kg-1]

    real(kind_phys),     intent(in)    :: piln(:, :)                   ! Log of interface pressures [ln(Pa)]
    real(kind_phys),     intent(in)    :: zm(:, :)                     ! Midpoint altitudes above ground [m]
    real(kind_phys),     intent(in)    :: zi(:, :)                     ! Interface altitudes above ground [m]

    real(kind_phys),     intent(in)    :: nm(:, :)                     ! Midpoint Brunt-Vaisalla frequencies [s-1]
    real(kind_phys),     intent(in)    :: ni(:, :)                     ! Interface Brunt-Vaisalla frequencies [s-1]
    real(kind_phys),     intent(in)    :: rhoi(:, :)                   ! Density at interfaces [kg m-3]
    real(kind_phys),     intent(in)    :: kvt_gw(:, :)                 ! Molecular thermal diffusivity [m2 s-1]

    ! Ridge scheme parameters
    logical,             intent(in)    :: use_gw_rdg_resid             ! Enable residual (non-ridge) orography GW. [flag]
    real(kind_phys),     intent(in)    :: effgw_rdg_resid              ! Efficiency scaling factor associated with residual non-ridge topo. [1]
    real(kind_phys),     intent(in)    :: effgw_rdg_gamma              ! Efficiency scaling factor associated iwth anisotropic OGW [1]
    real(kind_phys),     intent(in)    :: effgw_rdg_gamma_max          ! Max efficiency associated with anisotropic OGW [1]
    real(kind_phys),     intent(in)    :: rdg_gamma_cd_llb             ! Drag coefficient for obstacles in low-level flow [1]
    logical,             intent(in)    :: trpd_leewv_rdg_gamma         ! Allow trapping for meso-gamma Ridges? [flag]

    ! Gamma ridge input data
    real(kind_phys),     intent(in), pointer :: gbxarg (:)             ! Grid box area for gamma ridges [km2]
    real(kind_phys),     intent(in), pointer :: hwdthg (:,:)           ! Gamma ridge half-width [km]
    real(kind_phys),     intent(in), pointer :: clngtg (:,:)           ! Gamma ridge length [km]
    real(kind_phys),     intent(in), pointer :: mxdisg (:,:)           ! Gamma ridge/obstacle height [m]
    real(kind_phys),     intent(in), pointer :: anixyg (:,:)           ! Gamma ridge anisotropy [1]
    real(kind_phys),     intent(in), pointer :: angllg (:,:)           ! Gamma ridge clockwise angle wrt north-south [degrees]

    ! Output tendencies
    real(kind_phys),     intent(inout) :: q_tend(:, :, :)              ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys),     intent(inout) :: s_tend(:, :)                 ! Dry static energy tendency [J kg-1 s-1]
    real(kind_phys),     intent(inout) :: u_tend(:, :)                 ! Zonal wind tendency [m s-2]
    real(kind_phys),     intent(inout) :: v_tend(:, :)                 ! Meridional wind tendency [m s-2]
    real(kind_phys),     intent(out)   :: flx_heat(:)                  ! Surface heat flux for check energy [W m-2]

    ! Diagnostic outputs
    real(kind_phys),     intent(out)   :: taurx(:,:)                   ! wave stress in zonal direction, interfaces [N m-2]
    real(kind_phys),     intent(out)   :: taury(:,:)                   ! wave stress in meridional direction, interfaces [N m-2]
    real(kind_phys),     intent(out)   :: tauardgx(:,:)                ! ridge based momentum flux profile, zonal, interfaces [N m-2]
    real(kind_phys),     intent(out)   :: tauardgy(:,:)                ! ridge based momentum flux profile, meridional, interfaces [N m-2]
    real(kind_phys),     intent(out)   :: utrdg(:,:)                   ! U tendency accummulator [m s-1]
    real(kind_phys),     intent(out)   :: vtrdg(:,:)                   ! V tendency accummulator [m s-1]
    real(kind_phys),     intent(out)   :: ttrdg(:,:)                   ! T tendency from orographic gravity wave drag [K s-1]

    ! Error handling
    character(len=512),  intent(out)   :: errmsg
    integer,             intent(out)   :: errflg

    ! Local variables
    real(kind_phys) :: isovar_zero(ncol)                               ! Zero (dummy) isotropic variance for gamma [m]
    real(kind_phys) :: isowgt_zero(ncol)                               ! Zero (dummy) isotropic weight for gamma [1]

    errmsg = ''
    errflg = 0

    ! Allocate and initialize zero arrays for gamma ridge scheme
    ! FIXME: only beta has isovar/isowgt support; in gamma this was omitted in the
    ! original code.
    isovar_zero(:) = 0._kind_phys
    isowgt_zero(:) = 0._kind_phys

    ! Apply negative value correction for gamma ridge maximum displacement
    where (mxdisg(:, :) < 0._kind_phys)
      mxdisg(:, :) = 0._kind_phys
    end where

    call gw_rdg_calc(band_oro          = band_oro, &
                     vramp             = vramp(:), &
                     type              = 'GAMMA', &
                     pver              = pver, &
                     ncol              = ncol, &
                     pcnst             = pcnst, &
                     prdg              = prdg, &
                     n_rdg             = n_rdg_gamma, &
                     dt                = dt, &
                     pi                = pi, &
                     cpair             = cpair, &
                     rair              = rair, &
                     u                 = u(:, :), &
                     v                 = v(:, :), &
                     t                 = t(:, :), &
                     p                 = p, &
                     piln              = piln(:, :), &
                     zm                = zm(:, :), &
                     zi                = zi(:, :), &
                     nm                = nm(:, :), &
                     ni                = ni(:, :), &
                     rhoi              = rhoi(:, :), &
                     kvtt              = kvt_gw(:, :), &
                     q                 = q(:, :, :), &
                     dse               = dse(:, :), &
                     effgw_rdg         = effgw_rdg_gamma, &
                     effgw_rdg_max     = effgw_rdg_gamma_max, &
                     effgw_rdg_resid   = effgw_rdg_resid, &
                     luse_gw_rdg_resid = use_gw_rdg_resid, &
                     hwdth             = hwdthg(:, :), &
                     clngt             = clngtg(:, :), &
                     gbxar             = gbxarg(:), &
                     mxdis             = mxdisg(:, :), &
                     angll             = angllg(:, :), &
                     anixy             = anixyg(:, :), &
                     isovar            = isovar_zero(:), &    ! FIXME: check inconsistency in gw. hplin 8/7/25
                     isowgt            = isowgt_zero(:), &    ! FIXME: check inconsistency in gw. hplin 8/7/25
                     rdg_cd_llb        = rdg_gamma_cd_llb, &
                     trpd_leewv        = trpd_leewv_rdg_gamma, &
                     u_tend            = u_tend(:, :), &
                     v_tend            = v_tend(:, :), &
                     s_tend            = s_tend(:, :), &
                     q_tend            = q_tend(:, :, :), &
                     flx_heat          = flx_heat(:), &
                     taurx             = taurx(:,:), &
                     taury             = taury(:,:), &
                     tauardgx          = tauardgx(:,:), &
                     tauardgy          = tauardgy(:,:), &
                     utrdg             = utrdg(:,:), &
                     vtrdg             = vtrdg(:,:), &
                     ttrdg             = ttrdg(:,:), &
                     errmsg            = errmsg, &
                     errflg            = errflg)

  end subroutine gravity_wave_drag_ridge_gamma_run

  subroutine gw_rdg_resid_src(ncol, pver, band, p, &
                              u, v, t, mxdis, kwvrdg, zi, nm, &
                              src_level, tend_level, tau, ubm, ubi, xv, yv, &
                              ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, c, tauoro, errmsg, errflg)
    use gw_common, only: rair, GWBand
    use gw_utils, only: dot_2d, midpoint_interp, get_unit_vector
    !-----------------------------------------------------------------------
    ! Orographic source for multiple gravity wave drag parameterization.
    !
    ! The stress is returned for a single wave with c=0, over orography.
    ! For points where the orographic variance is small (including ocean),
    ! the returned stress is zero.
    !------------------------------Arguments--------------------------------
    ! Column dimension.
    integer, intent(in) :: ncol
    integer, intent(in) :: pver

    ! Band to emit orographic waves in.
    ! Regardless, we will only ever emit into l = 0.
    type(GWBand), intent(in) :: band
    ! Pressure coordinates.
    type(coords1d), intent(in) :: p

    ! Midpoint zonal/meridional winds. ( m s-1)
    real(kind_phys), intent(in) :: u(ncol, pver), v(ncol, pver)
    ! Midpoint temperatures. (K)
    real(kind_phys), intent(in) :: t(ncol, pver)
    ! Height estimate for ridge (m) [anisotropic orography].
    real(kind_phys), intent(in) :: mxdis(ncol)
    ! horiz wavenumber [anisotropic orography].
    real(kind_phys), intent(in) :: kwvrdg(ncol)
    ! Interface altitudes above ground (m).
    real(kind_phys), intent(in) :: zi(ncol, pver + 1)
    ! Midpoint Brunt-Vaisalla frequencies (s-1).
    real(kind_phys), intent(in) :: nm(ncol, pver)

    ! Indices of top gravity wave source level and lowest level where wind
    ! tendencies are allowed.
    integer, intent(out) :: src_level(ncol)
    integer, intent(out) :: tend_level(ncol)

    ! Averages over source region.
    real(kind_phys), intent(out) :: nsrc(ncol) ! B-V frequency.
    real(kind_phys), intent(out) :: rsrc(ncol) ! Density.
    real(kind_phys), intent(out) :: usrc(ncol) ! Zonal wind.
    real(kind_phys), intent(out) :: vsrc(ncol) ! Meridional wind.
    real(kind_phys), intent(out) :: ubmsrc(ncol) ! On-obstacle wind.
    ! normalized wavenumber
    real(kind_phys), intent(out) :: m2src(ncol)

    ! Wave Reynolds stress.
    real(kind_phys), intent(out) :: tau(ncol, -band%ngwv:band%ngwv, pver + 1)
    ! Projection of wind at midpoints and interfaces.
    real(kind_phys), intent(out) :: ubm(ncol, pver), ubi(ncol, pver + 1)
    ! Unit vectors of source wind (zonal and meridional components).
    real(kind_phys), intent(out) :: xv(ncol), yv(ncol)
    ! Phase speeds.
    real(kind_phys), intent(out) :: c(ncol, -band%ngwv:band%ngwv)
    ! source level mom. flux
    real(kind_phys), intent(out) :: tauoro(ncol)
    character(len=512), intent(out):: errmsg
    integer, intent(out)           :: errflg

    !---------------------------Local Storage-------------------------------
    ! Column and level indices.
    integer :: i, k

    ! Surface streamline displacement height (2*sgh).
    real(kind_phys) :: hdsp(ncol)

    ! Difference in interface pressure across source region.
    real(kind_phys) :: dpsrc(ncol)
    ! Wind speed in source region.
    real(kind_phys) :: wmsrc(ncol)

    real(kind_phys) :: Fcrit_res, sghmax

    errmsg = ''
    errflg = 0

    !--------------------------------------------------------------------------
    ! Check that ngwav is equal to zero, otherwise end the job
    !--------------------------------------------------------------------------
    if (band%ngwv /= 0) then
      errmsg = 'gw_rdg_src :: ERROR - band%ngwv must be zero and it is not'
      errflg = 1
      return
    end if

    !--------------------------------------------------------------------------
    ! Average the basic state variables for the wave source over the depth of
    ! the orographic standard deviation. Here we assume that the appropiate
    ! values of wind, stability, etc. for determining the wave source are
    ! averages over the depth of the atmosphere penterated by the typical
    ! mountain.
    ! Reduces to the bottom midpoint values when mxdis=0, such as over ocean.
    !--------------------------------------------------------------------------
    Fcrit_res = 1.0_kind_phys
    hdsp = mxdis ! no longer multipied by 2
    where (hdsp < 10._kind_phys)
      hdsp = 0._kind_phys
    end where

    src_level = pver + 1

    tau(:, 0, :) = 0.0_kind_phys

    ! Find depth of "source layer" for mountain waves
    ! i.e., between ground and mountain top
    do k = pver, 1, -1
      do i = 1, ncol
        ! Need to have h >= z(k+1) here or code will bomb when h=0.
        if ((hdsp(i) >= zi(i, k + 1)) .and. (hdsp(i) < zi(i, k))) then
          src_level(i) = k
        end if
      end do
    end do

    rsrc = 0._kind_phys
    usrc = 0._kind_phys
    vsrc = 0._kind_phys
    nsrc = 0._kind_phys
    do i = 1, ncol
      do k = pver, src_level(i), -1
        rsrc(i) = rsrc(i) + p%mid(i, k)/(rair*t(i, k))*p%del(i, k)
        usrc(i) = usrc(i) + u(i, k)*p%del(i, k)
        vsrc(i) = vsrc(i) + v(i, k)*p%del(i, k)
        nsrc(i) = nsrc(i) + nm(i, k)*p%del(i, k)
      end do
    end do

    do i = 1, ncol
      dpsrc(i) = p%ifc(i, pver + 1) - p%ifc(i, src_level(i))
    end do

    rsrc = rsrc/dpsrc
    usrc = usrc/dpsrc
    vsrc = vsrc/dpsrc
    nsrc = nsrc/dpsrc

    ! Get the unit vector components and magnitude at the surface.
    call get_unit_vector(usrc, vsrc, xv, yv, wmsrc)

    ubmsrc = wmsrc

    ! Project the local wind at midpoints onto the source wind.
    do k = 1, pver
      ubm(:, k) = dot_2d(u(:, k), v(:, k), xv, yv)
    end do

    ! Compute the interface wind projection by averaging the midpoint winds.
    ! Use the top level wind at the top interface.
    ubi(:, 1) = ubm(:, 1)

    ubi(:, 2:pver) = midpoint_interp(ubm)

    ! The minimum stratification allowing GW behavior
    ! should really depend on horizontal scale since
    !
    !      m^2 ~ (N/U)^2 - k^2
    !
    m2src = ((nsrc/(ubmsrc + 0.01_kind_phys))**2 - kwvrdg**2)/((nsrc/(ubmsrc + 0.01_kind_phys))**2)

    ! Compute the interface wind projection by averaging the midpoint winds.
    ! Use the top level wind at the top interface.
    ubi(:, 1) = ubm(:, 1)
    ubi(:, 2:pver) = midpoint_interp(ubm)
    ubi(:, pver + 1) = ubm(:, pver)

    ! Determine the orographic c=0 source term following McFarlane (1987).
    ! (DOI: https://doi.org/10.1175/1520-0469(1987)044<1775:TEOOEG>2.0.CO;2)
    ! Set the source top interface index to pver, if the orographic term is
    ! zero.
    do i = 1, ncol
      if ((src_level(i) > 0) .and. (m2src(i) > orom2min)) then
        sghmax = Fcrit_res*(ubmsrc(i)/nsrc(i))**2
        tauoro(i) = 0.5_kind_phys*kwvrdg(i)*min(hdsp(i)**2, sghmax)* &
                    rsrc(i)*nsrc(i)*ubmsrc(i)
      else
        tauoro(i) = 0._kind_phys
      end if
    end do

    do i = 1, ncol
      do k = src_level(i), pver + 1
        tau(i, 0, k) = tauoro(i)
      end do
    end do

    ! Allow wind tendencies all the way to the model bottom.
    tend_level = pver

    ! No spectrum; phase speed is just 0.
    c = 0._kind_phys

  end subroutine gw_rdg_resid_src

  subroutine gw_rdg_src(ncol, pver, pi, rair, band, p, &
                        u, v, t, mxdis, angxy, anixy, kwvrdg, iso, zi, nm, &
                        src_level, tend_level, bwv_level, tlb_level, tau, ubm, ubi, xv, yv, &
                        ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, c, &
                        errmsg, errflg)

    use gw_common, only: GWBand
    use gw_utils, only: dot_2d, midpoint_interp
    !-----------------------------------------------------------------------
    ! Orographic source for multiple gravity wave drag parameterization.
    !
    ! The stress is returned for a single wave with c=0, over orography.
    ! For points where the orographic variance is small (including ocean),
    ! the returned stress is zero.
    !------------------------------Arguments--------------------------------
    ! Column dimension.
    integer, intent(in) :: ncol
    integer, intent(in) :: pver
    real(kind_phys), intent(in) :: pi
    real(kind_phys), intent(in) :: rair

    ! Band to emit orographic waves in.
    ! Regardless, we will only ever emit into l = 0.
    type(GWBand), intent(in) :: band
    ! Pressure coordinates.
    type(coords1d), intent(in) :: p

    ! Midpoint zonal/meridional winds. ( m s-1)
    real(kind_phys), intent(in) :: u(:, :), v(:, :)
    ! Midpoint temperatures. (K)
    real(kind_phys), intent(in) :: t(:, :)
    ! Height estimate for ridge (m) [anisotropic orography].
    real(kind_phys), intent(in) :: mxdis(:)
    ! Angle of ridge axis w/resp to north (degrees) [anisotropic orography].
    real(kind_phys), intent(in) :: angxy(:)
    ! Anisotropy parameter [anisotropic orography].
    real(kind_phys), intent(in) :: anixy(:)
    ! horiz wavenumber [anisotropic orography].
    real(kind_phys), intent(in) :: kwvrdg(:)
    ! Isotropic source flag [anisotropic orography].
    integer, intent(in)  :: iso(:)
    ! Interface altitudes above ground (m).
    real(kind_phys), intent(in) :: zi(:, :)
    ! Midpoint Brunt-Vaisalla frequencies (s-1).
    real(kind_phys), intent(in) :: nm(:, :)

    ! Indices of top gravity wave source level and lowest level where wind
    ! tendencies are allowed.
    integer, intent(out) :: src_level(:)
    integer, intent(out) :: tend_level(:)
    integer, intent(out) :: bwv_level(:), tlb_level(:)

    ! Averages over source region.
    real(kind_phys), intent(out) :: nsrc(:) ! B-V frequency.
    real(kind_phys), intent(out) :: rsrc(:) ! Density.
    real(kind_phys), intent(out) :: usrc(:) ! Zonal wind.
    real(kind_phys), intent(out) :: vsrc(:) ! Meridional wind.
    real(kind_phys), intent(out) :: ubmsrc(:) ! On-ridge wind.
    ! Top of low-level flow layer.
    real(kind_phys), intent(out) :: tlb(:)
    ! Bottom of linear wave region.
    real(kind_phys), intent(out) :: bwv(:)
    ! normalized wavenumber
    real(kind_phys), intent(out) :: m2src(:)

    ! Wave Reynolds stress.
    real(kind_phys), intent(out) :: tau(ncol, -band%ngwv:band%ngwv, pver + 1)
    ! Projection of wind at midpoints and interfaces.
    real(kind_phys), intent(out) :: ubm(:, :), ubi(:, :)
    ! Unit vectors of source wind (zonal and meridional components).
    real(kind_phys), intent(out) :: xv(:), yv(:)
    ! Phase speeds.
    real(kind_phys), intent(out) :: c(ncol, -band%ngwv:band%ngwv)
    ! Froude numbers for flow/drag regimes
    real(kind_phys), intent(out) :: Fr1(:), Fr2(:), Frx(:)

    character(len=512), intent(out)               :: errmsg
    integer, intent(out)                          :: errflg

    ! Local variables
    integer :: i, k

    ! Surface streamline displacement height (2*sgh).
    real(kind_phys) :: hdsp(ncol)

    ! Difference in interface pressure across source region.
    real(kind_phys) :: dpsrc(ncol)
    ! Thickness of downslope wind region.
    real(kind_phys) :: ddw(ncol)
    ! Thickness of linear wave region.
    real(kind_phys) :: dwv(ncol)
    ! Wind speed in source region.
    real(kind_phys) :: wmsrc(ncol)

    real(kind_phys) :: ragl(ncol)

    errmsg = ''
    errflg = 0

    !--------------------------------------------------------------------------
    ! Check that ngwv is equal to zero, otherwise end the job
    !--------------------------------------------------------------------------
    if(band%ngwv /= 0) then
      errmsg = 'gw_rdg_src: band%ngwv must be zero but it is not'
      errflg = 1
      return
    endif

    !--------------------------------------------------------------------------
    ! Average the basic state variables for the wave source over the depth of
    ! the orographic standard deviation. Here we assume that the appropiate
    ! values of wind, stability, etc. for determining the wave source are
    ! averages over the depth of the atmosphere penterated by the typical
    ! mountain.
    ! Reduces to the bottom midpoint values when mxdis=0, such as over ocean.
    !--------------------------------------------------------------------------
    hdsp = mxdis ! no longer multipied by 2
    src_level = pver + 1
    bwv_level = -1
    tlb_level = -1

    tau(:, 0, :) = 0.0_kind_phys

    ! Find depth of "source layer" for mountain waves
    ! i.e., between ground and mountain top
    do k = pver, 1, -1
      do i = 1, ncol
        ! Need to have h >= z(k+1) here or code will bomb when h=0.
        if ((hdsp(i) >= zi(i, k + 1)) .and. (hdsp(i) < zi(i, k))) then
          src_level(i) = k
        end if
      end do
    end do

    rsrc = 0._kind_phys
    usrc = 0._kind_phys
    vsrc = 0._kind_phys
    nsrc = 0._kind_phys
    do i = 1, ncol
      do k = pver, src_level(i), -1
        rsrc(i) = rsrc(i) + p%mid(i, k)/(rair*t(i, k))*p%del(i, k)
        usrc(i) = usrc(i) + u(i, k)*p%del(i, k)
        vsrc(i) = vsrc(i) + v(i, k)*p%del(i, k)
        nsrc(i) = nsrc(i) + nm(i, k)*p%del(i, k)
      end do
    end do

    do i = 1, ncol
      dpsrc(i) = p%ifc(i, pver + 1) - p%ifc(i, src_level(i))
    end do

    rsrc = rsrc/dpsrc
    usrc = usrc/dpsrc
    vsrc = vsrc/dpsrc
    nsrc = nsrc/dpsrc

    wmsrc = sqrt(usrc**2 + vsrc**2)

    ! Get the unit vector components
    ! Want agl=0 with U>0 to give xv=1

    ragl = angxy*pi/180._kind_phys

    ! protect from wierd "bad" angles
    ! that may occur if hdsp is zero
    where (hdsp <= orohmin)
      ragl = 0._kind_phys
    end where

    yv = -sin(ragl)
    xv = cos(ragl)

    ! Kluge in possible "isotropic" obstacle.
    where ((iso == 1) .and. (wmsrc > orovmin))
      xv = usrc/wmsrc
      yv = vsrc/wmsrc
    end where

    ! Project the local wind at midpoints into the on-ridge direction
    do k = 1, pver
      ubm(:, k) = dot_2d(u(:, k), v(:, k), xv, yv)
    end do
    ubmsrc = dot_2d(usrc, vsrc, xv, yv)

    ! Ensure on-ridge wind is positive at source level
    do k = 1, pver
      ubm(:, k) = sign(ubmsrc*0._kind_phys + 1._kind_phys, ubmsrc)*ubm(:, k)
    end do

    ! Sean says just use 1._kind_phys as
    ! first argument
    xv = sign(ubmsrc*0._kind_phys + 1._kind_phys, ubmsrc)*xv
    yv = sign(ubmsrc*0._kind_phys + 1._kind_phys, ubmsrc)*yv

    ! Now make ubmsrc positive and protect
    ! against zero
    ubmsrc = abs(ubmsrc)
    ubmsrc = max(0.01_kind_phys, ubmsrc)

    ! The minimum stratification allowing GW behavior
    ! should really depend on horizontal scale since
    !
    !      m^2 ~ (N/U)^2 - k^2
    !
    ! Should also think about parameterizing
    ! trapped lee-waves.

    ! This needs to be made constistent with later
    ! treatment of nonhydrostatic effects.
    m2src = ((nsrc/(ubmsrc + 0.01_kind_phys))**2 - kwvrdg**2)/((nsrc/(ubmsrc + 0.01_kind_phys))**2)

    !-------------------------------------------------------------
    ! Calculate provisional limits (in Z [m]) for 3 regimes. This
    ! will modified later if wave breaking or trapping are
    ! diagnosed
    !
    !                                            ^
    !                                            | *** linear propagation ***
    !  (H) -------- mountain top -------------   | *** or wave breaking  ****
    !                                            | *** regimes  *************
    ! (BWV)------ bottom of linear waves ----    |
    !                    :                       |
    !                 *******                    |
    !                    :                       |
    ! (TLB)--- top of flow diversion layer---    '
    !                   :
    !        **** flow diversion *****
    !                    :
    !============================================

    !============================================
    ! For Dividing streamline para (DS2017)
    !--------------------------------------------
    ! High-drag downslope wind regime exists
    ! between bottom of linear waves and top of
    ! flow diversion. Linear waves can only
    ! attain vertical displacment of f1*U/N. So,
    ! bottom of linear waves is given by
    !
    !        BWV = H - Fr1*U/N
    !
    ! Downslope wind layer begins at BWV and
    ! extends below it until some maximum high
    ! drag obstacle height Fr2*U/N is attained
    ! (where Fr2 >= f1).  Below downslope wind
    ! there is flow diversion, so top of
    ! diversion layer (TLB) is equivalent to
    ! bottom of downslope wind layer and is;
    !
    !       TLB = H - Fr2*U/N
    !
    !-----------------------------------------

    ! Critical inverse Froude number
    !-----------------------------------------------
    Fr1(:) = Fr_c*1.00_kind_phys
    Frx(:) = hdsp(:)*nsrc(:)/abs(ubmsrc(:))/Fr_c

    if (do_divstream) then
      !------------------------------------------------
      ! Calculate Fr2(Frx) for DS2017
      !------------------------------------------------
      where (Frx <= Frx0)
        Fr2(:) = Fr1(:) + Fr1(:)*C_GammaMax*anixy(:)
      elsewhere((Frx > Frx0) .and. (Frx <= Frx1))
        Fr2(:) = Fr1(:) + Fr1(:)*C_GammaMax*anixy(:) &
                 *(Frx1 - Frx(:))/(Frx1 - Frx0)
      elsewhere(Frx > Frx1)
        Fr2(:) = Fr1(:)
      end where
    else
      !------------------------------------------
      ! Regime distinctions entirely carried by
      ! amplification of taudsw (next subr)
      !------------------------------------------
      Fr2(:) = Fr1(:)
    end if

    where (m2src > orom2min)
      ddw = Fr2*(abs(ubmsrc))/nsrc
    elsewhere
      ddw = 0._kind_phys
    end where

    ! If TLB is less than zero then obstacle is not
    ! high enough to produce an low-level diversion layer
    tlb = mxdis - ddw
    where (tlb < 0._kind_phys)
      tlb = 0._kind_phys
    end where
    do k = pver, pver/2, -1
      do i = 1, ncol
        if ((tlb(i) > zi(i, k + 1)) .and. (tlb(i) <= zi(i, k))) then
          tlb_level(i) = k
        end if
      end do
    end do

    ! Find *BOTTOM* of linear wave layer (BWV)
    !where ( nsrc > orostratmin )
    where (m2src > orom2min)
      dwv = Fr1*(abs(ubmsrc))/nsrc
    elsewhere
      dwv = -9.999e9_kind_phys ! if weak strat - no waves
    end where

    bwv = mxdis - dwv
    where ((bwv < 0._kind_phys) .or. (dwv < 0._kind_phys))
      bwv = 0._kind_phys
    end where
    do k = pver, 1, -1
      do i = 1, ncol
        if ((bwv(i) > zi(i, k + 1)) .and. (bwv(i) <= zi(i, k))) then
          bwv_level(i) = k + 1
        end if
      end do
    end do

    ! Compute the interface wind projection by averaging the midpoint winds.
    ! Use the top level wind at the top interface.
    ubi(:, 1) = ubm(:, 1)
    ubi(:, 2:pver) = midpoint_interp(ubm)
    ubi(:, pver + 1) = ubm(:, pver)

    ! Allow wind tendencies all the way to the model bottom.
    tend_level = pver

    ! No spectrum; phase speed is just 0.
    c = 0._kind_phys

    where (m2src < orom2min)
      tlb = mxdis
      tlb_level = src_level
    end where

  end subroutine gw_rdg_src

  subroutine gw_rdg_belowpeak(ncol, pver, band, rdg_cd_llb, &
                              t, mxdis, anixy, kwvrdg, zi, nm, ni, rhoi, &
                              src_level, tau, &
                              ubmsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, &
                              tauoro, taudsw, hdspwv, hdspdw)

    use gw_common, only: GWBand
    !-----------------------------------------------------------------------
    ! Orographic source for multiple gravity wave drag parameterization.
    !
    ! The stress is returned for a single wave with c=0, over orography.
    ! For points where the orographic variance is small (including ocean),
    ! the returned stress is zero.
    !------------------------------Arguments--------------------------------
    ! Column dimension.
    integer, intent(in) :: ncol

    integer, intent(in) :: pver
    ! Band to emit orographic waves in.
    ! Regardless, we will only ever emit into l = 0.
    type(GWBand), intent(in) :: band
    ! Drag coefficient for low-level flow
    real(kind_phys), intent(in) :: rdg_cd_llb

    ! Midpoint temperatures. (K)
    real(kind_phys), intent(in) :: t(:, :)
    ! Height estimate for ridge (m) [anisotropic orography].
    real(kind_phys), intent(in) :: mxdis(ncol)
    ! Anisotropy parameter [0-1] [anisotropic orography].
    real(kind_phys), intent(in) :: anixy(ncol)
    ! Inverse cross-ridge lengthscale (m-1) [anisotropic orography].
    real(kind_phys), intent(inout) :: kwvrdg(ncol)
    ! Interface altitudes above ground (m).
    real(kind_phys), intent(in) :: zi(:, :)
    ! Midpoint Brunt-Vaisalla frequencies (s-1).
    real(kind_phys), intent(in) :: nm(:, :)
    ! Interface Brunt-Vaisalla frequencies (s-1).
    real(kind_phys), intent(in) :: ni(:, :)
    ! Interface density (kg m-3).
    real(kind_phys), intent(in) :: rhoi(:, :)

    ! Indices of top gravity wave source level
    integer, intent(inout) :: src_level(ncol)

    ! Wave Reynolds stress.
    real(kind_phys), intent(inout) :: tau(ncol, -band%ngwv:band%ngwv, pver + 1)
    ! Top of low-level flow layer.
    real(kind_phys), intent(in)    :: tlb(ncol)
    ! Bottom of linear wave region.
    real(kind_phys), intent(in)    :: bwv(ncol)
    ! surface stress from linear waves.
    real(kind_phys), intent(out) :: tauoro(ncol)
    ! surface stress for downslope wind regime.
    real(kind_phys), intent(out) :: taudsw(ncol)

    ! Surface streamline displacement height for linear waves.
    real(kind_phys), intent(out) :: hdspwv(ncol)
    ! Surface streamline displacement height for downslope wind regime.
    real(kind_phys), intent(out) :: hdspdw(ncol)

    ! Froude numbers for flow/drag regimes
    real(kind_phys), intent(in) :: Fr1(ncol), Fr2(ncol), Frx(ncol)

    ! Averages over source region.
    real(kind_phys), intent(in) :: m2src(ncol) ! normalized non-hydro wavenumber
    real(kind_phys), intent(in) :: nsrc(ncol)  ! B-V frequency.
    real(kind_phys), intent(in) :: rsrc(ncol)  ! Density.
    real(kind_phys), intent(in) :: ubmsrc(ncol) ! On-ridge wind.

    !logical, intent(in), optional :: forcetlb

    !---------------------------Local Storage-------------------------------
    ! Column and level indices.
    integer :: i, k

    real(kind_phys) :: Coeff_LB(ncol), tausat, ubsrcx(ncol), dswamp
    real(kind_phys) :: taulin(ncol), BetaMax

    ! ubsrcx introduced to account for situations with high shear, strong strat.
    do i = 1, ncol
      ubsrcx(i) = max(ubmsrc(i), 0._kind_phys)
    end do

    do i = 1, ncol
      if (m2src(i) > orom2min) then
        hdspwv(i) = min(mxdis(i), Fr1(i)*ubsrcx(i)/nsrc(i))
      else
        hdspwv(i) = 0._kind_phys
      end if
    end do

    if (do_divstream) then
      do i = 1, ncol
        if (m2src(i) > orom2min) then
          hdspdw(i) = min(mxdis(i), Fr2(i)*ubsrcx(i)/nsrc(i))
        else
          hdspdw(i) = 0._kind_phys
        end if
      end do
    else
      do i = 1, ncol
        ! Needed only to mark where a DSW occurs
        if (m2src(i) > orom2min) then
          hdspdw(i) = mxdis(i)
        else
          hdspdw(i) = 0._kind_phys
        end if
      end do
    end if

    ! Calculate form drag coefficient ("CD")
    !--------------------------------------
    Coeff_LB = rdg_cd_llb*anixy

    ! Determine the orographic c=0 source term following McFarlane (1987).
    ! Set the source top interface index to pver, if the orographic term is
    ! zero.
    !
    ! This formula is basically from
    !
    !      tau(src) = rho * u' * w'
    ! where
    !      u' ~ N*h'  and w' ~ U*h'/b  (b="breite")
    !
    ! and 1/b has been replaced with k (kwvrdg)
    !
    do i = 1, ncol
      if ((src_level(i) > 0) .and. (m2src(i) > orom2min)) then
        tauoro(i) = kwvrdg(i)*(hdspwv(i)**2)*rsrc(i)*nsrc(i) &
                    *ubsrcx(i)
        taudsw(i) = kwvrdg(i)*(hdspdw(i)**2)*rsrc(i)*nsrc(i) &
                    *ubsrcx(i)
      else
        tauoro(i) = 0._kind_phys
        taudsw(i) = 0._kind_phys
      end if
    end do

    if (do_divstream) then
      do i = 1, ncol
        taulin(i) = 0._kind_phys
      end do
      !---------------------------------------
      ! Need linear drag when divstream is not used is used
      !---------------------------------------
    else
      do i = 1, ncol
        if ((src_level(i) > 0) .and. (m2src(i) > orom2min)) then
          taulin(i) = kwvrdg(i)*(mxdis(i)**2)*rsrc(i)*nsrc(i) &
                      *ubsrcx(i)
        else
          taulin(i) = 0._kind_phys
        end if
      end do
    end if

    if (do_divstream) then
      ! Amplify DSW between Frx=1. and Frx=Frx1
      do i = 1, ncol
        dswamp = 0._kind_phys
        BetaMax = C_BetaMax_DS*anixy(i)
        if ((Frx(i) > 1._kind_phys) .and. (Frx(i) <= Frx1)) then
          dswamp = (Frx(i) - 1._kind_phys)*(Frx1 - Frx(i))/(0.25_kind_phys*(Frx1 - 1._kind_phys)**2)
        end if
        taudsw(i) = (1._kind_phys + BetaMax*dswamp)*taudsw(i)
      end do
    else
      !-------------------
      ! Scinocca&McFarlane
      !--------------------
      do i = 1, ncol
        BetaMax = C_BetaMax_SM*anixy(i)
        if ((Frx(i) >= 1._kind_phys) .and. (Frx(i) < 1.5_kind_phys)) then
          dswamp = 2._kind_phys*BetaMax*(Frx(i) - 1._kind_phys)
        else if ((Frx(i) >= 1.5_kind_phys) .and. (Frx(i) < 3._kind_phys)) then
          dswamp = (1._kind_phys + BetaMax - (0.666_kind_phys**2))*(0.666_kind_phys*(3._kind_phys - Frx(i)))**2 &
                   + (1._kind_phys/Frx(i))**2 - 1._kind_phys
        else
          dswamp = 0._kind_phys
        end if
        if ((Frx(i) >= 1._kind_phys) .and. (Frx(i) < 3._kind_phys)) then
          taudsw(i) = (1._kind_phys + dswamp)*taulin(i) - tauoro(i)
        else
          taudsw(i) = 0._kind_phys
        end if
        ! This code defines "taudsw" as SUM of freely-propagating
        ! DSW enhancement. Different than in SM2000
        taudsw(i) = taudsw(i) + tauoro(i)
      end do
      !----------------------------------------------------
    end if

    do i = 1, ncol
      if (m2src(i) > orom2min) then
        where ((zi(i, :) < mxdis(i)) .and. (zi(i, :) >= bwv(i)))
          tau(i, 0, :) = tauoro(i)
          else where ((zi(i, :) < bwv(i)) .and. (zi(i, :) >= tlb(i)))
          tau(i, 0, :) = tauoro(i) + (taudsw(i) - tauoro(i))* &
                         (bwv(i) - zi(i, :))/ &
                         (bwv(i) - tlb(i))
        end where

        ! low-level form drag on obstacle. Quantity kwvrdg (~1/b) appears for consistency
        ! with tauoro and taudsw forms. Should be weighted by L*b/A_g before applied to flow.
        where ((zi(i, :) < tlb(i)) .and. (zi(i, :) >= 0._kind_phys))
          tau(i, 0, :) = taudsw(i) + &
                         Coeff_LB(i)*kwvrdg(i)*rsrc(i)*0.5_kind_phys*(ubsrcx(i)**2)*(tlb(i) - zi(i, :))
        end where

        if (do_smooth_regimes) then
          ! This blocks accounts for case where both mxdis and tlb fall
          ! between adjacent edges
          do k = 1, pver
            if ((zi(i, k) >= tlb(i)) .and. (zi(i, k + 1) < tlb(i)) .and. &
                (zi(i, k) >= mxdis(i)) .and. (zi(i, k + 1) < mxdis(i))) then
              src_level(i) = src_level(i) - 1
              tau(i, 0, k) = tauoro(i)
            end if
          end do
        end if

      else
        ! This block allows low-level dynamics to occur
        ! even if m2 is less than orom2min
        where ((zi(i, :) < tlb(i)) .and. (zi(i, :) >= 0._kind_phys))
          tau(i, 0, :) = taudsw(i) + &
                         Coeff_LB(i)*kwvrdg(i)*rsrc(i)*0.5_kind_phys* &
                         (ubsrcx(i)**2)*(tlb(i) - zi(i, :))
        end where
      end if
    end do

    ! This may be redundant with newest version of gw_drag_prof.
    ! That code reaches down to level k=src_level+1. (jtb 1/5/16)
    do i = 1, ncol
      k = src_level(i)
      if (ni(i, k) > orostratmin) then
        tausat = (Fr_c**2)*kwvrdg(i)*rhoi(i, k)*ubsrcx(i)**3/ &
                 (1._kind_phys*ni(i, k))
      else
        tausat = 0._kind_phys
      end if
      tau(i, 0, src_level(i)) = min(tauoro(i), tausat)
    end do

    ! Final clean-up. Do nothing if obstacle less than orohmin
    do i = 1, ncol
      if (mxdis(i) < orohmin) then
        tau(i, 0, :) = 0._kind_phys
        tauoro(i) = 0._kind_phys
        taudsw(i) = 0._kind_phys
      end if
    end do

    ! Disable vertical propagation if Scorer param is
    ! too small.
    do i = 1, ncol
      if (m2src(i) <= orom2min) then
        src_level(i) = 1
      end if
    end do

  end subroutine gw_rdg_belowpeak

  ! Parameterization of high-drag regimes and trapped lee-waves for CAM
  subroutine gw_rdg_break_trap(ncol, pver, pi, band, &
                               zi, nm, ni, ubm, ubi, rhoi, kwvrdg, bwv, tlb, wbr, &
                               src_level, tlb_level, &
                               hdspwv, hdspdw, mxdis, &
                               tauoro, taudsw, tau, &
                               ldo_trapped_waves, wdth_kwv_scale_in)
    use gw_common, only: GWBand

    ! Input arguments
    integer, intent(in) :: ncol
    integer, intent(in) :: pver
    real(kind_phys), intent(in) :: pi
    ! Band to emit orographic waves in.
    ! Regardless, we will only ever emit into l = 0.
    type(GWBand), intent(in) :: band

    ! Height estimate for ridge (m) [anisotropic orography].
    !real(kind_phys), intent(in) :: mxdis(:)
    ! Horz wavenumber for ridge (1/m) [anisotropic orography].
    real(kind_phys), intent(in) :: kwvrdg(:)
    ! Interface altitudes above ground (m).
    real(kind_phys), intent(in) :: zi(:, :)
    ! Midpoint Brunt-Vaisalla frequencies (s-1).
    real(kind_phys), intent(in) :: nm(:, :)
    ! Interface Brunt-Vaisalla frequencies (s-1).
    real(kind_phys), intent(in) :: ni(:, :)

    ! Indices of gravity wave sources.
    integer, intent(inout) :: src_level(:), tlb_level(:)

    ! Wave Reynolds stress.
    real(kind_phys), intent(inout) :: tau(ncol, -band%ngwv:band%ngwv, pver + 1)
    ! Wave Reynolds stresses at source.
    real(kind_phys), intent(in)    :: taudsw(:)
    real(kind_phys), intent(inout) :: tauoro(:)
    ! Projection of wind at midpoints and interfaces.
    real(kind_phys), intent(in) :: ubm(:, :)
    real(kind_phys), intent(in) :: ubi(:, :)
    ! Interface density (kg m-3).
    real(kind_phys), intent(in) :: rhoi(:, :)

    ! Top of low-level flow layer.
    real(kind_phys), intent(in) :: tlb(:)
    ! Bottom of linear wave region.
    real(kind_phys), intent(in) :: bwv(:)

    ! Surface streamline displacement height for linear waves.
    real(kind_phys), intent(in) :: hdspwv(:)
    ! Surface streamline displacement height for downslope wind regime.
    real(kind_phys), intent(in) :: hdspdw(:)
    ! Ridge height.
    real(kind_phys), intent(in) :: mxdis(:)

    ! Wave breaking level
    real(kind_phys), intent(out) :: wbr(:)

    logical, intent(in), optional :: ldo_trapped_waves
    real(kind_phys), intent(in), optional :: wdth_kwv_scale_in

    !---------------------------Local Storage-------------------------------
    ! Column and level indices.
    integer :: i, k, kp1, non_hydro
    real(kind_phys):: m2(ncol, pver), delz(ncol), tausat(ncol), trn(ncol)
    real(kind_phys):: wbrx(ncol)
    real(kind_phys):: phswkb(ncol, pver + 1)
    logical :: lldo_trapped_waves
    real(kind_phys):: wdth_kwv_scale
    ! Indices of important levels.
    integer :: trn_level(ncol)

    if (present(ldo_trapped_waves)) then
      lldo_trapped_waves = ldo_trapped_waves
      if (lldo_trapped_waves) then
        non_hydro = 1
      else
        non_hydro = 0
      end if
    else
      lldo_trapped_waves = .false.
      non_hydro = 0
    end if

    if (present(wdth_kwv_scale_in)) then
      wdth_kwv_scale = wdth_kwv_scale_in
    else
      wdth_kwv_scale = 1._kind_phys
    end if

    ! Calculate vertical wavenumber**2
    !---------------------------------
    m2 = (nm/(abs(ubm) + .01_kind_phys))**2
    do k = pver, 1, -1
      m2(:, k) = m2(:, k) - non_hydro*(wdth_kwv_scale*kwvrdg)**2
      ! sweeping up, zero out m2 above first occurence
      ! of m2(:,k)<=0
      kp1 = min(k + 1, pver)
      where ((m2(:, k) <= 0.0_kind_phys) .or. (m2(:, kp1) <= 0.0_kind_phys))
        m2(:, k) = 0._kind_phys
      end where
    end do

    ! Take square root of m**2 and
    ! do vertical integral to find
    ! WKB phase.
    !-----------------------------
    m2 = SQRT(m2)
    phswkb(:, :) = 0
    do k = pver, 1, -1
      where (zi(:, k) > tlb(:))
        delz(:) = min(zi(:, k) - zi(:, k + 1), zi(:, k) - tlb(:))
        phswkb(:, k) = phswkb(:, k + 1) + m2(:, k)*delz(:)
      end where
    end do

    ! Identify top edge of layer in which phswkb reaches 3*pi/2
    ! - approximately the "breaking level"
    !----------------------------------------------------------
    wbr(:) = 0._kind_phys
    wbrx(:) = 0._kind_phys
    if (do_smooth_regimes) then
      do k = pver, 1, -1
        where ((phswkb(:, k + 1) < 1.5_kind_phys*pi) .and. (phswkb(:, k) >= 1.5_kind_phys*pi) &
               .and. (hdspdw(:) > hdspwv(:)))
          wbr(:) = zi(:, k)
          ! Extrapolation to make regime
          ! transitions smoother
          wbrx(:) = zi(:, k) - (phswkb(:, k) - 1.5_kind_phys*pi) &
                    /(m2(:, k) + 1.e-6_kind_phys)
          src_level(:) = k - 1
        end where
      end do
    else
      do k = pver, 1, -1
        where ((phswkb(:, k + 1) < 1.5_kind_phys*pi) .and. (phswkb(:, k) >= 1.5_kind_phys*pi) &
               .and. (hdspdw(:) > hdspwv(:)))
          wbr(:) = zi(:, k)
          src_level(:) = k
        end where
      end do
    end if

    ! Adjust tauoro at new source levels if needed.
    ! This is problematic if Fr_c<1.0. Not sure why.
    !----------------------------------------------------------
    if (do_adjust_tauoro) then
      do i = 1, ncol
        if (wbr(i) > 0._kind_phys) then
          tausat(i) = (Fr_c**2)*kwvrdg(i)*rhoi(i, src_level(i)) &
                      *abs(ubi(i, src_level(i)))**3 &
                      /ni(i, src_level(i))
          tauoro(i) = min(tauoro(i), tausat(i))
        end if
      end do
    end if

    if (do_smooth_regimes) then
      do i = 1, ncol
      do k = 1, pver + 1
        if ((zi(i, k) <= wbr(i)) .and. (zi(i, k) > tlb(i))) then
          tau(i, 0, k) = tauoro(i) + (taudsw(i) - tauoro(i))* &
                         (wbrx(i) - zi(i, k))/ &
                         (wbrx(i) - tlb(i))
          tau(i, 0, k) = max(tau(i, 0, k), tauoro(i))
        end if
      end do
      end do
    else
      ! Following is for backwards B4B compatibility with earlier versions
      ! ("N1" and "N5" -- Note: "N5" used do_backward_compat=.true.)
      if (.not. do_backward_compat) then
        do i = 1, ncol
        do k = 1, pver + 1
          if ((zi(i, k) < wbr(i)) .and. (zi(i, k) >= tlb(i))) then
            tau(i, 0, k) = tauoro(i) + (taudsw(i) - tauoro(i))* &
                           (wbr(i) - zi(i, k))/ &
                           (wbr(i) - tlb(i))
          end if
        end do
        end do
      else
        do i = 1, ncol
        do k = 1, pver + 1
          if ((zi(i, k) <= wbr(i)) .and. (zi(i, k) > tlb(i))) then
            tau(i, 0, k) = tauoro(i) + (taudsw(i) - tauoro(i))* &
                           (wbr(i) - zi(i, k))/ &
                           (wbr(i) - tlb(i))
          end if
        end do
        end do
      end if
    end if

    if (lldo_trapped_waves) then

      ! Identify top edge of layer in which Scorer param drops below 0
      ! - approximately the "turning level"
      !----------------------------------------------------------
      trn(:) = 1.e8_kind_phys
      trn_level(:) = 0 ! pver+1
      where (m2(:, pver) <= 0._kind_phys)
        trn(:) = zi(:, pver)
        trn_level(:) = pver
      end where
      do k = pver - 1, 1, -1
        where ((m2(:, k + 1) > 0._kind_phys) .and. (m2(:, k) <= 0._kind_phys))
          trn(:) = zi(:, k)
          trn_level(:) = k
        end where
      end do

      do i = 1, ncol
        ! Case: Turning below mountain top
        if ((trn(i) < mxdis(i)) .and. (trn_level(i) >= 1)) then
          tau(i, 0, :) = tau(i, 0, :) - max(tauoro(i), taudsw(i))
          tau(i, 0, :) = max(tau(i, 0, :), 0._kind_phys)
          tau(i, 0, 1:tlb_level(i)) = 0._kind_phys
          src_level(i) = 1 ! disable any more tau calculation
        end if
        ! Case: Turning but no breaking
        if ((wbr(i) == 0._kind_phys) .and. (trn(i) > mxdis(i)) .and. (trn_level(i) >= 1)) then
          where ((zi(i, :) <= trn(i)) .and. (zi(i, :) >= bwv(i)))
            tau(i, 0, :) = tauoro(i)* &
                           (trn(i) - zi(i, :))/ &
                           (trn(i) - bwv(i))
          end where
          src_level(i) = 1 ! disable any more tau calculation
        end if
        ! Case: Turning AND breaking. Turning ABOVE breaking
        if ((wbr(i) > 0._kind_phys) .and. (trn(i) >= wbr(i)) .and. (trn_level(i) >= 1)) then
          where ((zi(i, :) <= trn(i)) .and. (zi(i, :) >= wbr(i)))
            tau(i, 0, :) = tauoro(i)* &
                           (trn(i) - zi(i, :))/ &
                           (trn(i) - wbr(i))
          end where
          src_level(i) = 1 ! disable any more tau calculation
        end if
        ! Case: Turning AND breaking. Turning BELOW breaking
        if ((wbr(i) > 0._kind_phys) .and. (trn(i) < wbr(i)) .and. (trn_level(i) >= 1)) then
          tauoro(i) = 0._kind_phys
          where ((zi(i, :) < wbr(i)) .and. (zi(i, :) >= tlb(i)))
            tau(i, 0, :) = tauoro(i) + (taudsw(i) - tauoro(i))* &
                           (wbr(i) - zi(i, :))/ &
                           (wbr(i) - tlb(i))
          end where
          src_level(i) = 1 ! disable any more tau calculation
        end if
      end do
    end if

  end subroutine gw_rdg_break_trap

  subroutine gw_rdg_calc(band_oro, &
                         vramp, &
                         type, pver, ncol, pcnst, prdg, n_rdg, dt, pi, cpair, rair, &
                         u, v, t, p, piln, zm, zi, &
                         nm, ni, rhoi, kvtt, q, dse, &
                         effgw_rdg, effgw_rdg_max, &
                         effgw_rdg_resid, luse_gw_rdg_resid, &
                         hwdth, clngt, gbxar, &
                         mxdis, angll, anixy, &
                         isovar, isowgt, &
                         rdg_cd_llb, trpd_leewv, &
                         u_tend, v_tend, s_tend, q_tend, flx_heat, &
                         taurx, taury, &
                         tauardgx, tauardgy, &
                         utrdg, vtrdg, ttrdg, &
                         errmsg, errflg)

    use coords_1d, only: Coords1D
    use gw_common, only: GWBand, gw_drag_prof, energy_change

    type(GWBand), intent(in) :: band_oro
    real(kind_phys), intent(in) :: vramp(:)
    ! A mid-scale "band" with only stationary waves (l = 0).
    character(len=5), intent(in) :: type         ! BETA or GAMMA
    integer, intent(in) :: pver
    integer, intent(in) :: ncol         ! number of atmospheric columns
    integer, intent(in) :: pcnst        ! number of atmospheric constituents
    integer, intent(in) :: prdg         !
    integer, intent(in) :: n_rdg
    real(kind_phys), intent(in) :: dt           ! Time step.
    real(kind_phys), intent(in) :: pi
    real(kind_phys), intent(in) :: cpair
    real(kind_phys), intent(in) :: rair

    real(kind_phys), intent(in) :: u(:, :)    ! Midpoint zonal winds. ( m s-1)
    real(kind_phys), intent(in) :: v(:, :)    ! Midpoint meridional winds. ( m s-1)
    real(kind_phys), intent(in) :: t(:, :)    ! Midpoint temperatures. (K)
    type(coords1d), intent(in) :: p               ! Pressure coordinates.
    real(kind_phys), intent(in) :: piln(:, :)  ! Log of interface pressures.
    real(kind_phys), intent(in) :: zm(:, :)   ! Midpoint altitudes above ground (m).
    real(kind_phys), intent(in) :: zi(:, :) ! Interface altitudes above ground (m).
    real(kind_phys), intent(in) :: nm(:, :)   ! Midpoint Brunt-Vaisalla frequencies (s-1).
    real(kind_phys), intent(in) :: ni(:, :) ! Interface Brunt-Vaisalla frequencies (s-1).
    real(kind_phys), intent(in) :: rhoi(:, :) ! Interface density (kg m-3).
    real(kind_phys), intent(in) :: kvtt(:, :) ! Molecular thermal diffusivity.
    real(kind_phys), intent(in) :: q(:, :, :)        ! Constituent array.
    real(kind_phys), intent(in) :: dse(:, :)  ! Dry static energy.

    real(kind_phys), intent(in) :: effgw_rdg       ! Tendency efficiency.
    real(kind_phys), intent(in) :: effgw_rdg_max
    real(kind_phys), intent(in) :: effgw_rdg_resid  ! Tendency efficiency.
    logical, intent(in) :: luse_gw_rdg_resid ! On-Off switch
    real(kind_phys), intent(in) :: hwdth(ncol, prdg) ! width of ridges.
    real(kind_phys), intent(in) :: clngt(ncol, prdg) ! length of ridges.
    real(kind_phys), intent(in) :: gbxar(ncol)      ! gridbox area

    real(kind_phys), intent(in) :: mxdis(ncol, prdg) ! Height estimate for ridge (m).
    real(kind_phys), intent(in) :: angll(ncol, prdg) ! orientation of ridges.
    real(kind_phys), intent(in) :: anixy(ncol, prdg) ! Anisotropy parameter.

    real(kind_phys), intent(in) :: isovar(ncol)     ! sqrt of residual variance
    real(kind_phys), intent(in) :: isowgt(ncol)     ! area frac of residual variance

    real(kind_phys), intent(in) :: rdg_cd_llb      ! Drag coefficient for low-level flow
    logical, intent(in) :: trpd_leewv

    real(kind_phys), intent(inout) :: s_tend(:, :)   ! temperature (K)
    real(kind_phys), intent(inout) :: u_tend(:, :)   ! meridional wind
    real(kind_phys), intent(inout) :: v_tend(:, :)   ! zonal wind
    real(kind_phys), intent(inout) :: q_tend(:, :, :) ! constituent array

    real(kind_phys), intent(out) :: flx_heat(:)

    ! for diagnostic output
    real(kind_phys), intent(out) :: taurx(ncol, pver + 1) ! wave stress in zonal direction [N m-2]
    real(kind_phys), intent(out) :: taury(ncol, pver + 1) ! wave stress in meridional direction [N m-2]
    real(kind_phys), intent(out) :: tauardgx(ncol, pver + 1) ! ridge based momentum flux profile, zonal, interfaces [N m-2]
    real(kind_phys), intent(out) :: tauardgy(ncol, pver + 1) ! ridge based momentum flux profile, meridional, interfaces [N m-2]
    real(kind_phys), intent(out) :: utrdg(ncol, pver) ! tendency accummulators
    real(kind_phys), intent(out) :: vtrdg(ncol, pver)
    real(kind_phys), intent(out) :: ttrdg(ncol, pver)

    character(len=512), intent(out)               :: errmsg
    integer, intent(out)                          :: errflg
    character(len=*), parameter :: sub = 'gw_rdg_calc'

    !---------------------------Local storage-------------------------------

    integer :: k, m, nn

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
    real(kind_phys) :: ubm(ncol, pver)
    real(kind_phys) :: ubi(ncol, pver + 1)

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

    real(kind_phys) :: utgw(ncol, pver)       ! zonal wind tendency
    real(kind_phys) :: vtgw(ncol, pver)       ! meridional wind tendency
    real(kind_phys) :: ttgw(ncol, pver)       ! temperature tendency
    real(kind_phys) :: qtgw(ncol, pver, pcnst) ! constituents tendencies

    ! Effective gravity wave diffusivity at interfaces.
    real(kind_phys) :: egwdffi(ncol, pver + 1)

    ! Temperature tendencies from diffusion and kinetic energy.
    real(kind_phys) :: dttdf(ncol, pver)
    real(kind_phys) :: dttke(ncol, pver)

    ! Wave stress in zonal/meridional direction
    real(kind_phys) :: taurx0(ncol, pver + 1)
    real(kind_phys) :: taury0(ncol, pver + 1)
    ! Provisional absolute wave stress from gw_drag_prof
    real(kind_phys) :: tau_diag(ncol, pver + 1)

    ! Energy change used by fixer.
    real(kind_phys) :: de(ncol)

    ! Wavenumber fields
    real(kind_phys) :: tau(ncol, -band_oro%ngwv:band_oro%ngwv, pver+1) ! wave Reynolds stress
    real(kind_phys) :: gwut(ncol, pver, -band_oro%ngwv:band_oro%ngwv)  ! gravity wave wind tendency for each wave
    real(kind_phys) :: phase_speeds(ncol, -band_oro%ngwv:band_oro%ngwv)! Wave phase speeds for each column

    errmsg = ''
    errflg = 0

    ! initialize accumulated momentum fluxes and tendencies
    taurx = 0._kind_phys
    taury = 0._kind_phys
    ttrdg = 0._kind_phys
    utrdg = 0._kind_phys
    vtrdg = 0._kind_phys
    tau_diag = -9999._kind_phys

    do nn = 1, n_rdg
      kwvrdg = 0.001_kind_phys/(hwdth(:, nn) + 0.001_kind_phys) ! this cant be done every time step !!!
      isoflag = 0
      effgw = effgw_rdg*(hwdth(1:ncol, nn)*clngt(1:ncol, nn))/gbxar(1:ncol)
      effgw = min(effgw_rdg_max, effgw)

      call gw_rdg_src(ncol, pver, pi, rair, band_oro, p, &
                      u, v, t, mxdis(:, nn), angll(:, nn), anixy(:, nn), kwvrdg, isoflag, zi, nm, &
                      src_level, tend_level, bwv_level, tlb_level, tau, ubm, ubi, xv, yv, &
                      ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, phase_speeds, &
                      errmsg, errflg)

      if(errflg /= 0) return

      call gw_rdg_belowpeak(ncol, pver, band_oro, rdg_cd_llb, &
                            t, mxdis(:, nn), anixy(:, nn), kwvrdg, &
                            zi, nm, ni, rhoi, &
                            src_level, tau, &
                            ubmsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, &
                            tauoro, taudsw, hdspwv, hdspdw)

      call gw_rdg_break_trap(ncol, pver, pi, band_oro, &
                             zi, nm, ni, ubm, ubi, rhoi, kwvrdg, bwv, tlb, wbr, &
                             src_level, tlb_level, hdspwv, hdspdw, mxdis(:, nn), &
                             tauoro, taudsw, tau, &
                             ldo_trapped_waves=trpd_leewv)

      call gw_drag_prof(ncol, band_oro, p, src_level, tend_level, dt, &
                        t, vramp, &
                        piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                        effgw, phase_speeds, kvtt, q, dse, tau, utgw, vtgw, &
                        ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                        kwvrdg=kwvrdg, &
                        satfac_in=1._kind_phys, lapply_vdiff=gw_rdg_do_vdiff, tau_diag=tau_diag)

      ! Add the tendencies from each ridge to the totals.
      do k = 1, pver
        ! diagnostics
        utrdg(:, k) = utrdg(:, k) + utgw(:, k)
        vtrdg(:, k) = vtrdg(:, k) + vtgw(:, k)
        ttrdg(:, k) = ttrdg(:, k) + ttgw(:, k)
        ! physics tendencies
        u_tend(:ncol, k) = u_tend(:ncol, k) + utgw(:, k)
        v_tend(:ncol, k) = v_tend(:ncol, k) + vtgw(:, k)
        s_tend(:ncol, k) = s_tend(:ncol, k) + ttgw(:, k)
      end do

      do m = 1, pcnst
        do k = 1, pver
          q_tend(:ncol, k, m) = q_tend(:ncol, k, m) + qtgw(:, k, m)
        end do
      end do

      do k = 1, pver + 1
        taurx0(:, k) = tau(:, 0, k)*xv
        taury0(:, k) = tau(:, 0, k)*yv
        taurx(:, k) = taurx(:, k) + taurx0(:, k)
        taury(:, k) = taury(:, k) + taury0(:, k)
      end do

    end do ! end of loop over multiple ridges

    ! Save ridge-based momentum profile before residual is applied, for diagnostics.
    tauardgx = taurx
    tauardgy = taury

    if (luse_gw_rdg_resid) then
      ! Add additional GW from residual variance. Assumed isotropic
      kwvrdg = 0.001_kind_phys/(100._kind_phys)
      effgw = effgw_rdg_resid*isowgt
      tauoro = 0._kind_phys

      call gw_rdg_resid_src(ncol, pver, band_oro, p, &
                            u, v, t, isovar, kwvrdg, zi, nm, &
                            src_level, tend_level, tau, ubm, ubi, xv, yv, &
                            ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, phase_speeds, tauoro, errmsg, errflg)

      call gw_drag_prof(ncol, band_oro, p, src_level, tend_level, dt, &
                        t, vramp, &
                        piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                        effgw, phase_speeds, kvtt, q, dse, tau, utgw, vtgw, &
                        ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                        kwvrdg=kwvrdg, &
                        satfac_in=1._kind_phys, lapply_vdiff=gw_rdg_do_vdiff, tau_diag=tau_diag)

      ! Add the tendencies from isotropic residual to the totals.
      do k = 1, pver
        ! diagnostics
        utrdg(:, k) = utrdg(:, k) + utgw(:, k)
        vtrdg(:, k) = vtrdg(:, k) + vtgw(:, k)
        ttrdg(:, k) = ttrdg(:, k) + ttgw(:, k)
        ! physics tendencies
        u_tend(:ncol, k) = u_tend(:ncol, k) + utgw(:, k)
        v_tend(:ncol, k) = v_tend(:ncol, k) + vtgw(:, k)
        s_tend(:ncol, k) = s_tend(:ncol, k) + ttgw(:, k)
      end do

      do m = 1, pcnst
        do k = 1, pver
          q_tend(:ncol, k, m) = q_tend(:ncol, k, m) + qtgw(:, k, m)
        end do
      end do

      do k = 1, pver + 1
        taurx0(:, k) = tau(:, 0, k)*xv
        taury0(:, k) = tau(:, 0, k)*yv
        taurx(:, k) = taurx(:, k) + taurx0(:, k)
        taury(:, k) = taury(:, k) + taury0(:, k)
      end do
    end if ! end of residual variance calc

    ! Calculate energy change for output to CAM's energy checker.
    call energy_change(dt, p, u, v, u_tend(:ncol, :), &
                       v_tend(:ncol, :), s_tend(:ncol, :), de)
    flx_heat(:ncol) = de

    ! Convert T tendency to K s-1
    ! FIXME: some places use cpairv (e.g., orographic) but cpair is used here. hplin 8/14/25
    ttrdg = ttrdg / cpair

  end subroutine gw_rdg_calc

end module gravity_wave_drag_ridge
