! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013 and CCPPized in August 2025.
!
! Beres, J. H., M. J. Alexander, and J. R. Holton, 2004:
! A Method of Specifying the Gravity Wave Spectrum above Convection Based on Latent Heating Properties and Background Wind
! J. Atmos. Sci., 61, 324–337, https://doi.org/10.1175/1520-0469(2004)061<0324:AMOSTG>2.0.CO;2.
module gravity_wave_drag_convection
  use ccpp_kinds, only: kind_phys
  use gw_common, only: GWBand

  implicit none
  private

  ! Public CCPP-compliant interfaces.
  public :: gravity_wave_drag_convection_deep_init
  public :: gravity_wave_drag_convection_deep_run

  public :: gravity_wave_drag_convection_shallow_init
  public :: gravity_wave_drag_convection_shallow_run

  type :: BeresSourceDesc
    ! Whether wind speeds are shifted to be relative to storm cells.
    logical :: storm_shift
    ! Index for level where wind speed is used as the source speed.
    integer :: k
    ! Heating depths below this value [m] will be ignored.
    real(kind_phys) :: min_hdepth
    ! Table bounds, for convenience. (Could be inferred from shape(mfcc).)
    integer :: maxh
    integer :: maxuh
    ! Heating depths [m].
    real(kind_phys), allocatable :: hd(:)
    ! Table of source spectra.
    real(kind_phys), allocatable :: mfcc(:, :, :)
  end type BeresSourceDesc

  ! Beres settings and table.
  type(BeresSourceDesc), public :: beres_dp_desc
  type(BeresSourceDesc), public :: beres_sh_desc
  type(GWBand)          :: band_mid
  logical               :: is_band_initialized = .false.

contains

  ! Initialize gravity waves from convection and read in source spectra.
!> \section arg_table_gravity_wave_drag_convection_deep_init Argument Table
!! \htmlinclude gravity_wave_drag_convection_deep_init.html
  subroutine gravity_wave_drag_convection_deep_init(&
             pver, &
             masterproc, iulog, &
             gw_drag_file_dp, &
             pref_edge, &
             gw_delta_c, &
             pgwv, &
             errmsg, errflg)

    use gw_common, only: wavelength_mid

    integer, intent(in)                           :: pver
    logical, intent(in)                           :: masterproc
    integer, intent(in)                           :: iulog
    character(len=*),   intent(in)                :: gw_drag_file_dp
    real(kind_phys),    intent(in)                :: pref_edge(:)
    real(kind_phys),    intent(in)                :: gw_delta_c
    integer, intent(in)                           :: pgwv

    character(len=512), intent(out)               :: errmsg
    integer, intent(out)                          :: errflg

    integer :: k

    character(len=*), parameter :: sub = 'gravity_wave_drag_convection_deep_init'

    ! Initialize error variables
    errmsg = ''
    errflg = 0

    if (.not. is_band_initialized) then
      band_mid = GWBand(pgwv, gw_delta_c, 1.0_kind_phys, wavelength_mid)
      is_band_initialized = .true.
    endif

    ! Set the deep scheme specification components.
    beres_dp_desc%storm_shift = .false.

    do k = 0, pver
      ! 700 hPa index
      if (pref_edge(k + 1) < 70000._kind_phys) beres_dp_desc%k = k + 1
    end do

    if (masterproc) then
      write (iulog, *) sub // ': Beres deep level =', beres_dp_desc%k
    end if

    ! Don't use deep convection heating depths below this limit.
    ! This is important for QBO. Bad result if left at 2.5 km.
    beres_dp_desc%min_hdepth = 1000._kind_phys

    ! Check that deep gw file is set in namelist
    if (trim(gw_drag_file_dp) == "" .or. trim(gw_drag_file_dp) == "UNSET_PATH") then
      write (errmsg, '(a, a)') sub, "No gw_drag_file provided for Beres deep ", &
        "scheme. Set this via namelist."
      errflg = 1
      return
    end if

    call gw_init_beres_desc(masterproc, iulog, gw_drag_file_dp, band_mid, beres_dp_desc, errmsg, errflg)
  end subroutine gravity_wave_drag_convection_deep_init

  ! Convective gravity waves (Beres scheme, deep).
!> \section arg_table_gravity_wave_drag_convection_deep_run Argument Table
!! \htmlinclude gravity_wave_drag_convection_deep_run.html
  subroutine gravity_wave_drag_convection_deep_run(&
             ncol, pver, pverp, pcnst, &
             dt, &
             p, vramp, &
             pi, cpair, &
             effgw_beres_dp, &
             gw_apply_tndmax, &
             u, v, t, q, dse, &
             piln, &
             rhoi, nm, ni, &
             kvt_gw, &
             ttend_dp, &
             zm, &
             lat, &
             tend_q, tend_u, tend_v, tend_s, &
             flx_heat, &
             src_level, tend_level, ubm, ubi, xv, yv, &
             hdepth, maxq0, &
             utgw, vtgw, ttgw, qtgw, &
             egwdffi_tot, dttdf, dttke, &
             taucd_west, taucd_east, taucd_south, taucd_north, &
             errmsg, errflg)

    use coords_1d, only: coords1d
    use gw_common, only: energy_change, energy_fixer
    use gw_common, only: momentum_flux, momentum_fixer
    use gw_common, only: gw_drag_prof
    use gw_common, only: calc_taucd
    use gw_common, only: west, east, south, north

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: pverp
    integer,            intent(in)    :: pcnst
    real(kind_phys),    intent(in)    :: dt
    type(coords1d),     intent(in)    :: p                        ! Pressure coordinates [Pa]
    real(kind_phys),    intent(in)    :: vramp(:)                 ! Ramping profile for gravity wave drag [1]
    real(kind_phys),    intent(in)    :: pi                       ! Mathematical constant pi [1]
    real(kind_phys),    intent(in)    :: cpair                    ! Specific heat of dry air at constant pressure [J kg-1 K-1]
    real(kind_phys),    intent(in)    :: effgw_beres_dp           ! Efficiency factor for deep convective gravity waves [1]
    logical,            intent(in)    :: gw_apply_tndmax          ! Whether or not to apply tendency max [flag]
    real(kind_phys),    intent(in)    :: u(:,:)                   ! Zonal wind at midpoints [m s-1]
    real(kind_phys),    intent(in)    :: v(:,:)                   ! Meridional wind at midpoints [m s-1]
    real(kind_phys),    intent(in)    :: t(:,:)                   ! Temperature at midpoints [K]
    real(kind_phys),    intent(in)    :: q(:,:,:)                 ! Constituent mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)    :: dse(:,:)                 ! Dry static energy [J kg-1]
    real(kind_phys),    intent(in)    :: piln(:, :)               ! Natural logarithm of pressure at interfaces [ln(Pa)]
    real(kind_phys),    intent(in)    :: rhoi(:, :)               ! Density at interfaces [kg m-3]
    real(kind_phys),    intent(in)    :: nm(:, :)                 ! Brunt-Vaisalla frequency at midpoints [s-1]
    real(kind_phys),    intent(in)    :: ni(:, :)                 ! Brunt-Vaisalla frequency at interfaces [s-1]
    real(kind_phys),    intent(in)    :: kvt_gw(:, :)             ! Molecular thermal diffusivity at interfaces [m2 s-1]
    real(kind_phys),    intent(in)    :: ttend_dp(:,:)            ! Temperature tendency from deep convection [K s-1]
    real(kind_phys),    intent(in)    :: zm(:,:)                  ! Geopotential height at midpoints [m]
    real(kind_phys),    intent(in)    :: lat(:)                   ! Latitude [rad]

    real(kind_phys),    intent(inout) :: tend_q(:, :, :)          ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys),    intent(inout) :: tend_u(:, :)             ! Zonal wind tendency [m s-2]
    real(kind_phys),    intent(inout) :: tend_v(:, :)             ! Meridional wind tendency [m s-2]
    real(kind_phys),    intent(inout) :: tend_s(:, :)             ! Dry static energy tendency [J kg-1 s-1]
    real(kind_phys),    intent(inout) :: flx_heat(:)              ! Surface heat flux for energy conservation check [W m-2]

    integer,            intent(out)   :: src_level(:)             ! Vertical level index of gravity wave source [index]
    integer,            intent(out)   :: tend_level(:)            ! Lowest vertical level index where tendencies are applied [index]
    real(kind_phys),    intent(out)   :: ubm(:, :)                ! Wind projection at midpoints along source wind direction [m s-1]
    real(kind_phys),    intent(out)   :: ubi(:, :)                ! Wind projection at interfaces along source wind direction [m s-1]
    real(kind_phys),    intent(out)   :: xv(:)                    ! Zonal component of source wind unit vector [1]
    real(kind_phys),    intent(out)   :: yv(:)                    ! Meridional component of source wind unit vector [1]
    real(kind_phys),    intent(out)   :: hdepth(:)                ! Convective heating depth [m]
    real(kind_phys),    intent(out)   :: maxq0(:)                 ! Maximum daily heating rate [K day-1]
    real(kind_phys),    intent(out)   :: utgw(:, :)               ! Zonal wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)   :: vtgw(:, :)               ! Meridional wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)   :: ttgw(:, :)               ! Temperature tendency from gravity waves [K s-1]
    real(kind_phys),    intent(out)   :: qtgw(:, :, :)            ! Constituent tendencies from gravity waves [kg kg-1 s-1]
    real(kind_phys),    intent(inout) :: egwdffi_tot(:, :)        ! Effective diffusivity coefficient from gravity waves, interfaces [m2 s-1]
    real(kind_phys),    intent(out)   :: dttdf(:, :)              ! Temperature tendency from diffusion [K s-1]
    real(kind_phys),    intent(out)   :: dttke(:, :)              ! Temperature tendency from kinetic energy dissipation [K s-1]

    ! Copies of taucd in each direction for diagnostic.
    real(kind_phys),    intent(out)   :: taucd_west(:, :)         ! Reynolds stress for waves in W direction, interfaces [N m-2]
    real(kind_phys),    intent(out)   :: taucd_east(:, :)         ! Reynolds stress for waves in E direction, interfaces [N m-2]
    real(kind_phys),    intent(out)   :: taucd_south(:, :)        ! Reynolds stress for waves in S direction, interfaces [N m-2]
    real(kind_phys),    intent(out)   :: taucd_north(:, :)        ! Reynolds stress for waves in N direction, interfaces [N m-2]

    character(len=512), intent(out)   :: errmsg
    integer, intent(out)              :: errflg

    integer :: k, m

    ! Wavenumber fields
    real(kind_phys) :: tau(ncol, -band_mid%ngwv:band_mid%ngwv, pverp)
    real(kind_phys) :: gwut(ncol, pver, -band_mid%ngwv:band_mid%ngwv)
    real(kind_phys) :: phase_speeds(ncol, -band_mid%ngwv:band_mid%ngwv)

    ! Reynolds stress for waves propagating in each cardinal direction.
    real(kind_phys) :: taucd(ncol, pverp, 4)

    ! Momentum fluxes used by fixer.
    real(kind_phys) :: um_flux(ncol), vm_flux(ncol)

    ! Energy change used by fixer.
    real(kind_phys) :: de(ncol)

    real(kind_phys) :: effgw(ncol)

    real(kind_phys) :: egwdffi(ncol, pverp)

    errmsg = ''
    errflg = 0

    tau  = 0._kind_phys
    gwut = 0._kind_phys
    phase_speeds = 0._kind_phys

    egwdffi(:,:) = 0._kind_phys

    ! Efficiency of gravity wave momentum transfer.
    ! This is really only to remove the pole points.
    where (pi/2._kind_phys - abs(lat(:)) >= 4*epsilon(1._kind_phys))
      effgw(:) = effgw_beres_dp
    elsewhere
      effgw(:) = 0._kind_phys
    end where

    ! Determine wave sources for Beres deep scheme
    call gw_beres_src( &
      ncol        = ncol, &
      pver        = pver, &
      desc        = beres_dp_desc, &
      u           = u(:,:), &
      v           = v(:,:), &
      netdt       = ttend_dp(:,:), &
      zm          = zm(:,:), &
      ! Output arguments
      src_level   = src_level(:), &
      tend_level  = tend_level(:), &
      tau         = tau(:,-band_mid%ngwv:band_mid%ngwv,:pverp), &
      ubm         = ubm(:,:pver), &
      ubi         = ubi(:,:pverp), &
      xv          = xv(:), &
      yv          = yv(:), &
      c           = phase_speeds(:,-band_mid%ngwv:band_mid%ngwv), &
      hdepth      = hdepth(:), &
      maxq0       = maxq0(:))

    ! Solve for the drag profile with Beres source spectrum.
    call gw_drag_prof( &
      ncol                = ncol, &
      band                = band_mid, &
      p                   = p, &
      src_level           = src_level, &
      tend_level          = tend_level, &
      dt                  = dt, &
      t                   = t(:,:), &
      vramp               = vramp, &
      piln                = piln(:,:), &
      rhoi                = rhoi(:,:), &
      nm                  = nm(:,:), &
      ni                  = ni(:,:), &
      ubm                 = ubm(:,:), &
      ubi                 = ubi(:,:), &
      xv                  = xv(:), &
      yv                  = yv(:), &
      effgw               = effgw(:), &
      c                   = phase_speeds(:,-band_mid%ngwv:band_mid%ngwv), &
      kvtt                = kvt_gw(:,:), &
      q                   = q(:,:,:), &
      dse                 = dse(:,:), &
      lapply_effgw_in     = gw_apply_tndmax, &
      ! Input/output arguments
      tau                 = tau(:,-band_mid%ngwv:band_mid%ngwv,:pverp), &
      ! Output arguments
      utgw                = utgw(:,:pver), &
      vtgw                = vtgw(:,:pver), &
      ttgw                = ttgw(:,:pver), &
      qtgw                = qtgw(:,:pver,:pcnst), &
      egwdffi             = egwdffi(:,:pverp), &
      gwut                = gwut(:,:pver,-band_mid%ngwv:band_mid%ngwv), &
      dttdf               = dttdf(:,:pver), &
      dttke               = dttke(:,:pver))

    ! Project stress into directional components.
    taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

    ! Make copies for diagnostics.
    taucd_west(:,:pverp)  = taucd(:,:pverp,west)
    taucd_east(:,:pverp)  = taucd(:,:pverp,east)
    taucd_south(:,:pverp) = taucd(:,:pverp,south)
    taucd_north(:,:pverp) = taucd(:,:pverp,north)

    ! Add the diffusion coefficients
    do k = 1, pverp
      egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
    end do

    ! Store constituents tendencies
    do m = 1, pcnst
      do k = 1, pver
         tend_q(:,k,m) = tend_q(:,k,m) + qtgw(:,k,m)
      end do
    end do

    ! Find momentum flux, and use it to fix the wind tendencies below
    ! the gravity wave region.
    call momentum_flux(tend_level, taucd, um_flux, vm_flux)
    call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

    ! Add the momentum tendencies to the output tendency arrays.
    do k = 1, pver
      tend_u(:,k) = tend_u(:,k) + utgw(:,k)
      tend_v(:,k) = tend_v(:,k) + vtgw(:,k)
    end do

    ! Find energy change in the current state, and use fixer to apply
    ! the difference in lower levels.
    call energy_change(dt, p, u, v, tend_u(:,:), &
        tend_v(:,:), tend_s(:,:)+ttgw, de)
    call energy_fixer(tend_level, p, de-flx_heat(:), ttgw)

    do k = 1, pver
      tend_s(:,k) = tend_s(:,k) + ttgw(:,k)
    end do

    ! Change ttgw to a temperature tendency before outputting it.
    ! FIXME: some places use cpairv (e.g., orographic) but cpair is used here. hplin 8/14/25
    ttgw = ttgw / cpair

  end subroutine gravity_wave_drag_convection_deep_run


  ! Initialize gravity waves from convection and read in source spectra.
!> \section arg_table_gravity_wave_drag_convection_shallow_init Argument Table
!! \htmlinclude gravity_wave_drag_convection_shallow_init.html
  subroutine gravity_wave_drag_convection_shallow_init(&
             pver, &
             masterproc, iulog, &
             gw_drag_file_sh, &
             pref_edge, &
             gw_delta_c, &
             pgwv, &
             errmsg, errflg)

    use gw_common, only: wavelength_mid

    integer, intent(in)                           :: pver
    logical, intent(in)                           :: masterproc
    integer, intent(in)                           :: iulog
    character(len=*),   intent(in)                :: gw_drag_file_sh
    real(kind_phys),    intent(in)                :: pref_edge(:)
    real(kind_phys),    intent(in)                :: gw_delta_c
    integer, intent(in)                           :: pgwv

    character(len=512), intent(out)               :: errmsg
    integer, intent(out)                          :: errflg

    integer :: k

    character(len=*), parameter :: sub = 'gravity_wave_drag_convection_shallow_init'

    ! Initialize error variables
    errmsg = ''
    errflg = 0

    if (.not. is_band_initialized) then
      band_mid = GWBand(pgwv, gw_delta_c, 1.0_kind_phys, wavelength_mid)
      is_band_initialized = .true.
    endif

    ! Set the shallow scheme specification components.
    beres_sh_desc%storm_shift = .false.

    do k = 0, pver
      ! 900 hPa index
      if (pref_edge(k + 1) < 90000._kind_phys) beres_sh_desc%k = k + 1
    end do

    if (masterproc) then
      write (iulog, *) sub //': Beres shallow level =', beres_sh_desc%k
    end if

    ! Use all heating depths for shallow convection.
    beres_sh_desc%min_hdepth = 0._kind_phys

    ! Check that shallow gw file is set in namelist
    if (trim(gw_drag_file_sh) == "" .or. trim(gw_drag_file_sh) == "UNSET_PATH") then
      write (errmsg, '(a, a)') sub, "No gw_drag_file provided for Beres shallow ", &
        "scheme. Set this via namelist."
      errflg = 1
      return
    end if

    call gw_init_beres_desc(masterproc, iulog, gw_drag_file_sh, band_mid, beres_sh_desc, errmsg, errflg)
  end subroutine gravity_wave_drag_convection_shallow_init

  ! Convective gravity waves (Beres scheme, shallow).
!> \section arg_table_gravity_wave_drag_convection_shallow_run Argument Table
!! \htmlinclude gravity_wave_drag_convection_shallow_run.html
  subroutine gravity_wave_drag_convection_shallow_run(&
             ncol, pver, pverp, pcnst, &
             dt, &
             p, vramp, &
             pi, cpair, &
             effgw_beres_sh, &
             gw_apply_tndmax, &
             u, v, t, q, dse, &
             piln, &
             rhoi, nm, ni, &
             kvt_gw, &
             ttend_sh, &
             zm, &
             lat, &
             tend_q, tend_u, tend_v, tend_s, &
             flx_heat, &
             src_level, tend_level, ubm, ubi, xv, yv, &
             hdepth, maxq0, &
             utgw, vtgw, ttgw, qtgw, &
             egwdffi_tot, dttdf, dttke, &
             errmsg, errflg)

    use coords_1d, only: coords1d
    use gw_common, only: energy_change, energy_fixer
    use gw_common, only: momentum_flux, momentum_fixer
    use gw_common, only: gw_drag_prof
    use gw_common, only: calc_taucd

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: pverp
    integer,            intent(in)    :: pcnst
    real(kind_phys),    intent(in)    :: dt
    type(coords1d),     intent(in)    :: p                        ! Pressure coordinates [Pa]
    real(kind_phys),    intent(in)    :: vramp(:)                 ! Ramping profile for gravity wave drag [1]
    real(kind_phys),    intent(in)    :: pi                       ! Mathematical constant pi [1]
    real(kind_phys),    intent(in)    :: cpair                    ! Specific heat of dry air at constant pressure [J kg-1 K-1]
    real(kind_phys),    intent(in)    :: effgw_beres_sh           ! Efficiency factor for shallow convective gravity waves [1]
    logical,            intent(in)    :: gw_apply_tndmax          ! Whether or not to apply tendency max [flag]
    real(kind_phys),    intent(in)    :: u(:,:)                   ! Zonal wind at midpoints [m s-1]
    real(kind_phys),    intent(in)    :: v(:,:)                   ! Meridional wind at midpoints [m s-1]
    real(kind_phys),    intent(in)    :: t(:,:)                   ! Temperature at midpoints [K]
    real(kind_phys),    intent(in)    :: q(:,:,:)                 ! Constituent mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)    :: dse(:,:)                 ! Dry static energy [J kg-1]
    real(kind_phys),    intent(in)    :: piln(:, :)               ! Natural logarithm of pressure at interfaces [ln(Pa)]
    real(kind_phys),    intent(in)    :: rhoi(:, :)               ! Density at interfaces [kg m-3]
    real(kind_phys),    intent(in)    :: nm(:, :)                 ! Brunt-Vaisalla frequency at midpoints [s-1]
    real(kind_phys),    intent(in)    :: ni(:, :)                 ! Brunt-Vaisalla frequency at interfaces [s-1]
    real(kind_phys),    intent(in)    :: kvt_gw(:, :)             ! Molecular thermal diffusivity at interfaces [m2 s-1]
    real(kind_phys),    intent(in)    :: ttend_sh(:,:)            ! Temperature tendency from shallow convection [K s-1]
    real(kind_phys),    intent(in)    :: zm(:,:)                  ! Geopotential height at midpoints [m]
    real(kind_phys),    intent(in)    :: lat(:)                   ! Latitude [rad]

    real(kind_phys),    intent(inout) :: tend_q(:, :, :)          ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys),    intent(inout) :: tend_u(:, :)             ! Zonal wind tendency [m s-2]
    real(kind_phys),    intent(inout) :: tend_v(:, :)             ! Meridional wind tendency [m s-2]
    real(kind_phys),    intent(inout) :: tend_s(:, :)             ! Dry static energy tendency [J kg-1 s-1]
    real(kind_phys),    intent(inout) :: flx_heat(:)              ! Surface heat flux for energy conservation check [W m-2]

    integer,            intent(out)   :: src_level(:)             ! Vertical level index of gravity wave source [index]
    integer,            intent(out)   :: tend_level(:)            ! Lowest vertical level index where tendencies are applied [index]
    real(kind_phys),    intent(out)   :: ubm(:, :)                ! Wind projection at midpoints along source wind direction [m s-1]
    real(kind_phys),    intent(out)   :: ubi(:, :)                ! Wind projection at interfaces along source wind direction [m s-1]
    real(kind_phys),    intent(out)   :: xv(:)                    ! Zonal component of source wind unit vector [1]
    real(kind_phys),    intent(out)   :: yv(:)                    ! Meridional component of source wind unit vector [1]
    real(kind_phys),    intent(out)   :: hdepth(:)                ! Convective heating depth [m]
    real(kind_phys),    intent(out)   :: maxq0(:)                 ! Maximum daily heating rate [K day-1]
    real(kind_phys),    intent(out)   :: utgw(:, :)               ! Zonal wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)   :: vtgw(:, :)               ! Meridional wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)   :: ttgw(:, :)               ! Temperature tendency from gravity waves [K s-1]
    real(kind_phys),    intent(out)   :: qtgw(:, :, :)            ! Constituent tendencies from gravity waves [kg kg-1 s-1]
    real(kind_phys),    intent(inout) :: egwdffi_tot(:, :)        ! Effective diffusivity coefficient from gravity waves, interfaces [m2 s-1]
    real(kind_phys),    intent(out)   :: dttdf(:, :)              ! Temperature tendency from diffusion [K s-1]
    real(kind_phys),    intent(out)   :: dttke(:, :)              ! Temperature tendency from kinetic energy dissipation [K s-1]

    character(len=512), intent(out)   :: errmsg
    integer, intent(out)              :: errflg

    integer :: k, m

    ! Wavenumber fields
    real(kind_phys) :: tau(ncol, -band_mid%ngwv:band_mid%ngwv, pverp)
    real(kind_phys) :: gwut(ncol, pver, -band_mid%ngwv:band_mid%ngwv)
    real(kind_phys) :: phase_speeds(ncol, -band_mid%ngwv:band_mid%ngwv)

    ! Reynolds stress for waves propagating in each cardinal direction.
    real(kind_phys) :: taucd(ncol, pverp, 4)

    ! Momentum fluxes used by fixer.
    real(kind_phys) :: um_flux(ncol), vm_flux(ncol)

    ! Energy change used by fixer.
    real(kind_phys) :: de(ncol)

    real(kind_phys) :: effgw(ncol)

    real(kind_phys) :: egwdffi(ncol, pverp)

    errmsg = ''
    errflg = 0

    tau  = 0._kind_phys
    gwut = 0._kind_phys
    phase_speeds = 0._kind_phys

    egwdffi(:,:) = 0._kind_phys

    ! Efficiency of gravity wave momentum transfer.
    ! This is really only to remove the pole points.
    where (pi/2._kind_phys - abs(lat(:)) >= 4*epsilon(1._kind_phys))
      effgw(:) = effgw_beres_sh
    elsewhere
      effgw(:) = 0._kind_phys
    end where

    ! Determine wave sources for Beres deep scheme
    call gw_beres_src( &
      ncol        = ncol, &
      pver        = pver, &
      desc        = beres_sh_desc, &
      u           = u(:,:), &
      v           = v(:,:), &
      netdt       = ttend_sh(:,:), &
      zm          = zm(:,:), &
      ! Output arguments
      src_level   = src_level(:), &
      tend_level  = tend_level(:), &
      tau         = tau(:,-band_mid%ngwv:band_mid%ngwv,:pverp), &
      ubm         = ubm(:,:pver), &
      ubi         = ubi(:,:pverp), &
      xv          = xv(:), &
      yv          = yv(:), &
      c           = phase_speeds(:,-band_mid%ngwv:band_mid%ngwv), &
      hdepth      = hdepth(:), &
      maxq0       = maxq0(:))

    ! Solve for the drag profile with Beres source spectrum.
    call gw_drag_prof( &
      ncol                = ncol, &
      band                = band_mid, &
      p                   = p, &
      src_level           = src_level, &
      tend_level          = tend_level, &
      dt                  = dt, &
      t                   = t(:,:), &
      vramp               = vramp, &
      piln                = piln(:,:), &
      rhoi                = rhoi(:,:), &
      nm                  = nm(:,:), &
      ni                  = ni(:,:), &
      ubm                 = ubm(:,:), &
      ubi                 = ubi(:,:), &
      xv                  = xv(:), &
      yv                  = yv(:), &
      effgw               = effgw(:), &
      c                   = phase_speeds(:,-band_mid%ngwv:band_mid%ngwv), &
      kvtt                = kvt_gw(:,:), &
      q                   = q(:,:,:), &
      dse                 = dse(:,:), &
      lapply_effgw_in     = gw_apply_tndmax, &
      ! Input/output arguments
      tau                 = tau(:,-band_mid%ngwv:band_mid%ngwv,:pverp), &
      ! Output arguments
      utgw                = utgw(:,:pver), &
      vtgw                = vtgw(:,:pver), &
      ttgw                = ttgw(:,:pver), &
      qtgw                = qtgw(:,:pver,:pcnst), &
      egwdffi             = egwdffi(:,:pverp), &
      gwut                = gwut(:,:pver,-band_mid%ngwv:band_mid%ngwv), &
      dttdf               = dttdf(:,:pver), &
      dttke               = dttke(:,:pver))

    ! Project stress into directional components.
    taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

    ! Add the diffusion coefficients
    do k = 1, pverp
      egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
    end do

    ! Store constituents tendencies
    do m = 1, pcnst
      do k = 1, pver
         tend_q(:,k,m) = tend_q(:,k,m) + qtgw(:,k,m)
      end do
    end do

    ! Find momentum flux, and use it to fix the wind tendencies below
    ! the gravity wave region.
    call momentum_flux(tend_level, taucd, um_flux, vm_flux)
    call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

    ! Add the momentum tendencies to the output tendency arrays.
    do k = 1, pver
      tend_u(:,k) = tend_u(:,k) + utgw(:,k)
      tend_v(:,k) = tend_v(:,k) + vtgw(:,k)
    end do

    ! Find energy change in the current state, and use fixer to apply
    ! the difference in lower levels.
    call energy_change(dt, p, u, v, tend_u(:,:), &
        tend_v(:,:), tend_s(:,:)+ttgw, de)
    call energy_fixer(tend_level, p, de-flx_heat(:), ttgw)

    do k = 1, pver
      tend_s(:,k) = tend_s(:,k) + ttgw(:,k)
    end do

    ! Change ttgw to a temperature tendency before outputting it.
    ! FIXME: some places use cpairv (e.g., orographic) but cpair is used here. hplin 8/14/25
    ttgw = ttgw / cpair

  end subroutine gravity_wave_drag_convection_shallow_run

!==========================================================================

  ! Initialization / IO routine used by the two init routines.
  subroutine gw_init_beres_desc(amIRoot, iulog, file_path, band, desc, errmsg, errflg)
    use ccpp_io_reader, only: abstract_netcdf_reader_t, create_netcdf_reader_t

    logical,               intent(in)        :: amIRoot
    integer,               intent(in)        :: iulog
    type(GWBand),          intent(in)        :: band
    type(BeresSourceDesc), intent(inout)     :: desc
    character(len=*),      intent(in)        :: file_path
    character(len=512),    intent(out)       :: errmsg
    integer,               intent(out)       :: errflg

    character(len=*), parameter              :: sub = 'gw_init_beres_desc'
    real(kind_phys),  allocatable            :: tmp_var1d(:)
    real(kind_phys),  allocatable            :: file_mfcc(:, :, :) ! lookup table from the file f(depth, wind, phase speed)
    class(abstract_netcdf_reader_t), allocatable :: reader

    ! Number of wavenumbers in the input file.
    integer :: ngwv_file

    ! Read Beres file.
    reader = create_netcdf_reader_t()

    ! Open file
    call reader%open_file(file_path, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if

    ! Get HD (heating depth) dimension.
    call reader%get_var('HD', desc%hd, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if
    desc%maxh = size(desc%hd)

    ! Get MW (mean wind) dimension.
    call reader%get_var('MW', tmp_var1d, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if
    desc%maxuh = size(tmp_var1d)
    deallocate(tmp_var1d, stat=errflg)

    ! Get PS (phase speed) dimension.
    call reader%get_var('PS', tmp_var1d, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if
    ngwv_file = size(tmp_var1d)
    deallocate(tmp_var1d, stat=errflg)

    ! Number in each direction is half of total (and minus phase speed of 0).
    desc%maxuh = (desc%maxuh - 1)/2
    ngwv_file = (ngwv_file - 1)/2

    ! Check for inconsistency between file and model variable
    if (ngwv_file < band%ngwv) then
      write (errmsg, '(a, a, i4, a, i4)') sub, ": PhaseSpeed in lookup table file, ", &
        ngwv_file, ", does not cover the whole spectrum implied by the model ngwv, ", band%ngwv
      errflg = 1
      return
    end if

    ! While not currently documented in the file, it uses kilometers. Convert
    ! to meters.
    desc%hd = desc%hd*1000._kind_phys

    ! Allocate mfcc. "desc%maxh" and "desc%maxuh" are from the file, but the
    ! model determines wavenumber dimension.
    allocate (desc%mfcc(desc%maxh, -desc%maxuh:desc%maxuh, &
                        -band%ngwv:band%ngwv), stat=errflg, errmsg=errmsg)
    if(errflg /= 0) return

    ! Get mfcc data.
    call reader%get_var('mfcc', file_mfcc, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if

    desc%mfcc(:, -desc%maxuh:desc%maxuh, -band%ngwv:band%ngwv) = file_mfcc(:, :, ngwv_file-band%ngwv+1:ngwv_file+band%ngwv+1)
    ! 2*band%ngwv = (ngwv_file+band%ngwv+1-ngwv_file+band%ngwv-1)
    ! the file may cover more than what is requested of ngwv, so the bounds (lower and upper) to be specified explicitly.

    ! Close file
    call reader%close_file(errmsg, errflg)
    if (errflg /= 0) then
      return
    end if

    if (amIRoot) then
      write (iulog, *) "gravity_wave_drag_convection: Read in source spectra from file."
      write (iulog, *) "gravity_wave_drag_convection: mfcc max, min = ", &
        maxval(desc%mfcc), ", ", minval(desc%mfcc)
    end if

  end subroutine gw_init_beres_desc

  ! Driver for multiple gravity wave drag parameterizations.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  subroutine gw_beres_src(ncol, pver, &
                          desc, &
                          u, v, &
                          netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
                          c, hdepth, maxq0)
    use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp, index_of_nearest
    use gw_common, only: qbo_hdepth_scaling

!------------------------------Arguments--------------------------------
    ! Column dimension.
    integer, intent(in) :: ncol
    integer, intent(in) :: pver
    type(BeresSourceDesc) :: desc

    ! Midpoint zonal/meridional winds.
    real(kind_phys), intent(in) :: u(:, :), v(:, :)
    ! Heating rate due to convection.
    real(kind_phys), intent(in) :: netdt(:, :)
    ! Midpoint altitudes.
    real(kind_phys), intent(in) :: zm(:, :)

    ! Indices of top gravity wave source level and lowest level where wind
    ! tendencies are allowed.
    integer, intent(out) :: src_level(:)
    integer, intent(out) :: tend_level(:)

    ! Wave Reynolds stress.
    real(kind_phys), intent(out) :: tau(ncol, -band_mid%ngwv:band_mid%ngwv, pver + 1)
    ! Projection of wind at midpoints and interfaces.
    real(kind_phys), intent(out) :: ubm(:, :), ubi(:, :)
    ! Unit vectors of source wind (zonal and meridional components).
    real(kind_phys), intent(out) :: xv(:), yv(:)
    ! Phase speeds.
    real(kind_phys), intent(out) :: c(ncol, -band_mid%ngwv:band_mid%ngwv)

    ! Heating depth [m] and maximum heating in each column.
    real(kind_phys), intent(out) :: hdepth(:), maxq0(:)

!---------------------------Local Storage-------------------------------
    ! Column and level indices.
    integer :: i, k

    ! Zonal/meridional wind at roughly the level where the convection occurs.
    real(kind_phys) :: uconv(ncol), vconv(ncol)

    ! Maximum heating rate.
    real(kind_phys) :: q0(ncol)

    ! Bottom/top heating range index.
    integer  :: boti(ncol), topi(ncol)
    ! Index for looking up heating depth dimension in the table.
    integer  :: hd_idx(ncol)
    ! Mean wind in heating region.
    real(kind_phys) :: uh(ncol)
    ! Min/max wavenumber for critical level filtering.
    integer :: Umini(ncol), Umaxi(ncol)
    ! Source level tau for a column.
    real(kind_phys) :: tau0(-band_mid%ngwv:band_mid%ngwv)
    ! Speed of convective cells relative to storm.
    real(kind_phys) :: CS(ncol)
    ! Index to shift spectra relative to ground.
    integer :: shift

    ! Heating rate conversion factor.
    real(kind_phys), parameter :: CF = 20._kind_phys
    ! Averaging length.
    real(kind_phys), parameter :: AL = 1.0e5_kind_phys

    character(len=*), parameter :: sub = 'gw_beres_src'

    !----------------------------------------------------------------------
    ! Initialize tau array
    !----------------------------------------------------------------------

    tau = 0.0_kind_phys
    hdepth = 0.0_kind_phys
    q0 = 0.0_kind_phys
    tau0 = 0.0_kind_phys

    !------------------------------------------------------------------------
    ! Determine wind and unit vectors approximately at the source level, then
    ! project winds.
    !------------------------------------------------------------------------

    ! Source wind speed and direction.
    uconv = u(:, desc%k)
    vconv = v(:, desc%k)

    ! Get the unit vector components and magnitude at the source level.
    call get_unit_vector(uconv, vconv, xv, yv, ubi(:, desc%k + 1))

    ! Project the local wind at midpoints onto the source wind.
    do k = 1, pver
      ubm(:, k) = dot_2d(u(:, k), v(:, k), xv, yv)
    end do

    ! Compute the interface wind projection by averaging the midpoint winds.
    ! Use the top level wind at the top interface.
    ubi(:, 1) = ubm(:, 1)

    ubi(:, 2:pver) = midpoint_interp(ubm)

    !-----------------------------------------------------------------------
    ! Calculate heating depth.
    !
    ! Heating depth is defined as the first height range from the bottom in
    ! which heating rate is continuously positive.
    !-----------------------------------------------------------------------

    ! First find the indices for the top and bottom of the heating range.
    boti = 0
    topi = 0
    do k = pver, 1, -1
      do i = 1, ncol
        if (boti(i) == 0) then
          ! Detect if we are outside the top of range (where z = 20 km).
          if (zm(i, k) >= 20000._kind_phys) then
            boti(i) = k
            topi(i) = k
          else
            ! First spot where heating rate is positive.
            if (netdt(i, k) > 0.0_kind_phys) boti(i) = k
          end if
        end if
      end do
      ! When all done, exit
      if (all(boti /= 0)) exit
    end do

    do k = 1, pver
      do i = 1, ncol
        if (topi(i) == 0) then
          ! First spot where heating rate is positive.
          if ((netdt(i, k) > 0.0_kind_phys) .AND. (zm(i, k) <= 20000._kind_phys)) topi(i) = k - 1
        end if
      end do
      ! When all done, exit
      if (all(topi /= 0)) exit
    end do

    ! Heating depth in m.
    hdepth = [((zm(i, topi(i)) - zm(i, boti(i))), i=1, ncol)]

    ! J. Richter: this is an effective reduction of the GW phase speeds (needed to drive the QBO)
    hdepth = max(1000._kind_phys, hdepth*qbo_hdepth_scaling)

    hd_idx = index_of_nearest(hdepth, desc%hd)

    ! hd_idx=0 signals that a heating depth is too shallow, i.e. that it is
    ! either not big enough for the lowest table entry, or it is below the
    ! minimum allowed for this convection type.
    ! Values above the max in the table still get the highest value, though.
    where (hdepth < max(desc%min_hdepth, desc%hd(1))) hd_idx = 0

    ! Maximum heating rate.
    do k = minval(topi), maxval(boti)
      where (k >= topi .and. k <= boti)
        q0 = max(q0, netdt(:, k))
      end where
    end do

    !output max heating rate in K/day
    maxq0 = q0*24._kind_phys*3600._kind_phys

    ! Multipy by conversion factor
    q0 = q0*CF

    if (desc%storm_shift) then

      ! Find the cell speed where the storm speed is > 10 m/s.
      ! Storm speed is taken to be the source wind speed.
      CS = sign(max(abs(ubm(:, desc%k)) - 10._kind_phys, 0._kind_phys), ubm(:, desc%k))

      ! Average wind in heating region, relative to storm cells.
      uh = 0._kind_phys
      do k = minval(topi), maxval(boti)
        where (k >= topi .and. k <= boti)
          uh = uh + ubm(:, k)/(boti - topi + 1)
        end where
      end do

      uh = uh - CS

    else

      ! For shallow convection, wind is relative to ground, and "heating
      ! region" wind is just the source level wind.
      uh = ubm(:, desc%k)

    end if

    ! Limit uh to table range.
    uh = min(uh, real(desc%maxuh, kind_phys))
    uh = max(uh, -real(desc%maxuh, kind_phys))

    ! Speeds for critical level filtering.
    Umini = band_mid%ngwv
    Umaxi = -band_mid%ngwv
    do k = minval(topi), maxval(boti)
      where (k >= topi .and. k <= boti)
        Umini = min(Umini, nint(ubm(:, k)/band_mid%dc))
        Umaxi = max(Umaxi, nint(ubm(:, k)/band_mid%dc))
      end where
    end do

    Umini = max(Umini, -band_mid%ngwv)
    Umaxi = min(Umaxi, band_mid%ngwv)

    !-----------------------------------------------------------------------
    ! Gravity wave sources
    !-----------------------------------------------------------------------
    ! Start loop over all columns.
    !-----------------------------------------------------------------------
    do i = 1, ncol

      !---------------------------------------------------------------------
      ! Look up spectrum only if the heating depth is large enough, else set
      ! tau0 = 0.
      !---------------------------------------------------------------------

      if (hd_idx(i) > 0) then

        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------

        tau0 = desc%mfcc(hd_idx(i), nint(uh(i)), :)

        if (desc%storm_shift) then
          ! For deep convection, the wind was relative to storm cells, so
          ! shift the spectrum so that it is now relative to the ground.
          shift = -nint(CS(i)/band_mid%dc)
          tau0 = eoshift(tau0, shift)
        end if

        ! Adjust magnitude.
        tau0 = tau0*q0(i)*q0(i)/AL

        ! Adjust for critical level filtering.
        tau0(Umini(i):Umaxi(i)) = 0.0_kind_phys

        tau(i, :, topi(i) + 1) = tau0

      end if ! heating depth above min and not at the pole

    end do

    !-----------------------------------------------------------------------
    ! End loop over all columns.
    !-----------------------------------------------------------------------

    ! Output the source level.
    src_level = topi
    tend_level = topi

    ! Set phase speeds; just use reference speeds.
    c = spread(band_mid%cref, 1, ncol)

  end subroutine gw_beres_src

end module gravity_wave_drag_convection
