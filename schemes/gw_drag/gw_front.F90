! This module handles gravity waves from frontal sources, and was extracted
! from gw_drag in May 2013.
module gw_front
  use ccpp_kinds, only: kind_phys
  use gw_common, only: GWBand, unset_kind_phys

  implicit none
  private
  save

  public :: gravity_wave_frontogenesis_init
  public :: gravity_wave_frontogenesis_run
  public :: gravity_wave_frontogenesis_inertial_run

  ! Tuneable settings.
  type CMSourceDesc
    ! Source level.
    integer :: ksrc
    ! Level at which to check whether the frontogenesis function is above
    ! the critical threshold.
    integer :: kfront
    ! Frontogenesis function critical threshold.
    real(kind_phys) :: frontgfc = huge(1._kind_phys)
    ! The stress spectrum to produce at the source level.
    real(kind_phys), allocatable :: src_tau(:)
  end type CMSourceDesc

  ! Frontogenesis wave settings.
  type(CMSourceDesc) :: cm_desc
  type(CMSourceDesc) :: cm_igw_desc

  type(GWBand)       :: band_mid
  type(GWBand)       :: band_long

  ! Parameters for the IGW polar taper.
  real(kind_phys) :: degree2radian
  real(kind_phys) :: al0
  real(kind_phys) :: dlat0

contains

  subroutine gravity_wave_frontogenesis_init(&
             pver, pi, &
             masterproc, iulog, &
             pref_edge, frontgfc, &
             gw_delta_c, pgwv, pgwv_long, &
             taubgnd, taubgnd_igw, &
             effgw_cm, effgw_cm_igw, use_gw_front, use_gw_front_igw, &
             front_gaussian_width, &
             errmsg, errflg)

    use gw_common, only: wavelength_mid, wavelength_long

    integer,            intent(in)                :: pver
    real(kind_phys),    intent(in)                :: pi                       ! pi_constant [1]
    logical,            intent(in)                :: masterproc
    integer,            intent(in)                :: iulog
    real(kind_phys),    intent(in)                :: pref_edge(:)             ! Reference pressure at interfaces [Pa]
    real(kind_phys),    intent(in)                :: frontgfc                 ! Frontogenesis function critical threshold [1]
    real(kind_phys),    intent(in)                :: gw_delta_c               ! Gravity wave phase speed interval [m s-1]
    integer,            intent(in)                :: pgwv                     ! Gravity wave spectrum dimension (wave numbers are from -pgwv to pgwv).
    integer,            intent(in)                :: pgwv_long
    real(kind_phys),    intent(in)                :: taubgnd                  ! Background source strength (used for waves from frontogenesis). waves [N m-2]
    real(kind_phys),    intent(in)                :: taubgnd_igw              ! Background source strength (used for inertial waves from frontogenesis). [N m-2]
    real(kind_phys),    intent(in)                :: effgw_cm                 ! Efficiency associated with gravity waves from frontogenesis. [1]
    real(kind_phys),    intent(in)                :: effgw_cm_igw             ! Efficiency associated with inertial gravity waves from frontogenesis. [1]
    logical,            intent(in)                :: use_gw_front             ! Use frontogenesis gravity waves [flag]
    logical,            intent(in)                :: use_gw_front_igw         ! Use inertial frontogenesis gravity waves [flag]
    real(kind_phys),    intent(in)                :: front_gaussian_width     ! Width of gaussian used to create frontogenesis tau profile [m s-1]
    character(len=512), intent(out)               :: errmsg
    integer,            intent(out)               :: errflg

    integer :: k, l

    ! Bottom level for frontal waves.
    integer :: kbot_front

    ! Index for levels at specific pressures.
    integer :: kfront

    character(len=*), parameter :: sub = 'gravity_wave_frontogenesis_init'

    ! Initialize error variables
    errmsg = ''
    errflg = 0

    if (use_gw_front .or. use_gw_front_igw) then

      ! Check that deep gw file is set in namelist
      if (frontgfc == unset_kind_phys) then
        errmsg = "Frontogenesis enabled, but frontgfc was not set!"
        errflg = 1
        return
      end if

      do k = 0, pver
        ! Check frontogenesis at 600 hPa.
        if (pref_edge(k + 1) < 60000._kind_phys) kfront = k + 1
      end do

      ! Source waves from 500 hPa.
      kbot_front = maxloc(pref_edge, 1, (pref_edge < 50000._kind_phys)) - 1

      if (masterproc) then
        write (iulog, *) 'KFRONT      =', kfront
        write (iulog, *) 'KBOT_FRONT  =', kbot_front
        write (iulog, *) ' '
      end if
    end if

    if (use_gw_front) then
      band_mid = GWBand(pgwv, gw_delta_c, 1.0_kind_phys, wavelength_mid)
      if (masterproc) then
        write (iulog, *) ' '
        write (iulog, *) "gravity_wave_frontogenesis_init: band_mid%ngwv = ", band_mid%ngwv
        do l = -band_mid%ngwv, band_mid%ngwv
          write (iulog, '(A,I0,A,F7.2)') &
            "gravity_wave_frontogenesis_init: band_mid%cref(", l, ") = ", band_mid%cref(l)
        end do
        write (iulog, *) 'gravity_wave_frontogenesis_init: band_mid%kwv = ', band_mid%kwv
        write (iulog, *) 'gravity_wave_frontogenesis_init: band_mid%fcrit2 = ', band_mid%fcrit2
      end if
    end if

    if (use_gw_front_igw) then
      band_long = GWBand(pgwv_long, gw_delta_c, 1.0_kind_phys, wavelength_long)
      if (masterproc) then
        write (iulog, *) ' '
        write (iulog, *) "gravity_wave_frontogenesis_init: band_long%ngwv = ", band_long%ngwv
        do l = -band_long%ngwv, band_long%ngwv
          write (iulog, '(A,I2,A,F7.2)') &
            "gravity_wave_frontogenesis_init: band_long%cref(", l, ") = ", band_long%cref(l)
        end do
        write (iulog, *) 'gravity_wave_frontogenesis_init: band_long%kwv = ', band_long%kwv
        write (iulog, *) 'gravity_wave_frontogenesis_init: band_long%fcrit2 = ', band_long%fcrit2
        write (iulog, *) ' '
      end if
    end if

    if (use_gw_front) then
      ! Check that deep gw file is set in namelist
      if (taubgnd == unset_kind_phys .or. effgw_cm == unset_kind_phys) then
        errmsg = "Frontogenesis mid-scale waves enabled, but not all required namelist variables were set!"
        errflg = 1
        return
      end if

      if (masterproc) then
        write (iulog, *) 'gw_init: gw spectrum taubgnd, effgw_cm = ', taubgnd, effgw_cm
      end if

      cm_desc = gaussian_cm_desc(band_mid, pi, kbot_front, kfront, frontgfc, &
                                 taubgnd, front_gaussian_width)

    end if

    if (use_gw_front_igw) then
      ! Check that deep gw file is set in namelist
      if (effgw_cm_igw == unset_kind_phys .or. taubgnd_igw == unset_kind_phys) then
        write (errmsg, '(a, a)') sub, &
          " Frontogenesis inertial waves enabled, but not all required namelist variables were set!"
        errflg = 1
        return
      end if

      if (masterproc) then
        write (iulog, *) 'gw_init: gw spectrum taubgnd_igw, effgw_cm_igw = ', taubgnd_igw, effgw_cm_igw
      end if

      cm_igw_desc = gaussian_cm_desc(band_long, pi, kbot_front, kfront, frontgfc, &
                                     taubgnd_igw, front_gaussian_width)

      ! Parameters for the IGW polar taper.
      degree2radian = pi/180._kind_phys
      al0 = 82.5_kind_phys * degree2radian
      dlat0 = 5.0_kind_phys * degree2radian

    end if
  end subroutine gravity_wave_frontogenesis_init

  ! Frontally generated gravity waves
!> \section arg_table_gravity_wave_frontogenesis_run Argument Table
!! \htmlinclude gravity_wave_frontogenesis_run.html
  subroutine gravity_wave_frontogenesis_run( &
             ncol, pver, pcnst, &
             dt, &
             cpair, &
             p, vramp, &
             piln, rhoi, nm, ni, &
             effgw_cm, gw_polar_taper, gw_apply_tndmax, &
             lat, &
             u, v, t, q, dse, &
             frontgf, &
             kvt_gw, &
             tend_q, tend_u, tend_v, tend_s, egwdffi_tot, &
             src_level, tend_level, ubm, ubi, xv, yv, &
             utgw, vtgw, ttgw, qtgw, &
             dttdf, dttke, &
             taucd_west, taucd_east, taucd_south, taucd_north, &
             utend1, utend2, utend3, utend4, utend5, &
             flx_heat, &
             errmsg, errflg)

    use coords_1d, only: Coords1D
    use gw_common, only: energy_change, energy_fixer
    use gw_common, only: momentum_flux, momentum_fixer
    use gw_common, only: gw_drag_prof
    use gw_common, only: calc_taucd
    use gw_common, only: west, east, south, north
    use gw_common, only: find_bin

    integer,            intent(in)                :: ncol
    integer,            intent(in)                :: pver
    integer,            intent(in)                :: pcnst
    real(kind_phys),    intent(in)                :: dt                       ! Physics timestep [s]
    real(kind_phys),    intent(in)                :: cpair                    ! Specific heat of dry air at constant pressure [J kg-1 K-1]
    type(Coords1D),     intent(in)                :: p                        ! Pressure coordinates [Pa]
    real(kind_phys),    pointer, intent(in)       :: vramp(:)                 ! Gravity wave drag tapering coefficients [1]
    real(kind_phys),    intent(in)                :: piln(:, :)               ! Natural log of pressure at interfaces [ln(Pa)]
    real(kind_phys),    intent(in)                :: rhoi(:, :)               ! Density at interfaces [kg m-3]
    real(kind_phys),    intent(in)                :: nm(:, :)                 ! Brunt-Vaisalla frequency at midpoints [s-1]
    real(kind_phys),    intent(in)                :: ni(:, :)                 ! Brunt-Vaisalla frequency at interfaces [s-1]
    real(kind_phys),    intent(in)                :: effgw_cm                 ! Efficiency factor for frontogenesis gravity waves [1]
    logical,            intent(in)                :: gw_polar_taper           ! Apply polar tapering to frontogenesis efficiency [flag]
    logical,            intent(in)                :: gw_apply_tndmax          ! Maximum wind tendency limiter from stress divergence [flag]
    real(kind_phys),    intent(in)                :: lat(:)                   ! Latitude [rad]
    real(kind_phys),    intent(in)                :: u(:,:)                   ! Zonal wind at midpoints [m s-1]
    real(kind_phys),    intent(in)                :: v(:,:)                   ! Meridional wind at midpoints [m s-1]
    real(kind_phys),    intent(in)                :: t(:,:)                   ! Temperature at midpoints [K]
    real(kind_phys),    intent(in)                :: q(:,:,:)                 ! Constituent mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)                :: dse(:,:)                 ! Dry static energy [J kg-1]
    real(kind_phys),    intent(in)                :: frontgf(:, :)            ! Frontogenesis function [1]
    real(kind_phys),    intent(in)                :: kvt_gw(:,:)              ! Eddy diffusion coefficient for heat [m2 s-1]

    real(kind_phys),    intent(inout)             :: tend_q(:, :, :)          ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys),    intent(inout)             :: tend_u(:, :)             ! Zonal wind tendency [m s-2]
    real(kind_phys),    intent(inout)             :: tend_v(:, :)             ! Meridional wind tendency [m s-2]
    real(kind_phys),    intent(inout)             :: tend_s(:, :)             ! Dry static energy tendency [J kg-1 s-1]
    real(kind_phys),    intent(inout)             :: egwdffi_tot(:, :)        ! Total eddy diffusion coefficient from gravity waves [m2 s-1]

    integer,            intent(out)               :: src_level(:)             ! Vertical level index of gravity wave source [index]
    integer,            intent(out)               :: tend_level(:)            ! Lowest vertical level index where tendencies are applied [index]
    real(kind_phys),    intent(out)               :: ubm(:, :)                ! Wind projection at midpoints along source wind direction [m s-1]
    real(kind_phys),    intent(out)               :: ubi(:, :)                ! Wind projection at interfaces along source wind direction [m s-1]
    real(kind_phys),    intent(out)               :: xv(:)                    ! Zonal component of source wind unit vector [1]
    real(kind_phys),    intent(out)               :: yv(:)                    ! Meridional component of source wind unit vector [1]
    real(kind_phys),    intent(out)               :: utgw(:, :)               ! Zonal wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)               :: vtgw(:, :)               ! Meridional wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)               :: ttgw(:, :)               ! Temperature tendency from gravity waves [K s-1]
    real(kind_phys),    intent(out)               :: qtgw(:, :, :)            ! Constituent tendencies from gravity waves [kg kg-1 s-1]
    real(kind_phys),    intent(out)               :: dttdf(:, :)              ! Temperature tendency from diffusion [K s-1]
    real(kind_phys),    intent(out)               :: dttke(:, :)              ! Temperature tendency from kinetic energy dissipation [K s-1]
    real(kind_phys),    intent(out)               :: flx_heat(:)              ! Surface heat flux for energy conservation check [W m-2]

    ! Copies of taucd in each direction for diagnostic.
    real(kind_phys),    intent(out)               :: taucd_west(:, :)         ! Reynolds stress for waves in W direction, interfaces [N m-2]
    real(kind_phys),    intent(out)               :: taucd_east(:, :)         ! Reynolds stress for waves in E direction, interfaces [N m-2]
    real(kind_phys),    intent(out)               :: taucd_south(:, :)        ! Reynolds stress for waves in S direction, interfaces [N m-2]
    real(kind_phys),    intent(out)               :: taucd_north(:, :)        ! Reynolds stress for waves in N direction, interfaces [N m-2]

    ! Wind tendencies broken across five spectral bins for diagnostic.
    real(kind_phys),    intent(out)               :: utend1(:, :)             ! U tendency c < -40 [m s-2]
    real(kind_phys),    intent(out)               :: utend2(:, :)             ! U tendency -40 < c < -15 [m s-2]
    real(kind_phys),    intent(out)               :: utend3(:, :)             ! U tendency -15 < c < 15 [m s-2]
    real(kind_phys),    intent(out)               :: utend4(:, :)             ! U tendency 15 < c < 40 [m s-2]
    real(kind_phys),    intent(out)               :: utend5(:, :)             ! U tendency c > 40 [m s-2]

    character(len=512), intent(out)               :: errmsg
    integer,            intent(out)               :: errflg

    ! Local variables
    integer :: i, k, m, l, stat

    ! Wavenumber fields
    real(kind_phys) :: tau(ncol, -band_mid%ngwv:band_mid%ngwv, pver + 1)
    real(kind_phys) :: gwut(ncol, pver, -band_mid%ngwv:band_mid%ngwv)
    real(kind_phys) :: phase_speeds(ncol, -band_mid%ngwv:band_mid%ngwv)

    real(kind_phys) :: utend(ncol, pver, 5)
    real(kind_phys) :: ix(ncol, -band_mid%ngwv:band_mid%ngwv)
    real(kind_phys) :: iy(ncol, -band_mid%ngwv:band_mid%ngwv)

    real(kind_phys) :: egwdffi(ncol, pver+1)
    real(kind_phys) :: effgw(ncol)
    real(kind_phys) :: taucd(ncol, pver+1, 4)
    real(kind_phys) :: um_flux(ncol)
    real(kind_phys) :: vm_flux(ncol)
    real(kind_phys) :: de(ncol)

    errmsg = ''
    errflg = 0

    ! Efficiency of gravity wave momentum transfer
    effgw = effgw_cm

    ! Frontogenesis is too high at the poles (at least for the FV dycore), so introduce a polar taper
    if (gw_polar_taper) then
      effgw = effgw*cos(lat(:ncol))
    endif

    call gw_cm_src(ncol, pver, band_mid, &
                   cm_desc, &
                   u, v, frontgf(:ncol, :), &
                   src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds)

    ! Solve for the drag profile with C&M source spectrum
    call gw_drag_prof(ncol, band_mid, p, src_level, tend_level, dt, &
                      t, vramp, &
                      piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                      effgw, phase_speeds, kvt_gw, q, dse, tau, utgw, vtgw, &
                      ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                      lapply_effgw_in=gw_apply_tndmax)

    ! Project stress into directional components
    taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

    ! Make copies for diagnostics.
    taucd_west(:,:)  = taucd(:,:,west)
    taucd_east(:,:)  = taucd(:,:,east)
    taucd_south(:,:) = taucd(:,:,south)
    taucd_north(:,:) = taucd(:,:,north)

    ! Add the diffusion coefficients
    do k = 1, pver + 1
      egwdffi_tot(:, k) = egwdffi_tot(:, k) + egwdffi(:, k)
    end do

    ! Add the constituent tendencies
    do m = 1, pcnst
      do k = 1, pver
        tend_q(:ncol, k, m) = tend_q(:ncol, k, m) + qtgw(:, k, m)
      end do
    end do

    ! Find momentum flux, and use it to fix the wind tendencies below
    ! the gravity wave region
    call momentum_flux(tend_level, taucd, um_flux, vm_flux)
    call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

    ! Add the momentum tendencies to the output tendency arrays
    do k = 1, pver
      tend_u(:ncol, k) = tend_u(:ncol, k) + utgw(:, k)
      tend_v(:ncol, k) = tend_v(:ncol, k) + vtgw(:, k)
    end do

    ! Diagnostic: accumulate wind tendencies binned according to phase speed.
    ix = find_bin(phase_speeds)
    do l = -band_mid%ngwv, band_mid%ngwv
      do k = 1, pver
        do i = 1, ncol
          utend(i, k, ix(i, l)) = utend(i, k, ix(i, l)) + gwut(i, k, l)
        end do
      end do
    end do

    ! Get the zonal part
    do l = 1, 5
      do k = 1, pver
        utend(:ncol, k, l) = utend(:ncol, k, l) * xv(:ncol)
      enddo
    enddo

    ! Copy into sub-tendency arrays for export
    utend1(:ncol,:pver) = utend(:ncol, :pver, 1)
    utend2(:ncol,:pver) = utend(:ncol, :pver, 2)
    utend3(:ncol,:pver) = utend(:ncol, :pver, 3)
    utend4(:ncol,:pver) = utend(:ncol, :pver, 4)
    utend5(:ncol,:pver) = utend(:ncol, :pver, 5)

    ! Find energy change in the current state, and use fixer to apply
    ! the difference in lower levels
    call energy_change(dt, p, u, v, tend_u(:ncol, :), &
                       tend_v(:ncol, :), tend_s(:ncol, :) + ttgw, de)
    call energy_fixer(tend_level, p, de - flx_heat(:ncol), ttgw)

    do k = 1, pver
      tend_s(:ncol, k) = tend_s(:ncol, k) + ttgw(:, k)
    end do

    ! Convert dry static energy tendency to temperature tendency for output
    ! FIXME: some places use cpairv (e.g., orographic) but cpair is used here. hplin 8/14/25
    ttgw = ttgw / cpair

  end subroutine gravity_wave_frontogenesis_run

  ! Frontaly generated intertial gravity waves
!> \section arg_table_gravity_wave_frontogenesis_inertial_run Argument Table
!! \htmlinclude gravity_wave_frontogenesis_inertial_run.html
  subroutine gravity_wave_frontogenesis_inertial_run( &
             ncol, pver, pcnst, &
             dt, &
             cpair, &
             p, vramp, &
             piln, rhoi, nm, ni, &
             effgw_cm_igw, gw_polar_taper, gw_apply_tndmax, &
             lat, &
             u, v, t, q, dse, &
             frontgf, &
             kvt_gw, &
             tend_q, tend_u, tend_v, tend_s, egwdffi_tot, &
             src_level, tend_level, ubm, ubi, xv, yv, &
             utgw, vtgw, ttgw, qtgw, &
             dttdf, dttke, &
             flx_heat, &
             errmsg, errflg)

    use coords_1d, only: Coords1D
    use gw_common, only: energy_change, energy_fixer
    use gw_common, only: momentum_flux, momentum_fixer
    use gw_common, only: gw_drag_prof
    use gw_common, only: calc_taucd
    use gw_common, only: coriolis_speed, adjust_inertial

    integer,            intent(in)                :: ncol
    integer,            intent(in)                :: pver
    integer,            intent(in)                :: pcnst
    real(kind_phys),    intent(in)                :: dt                       ! Physics timestep [s]
    real(kind_phys),    intent(in)                :: cpair                    ! Specific heat of dry air at constant pressure [J kg-1 K-1]
    type(Coords1D),     intent(in)                :: p                        ! Pressure coordinates [Pa]
    real(kind_phys),    pointer, intent(in)       :: vramp(:)                 ! Gravity wave drag tapering coefficients [1]
    real(kind_phys),    intent(in)                :: piln(:, :)               ! Natural log of pressure at interfaces [ln(Pa)]
    real(kind_phys),    intent(in)                :: rhoi(:, :)               ! Density at interfaces [kg m-3]
    real(kind_phys),    intent(in)                :: nm(:, :)                 ! Brunt-Vaisalla frequency at midpoints [s-1]
    real(kind_phys),    intent(in)                :: ni(:, :)                 ! Brunt-Vaisalla frequency at interfaces [s-1]
    real(kind_phys),    intent(in)                :: effgw_cm_igw             ! Efficiency factor for inertial frontogenesis gravity waves [1]
    logical,            intent(in)                :: gw_polar_taper           ! Apply polar tapering to frontogenesis efficiency [flag]
    logical,            intent(in)                :: gw_apply_tndmax          ! Maximum wind tendency limiter from stress divergence [flag]

    real(kind_phys),    intent(in)                :: lat(:)                   ! Latitude [rad]
    real(kind_phys),    intent(in)                :: u(:,:)                   ! Zonal wind at midpoints [m s-1]
    real(kind_phys),    intent(in)                :: v(:,:)                   ! Meridional wind at midpoints [m s-1]
    real(kind_phys),    intent(in)                :: t(:,:)                   ! Temperature at midpoints [K]
    real(kind_phys),    intent(in)                :: q(:,:,:)                 ! Constituent mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)                :: dse(:,:)                 ! Dry static energy [J kg-1]
    real(kind_phys),    intent(in)                :: frontgf(:, :)            ! Frontogenesis function [1]
    real(kind_phys),    intent(in)                :: kvt_gw(:,:)              ! Eddy diffusion coefficient for heat [m2 s-1]

    real(kind_phys),    intent(inout)             :: tend_q(:, :, :)          ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys),    intent(inout)             :: tend_u(:, :)             ! Zonal wind tendency [m s-2]
    real(kind_phys),    intent(inout)             :: tend_v(:, :)             ! Meridional wind tendency [m s-2]
    real(kind_phys),    intent(inout)             :: tend_s(:, :)             ! Dry static energy tendency [J kg-1 s-1]
    real(kind_phys),    intent(inout)             :: egwdffi_tot(:, :)        ! Total eddy diffusion coefficient from gravity waves [m2 s-1]

    integer,            intent(out)               :: src_level(:)             ! Vertical level index of gravity wave source [index]
    integer,            intent(out)               :: tend_level(:)            ! Lowest vertical level index where tendencies are applied [index]
    real(kind_phys),    intent(out)               :: ubm(:, :)                ! Wind projection at midpoints along source wind direction [m s-1]
    real(kind_phys),    intent(out)               :: ubi(:, :)                ! Wind projection at interfaces along source wind direction [m s-1]
    real(kind_phys),    intent(out)               :: xv(:)                    ! Zonal component of source wind unit vector [1]
    real(kind_phys),    intent(out)               :: yv(:)                    ! Meridional component of source wind unit vector [1]
    real(kind_phys),    intent(out)               :: utgw(:, :)               ! Zonal wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)               :: vtgw(:, :)               ! Meridional wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)               :: ttgw(:, :)               ! Temperature tendency from gravity waves [K s-1]
    real(kind_phys),    intent(out)               :: qtgw(:, :, :)            ! Constituent tendencies from gravity waves [kg kg-1 s-1]
    real(kind_phys),    intent(out)               :: dttdf(:, :)              ! Temperature tendency from diffusion [K s-1]
    real(kind_phys),    intent(out)               :: dttke(:, :)              ! Temperature tendency from kinetic energy dissipation [K s-1]
    real(kind_phys),    intent(out)               :: flx_heat(:)              ! Surface heat flux for energy conservation check [W m-2]
    character(len=512), intent(out)               :: errmsg
    integer,            intent(out)               :: errflg

    ! Local variables
    integer :: k, m, stat

    ! Wavenumber fields
    real(kind_phys) :: tau(ncol, -band_long%ngwv:band_long%ngwv, pver + 1)
    real(kind_phys) :: gwut(ncol, pver, -band_long%ngwv:band_long%ngwv)
    real(kind_phys) :: phase_speeds(ncol, -band_long%ngwv:band_long%ngwv)
    real(kind_phys) :: ro_adjust(ncol, -band_long%ngwv:band_long%ngwv, pver + 1) ! Adjustment for inertial gravity waves.

    ! Coriolis characteristic speed.
    real(kind_phys) :: u_coriolis(ncol)

    real(kind_phys) :: egwdffi(ncol, pver+1)
    real(kind_phys) :: effgw(ncol)                              ! Efficiency factor for gravity wave source.
    real(kind_phys) :: taucd(ncol, pver+1, 4)
    real(kind_phys) :: um_flux(ncol)
    real(kind_phys) :: vm_flux(ncol)
    real(kind_phys) :: de(ncol)
    real(kind_phys) :: al0, dlat0

    character(len=*), parameter :: sub = 'gravity_wave_frontogenesis_inertial_run'

    errmsg = ''
    errflg = 0

    ! Efficiency of gravity wave momentum transfer
    effgw(:) = effgw_cm_igw

    ! Frontogenesis is too high at the poles (at least for the FV dycore), so introduce a polar taper
    if (gw_polar_taper) then
      ! Polar taper parameters
      al0 = 70._kind_phys * degree2radian
      dlat0 = 10._kind_phys * degree2radian

      where (abs(lat(:ncol)) <= 89._kind_phys * degree2radian)
        effgw = effgw * 0.25_kind_phys * &
                (1._kind_phys + tanh((lat(:ncol) + al0) / dlat0)) * &
                (1._kind_phys - tanh((lat(:ncol) - al0) / dlat0))
      elsewhere
        effgw = 0._kind_phys
      end where
    end if

    call gw_cm_src(ncol, pver, band_long, &
                   cm_igw_desc, &
                   u, v, frontgf(:ncol, :), &
                   src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds)

    u_coriolis = coriolis_speed(band_long, lat(:ncol))
    call adjust_inertial(band_long, tend_level, u_coriolis, phase_speeds, ubi, &
                         tau, ro_adjust)

    ! Solve for the drag profile with C&M source spectrum
    call gw_drag_prof(ncol, band_long, p, src_level, tend_level, dt, &
                      t, vramp, &
                      piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                      effgw, phase_speeds, kvt_gw, q, dse, tau, utgw, vtgw, &
                      ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                      ro_adjust=ro_adjust, lapply_effgw_in=gw_apply_tndmax)

    ! Project stress into directional components
    taucd = calc_taucd(ncol, band_long%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

    ! Add the diffusion coefficients
    do k = 1, pver + 1
      egwdffi_tot(:, k) = egwdffi_tot(:, k) + egwdffi(:, k)
    end do

    ! Add the constituent tendencies
    do m = 1, pcnst
      do k = 1, pver
        tend_q(:ncol, k, m) = tend_q(:ncol, k, m) + qtgw(:, k, m)
      end do
    end do

    ! Find momentum flux, and use it to fix the wind tendencies below
    ! the gravity wave region
    call momentum_flux(tend_level, taucd, um_flux, vm_flux)
    call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

    ! Add the momentum tendencies to the output tendency arrays
    do k = 1, pver
      tend_u(:ncol, k) = tend_u(:ncol, k) + utgw(:, k)
      tend_v(:ncol, k) = tend_v(:ncol, k) + vtgw(:, k)
    end do

    ! Find energy change in the current state, and use fixer to apply
    ! the difference in lower levels
    call energy_change(dt, p, u, v, tend_u(:ncol, :), &
                       tend_v(:ncol, :), tend_s(:ncol, :) + ttgw, de)
    call energy_fixer(tend_level, p, de - flx_heat(:ncol), ttgw)

    do k = 1, pver
      tend_s(:ncol, k) = tend_s(:ncol, k) + ttgw(:, k)
    end do

    ! Convert dry static energy tendency to temperature tendency for output
    ! FIXME: some places use cpairv (e.g., orographic) but cpair is used here. hplin 8/14/25
    ttgw = ttgw / cpair

  end subroutine gravity_wave_frontogenesis_inertial_run


!==========================================================================

  ! Create a flat profile to be launched (all wavenumbers have the same
  ! source strength, except that l=0 is excluded).
  function flat_cm_desc(band, ksrc, kfront, frontgfc, taubgnd) result(desc)
    ! Wavelengths triggered by frontogenesis.
    type(GWBand), intent(in) :: band
    ! The following are used to set the corresponding object components.
    integer, intent(in) :: ksrc
    integer, intent(in) :: kfront
    real(kind_phys), intent(in) :: frontgfc
    ! Amount of stress to launch from each wavelength.
    real(kind_phys), intent(in) :: taubgnd

    type(CMSourceDesc) :: desc

    desc%ksrc = ksrc
    desc%kfront = kfront
    desc%frontgfc = frontgfc

    allocate (desc%src_tau(-band%ngwv:band%ngwv))
    desc%src_tau = taubgnd

    ! Prohibit wavenumber 0.
    desc%src_tau(0) = 0._kind_phys

  end function flat_cm_desc

  ! Create a source tau profile that is a gaussian over wavenumbers (l=0 is
  ! excluded).
  function gaussian_cm_desc(band, pi, ksrc, kfront, frontgfc, height, width) &
    result(desc)

    use shr_spfn_mod, only: erfc => shr_spfn_erfc

    ! Wavelengths triggered by frontogenesis.
    type(GWBand), intent(in) :: band
    real(kind_phys), intent(in) :: pi
    ! The following are used to set the corresponding object components.
    integer, intent(in) :: ksrc
    integer, intent(in) :: kfront
    real(kind_phys), intent(in) :: frontgfc
    ! Parameters of gaussian.
    real(kind_phys), intent(in) :: height
    real(kind_phys), intent(in) :: width

    type(CMSourceDesc) :: desc

    ! Bounds used to average bins of the gaussian.
    real(kind_phys) :: gaussian_bounds(2*band%ngwv + 2)

    ! Wavenumber index.
    integer :: l

    desc%ksrc = ksrc
    desc%kfront = kfront
    desc%frontgfc = frontgfc

    allocate (desc%src_tau(-band%ngwv:band%ngwv))

    ! Find the boundaries of each bin.
    gaussian_bounds(:2*band%ngwv + 1) = band%cref - 0.5_kind_phys*band%dc
    gaussian_bounds(2*band%ngwv + 2) = band%cref(band%ngwv) + 0.5_kind_phys*band%dc

    ! Integral of the gaussian at bin interfaces (from the point to
    ! positive infinity).
    gaussian_bounds = &
      [(erfc(gaussian_bounds(l)/width)*height*width*sqrt(pi)/2._kind_phys, &
        l=1, 2*band%ngwv + 2)]

    ! Get average in each bin using integral from right to left side.
    desc%src_tau = &
      (gaussian_bounds(:2*band%ngwv + 1) - gaussian_bounds(2:))/band%dc

    ! Prohibit wavenumber 0.
    desc%src_tau(0) = 0._kind_phys

  end function gaussian_cm_desc

  subroutine gw_cm_src(ncol, pver, band, &
                       desc, &
                       u, v, frontgf, &
                       src_level, tend_level, tau, ubm, ubi, xv, yv, c)
    use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
    !-----------------------------------------------------------------------
    ! Driver for multiple gravity wave drag parameterization.
    !
    ! The parameterization is assumed to operate only where water vapor
    ! concentrations are negligible in determining the density.
    !-----------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    ! Column dimension.
    integer, intent(in) :: ncol
    integer, intent(in) :: pver

    ! Wavelengths triggered by frontogenesis.
    type(CMSourceDesc), intent(in)  :: desc
    type(GWBand), intent(in)        :: band

    ! Midpoint zonal/meridional winds.
    real(kind_phys), intent(in) :: u(ncol, pver), v(ncol, pver)
    ! Frontogenesis function.
    real(kind_phys), intent(in) :: frontgf(:, :)

    ! Indices of top gravity wave source level and lowest level where wind
    ! tendencies are allowed.
    integer, intent(out) :: src_level(ncol)
    integer, intent(out) :: tend_level(ncol)

    ! Wave Reynolds stress.
    real(kind_phys), intent(out) :: tau(ncol, -band%ngwv:band%ngwv, pver + 1)
    ! Projection of wind at midpoints and interfaces.
    real(kind_phys), intent(out) :: ubm(ncol, pver), ubi(ncol, pver + 1)
    ! Unit vectors of source wind (zonal and meridional components).
    real(kind_phys), intent(out) :: xv(ncol), yv(ncol)
    ! Phase speeds.
    real(kind_phys), intent(out) :: c(ncol, -band%ngwv:band%ngwv)

    !---------------------------Local Storage-------------------------------
    ! Column and wavenumber indices.
    integer :: k, l

    ! Whether or not to launch waves in this column.
    logical :: launch_wave(ncol)

    ! Zonal/meridional wind averaged over source region.
    real(kind_phys) :: usrc(ncol), vsrc(ncol)

    !------------------------------------------------------------------------
    ! Determine the source layer wind and unit vectors, then project winds.
    !------------------------------------------------------------------------

    ! Just use the source level interface values for the source wind speed
    ! and direction (unit vector).
    src_level = desc%ksrc
    tend_level = desc%ksrc
    usrc = 0.5_kind_phys*(u(:, desc%ksrc + 1) + u(:, desc%ksrc))
    vsrc = 0.5_kind_phys*(v(:, desc%ksrc + 1) + v(:, desc%ksrc))

    ! Get the unit vector components and magnitude at the surface.
    call get_unit_vector(usrc, vsrc, xv, yv, ubi(:, desc%ksrc + 1))

    ! Project the local wind at midpoints onto the source wind.
    do k = 1, desc%ksrc
      ubm(:, k) = dot_2d(u(:, k), v(:, k), xv, yv)
    end do

    ! Compute the interface wind projection by averaging the midpoint winds.
    ! Use the top level wind at the top interface.
    ubi(:, 1) = ubm(:, 1)

    ubi(:, 2:desc%ksrc) = midpoint_interp(ubm(:, 1:desc%ksrc))

    !-----------------------------------------------------------------------
    ! Gravity wave sources
    !-----------------------------------------------------------------------

    tau = 0._kind_phys

    ! GW generation depends on frontogenesis at specified level (may be below
    ! actual source level).
    launch_wave = (frontgf(:, desc%kfront) > desc%frontgfc)

    do l = -band%ngwv, band%ngwv
      where (launch_wave)
        tau(:, l, desc%ksrc + 1) = desc%src_tau(l)
      end where
    end do

    ! Set phase speeds as reference speeds plus the wind speed at the source
    ! level.
    c = spread(band%cref, 1, ncol) + &
        spread(ubi(:, desc%ksrc + 1), 2, 2*band%ngwv + 1)

  end subroutine gw_cm_src

end module gw_front
