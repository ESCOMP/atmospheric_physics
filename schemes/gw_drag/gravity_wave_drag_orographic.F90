! Stationary gravity waves from orographic sources
! Extracted: May 2013
! CCPPized:  August 2025
module gravity_wave_drag_orographic

  use ccpp_kinds, only: kind_phys
  use gw_common, only: GWBand

  implicit none
  private
  save

  public :: gravity_wave_drag_orographic_init
  public :: gravity_wave_drag_orographic_run

  type(GWBand) :: band_oro

  ! Limiters (min/max values)
  ! min surface displacement height for orographic waves
  real(kind_phys), parameter :: orohmin = 10._kind_phys
  ! min wind speed for orographic waves
  real(kind_phys), parameter :: orovmin = 2._kind_phys

contains

  subroutine gravity_wave_drag_orographic_init(&
             gw_delta_c, &
             effgw_oro, &
             errmsg, errflg)

    use gw_common, only: wavelength_mid, unset_kind_phys

    real(kind_phys),    intent(in)                :: gw_delta_c
    real(kind_phys),    intent(in)                :: effgw_oro
    character(len=512), intent(out)               :: errmsg
    integer, intent(out)                          :: errflg

    errmsg = ''
    errflg = 0

    if(effgw_oro == unset_kind_phys) then
      errflg = 1
      errmsg = "gravity_wave_drag_orographic_init: orographic gravity wave drag is on but effgw_oro is unset."
      return
    endif

    ! Band to emit orographic waves in.
    ! Regardless, we will only ever emit into l = 0.
    band_oro    = GWBand(0,         gw_delta_c,      1.0_kind_phys, wavelength_mid)

  end subroutine gravity_wave_drag_orographic_init

  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  subroutine gravity_wave_drag_orographic_run( &
             ncol, pver, pcnst, &
             dt, &
             rair, &
             cpairv, &
             p, vramp, &
             piln, rhoi, nm, ni, &
             effgw_oro, gw_lndscl_sgh, gw_oro_south_fac, gw_apply_tndmax, &
             landfrac, &
             lat, &
             u, v, t, q, dse, &
             sgh, zm, &
             kvt_gw, &
             tend_q, tend_u, tend_v, tend_s, &
             src_level, tend_level, ubm, ubi, xv, yv, &
             utgw, vtgw, ttgw, qtgw, &
             dttdf, dttke, egwdffi_tot, &
             flx_heat, &
             tau0x, tau0y, taua, &
             errmsg, errflg)

    use coords_1d, only: Coords1D
    use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
    use gw_common, only: gw_drag_prof, energy_change

    integer,            intent(in)                :: ncol
    integer,            intent(in)                :: pver
    integer,            intent(in)                :: pcnst
    real(kind_phys),    intent(in)                :: dt                       ! Physics timestep [s]
    real(kind_phys),    intent(in)                :: rair                     ! Gas constant for dry air [J kg-1 K-1]
    real(kind_phys),    intent(in)                :: cpairv(:,:)              ! Specific heat of dry air at constant pressure [J kg-1 K-1]
    type(Coords1D),     intent(in)                :: p                        ! Pressure coordinates [Pa]
    real(kind_phys),    pointer, intent(in)       :: vramp(:)                 ! Gravity wave drag tapering coefficients [1]
    real(kind_phys),    intent(in)                :: piln(:, :)               ! Natural log of pressure at interfaces [ln(Pa)]
    real(kind_phys),    intent(in)                :: rhoi(:, :)               ! Density at interfaces [kg m-3]
    real(kind_phys),    intent(in)                :: nm(:, :)                 ! Brunt-Vaisalla frequency at midpoints [s-1]
    real(kind_phys),    intent(in)                :: ni(:, :)                 ! Brunt-Vaisalla frequency at interfaces [s-1]
    real(kind_phys),    intent(in)                :: effgw_oro                ! Efficiency factor for orographic gravity waves [1]
    logical,            intent(in)                :: gw_lndscl_sgh            ! Scale SGH by land fraction? [flag]
    real(kind_phys),    intent(in)                :: gw_oro_south_fac         ! Southern hemisphere scaling factor for orographic waves [1]
    logical,            intent(in)                :: gw_apply_tndmax          ! Maximum wind tendency limiter from stress divergence [flag]
    real(kind_phys),    intent(in)                :: landfrac(:)              ! Land fraction [fraction]
    real(kind_phys),    intent(in)                :: lat(:)                   ! Latitude [rad]
    real(kind_phys),    intent(in)                :: u(:,:)                   ! Zonal wind at midpoints [m s-1]
    real(kind_phys),    intent(in)                :: v(:,:)                   ! Meridional wind at midpoints [m s-1]
    real(kind_phys),    intent(in)                :: t(:,:)                   ! Temperature at midpoints [K]
    real(kind_phys),    intent(in)                :: q(:,:,:)                 ! Constituent mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)                :: dse(:,:)                 ! Dry static energy [J kg-1]
    real(kind_phys),    intent(in)                :: sgh(:)                   ! Standard deviation of orography [m]
    real(kind_phys),    intent(in)                :: zm(:,:)                  ! Geopotential height at midpoints [m]
    real(kind_phys),    intent(in)                :: kvt_gw(:,:)              ! Eddy diffusion coefficient for heat [m2 s-1]

    real(kind_phys),    intent(inout)             :: tend_q(:, :, :)          ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys),    intent(inout)             :: tend_u(:, :)             ! Zonal wind tendency [m s-2]
    real(kind_phys),    intent(inout)             :: tend_v(:, :)             ! Meridional wind tendency [m s-2]
    real(kind_phys),    intent(inout)             :: tend_s(:, :)             ! Dry static energy tendency [J kg-1 s-1]

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
    real(kind_phys),    intent(out)               :: egwdffi_tot(:, :)        ! Total eddy diffusion coefficient from gravity waves [m2 s-1]
    real(kind_phys),    intent(out)               :: flx_heat(:)              ! Surface heat flux for energy conservation check [W m-2]
    real(kind_phys),    intent(out)               :: tau0x(:)                 ! Zonal gravity wave surface stress [N m-2]
    real(kind_phys),    intent(out)               :: tau0y(:)                 ! Meridional gravity wave surface stress [N m-2]
    real(kind_phys),    intent(out)               :: taua(:,:)                ! Total stress from orographic scheme [N m-2]
    character(len=512), intent(out)               :: errmsg
    integer,            intent(out)               :: errflg

    integer :: i, k, m

    ! Wavenumber fields
    real(kind_phys) :: tau(ncol, -band_oro%ngwv:band_oro%ngwv, pver+1)
    real(kind_phys) :: gwut(ncol, pver, -band_oro%ngwv:band_oro%ngwv)
    real(kind_phys) :: phase_speeds(ncol, -band_oro%ngwv:band_oro%ngwv)

    ! Scaled sgh when gw_lndscl_sgh
    real(kind_phys) :: sgh_local(ncol)

    real(kind_phys) :: effgw(ncol)

    real(kind_phys) :: egwdffi(ncol, pver+1)

    ! Surface streamline displacement height (2*sgh).
    real(kind_phys) :: hdsp(ncol)

    ! Max orographic standard deviation to use.
    real(kind_phys) :: sghmax

    ! c=0 stress from orography.
    real(kind_phys) :: tauoro(ncol)

    ! Averages over source region.
    real(kind_phys) :: nsrc(ncol) ! B-V frequency.
    real(kind_phys) :: rsrc(ncol) ! Density.
    real(kind_phys) :: usrc(ncol) ! Zonal wind.
    real(kind_phys) :: vsrc(ncol) ! Meridional wind.

    ! Difference in interface pressure across source region.
    real(kind_phys) :: dpsrc(ncol)

    ! Energy change used by fixer.
    real(kind_phys) :: de(ncol)

    ! Scale sgh by land fraction?
    if (gw_lndscl_sgh) then
      where (landfrac(:ncol) >= epsilon(1._kind_phys))
        effgw = effgw_oro * landfrac(:ncol)
        sgh_local = sgh(:ncol) / sqrt(landfrac(:ncol))
      elsewhere
        effgw = 0._kind_phys
        sgh_local = 0._kind_phys
      end where
    else
      effgw(:ncol) = effgw_oro
      sgh_local(:ncol) = sgh(:ncol)
    endif

    !--------------------------------------------------------------------------
    ! Average the basic state variables for the wave source over the depth of
    ! the orographic standard deviation. Here we assume that the appropiate
    ! values of wind, stability, etc. for determining the wave source are
    ! averages over the depth of the atmosphere penterated by the typical
    ! mountain.
    ! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
    !--------------------------------------------------------------------------
    hdsp(:) = 2.0_kind_phys*sgh_local(:)

    k = pver
    src_level = k - 1
    rsrc = p%mid(:, k)/(rair*t(:, k))*p%del(:, k)
    usrc = u(:, k)*p%del(:, k)
    vsrc = v(:, k)*p%del(:, k)
    nsrc = nm(:, k)*p%del(:, k)

    do k = pver - 1, 1, -1
      do i = 1, ncol
        if (hdsp(i) > sqrt(zm(i, k)*zm(i, k + 1))) then
          src_level(i) = k - 1
          rsrc(i) = rsrc(i) + &
                    p%mid(i, k)/(rair*t(i, k))*p%del(i, k)
          usrc(i) = usrc(i) + u(i, k)*p%del(i, k)
          vsrc(i) = vsrc(i) + v(i, k)*p%del(i, k)
          nsrc(i) = nsrc(i) + nm(i, k)*p%del(i, k)
        end if
      end do
      ! Break the loop when all source levels found.
      if (all(src_level >= k)) exit
    end do

    do i = 1, ncol
      dpsrc(i) = p%ifc(i, pver + 1) - p%ifc(i, src_level(i) + 1)
    end do

    rsrc = rsrc/dpsrc
    usrc = usrc/dpsrc
    vsrc = vsrc/dpsrc
    nsrc = nsrc/dpsrc

    ! Get the unit vector components and magnitude at the surface.
    call get_unit_vector(usrc, vsrc, xv, yv, ubi(:, pver + 1))

    ! Project the local wind at midpoints onto the source wind.
    do k = 1, pver
      ubm(:, k) = dot_2d(u(:, k), v(:, k), xv, yv)
    end do

    ! Compute the interface wind projection by averaging the midpoint winds.
    ! Use the top level wind at the top interface.
    ubi(:, 1) = ubm(:, 1)

    ubi(:, 2:pver) = midpoint_interp(ubm)

    ! Determine the orographic c=0 source term following McFarlane (1987).
    ! Set the source top interface index to pver, if the orographic term is
    ! zero.
    do i = 1, ncol
      if ((ubi(i, pver + 1) > orovmin) .and. (hdsp(i) > orohmin)) then
        sghmax = band_oro%fcrit2*(ubi(i, pver + 1)/nsrc(i))**2
        tauoro(i) = 0.5_kind_phys*band_oro%kwv*min(hdsp(i)**2, sghmax)* &
                    rsrc(i)*nsrc(i)*ubi(i, pver + 1)
      else
        tauoro(i) = 0._kind_phys
        src_level(i) = pver
      end if
    end do

    ! Set the phase speeds and wave numbers in the direction of the source
    ! wind. Set the source stress magnitude (positive only, note that the
    ! sign of the stress is the same as (c-u).
    tau = 0._kind_phys
    do k = pver, minval(src_level), -1
      where (src_level <= k) tau(:, 0, k + 1) = tauoro
    end do

    ! Allow wind tendencies all the way to the model bottom.
    tend_level = pver

    ! No spectrum; phase speed is just 0.
    phase_speeds(:,:) = 0._kind_phys

    ! Scale tau in the southern hemisphere.
    do i = 1, ncol
      if (lat(i) < 0._kind_phys) then
        tau(i, :, :) = tau(i, :, :)*gw_oro_south_fac
      end if
    end do

    ! Solve for the drag profile with orographic sources.
    call gw_drag_prof(ncol, band_oro, p, src_level, tend_level, dt, &
                      t, vramp, &
                      piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                      effgw, phase_speeds, kvt_gw, q, dse, tau, utgw, vtgw, &
                      ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                      lapply_effgw_in=gw_apply_tndmax)

    ! For orographic waves, don't bother with taucd, since there are no
    ! momentum conservation routines or directional diagnostics.

    ! Add the diffusion coefficients
    do k = 1, pver + 1
      egwdffi_tot(:, k) = egwdffi_tot(:, k) + egwdffi(:, k)
    end do

    ! Add the orographic tendencies to the spectrum tendencies.
    ! Don't calculate fixers, since we are too close to the ground to
    ! spread momentum/energy differences across low layers.
    do k = 1, pver
      tend_u(:ncol, k) = tend_u(:ncol, k) + utgw(:, k)
      tend_v(:ncol, k) = tend_v(:ncol, k) + vtgw(:, k)
      tend_s(:ncol, k) = tend_s(:ncol, k) + ttgw(:, k)
      ! Convert to temperature tendency for output.
      ttgw(:, k) = ttgw(:, k)/cpairv(:ncol, k)
    end do

    ! Calculate energy change for output to CAM's energy checker.
    ! This is sort of cheating; we don't have a good a priori idea of the
    ! energy coming from surface stress, so we just integrate what we and
    ! actually have so far and OVERWRITE flx_heat with that.
    call energy_change(dt, p, u, v, tend_u(:ncol, :), &
                       tend_v(:ncol, :), tend_s(:ncol, :), de)
    flx_heat(:ncol) = de

    ! Store constituents tendencies
    do m = 1, pcnst
      do k = 1, pver
         tend_q(:ncol,k,m) = tend_q(:ncol,k,m) + qtgw(:ncol,k,m)
      end do
    end do

    ! Compute surface stresses for diagnostic
    tau0x(:ncol) = tau(:ncol,0,pver+1) * xv(:ncol)
    tau0y(:ncol) = tau(:ncol,0,pver+1) * yv(:ncol)
    taua(:ncol, :pver+1) = tau(:ncol,0,:pver+1)

  end subroutine gravity_wave_drag_orographic_run

end module gravity_wave_drag_orographic
