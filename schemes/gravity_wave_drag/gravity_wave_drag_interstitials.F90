! Interstitial schemes for preparing fields needed for various gravity wave drag schemes
module gravity_wave_drag_interstitials
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: gravity_wave_drag_prepare_profiles_run
  public :: gravity_wave_drag_prepare_profiles_timestep_final

contains

  ! (originally gw_prof in gw_common)
  ! Compute profiles of background state quantities for multiple
  ! gravity wave drag parameterizations, and prepare vertical pressure
  ! coordinates.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
!> \section arg_table_gravity_wave_drag_prepare_profiles_run Argument Table
!! \htmlinclude gravity_wave_drag_prepare_profiles_run.html
  subroutine gravity_wave_drag_prepare_profiles_run( &
    ncol, pver, &
    cpair, rair, gravit, &
    pint, &
    t, &
    do_molec_diff, nbot_molec, &
    cpairv, kvt, &
    p, &
    rhoi, &
    nm, ni, &
    egwdffi, &
    tend_q, tend_u, tend_v, tend_s, &
    kvt_gw, &
    flx_heat, &
    scheme_name, &
    errmsg, errflg)

    use gw_utils, only: midpoint_interp
    use coords_1d, only: Coords1D

    ! Input arguments
    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver

    real(kind_phys), intent(in)  :: cpair      ! Specific heat of dry air, constant pressure [J K-1 kg-1]
    real(kind_phys), intent(in)  :: rair       ! Dry air gas constant [J K-1 kg-1]
    real(kind_phys), intent(in)  :: gravit     ! standard_gravitational_acceleration [m s-2]
    real(kind_phys), intent(in)  :: pint(:, :) ! Interface pressures [Pa]
    real(kind_phys), intent(in)  :: t(:, :)    ! Midpoint temperatures [K]

    logical,         intent(in)  :: do_molec_diff ! do_molecular_diffusion [flag]
    integer,         intent(in)  :: nbot_molec ! index_of_pressure_at_bottom_of_molecular_diffusion [index]
    real(kind_phys), intent(in)  :: cpairv(:,:)! composition_dependent_specific_heat_of_dry_air_at_constant_pressure [J kg-1 K-1]
    real(kind_phys), intent(in)  :: kvt(:,:)   ! molecular_kinematic_temperature_conductivity_at_interfaces

    ! Output arguments
    type(coords1d),  intent(out) :: p          ! Vertical coordinate for gravity waves [Coords1D, Pa]

    real(kind_phys), intent(out) :: rhoi(:, :) ! Density at layer interfaces [kg m-3]
    real(kind_phys), intent(out) :: nm(:, :)   ! Midpoint Brunt-Vaisalla frequencies [s-1]
    real(kind_phys), intent(out) :: ni(:, :)   ! Interface Brunt-Vaisalla frequencies [s-1]
    real(kind_phys), intent(out) :: egwdffi(:, :)   ! Effective diffusion coefficient from gravity waves [m2 s-1]

    ! Initialize tendency accumulators to make framework happy about intent(inout) in the individual gw schemes.
    real(kind_phys), intent(out) :: tend_q(:, :, :) ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys), intent(out) :: tend_u(:, :)    ! Zonal wind tendency [m s-2]
    real(kind_phys), intent(out) :: tend_v(:, :)    ! Meridional wind tendency [m s-2]
    real(kind_phys), intent(out) :: tend_s(:, :)    ! Dry static energy tendency [J kg-1 s-1]

    real(kind_phys), intent(out) :: kvt_gw(:,:)! scaled_molecular_kinematic_temperature_conductivity_at_interfaces_for_gravity_wave_drag [m2 s-1]

    ! for energy checker:
    real(kind_phys), intent(out) :: flx_heat(:) ! net_sensible_heat_flux_through_top_and_bottom_of_atmosphere_column [W m-2]
    character(len=*), intent(out):: scheme_name

    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errflg

    ! Local variables
    integer :: i, k
    real(kind_phys) :: dtdp                    ! dt/dp
    real(kind_phys) :: n2                      ! Brunt-Vaisalla frequency squared [s-2]
    real(kind_phys) :: ti(ncol, pver + 1)      ! Interface temperature [K]

    ! Minimum value of Brunt-Vaisalla frequency squared.
    real(kind_phys), parameter :: n2min = 5.e-5_kind_phys

    errmsg = ''
    errflg = 0

    !------------------------------------------------------------------------
    ! Create Coords1D coordinate for gravity waves.
    !------------------------------------------------------------------------
    p = Coords1D(pint(:ncol, :pver+1))

    !------------------------------------------------------------------------
    ! Initialization of effective diffusion coefficients,
    ! accummulated tendencies from all gravity wave schemes,
    ! and quantities for energy checker.
    !------------------------------------------------------------------------
    egwdffi(:, :) = 0._kind_phys
    flx_heat(:) = 0._kind_phys
    scheme_name = 'gravity_wave_drag'

    tend_q(:,:,:) = 0._kind_phys
    tend_u(:,:)   = 0._kind_phys
    tend_v(:,:)   = 0._kind_phys
    tend_s(:,:)   = 0._kind_phys

    !------------------------------------------------------------------------
    ! Calculate local molecular diffusivity
    !------------------------------------------------------------------------
    if(do_molec_diff) then
      kvt_gw(:,:) = kvt(:,:)

      ! Use linear extrapolation of cpairv to top interface.
      kvt_gw(:, 1) = kvt_gw(:, 1)/ &
                   (1.5_kind_phys*cpairv(:ncol, 1) - &
                    0.5_kind_phys*cpairv(:ncol, 2))

      ! Interpolate cpairv to other interfaces.
      do k = 2, nbot_molec
        kvt_gw(:, k) = kvt_gw(:, k)/ &
                     (cpairv(:ncol, k + 1) + cpairv(:ncol, k))*2._kind_phys
      end do

    else
      kvt_gw(:,:) = 0._kind_phys
    endif

    !------------------------------------------------------------------------
    ! Determine the interface densities and Brunt-Vaisala frequencies.
    !------------------------------------------------------------------------

    ! The top interface values are calculated assuming an isothermal
    ! atmosphere above the top level.
    k = 1
    do i = 1, ncol
      ti(i, k) = t(i, k)
      rhoi(i, k) = p%ifc(i, k)/(rair*ti(i, k))
      ni(i, k) = sqrt(gravit*gravit/(cpair*ti(i, k)))
    end do

    ! Interior points use centered differences.
    ti(:, 2:pver) = midpoint_interp(t)
    do k = 2, pver
      do i = 1, ncol
        rhoi(i, k) = p%ifc(i, k)/(rair*ti(i, k))
        dtdp = (t(i, k) - t(i, k - 1))*p%rdst(i, k - 1)
        n2 = gravit*gravit/ti(i, k)*(1._kind_phys/cpair - rhoi(i, k)*dtdp)
        ni(i, k) = sqrt(max(n2min, n2))
      end do
    end do

    ! Bottom interface uses bottom level temperature, density; next interface
    ! B-V frequency.
    k = pver + 1
    do i = 1, ncol
      ti(i, k) = t(i, k - 1)
      rhoi(i, k) = p%ifc(i, k)/(rair*ti(i, k))
      ni(i, k) = ni(i, k - 1)
    end do

    !------------------------------------------------------------------------
    ! Determine the midpoint Brunt-Vaisala frequencies.
    !------------------------------------------------------------------------
    nm = midpoint_interp(ni)

  end subroutine gravity_wave_drag_prepare_profiles_run

!> \section arg_table_gravity_wave_drag_prepare_profiles_timestep_final Argument Table
!! \htmlinclude gravity_wave_drag_prepare_profiles_timestep_final.html
  subroutine gravity_wave_drag_prepare_profiles_timestep_final(p, errmsg, errflg)

    use coords_1d, only: Coords1D

    type(coords1d),     intent(inout)  :: p          ! Vertical coordinate for gravity waves [Coords1D, Pa]
    character(len=512), intent(out)    :: errmsg
    integer, intent(out)               :: errflg

    errmsg = ''
    errflg = 0

    call p%finalize()

  end subroutine gravity_wave_drag_prepare_profiles_timestep_final

end module gravity_wave_drag_interstitials
