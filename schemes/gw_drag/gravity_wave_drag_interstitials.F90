! Interstitial schemes for preparing fields needed for various gravity wave drag schemes
module gravity_wave_drag_interstitials
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: gravity_wave_drag_prepare_profiles_run
  public :: gravity_wave_drag_prepare_profiles_timestep_finalize

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
    p, &
    rhoi, &
    nm, ni, &
    errmsg, errflg)

    use gw_utils, only: midpoint_interp
    use coords_1d, only: Coords1D

    ! Input arguments
    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver

    real(kind_phys), intent(in)  :: cpair      ! Specific heat of dry air, constant pressure [J K-1 kg-1]
    real(kind_phys), intent(in)  :: rair       ! Dry air gas constant [J K-1 kg-1]
    real(kind_phys), intent(in)  :: gravit     ! [m s-2]
    real(kind_phys), intent(in)  :: pint(:, :) ! Interface pressures [Pa]
    real(kind_phys), intent(in)  :: t(:, :)    ! Midpoint temperatures [K]

    ! Output arguments
    type(Coords1D),  intent(out) :: p          ! Vertical coordinate for gravity waves [Coords1D, Pa]

    real(kind_phys), intent(out) :: rhoi(:, :) ! Density at layer interfaces [kg m-3]
    real(kind_phys), intent(out) :: nm(:, :)   ! Midpoint Brunt-Vaisalla frequencies [s-1]
    real(kind_phys), intent(out) :: ni(:, :)   ! Interface Brunt-Vaisalla frequencies [s-1]

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

!> \section arg_table_gravity_wave_drag_prepare_profiles_timestep_finalize Argument Table
!! \htmlinclude gravity_wave_drag_prepare_profiles_timestep_finalize.html
  subroutine gravity_wave_drag_prepare_profiles_timestep_finalize(p, errmsg, errflg)

    use coords_1d, only: Coords1D

    type(Coords1D),     intent(inout)  :: p          ! Vertical coordinate for gravity waves [Coords1D, Pa]
    character(len=512), intent(out)    :: errmsg
    integer, intent(out)               :: errflg

    errmsg = ''
    errflg = 0

    call p%finalize()

  end subroutine gravity_wave_drag_prepare_profiles_timestep_finalize

end module gravity_wave_drag_interstitials
