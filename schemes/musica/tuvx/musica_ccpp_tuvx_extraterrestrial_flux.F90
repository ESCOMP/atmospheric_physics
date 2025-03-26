module musica_ccpp_tuvx_extraterrestrial_flux
  use ccpp_kinds, only: kind_phys

  implicit none

  private
  public :: create_extraterrestrial_flux_profile, set_extraterrestrial_flux_values

  !> Label for extraterrestrial_flux in TUV-x
  character(len=*), parameter, public :: extraterrestrial_flux_label = "extraterrestrial flux"
  !> Unit for extraterrestrial_flux in TUV-x
  character(len=*), parameter, public :: extraterrestrial_flux_unit = "photon cm-2 s-1"
  !> Wavelength grid interface values
  real(kind_phys), protected, allocatable :: wavelength_grid_interfaces_(:) ! nm
  !> Default value of number of wavelength grid bins
  integer, parameter :: DEFAULT_NUM_WAVELENGTH_BINS = 0
  !> Number of wavelength grid bins
  integer, protected :: num_wavelength_bins_ = DEFAULT_NUM_WAVELENGTH_BINS

contains

  !> Creates a TUV-x extraterrestrial flux profile from the host-model wavelength grid
  function create_extraterrestrial_flux_profile(wavelength_grid, &
      wavelength_grid_interfaces, errmsg, errcode) result( profile )
    use musica_util,                      only: error_t
    use musica_ccpp_util,                 only: has_error_occurred
    use musica_ccpp_tuvx_wavelength_grid, only: m_to_nm
    use musica_tuvx_grid,                 only: grid_t
    use musica_tuvx_profile,              only: profile_t

    type(grid_t),     intent(inout) :: wavelength_grid
    real(kind_phys),  intent(in)    :: wavelength_grid_interfaces(:) ! m
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode
    type(profile_t),  pointer       :: profile

    ! local variables
    type(error_t) :: error

    profile => profile_t( extraterrestrial_flux_label, extraterrestrial_flux_unit, &
                          wavelength_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    num_wavelength_bins_ = wavelength_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    allocate(wavelength_grid_interfaces_( size( wavelength_grid_interfaces ) ))

    wavelength_grid_interfaces_(:) = wavelength_grid_interfaces(:) * m_to_nm

  end function create_extraterrestrial_flux_profile

  !> Sets TUV-x extraterrestrial flux midpoints
  !
  ! Extraterrestrial flux is read from data files and interpolated to the
  ! TUV-x wavelength grid. CAM extraterrestrial flux values are multiplied by the
  ! width of the wavelength bins to get the TUV-x units of photon cm-2 s-1
  !
  ! TUV-x only uses mid-point values for extraterrestrial flux
  subroutine set_extraterrestrial_flux_values(profile, photolysis_wavelength_grid_interfaces, &
      extraterrestrial_flux, errmsg, errcode)
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t
    use ccpp_kinds,          only: kind_phys
    use ccpp_tuvx_utils,     only: rebin

    type(profile_t),  intent(inout) :: profile
    real(kind_phys),  intent(in)    :: photolysis_wavelength_grid_interfaces(:) ! nm
    real(kind_phys),  intent(in)    :: extraterrestrial_flux(:)                 ! photons cm-2 s-1 nm-1
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: midpoints(num_wavelength_bins_)

    if (.not. allocated(wavelength_grid_interfaces_)) then
      errmsg = "[MUSICA Error] Failed to allocate the TUV-x wavelength grid interface array"
      errcode = 1
      return
    end if

    if (num_wavelength_bins_ <= DEFAULT_NUM_WAVELENGTH_BINS) then
      errmsg = "[MUSICA Error] Invalid size of TUV-x wavelength bins."
      errcode = 1
      deallocate( wavelength_grid_interfaces_ )
      return
    end if

    ! Regrid normalized flux to TUV-x wavelength grid
    call rebin( size(photolysis_wavelength_grid_interfaces) - 1, num_wavelength_bins_,     &
                photolysis_wavelength_grid_interfaces, wavelength_grid_interfaces_, &
                extraterrestrial_flux, midpoints )

    ! Convert normalized flux to flux on TUV-x wavelength grid
    midpoints = midpoints * ( wavelength_grid_interfaces_(2 : num_wavelength_bins_ + 1) &
                 - wavelength_grid_interfaces_(1 :num_wavelength_bins_) )

    call profile%set_midpoint_values( midpoints, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) then
      deallocate( wavelength_grid_interfaces_ )
      return
    end if

  end subroutine set_extraterrestrial_flux_values

end module musica_ccpp_tuvx_extraterrestrial_flux