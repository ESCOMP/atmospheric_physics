module musica_ccpp_tuvx_surface_albedo
  implicit none

  private
  public :: create_surface_albedo_profile, set_surface_albedo_values, &
            surface_albedo_label, surface_albedo_unit

  !> Label for surface albedo in TUV-x
  character(len=*), parameter :: surface_albedo_label = "surface albedo"
  !> Unit for surface albedo in TUV-x
  character(len=*), parameter :: surface_albedo_unit = "none"
  !> Default value of number of wavelength bins
  integer, parameter :: DEFAULT_NUM_WAVELENGTH_BINS = 0
  !> Number of wavelength bins
  integer, protected :: num_wavelength_bins = DEFAULT_NUM_WAVELENGTH_BINS

contains

  !> Creates a TUV-x surface albedo profile from the host-model wavelength grid
  function create_surface_albedo_profile( wavelength_grid, errmsg, errcode ) & 
      result( profile )
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_grid,    only: grid_t
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(grid_t),     intent(inout) :: wavelength_grid
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode
    type(profile_t),  pointer       :: profile

    ! local variables
    type(error_t) :: error

    num_wavelength_bins = wavelength_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    profile => profile_t( surface_albedo_label, surface_albedo_unit, &
                          wavelength_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_surface_albedo_profile

  !> Sets TUV-x surface albedo values
  !!
  !! CAM uses a single value for surface albedo at all wavelengths
  subroutine set_surface_albedo_values( profile, host_surface_albedo, &
                                        errmsg, errcode )
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(profile_t),  intent(inout) :: profile
    real(kind_phys),  intent(in)    :: host_surface_albedo ! unitless
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: surface_albedo_interfaces(num_wavelength_bins + 1)

    if (size(surface_albedo_interfaces) <= DEFAULT_NUM_WAVELENGTH_BINS + 1) then
      errmsg = "[MUSICA Error] Invalid size of TUV-x wavelength interfaces."
      errcode = 1
      return
    end if

    surface_albedo_interfaces(:) = host_surface_albedo

    call profile%set_edge_values( surface_albedo_interfaces, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_surface_albedo_values

end module musica_ccpp_tuvx_surface_albedo