module musica_ccpp_tuvx_surface_albedo
  implicit none

  private
  public :: create_surface_albedo_profile, set_surface_albedo_values

  !> Label for surface albedo in TUV-x
  character(len=*), parameter :: surface_albedo_label = "surface_albedo"
  !> Units for surface albedo in TUV-x
  character(len=*), parameter :: surface_albedo_units = ""

contains

  !> Creates a TUV-x surface albedo profile from the host-model wavelength grid
  function create_surface_albedo_profile( wavelength_grid, errmsg, errcode ) & 
      result( profile )

    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_grid,    only: grid_t
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    ! Arguments
    type(grid_t),     intent(in)  :: wavelength_grid
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode

    ! Return value
    type(profile_t),  pointer :: profile

    ! Local variables
    type(error_t) :: error

    profile => profile_t( surface_albedo_label, surface_albedo_units, &
                          wavelength_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_surface_albedo_profile

  !> Sets TUV-x temperature edges from the host-model temperature midpoints
  subroutine set_surface_albedo_values( profile, host_surface_albedo, &
      errmsg, errcode )

    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    ! Arguments
    type(profile_t),  intent(inout) :: profile
    real(kind_phys),  intent(in)    :: host_surface_albedo ! unitless
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! Local variables
    type(error_t)   :: error

    ! TODO(jiwon) - is this correct?
    ! real(r8) :: albedos(this%n_wavelength_bins_ + 1)
    real(kind_phys) :: surface_albedo(size(host_surface_albedo))

    surface_albedo = host_surface_albedo

    call profile%set_edge_values( surface_albedo, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_surface_albedo

end module musica_ccpp_tuvx_surface_albedo