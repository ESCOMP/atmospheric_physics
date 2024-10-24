module musica_ccpp_tuvx_temperature
  implicit none

  private
  public :: create_temperature_profile, set_temperature_values, &
            temperature_label, temperature_unit

  !> Label for temperature in TUV-x
  character(len=*), parameter :: temperature_label = "temperature"
  !> Unit for temperature in TUV-x
  character(len=*), parameter :: temperature_unit = "K"

contains

  !> Creates a TUV-x temperature profile
  function create_temperature_profile(height_grid, errmsg, errcode) &
      result(profile)
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_grid,    only: grid_t
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(grid_t),     intent(in)  :: height_grid
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode
    type(profile_t),  pointer     :: profile

    ! local variables
    type(error_t) :: error

    profile => profile_t( temperature_label, temperature_unit, &
                          height_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_temperature_profile

  !> Sets TUV-x temperature values from host-model temperatures
  !!
  !! See description of `musica_ccpp_tuvx_height_grid.F90` for
  !! CAM-SIMA <-> TUV-x height grid mapping
  subroutine set_temperature_values(profile, host_midpoint_temperatures, &
                                    host_surface_temperature, errmsg, errcode)
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t
    use ccpp_kinds,          only: kind_phys

    type(profile_t),  intent(inout) :: profile
    real(kind_phys),  intent(in)    :: host_midpoint_temperatures(:) ! K
    real(kind_phys),  intent(in)    :: host_surface_temperature      ! K
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: interfaces(size(host_midpoint_temperatures)+2)
    integer         :: n_host_midpoint_temperatures

    n_host_midpoint_temperatures = size(host_midpoint_temperatures)

    interfaces(1) = host_surface_temperature
    interfaces(2:n_host_midpoint_temperatures+1) = host_midpoint_temperatures(n_host_midpoint_temperatures:1:-1)
    interfaces(n_host_midpoint_temperatures+2) = host_midpoint_temperatures(1)

    call profile%set_edge_values( interfaces, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_temperature_values

end module musica_ccpp_tuvx_temperature