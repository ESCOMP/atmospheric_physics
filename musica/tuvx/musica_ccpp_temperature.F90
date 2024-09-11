module musica_ccpp_temperature
  implicit none

  private
  public :: create_temperature_profile, set_temperatures

  !> Label for temperature in TUV-x
  character(len=*), parameter :: temperature_label = "temperature"
  !> Units for temperature in TUV-x
  character(len=*), parameter :: temperature_units = "K"

contains

  !> Creates a TUVX height grid from the host-model height grid
  function create_temperature_profile( height_grid, errmsg, errcode ) & 
      result( profile )

    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_grid,    only: grid_t
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(grid_t),     intent(in)  :: height_grid
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode
    type(profile_t),  pointer     :: profile

    type(error_t) :: error
    
    profile => profile_t( temperature_label, temperature_units, &
                          height_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_temperature_profile

  !> Sets TUVX temperatures from the host-model temperatures
  subroutine set_temperatures( profile, host_temperature_mid, &
      host_surface_temperature, errmsg, errcode )

    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t
    
    type(profile_t),  intent(inout) :: profile
    real(kind_phys),  intent(in)    :: host_temperature_mid(:)  ! K
    real(kind_phys),  intent(in)    :: host_surface_temperature ! K
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    type(error_t)   :: error
    real(kind_phys) :: midpoints(size(host_temperature_mid)+2)
    integer         :: n_host_midpoints

    n_host_midpoints = size(host_temperature_mid)

    midpoints(1) = host_surface_temperature
    midpoints(2:n_host_midpoints+1) = host_temperature_mid(n_host_midpoints:1:-1)
    midpoints(n_host_midpoints+2) = host_temperature_mid(1)

    call profile%set_midpoint_values(midpoints, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) return

 end subroutine set_temperatures
end module