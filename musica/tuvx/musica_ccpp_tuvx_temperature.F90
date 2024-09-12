module musica_ccpp_tuvx_temperature
  implicit none

  private
  public :: create_temperature_profile, set_temperature_values

  !> Label for temperature in TUV-x
  character(len=*), parameter :: temperature_label = "temperature"
  !> Units for temperature in TUV-x
  character(len=*), parameter :: temperature_units = "K"

contains

  !> Creates a TUVX temperature profile from the host-model height grid
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

  !> Sets TUVX temperature edges from the host-model temperature midpoints
  subroutine set_temperature_values( profile, host_temperature_mid, &
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
    real(kind_phys) :: edges(size(host_temperature_mid)+2)
    integer         :: n_host_temperature_mid

    n_host_temperature_mid = size(host_temperature_mid)

    edges(1) = host_surface_temperature
    edges(2:n_host_temperature_mid+1) = host_temperature_mid(n_host_temperature_mid:1:-1)
    edges(n_host_temperature_mid+2) = host_temperature_mid(1)

    call profile%set_edge_values( edges, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_temperature_values

end module musica_ccpp_tuvx_temperature