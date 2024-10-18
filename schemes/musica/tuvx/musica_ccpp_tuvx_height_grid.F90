module musica_ccpp_tuvx_height_grid

  implicit none

  private
  public :: create_height_grid, set_height_grid_values, calculate_heights

  ! Conversions between the CAM-SIMA height grid and the TUVX height grid
  !
  !-----------------------------------------------------------------------
  ! Notes on the conversion between the host-model height grid and the TUVX
  !
  !  TUV-x heights are "bottom-up" and require atmospheric constituent
  !  concentrations at interfaces. Therefore, CAM-SIMA mid-points are used
  !  as TUV-x grid interfaces, with an additional layer introduced between
  !  the surface and the lowest CAM-SIMA mid-point, and a layer at the
  !  top of the TUV-x grid to hold species densities above the top CAM-SIMA
  !  mid-point.
  !
  !  Here,
  !    - i_int is the index of an interface
  !    - i_mid is the index of a mid-point
  !    - pver is the CCPP vertical_layer_dimension
  !    - pver+1 is the CCPP vertical_interface_dimension

  !  ---- (interface)  ===== (mid-point)
  !
  !        CAM                                  TUV-x
  ! ************************ (exo values) *****************************
  ! ------(top)------ i_int = 1           -------(top)------ i_int = pver + 2
  !                                       ================== i_mid = pver + 1
  ! ================= i_mid = 1           ------------------ i_int = pver + 1
  ! ----------------- i_int = 2           ================== i_mid = pver
  !                                       ------------------ i_int = pver
  !        ||
  !        ||                                     ||
  !                                               ||
  ! ----------------- i_int = pver
  ! ================= i_imd = pver        ------------------ i_int = 2
  !                                       ================== i_mid = 1
  ! -----(ground)---- i_int = pver+1      -----(ground)----- i_int = 1
  !
  !-----------------------------------------------------------------------

  !> Label for height grid in TUV-x
  character(len=*), parameter, public :: height_grid_label = "height"
  !> Unit for height grid in TUV-x
  character(len=*), parameter, public :: height_grid_unit = "km"

contains

  !> Creates a TUV-x height grid
  function create_height_grid( vertical_layer_dimension, &
      vertical_interface_dimension, errmsg, errcode ) result( height_grid )

    use musica_ccpp_util, only: has_error_occurred
    use musica_tuvx_grid, only: grid_t
    use musica_util,      only: error_t

    integer,          intent(in)  :: vertical_layer_dimension     ! (count)
    integer,          intent(in)  :: vertical_interface_dimension ! (count)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode
    type(grid_t),     pointer     :: height_grid

    ! local variable
    type(error_t) :: error

    height_grid => null()
    if ( vertical_layer_dimension < 1 ) then
      errmsg = "[MUSICA Error] Invalid vertical_layer_dimension."
      errcode = 1
      return
    end if
    if ( vertical_interface_dimension - vertical_layer_dimension /= 1 ) then
      errmsg = "[MUSICA Error] Invalid vertical_interface_dimension."
      errcode = 1
      return
    end if
    height_grid => grid_t( height_grid_label, height_grid_unit, &
                           vertical_interface_dimension, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_height_grid

  !> Sets TUV-x height grid values from the host-model height grid
  subroutine set_height_grid_values( height_grid, host_midpoints, &
      host_interfaces, errmsg, errcode )

    use ccpp_kinds,       only: kind_phys
    use musica_ccpp_util, only: has_error_occurred
    use musica_tuvx_grid, only: grid_t
    use musica_util,      only: error_t

    type(grid_t),     intent(inout) :: height_grid
    real(kind_phys),  intent(in)    :: host_midpoints(:)  ! km
    real(kind_phys),  intent(in)    :: host_interfaces(:) ! km
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: midpoints(size(host_midpoints)+1)
    real(kind_phys) :: interfaces(size(host_interfaces)+1)
    integer         :: n_host_midpoints, n_host_interfaces

    if ( size(midpoints) /= height_grid%number_of_sections( error ) ) then
      errmsg = "[MUSICA Error] Invalid size of TUV-x mid-point heights."
      errcode = 1
      return
    end if

    if ( size(interfaces) /= height_grid%number_of_sections( error ) + 1 ) then
      errmsg = "[MUSICA Error] Invalid size of TUV-x interface heights."
      errcode = 1
      return
    end if

    n_host_midpoints = size(host_midpoints)
    n_host_interfaces = size(host_interfaces)

    interfaces(1) = host_interfaces(n_host_interfaces)
    interfaces(2:n_host_interfaces) = host_midpoints(n_host_midpoints:1:-1)
    interfaces(n_host_interfaces+1) = host_interfaces(1)

    midpoints(1) = 0.5_kind_phys * &
              ( host_midpoints(n_host_midpoints) + host_interfaces(n_host_interfaces) )
    midpoints(2:n_host_midpoints) = host_interfaces(n_host_midpoints:2:-1)
    midpoints(n_host_midpoints+1) = 0.5_kind_phys * &
              ( interfaces(n_host_interfaces) + interfaces(n_host_interfaces+1) )

    call height_grid%set_edges( interfaces, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call height_grid%set_midpoints( midpoints, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_height_grid_values

  !> Calculates the heights needed for the TUV-x grid based on available data
  !!
  !! Uses the reciprocal of gravitational acceleration, the surface geopotential,
  !! and the geopotential height wrt the surface and midpoints/interfaces to calculate
  !! the midpoint/interface height above sea level in km needed for the TUV-x grid.
  !!
  !! The equation used is taked from CAMChem
  !! (see https://github.com/ESCOMP/CAM/blob/f0e489e9708ce7b91635f6d4997fbf1e390b0dbb/src/chemistry/mozart/mo_gas_phase_chemdr.F90#L514-L526)
  subroutine calculate_heights( geopotential_height_wrt_surface_at_midpoint, &
      geopotential_height_wrt_surface_at_interface, &
      surface_geopotential, reciprocal_of_gravitational_acceleration, &
      height_midpoints, height_interfaces )

    use ccpp_kinds,       only: kind_phys

    real(kind_phys), intent(in)  :: geopotential_height_wrt_surface_at_midpoint(:)  ! m
    real(kind_phys), intent(in)  :: geopotential_height_wrt_surface_at_interface(:) ! m
    real(kind_phys), intent(in)  :: surface_geopotential ! m2 s-2
    real(kind_phys), intent(in)  :: reciprocal_of_gravitational_acceleration ! s2 m-1
    real(kind_phys), intent(out) :: height_midpoints(:)  ! km
    real(kind_phys), intent(out) :: height_interfaces(:) ! km

    ! local variable
    real(kind_phys) :: surface_height ! m

    surface_height = surface_geopotential * reciprocal_of_gravitational_acceleration
    height_midpoints(:) = 0.001_kind_phys * &
                          ( geopotential_height_wrt_surface_at_midpoint(:) &
                            + surface_height )
    height_interfaces(:) = 0.001_kind_phys * &
                           ( geopotential_height_wrt_surface_at_interface(:) &
                             + surface_height )

  end subroutine calculate_heights

end module musica_ccpp_tuvx_height_grid
