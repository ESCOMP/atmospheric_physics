module musica_ccpp_tuvx_height_grid

  implicit none

  private
  public :: create_height_grid, set_height_grid_values

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
  !> Units for height grid in TUV-x
  character(len=*), parameter, public :: height_grid_units = "km"

contains

  !> Creates a TUVX height grid from the host-model height grid
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
    height_grid => grid_t( height_grid_label, height_grid_units, &
                           vertical_interface_dimension, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_height_grid

  !> Sets TUVX height grid values from the host-model height grid
  subroutine set_height_grid_values( height_grid, host_midpoints, &
      host_edges, errmsg, errcode )

    use ccpp_kinds,       only: kind_phys
    use musica_ccpp_util, only: has_error_occurred
    use musica_tuvx_grid, only: grid_t
    use musica_util,      only: error_t

    type(grid_t),     intent(inout) :: height_grid
    real(kind_phys),  intent(in)    :: host_midpoints(:) ! km
    real(kind_phys),  intent(in)    :: host_edges(:)     ! km
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    type(error_t) :: error
    real(kind_phys) :: midpoints(size(host_midpoints)+1)
    real(kind_phys) :: edges(size(host_edges)+1)
    integer :: n_host_midpoints, n_host_edges

    if ( size(midpoints) /= height_grid%number_of_sections( error ) ) then
      errmsg = "[MUSICA Error] Invalid size of TUV-x mid-point heights."
      errcode = 1
      return
    end if
    if ( has_error_occurred( error, errmsg, errcode ) ) return
    if ( size(edges) /= height_grid%number_of_sections( error ) + 1 ) then
      errmsg = "[MUSICA Error] Invalid size of TUV-x interface heights."
      errcode = 1
      return
    end if
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    n_host_midpoints = size(host_midpoints)
    n_host_edges = size(host_edges)

    edges(1) = host_edges(n_host_edges)
    edges(2:n_host_edges) = host_midpoints(n_host_midpoints:1:-1)
    edges(n_host_edges+1) = host_edges(1)

    midpoints(1) = 0.5_kind_phys * &
              ( host_midpoints(n_host_midpoints) + host_edges(n_host_edges) )
    midpoints(2:n_host_midpoints) = host_edges(n_host_midpoints:2:-1)
    midpoints(n_host_midpoints+1) = 0.5_kind_phys * &
              ( edges(n_host_edges) + edges(n_host_edges+1) )

    call height_grid%set_edges( edges, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return
    call height_grid%set_midpoints( midpoints, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_height_grid_values

end module musica_ccpp_tuvx_height_grid
