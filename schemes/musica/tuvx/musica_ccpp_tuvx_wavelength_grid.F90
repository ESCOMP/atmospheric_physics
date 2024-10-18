module musica_ccpp_tuvx_wavelength_grid

  implicit none

  private
  public :: create_wavelength_grid

  !> Label for height grid in TUV-x
  character(len=*), parameter, public :: wavelength_grid_label = "wavelength"
  !> Units for height grid in TUV-x
  character(len=*), parameter, public :: wavelength_grid_unit = "nm"

contains

  !> Creates a TUV-x wavelength grid
  function create_wavelength_grid( num_wavelength_bin, errmsg, errcode ) &
      result( wavelength_grid )

    use musica_ccpp_util, only: has_error_occurred
    use musica_tuvx_grid, only: grid_t
    use musica_util,      only: error_t

    integer,          intent(in)  :: num_wavelength_bin
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode
    type(grid_t),     pointer     :: wavelength_grid

    ! local variable
    type(error_t) :: error

    wavelength_grid => grid_t( wavelength_grid_label, wavelength_grid_unit, &
                               num_wavelength_bin, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_wavelength_grid

end module musica_ccpp_tuvx_wavelength_grid