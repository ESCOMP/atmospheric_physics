module musica_ccpp_tuvx_aerosol_optics
  implicit none

  private
  public :: set_aerosol_optics_values

  !> Label
  character(len=*), parameter, public :: aerosol_extinction_label = "ext"
  !> Unit
  character(len=*), parameter, public :: aerosol_extinction_unit = "units"
  !> Default value of number of wavelength bins
  integer, parameter :: DEFAULT_NUM_WAVELENGTH_BINS = 0
  !> Number of wavelength bins
  integer, protected :: num_wavelength_bins_ = DEFAULT_NUM_WAVELENGTH_BINS

contains

  subroutine set_aerosol_optics_values( profile, host_aerosol_extinction, &
                                        errmsg, errcode )
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_util,    only: has_error_occurred
    use musica_util,         only: error_t

    real(kind_phys),  intent(in)    :: host_aerosol_extinction ! units
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables

  end subroutine set_aerosol_optics_values

end module musica_ccpp_tuvx_aerosol_optics
