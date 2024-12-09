module musica_ccpp_tuvx_aerosol_optics
  implicit none

  private
  public :: set_aerosol_optics_values

  !> Label
  character(len=*), parameter, public :: aerosol_optical_depths_label = "optical depths"
  character(len=*), parameter, public :: aerosol_single_scattering_albedo_label = "single scattering albedos"
  character(len=*), parameter, public :: aerosol_asymmetry_factor_label = "asymmetry factor"
  character(len=*), parameter, public :: aerosol_visible_optical_depth_label = "550 nm optical depth"
  !> Unit
  character(len=*), parameter, public :: aerosol_optical_depths_unit = "none"
  character(len=*), parameter, public :: aerosol_single_scattering_albedo_unit = "none"
  character(len=*), parameter, public :: aerosol_asymmetry_factor_unit = "none"
  character(len=*), parameter, public :: aerosol_visible_optical_depth_unit = "none"
  !> Default value of number of wavelength bins
  integer, parameter :: DEFAULT_NUM_WAVELENGTH_BINS = 0
  !> Number of wavelength bins
  integer, protected :: num_wavelength_bins_ = DEFAULT_NUM_WAVELENGTH_BINS

contains

  subroutine set_aerosol_optics_values( profile, &
                                        host_aerosol_optical_depths, &
                                        host_aerosol_single_scattering_albedo, &
                                        host_aerosol_asymmetry_factor, &
                                        host_aerosol_visible_optical_depth, &
                                        errmsg, errcode )
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(profile_t),  intent(inout) :: profile
    real(kind_phys),  intent(in)    :: host_aerosol_optical_depths
    real(kind_phys),  intent(in)    :: host_aerosol_single_scattering_albedo
    real(kind_phys),  intent(in)    :: host_aerosol_asymmetry_factor
    real(kind_phys),  intent(in)    :: host_aerosol_visible_optical_depth
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables

  end subroutine set_aerosol_optics_values

end module musica_ccpp_tuvx_aerosol_optics
