module musica_ccpp_tuvx_aerosol_optics
  implicit none

  private
  public :: set_aerosol_optics_values

  !> Label
  character(len=*), parameter, public :: aerosol_optical_depth_label = "optical depths"
  character(len=*), parameter, public :: aerosol_single_scattering_albedo_label = "single scattering albedos"
  character(len=*), parameter, public :: aerosol_asymmetry_factor_label = "asymmetry factor"
  character(len=*), parameter, public :: aerosol_visible_optical_depth_label = "550 nm optical depth"
  !> Unit
  character(len=*), parameter, public :: aerosol_optical_depth_unit = "none"
  character(len=*), parameter, public :: aerosol_single_scattering_albedo_unit = "none"
  character(len=*), parameter, public :: aerosol_asymmetry_factor_unit = "none"
  character(len=*), parameter, public :: aerosol_visible_optical_depth_unit = "none"
  !> Default value of number of wavelength bins
  integer, parameter :: DEFAULT_NUM_WAVELENGTH_BINS = 0
  !> Number of wavelength bins
  integer, protected :: num_wavelength_bins_ = DEFAULT_NUM_WAVELENGTH_BINS

contains

  !> Creates a TUV-x aerosol optical depth profile from the host-model wavelength grid
  function create_aerosol_otpical_depth_profile( wavelength_grid, errmsg, errcode ) & 
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

    profile => profile_t( aerosol_optical_depth_label, aerosol_optical_depth_unit, &
                          wavelength_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    num_wavelength_bins_ = wavelength_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_aerosol_optical_depth_profile

  subroutine set_aerosol_optics_values( profile, &
                                        host_aerosol_optical_depth, &
                                        host_aerosol_single_scattering_albedo, &
                                        host_aerosol_asymmetry_factor, &
                                        host_aerosol_visible_optical_depth, &
                                        errmsg, errcode )
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(profile_t),  intent(inout) :: profile
    real(kind_phys),  intent(in)    :: host_aerosol_optical_depth
    real(kind_phys),  intent(in)    :: host_aerosol_single_scattering_albedo
    real(kind_phys),  intent(in)    :: host_aerosol_asymmetry_factor
    real(kind_phys),  intent(in)    :: host_aerosol_visible_optical_depth
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: aerosol_optical_depth_interfaces(num_wavelength_bins_ + 1)
    real(kind_phys) :: aerosol_single_scattering_albedo
    real(kind_phys) :: aerosol_asymmetry_factor
    real(kind_phys) :: aerosol_visible_optical_depth

    if (num_wavelength_bins_ <= DEFAULT_NUM_WAVELENGTH_BINS) then
      errmsg = "[MUSICA Error] Invalid size of TUV-x wavelength bins."
      errcode = 1
      return
    end if

    aerosol_optical_depth_interfaces(:) = host_aerosol_optical_depth
    aerosol_single_scattering_albedo = host_aerosol_single_scattering_albedo
    aerosol_asymmetry_factor = host_aerosol_asymmetry_factor
    aerosol_visible_optical_depth = host_aerosol_visible_optical_depth

    call profile%set_edge_values( aerosol_optical_depth_interfaces, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_aerosol_optics_values

end module musica_ccpp_tuvx_aerosol_optics
