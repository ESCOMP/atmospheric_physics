! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_tuvx_aerosol_optics
  implicit none

  private
  public :: create_aerosol_optics_radiator, set_aerosol_optics_values

  !> Label for aerosol optical properties in TUV-x
  character(len=*), parameter, public :: aerosol_optics_label = "aerosols"
  !> Label
  character(len=*), parameter, public :: \
    aerosol_optical_depth_label = "optical depths"
  character(len=*), parameter, public :: \
    aerosol_single_scattering_albedo_label = "single scattering albedos"
  character(len=*), parameter, public :: \
    aerosol_asymmetry_factor_label = "asymmetry factor"
  !> Unit
  character(len=*), parameter, public :: aerosol_optical_depth_unit = "none"
  character(len=*), parameter, public :: aerosol_single_scattering_albedo_unit = "none"
  character(len=*), parameter, public :: aerosol_asymmetry_factor_unit = "none"
  !> Default value of number of vertical levels
  integer, parameter :: DEFAULT_NUM_VERTICAL_LEVELS = 0
  !> Number of vertical levels
  integer, protected :: num_vertical_levels = DEFAULT_NUM_VERTICAL_LEVELS
  !> Default value of number of wavelength bins
  integer, parameter :: DEFAULT_NUM_WAVELENGTH_BINS = 0
  !> Number of wavelength bins
  integer, protected :: num_wavelength_bins = DEFAULT_NUM_WAVELENGTH_BINS
  !> Default value of number of streams
  integer, parameter :: DEFAULT_NUM_STREAMS = 1
  !> Number of streams
  integer, protected :: num_streams = DEFAULT_NUM_STREAMS

contains

  !> Creates a TUV-x aerosol optics radiator
  function create_aerosol_optics_radiator( height_grid, wavelength_grid, &
      errmsg, errcode ) result( radiator )
    use musica_ccpp_util,     only: has_error_occurred
    use musica_tuvx_grid,     only: grid_t
    use musica_tuvx_radiator, only: radiator_t
    use musica_util,          only: error_t

    type(grid_t),     intent(inout) :: height_grid
    type(grid_t),     intent(inout) :: wavelength_grid
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode
    type(radiator_t), pointer       :: radiator

    ! local variables
    type(error_t) :: error

    num_vertical_levels = height_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    num_wavelength_bins = wavelength_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    radiator => radiator_t( aerosol_optics_label, height_grid, wavelength_grid, &
                            error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_aerosol_optics_radiator

  !> Sets TUV-x aerosol optics values
  !! Without the optional arguments the optical properties are zero (no aerosol model).
  !! With them, the host-grid properties of an aerosol model (aerosol_model_t) are mapped
  !! to the TUV-x height grid:
  !!   optical depth = specific extinction (m2 per kg of air) * dp / g
  !! with the vertical order reversed (TUV-x counts from the surface) and the extra TUV-x
  !! top layer left without aerosol, as for the cloud optics.
  !! @param extinction Specific extinction per kg of air (m2 kg-1) [host layer, wavelength]
  !! @param single_scattering_albedo Single scattering albedo (unitless) [host layer, wavelength]
  !! @param asymmetry_factor Asymmetry factor (unitless) [host layer, wavelength]
  !! @param delta_pressure Pressure thickness of the host layers (Pa)
  !! @param reciprocal_of_gravitational_acceleration (s2 m-1)
  subroutine set_aerosol_optics_values( radiator, errmsg, errcode, extinction, &
      single_scattering_albedo, asymmetry_factor, delta_pressure, &
      reciprocal_of_gravitational_acceleration )
    use ccpp_kinds,           only: kind_phys
    use musica_ccpp_util,     only: has_error_occurred
    use musica_tuvx_radiator, only: radiator_t
    use musica_util,          only: error_t

    type(radiator_t), intent(inout) :: radiator
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode
    real(kind_phys), optional, intent(in) :: extinction(:,:)
    real(kind_phys), optional, intent(in) :: single_scattering_albedo(:,:)
    real(kind_phys), optional, intent(in) :: asymmetry_factor(:,:)
    real(kind_phys), optional, intent(in) :: delta_pressure(:)
    real(kind_phys), optional, intent(in) :: reciprocal_of_gravitational_acceleration

    ! local variables
    type(error_t)   :: error
    integer         :: i_level, i_host, num_host_levels
    real(kind_phys) :: \
      aerosol_optical_depth(num_vertical_levels, num_wavelength_bins)
    real(kind_phys) :: \
      aerosol_single_scattering_albedo(num_vertical_levels, num_wavelength_bins)
    real(kind_phys) :: \
      aerosol_asymmetry_factor(num_vertical_levels, num_wavelength_bins, num_streams)

    errmsg = ''
    errcode = 0
    aerosol_optical_depth(:,:) = 0.0_kind_phys
    aerosol_single_scattering_albedo(:,:) = 0.0_kind_phys
    aerosol_asymmetry_factor(:,:,:) = 0.0_kind_phys

    if ( present( extinction ) ) then
      if ( .not. ( present( single_scattering_albedo ) .and. present( asymmetry_factor ) .and. &
                   present( delta_pressure ) .and. &
                   present( reciprocal_of_gravitational_acceleration ) ) ) then
        errmsg = "[MUSICA Error] Aerosol extinction needs single scattering albedo, " // &
                 "asymmetry factor, pressure thickness and 1/g for TUV-x."
        errcode = 1
        return
      end if
      num_host_levels = size( extinction, dim=1 )
      if ( num_host_levels + 1 /= num_vertical_levels .or. &
           size( extinction, dim=2 ) /= num_wavelength_bins .or. &
           any( shape( single_scattering_albedo ) /= shape( extinction ) ) .or. &
           any( shape( asymmetry_factor ) /= shape( extinction ) ) .or. &
           size( delta_pressure ) /= num_host_levels ) then
        errmsg = "[MUSICA Error] Invalid size of aerosol optical properties for TUV-x."
        errcode = 1
        return
      end if
      do i_level = 1, num_host_levels
        i_host = num_host_levels - i_level + 1
        aerosol_optical_depth(i_level,:) = extinction(i_host,:) * delta_pressure(i_host) &
                                         * reciprocal_of_gravitational_acceleration
        aerosol_single_scattering_albedo(i_level,:) = single_scattering_albedo(i_host,:)
        aerosol_asymmetry_factor(i_level,:,1) = asymmetry_factor(i_host,:)
      end do
    end if

    call radiator%set_optical_depths( aerosol_optical_depth, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call radiator%set_single_scattering_albedos( aerosol_single_scattering_albedo, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call radiator%set_asymmetry_factors( aerosol_asymmetry_factor, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_aerosol_optics_values

end module musica_ccpp_tuvx_aerosol_optics
