! Copyright (C) 2024-2025 National Science Foundation-National Center for Atmospheric Research
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
  ! Temporarily setting optical properties to zero until aerosol optical
  ! property calculations are ported to CAM-SIMA.
  subroutine set_aerosol_optics_values( radiator, errmsg, errcode )
    use ccpp_kinds,           only: kind_phys
    use musica_ccpp_util,     only: has_error_occurred
    use musica_tuvx_radiator, only: radiator_t
    use musica_util,          only: error_t

    type(radiator_t), intent(inout) :: radiator
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: \
      aerosol_optical_depth(num_vertical_levels, num_wavelength_bins)
    real(kind_phys) :: \
      aerosol_single_scattering_albedo(num_vertical_levels, num_wavelength_bins)
    real(kind_phys) :: \
      aerosol_asymmetry_factor(num_vertical_levels, num_wavelength_bins, num_streams)

    aerosol_optical_depth(:,:) = 0.0_kind_phys
    aerosol_single_scattering_albedo(:,:) = 0.0_kind_phys
    aerosol_asymmetry_factor(:,:,:) = 0.0_kind_phys

    call radiator%set_optical_depths( aerosol_optical_depth, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call radiator%set_single_scattering_albedos( aerosol_single_scattering_albedo, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call radiator%set_asymmetry_factors( aerosol_asymmetry_factor, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_aerosol_optics_values

end module musica_ccpp_tuvx_aerosol_optics
