! Copyright (C) 2024 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_tuvx_cloud_optics
  implicit none

  private
  public :: create_cloud_optics_radiator, set_cloud_optics_values, &
            cloud_optics_label

  ! This module is used to set the optical properties of clouds in TUV-x.
  ! Optical properties are defined as a function of wavelength and height,
  ! and include the cloud optical depth, single scattering albedo,
  ! and asymmetry parameter.
  !
  ! See module_ccpp_tuvx_height_grid for the definition of the height grid
  ! and its mapping to the CAM-SIMA vertical grid.

  !> Label for cloud optical properties in TUV-x
  character(len=*), parameter :: cloud_optics_label = "clouds"
  !> Default value of number of vertical levels
  integer, parameter :: DEFAULT_NUM_VERTICAL_LEVELS = 0
  !> Number of vertical levels
  integer, protected :: num_vertical_levels = DEFAULT_NUM_VERTICAL_LEVELS
  !> Default value of number of wavelength bins
  integer, parameter :: DEFAULT_NUM_WAVELENGTH_BINS = 0
  !> Number of wavelength bins
  integer, protected :: num_wavelength_bins = DEFAULT_NUM_WAVELENGTH_BINS

contains

  !> Creates a TUV-x cloud optics radiator from the host-model wavelength grid
  function create_cloud_optics_radiator( height_grid, wavelength_grid, &
      errmsg, errcode ) result( radiator )
    use musica_ccpp_util,     only: has_error_occurred
    use musica_tuvx_grid,     only: grid_t
    use musica_tuvx_radiator, only: radiator_t
    use musica_util,          only: error_t

    type(grid_t),                        intent(inout) :: height_grid
    type(grid_t),                        intent(inout) :: wavelength_grid
    character(len=*),                    intent(out)   :: errmsg
    integer,                             intent(out)   :: errcode
    type(radiator_t), pointer                          :: radiator

    ! local variables
    type(error_t) :: error

    num_vertical_levels = height_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    num_wavelength_bins = wavelength_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    radiator => radiator_t( cloud_optics_label, height_grid, wavelength_grid, &
                            error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_cloud_optics_radiator

  !> Sets TUV-x cloud optics values
  subroutine set_cloud_optics_values( radiator, cloud_fraction, delta_pressure, &
                                      cloud_liquid_water_content, &
                                      reciprocal_of_gravitational_acceleration, &
                                      errmsg, errcode )
    use ccpp_kinds,           only: kind_phys
    use musica_ccpp_util,     only: has_error_occurred
    use musica_tuvx_radiator, only: radiator_t
    use musica_util,          only: error_t

    type(radiator_t), intent(inout) :: radiator
    real(kind_phys),  intent(in)    :: delta_pressure(:) ! pressure delta about vertical level midpoints (Pa)
    real(kind_phys),  intent(in)    :: cloud_fraction(:) ! (unitless)
    real(kind_phys),  intent(in)    :: cloud_liquid_water_content(:) ! (kg/kg)
    real(kind_phys),  intent(in)    :: reciprocal_of_gravitational_acceleration ! (s^2/m)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: od(num_vertical_levels)
    real(kind_phys) :: cloud_optical_depth(num_vertical_levels, num_wavelength_bins)
    integer         :: n_host_midpoints
    integer         :: i_level
    
    n_host_midpoints = size(cloud_fraction)
    if ( size(delta_pressure) /= n_host_midpoints ) then
      errmsg = "[MUSICA Error] Invalid size of cloud pressure delta for TUV-x."
      errcode = 1
      return
    end if
    if ( size(cloud_liquid_water_content) /= n_host_midpoints ) then
      errmsg = "[MUSICA Error] Invalid size of cloud liquid water content for TUV-x."
      errcode = 1
      return
    end if
    if ( n_host_midpoints + 1 /= num_vertical_levels ) then
      errmsg = "[MUSICA Error] Invalid size of cloud fraction for TUV-x."
      errcode = 1
      return
    end if

    ! Estimate cloud optical depth (od) [unitless] from cloud fraction (cf)
    ! [unitless] and liquid water content (lwc) [kg kg-1] by first calculating
    ! the cloud liquid water path (lwp) [kg m-2]:
    !   lwp = 1/g * lwc * dP / cf
    ! where g is the gravitational acceleration [m s-2] and dP is the change in
    ! pressure across the vertical level [Pa]. 
    ! The cloud optical depth is then estimated as:
    !   od = lwp * 155 * cf^1.5
    ! A constant cloud optical depth is used for all wavelengths.
    do i_level = 1, n_host_midpoints
      if ( cloud_fraction(i_level) > 0.0_kind_phys ) then
        od(i_level) = ( reciprocal_of_gravitational_acceleration &
                      * cloud_liquid_water_content(i_level) * delta_pressure(i_level) &
                      / cloud_fraction(i_level) ) * 155.0_kind_phys * cloud_fraction(i_level)**1.5_kind_phys
      else
        od(i_level) = 0.0_kind_phys
      end if
    end do
    do i_level = 1, n_host_midpoints
      cloud_optical_depth(i_level, :) = od(n_host_midpoints-i_level+1)
    end do
    cloud_optical_depth(num_vertical_levels, :) = 0.0_kind_phys
    call radiator%set_optical_depths( cloud_optical_depth, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_cloud_optics_values

end module musica_ccpp_tuvx_cloud_optics
