! Copyright (C) 2024 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_tuvx_cloud_optics

  use musica_ccpp_tuvx_cloud_optics

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  real, parameter :: ABS_ERROR = 1e-5

  call test_create_cloud_optics_radiator()

contains

  subroutine test_create_cloud_optics_radiator()

    use musica_util,                      only: error_t
    use musica_ccpp_tuvx_height_grid,     only: create_height_grid
    use musica_ccpp_tuvx_wavelength_grid, only: create_wavelength_grid
    use musica_tuvx_grid,                 only: grid_t
    use musica_tuvx_radiator,             only: radiator_t
    use ccpp_kinds,                       only: kind_phys

    integer, parameter        :: NUM_HOST_HEIGHT_MIDPOINTS = 2
    integer, parameter        :: NUM_HOST_HEIGHT_INTERFACES = 3
    integer, parameter        :: NUM_WAVELENGTH_MIDPOINTS = 3
    integer, parameter        :: NUM_WAVELENGTH_INTERFACES = 4
    real(kind_phys)           :: host_wavelength_interfaces(NUM_WAVELENGTH_INTERFACES) = [180.0e-9_kind_phys, 200.0e-9_kind_phys, 240.0e-9_kind_phys, 300.0e-9_kind_phys]
    real(kind_phys)           :: delta_pressure(NUM_HOST_HEIGHT_MIDPOINTS) = [100.0_kind_phys, 200.0_kind_phys]
    real(kind_phys)           :: cloud_fraction(NUM_HOST_HEIGHT_MIDPOINTS) = [0.1_kind_phys, 0.0_kind_phys]
    real(kind_phys)           :: liquid_water_content(NUM_HOST_HEIGHT_MIDPOINTS) = [0.0003_kind_phys, 0.0004_kind_phys]
    real(kind_phys)           :: reciprocal_of_gravitational_acceleration = 0.1_kind_phys
    real(kind_phys)           :: cloud_optical_depth(NUM_HOST_HEIGHT_MIDPOINTS+1, NUM_WAVELENGTH_MIDPOINTS)
    real(kind_phys)           :: expected_cloud_optical_depth(NUM_HOST_HEIGHT_MIDPOINTS+1, NUM_WAVELENGTH_MIDPOINTS) = &
      reshape([ 0.0_kind_phys, 0.14704591_kind_phys, 0.0_kind_phys, 0.0_kind_phys, 0.14704591_kind_phys, 0.0_kind_phys, 0.0_kind_phys, 0.14704591_kind_phys, 0.0_kind_phys, 0.0_kind_phys, 0.14704591_kind_phys, 0.0_kind_phys ], &
              [ NUM_HOST_HEIGHT_MIDPOINTS+1, NUM_WAVELENGTH_MIDPOINTS ])
    real(kind_phys)           :: single_scattering_albedo(NUM_HOST_HEIGHT_MIDPOINTS+1, NUM_WAVELENGTH_MIDPOINTS)
    real(kind_phys)           :: asymmetry_parameter(NUM_HOST_HEIGHT_MIDPOINTS+1, NUM_WAVELENGTH_MIDPOINTS,1)
    type(grid_t),     pointer :: height_grid => null()
    type(grid_t),     pointer :: wavelength_grid => null()
    type(radiator_t), pointer :: clouds => null()
    type(error_t)             :: error
    character(len=512)        :: errmsg
    integer                   :: errcode
    integer                   :: i

    height_grid => create_height_grid(NUM_HOST_HEIGHT_MIDPOINTS, NUM_HOST_HEIGHT_INTERFACES, &
                                      errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(height_grid))

    wavelength_grid => create_wavelength_grid(host_wavelength_interfaces, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(wavelength_grid))

    clouds => create_cloud_optics_radiator(height_grid, wavelength_grid, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(clouds))

    call set_cloud_optics_values(clouds, cloud_fraction, [ 0.0_kind_phys ], &
                                 liquid_water_content, &
                                 reciprocal_of_gravitational_acceleration, &
                                 errmsg, errcode)
    ASSERT(errcode == 1)

    call set_cloud_optics_values(clouds, cloud_fraction, delta_pressure, &
                                 [ 1.0_kind_phys ], &
                                 reciprocal_of_gravitational_acceleration, &
                                 errmsg, errcode)
    ASSERT(errcode == 1)

    call set_cloud_optics_values(clouds, [ 1.0_kind_phys ], delta_pressure, &
                                 liquid_water_content, &
                                 reciprocal_of_gravitational_acceleration, &
                                 errmsg, errcode)
    ASSERT(errcode == 1)

    call set_cloud_optics_values(clouds, cloud_fraction, delta_pressure, &
                                 liquid_water_content, &
                                 reciprocal_of_gravitational_acceleration, &
                                 errmsg, errcode)
    ASSERT(errcode == 0)

    call clouds%get_optical_depths(cloud_optical_depth, error)
    ASSERT(error%is_success())
    do i = 1, size(cloud_optical_depth, dim=1)
      do j = 1, size(cloud_optical_depth, dim=2)
        ASSERT_NEAR(cloud_optical_depth(i,j), expected_cloud_optical_depth(i,j), ABS_ERROR)
      end do
    end do

    call clouds%get_single_scattering_albedos(single_scattering_albedo, error)
    ASSERT(error%is_success())
    do i = 1, size(single_scattering_albedo, dim=1)
      do j = 1, size(single_scattering_albedo, dim=2)
        ASSERT_NEAR(single_scattering_albedo(i,j), 0.0_kind_phys, ABS_ERROR)
      end do
    end do

    call clouds%get_asymmetry_factors(asymmetry_parameter, error)
    ASSERT(error%is_success())
    do i = 1, size(asymmetry_parameter, dim=1)
      do j = 1, size(asymmetry_parameter, dim=2)
        ASSERT_NEAR(asymmetry_parameter(i,j,1), 0.0_kind_phys, ABS_ERROR)
      end do
    end do

    deallocate( height_grid )
    deallocate( wavelength_grid )
    deallocate( clouds )

  end subroutine test_create_cloud_optics_radiator

end program test_tuvx_cloud_optics
