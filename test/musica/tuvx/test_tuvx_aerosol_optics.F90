! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_tuvx_aerosol_optics

  use musica_ccpp_tuvx_aerosol_optics

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_create_aerosol_optics_radiator()
  call test_set_aerosol_optics_values_with_properties()

contains

  subroutine test_create_aerosol_optics_radiator()

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
    real(kind_phys)           :: host_wavelength_interfaces(NUM_WAVELENGTH_INTERFACES) = &
      [180.0_kind_phys, 200.0_kind_phys, 240.0_kind_phys, 300.0_kind_phys] ! nm
    real(kind_phys)           :: aerosol_optical_depth(NUM_HOST_HEIGHT_MIDPOINTS+1, NUM_WAVELENGTH_MIDPOINTS)
    real(kind_phys)           :: single_scattering_albedo(NUM_HOST_HEIGHT_MIDPOINTS+1, NUM_WAVELENGTH_MIDPOINTS)
    real(kind_phys)           :: asymmetry_parameter(NUM_HOST_HEIGHT_MIDPOINTS+1, NUM_WAVELENGTH_MIDPOINTS,1)
    type(grid_t),     pointer :: height_grid => null()
    type(grid_t),     pointer :: wavelength_grid => null()
    type(radiator_t), pointer :: aerosols => null()
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

    aerosols => create_aerosol_optics_radiator(height_grid, wavelength_grid, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(aerosols))

    call set_aerosol_optics_values(aerosols, errmsg, errcode)
    ASSERT(errcode == 0)

    call aerosols%get_optical_depths(aerosol_optical_depth, error)
    ASSERT(error%is_success())
    do i = 1, size(aerosol_optical_depth, dim=1)
      do j = 1, size(aerosol_optical_depth, dim=2)
        ASSERT_NEAR(aerosol_optical_depth(i,j), 0.0_kind_phys, ABS_ERROR)
      end do
    end do

    call aerosols%get_single_scattering_albedos(single_scattering_albedo, error)
    ASSERT(error%is_success())
    do i = 1, size(single_scattering_albedo, dim=1)
      do j = 1, size(single_scattering_albedo, dim=2)
        ASSERT_NEAR(single_scattering_albedo(i,j), 0.0_kind_phys, ABS_ERROR)
      end do
    end do

    call aerosols%get_asymmetry_factors(asymmetry_parameter, error)
    ASSERT(error%is_success())
    do i = 1, size(asymmetry_parameter, dim=1)
      do j = 1, size(asymmetry_parameter, dim=2)
        ASSERT_NEAR(asymmetry_parameter(i,j,1), 0.0_kind_phys, ABS_ERROR)
      end do
    end do

    deallocate( height_grid )
    deallocate( wavelength_grid )
    deallocate( aerosols )

  end subroutine test_create_aerosol_optics_radiator


  subroutine test_set_aerosol_optics_values_with_properties()

    use musica_util,                      only: error_t
    use musica_ccpp_tuvx_height_grid,     only: create_height_grid
    use musica_ccpp_tuvx_wavelength_grid, only: create_wavelength_grid
    use musica_tuvx_grid,                 only: grid_t
    use musica_tuvx_radiator,             only: radiator_t
    use ccpp_kinds,                       only: kind_phys

    integer, parameter        :: NUM_HOST_LAYERS = 2
    integer, parameter        :: NUM_WAVELENGTHS = 3
    real(kind_phys), parameter :: RG = 0.1_kind_phys        ! 1/g (s2 m-1)
    real(kind_phys)           :: host_wavelength_interfaces(NUM_WAVELENGTHS+1) ! nm
    real(kind_phys)           :: extinction(NUM_HOST_LAYERS, NUM_WAVELENGTHS)  ! m2 per kg air (layer 1 = top)
    real(kind_phys)           :: ssa(NUM_HOST_LAYERS, NUM_WAVELENGTHS)
    real(kind_phys)           :: asym(NUM_HOST_LAYERS, NUM_WAVELENGTHS)
    real(kind_phys)           :: delta_pressure(NUM_HOST_LAYERS)               ! Pa
    real(kind_phys)           :: od(NUM_HOST_LAYERS+1, NUM_WAVELENGTHS)
    real(kind_phys)           :: om(NUM_HOST_LAYERS+1, NUM_WAVELENGTHS)
    real(kind_phys)           :: g(NUM_HOST_LAYERS+1, NUM_WAVELENGTHS, 1)
    type(grid_t),     pointer :: height_grid => null()
    type(grid_t),     pointer :: wavelength_grid => null()
    type(radiator_t), pointer :: aerosols => null()
    type(error_t)             :: error
    character(len=512)        :: errmsg
    integer                   :: errcode, j

    host_wavelength_interfaces = [180.0_kind_phys, 200.0_kind_phys, 240.0_kind_phys, 300.0_kind_phys]
    extinction = reshape([1.0_kind_phys, 2.0_kind_phys, 3.0_kind_phys, 4.0_kind_phys, 5.0_kind_phys, &
                          6.0_kind_phys], [NUM_HOST_LAYERS, NUM_WAVELENGTHS])
    ssa(1,:) = 0.8_kind_phys
    ssa(2,:) = 0.9_kind_phys
    asym(1,:) = 0.6_kind_phys
    asym(2,:) = 0.7_kind_phys
    delta_pressure = [100.0_kind_phys, 200.0_kind_phys]
    height_grid => create_height_grid(NUM_HOST_LAYERS, NUM_HOST_LAYERS+1, errmsg, errcode)
    ASSERT(errcode == 0)
    wavelength_grid => create_wavelength_grid(host_wavelength_interfaces, errmsg, errcode)
    ASSERT(errcode == 0)
    aerosols => create_aerosol_optics_radiator(height_grid, wavelength_grid, errmsg, errcode)
    ASSERT(errcode == 0)

    call set_aerosol_optics_values(aerosols, errmsg, errcode, extinction=extinction, &
        single_scattering_albedo=ssa, asymmetry_factor=asym, delta_pressure=delta_pressure, &
        reciprocal_of_gravitational_acceleration=RG)
    ASSERT(errcode == 0)
    call aerosols%get_optical_depths(od, error)
    ASSERT(error%is_success())
    call aerosols%get_single_scattering_albedos(om, error)
    ASSERT(error%is_success())
    call aerosols%get_asymmetry_factors(g, error)
    ASSERT(error%is_success())
    ! TUV-x layer 1 (bottom) = host layer 2; extra top layer without aerosol
    do j = 1, NUM_WAVELENGTHS
      ASSERT_NEAR(od(1,j), extinction(2,j) * delta_pressure(2) * RG, 1.0e-12_kind_phys)
      ASSERT_NEAR(od(2,j), extinction(1,j) * delta_pressure(1) * RG, 1.0e-12_kind_phys)
      ASSERT(od(3,j) == 0.0_kind_phys)
      ASSERT_NEAR(om(1,j), 0.9_kind_phys, 1.0e-12_kind_phys)
      ASSERT_NEAR(om(2,j), 0.8_kind_phys, 1.0e-12_kind_phys)
      ASSERT_NEAR(g(1,j,1), 0.7_kind_phys, 1.0e-12_kind_phys)
      ASSERT_NEAR(g(2,j,1), 0.6_kind_phys, 1.0e-12_kind_phys)
    end do

    ! incomplete or mis-sized properties are errors
    call set_aerosol_optics_values(aerosols, errmsg, errcode, extinction=extinction)
    ASSERT(errcode /= 0)
    call set_aerosol_optics_values(aerosols, errmsg, errcode, extinction=extinction(:,1:2), &
        single_scattering_albedo=ssa(:,1:2), asymmetry_factor=asym(:,1:2), &
        delta_pressure=delta_pressure, reciprocal_of_gravitational_acceleration=RG)
    ASSERT(errcode /= 0)

    deallocate( height_grid )
    deallocate( wavelength_grid )
    deallocate( aerosols )

  end subroutine test_set_aerosol_optics_values_with_properties

end program test_tuvx_aerosol_optics
