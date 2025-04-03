! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_tuvx_extraterrestrial_flux

  use musica_ccpp_tuvx_extraterrestrial_flux

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_update_extraterrestrial_flux()

contains

  subroutine test_update_extraterrestrial_flux()
    use musica_util,                      only: error_t
    use musica_tuvx_grid,                 only: grid_t
    use musica_tuvx_profile,              only: profile_t
    use musica_ccpp_tuvx_wavelength_grid, only: create_wavelength_grid
    use ccpp_kinds,                       only: kind_phys

    integer, parameter       :: NUM_WAVELENGTH_BINS = 4
    integer, parameter       :: NUM_WAVELENGTH_GRID_INTERFACES = 5
    real, parameter          :: ABS_ERROR = 1e-6

    real(kind_phys)          :: host_interfaces(NUM_WAVELENGTH_GRID_INTERFACES) = &
      [140.0_kind_phys, 160.0_kind_phys, 180.0_kind_phys, 200.0_kind_phys, 220.0_kind_phys] ! nm
    real(kind_phys)          :: extraterrestrial_flux(NUM_WAVELENGTH_BINS) = &
      [1.6e13_kind_phys, 1.4e13_kind_phys, 1.2e13_kind_phys, 1.0e13_kind_phys] ! photons cm-2 s-1 nm-1
    real(kind_phys)          :: expected_extraterrestrial_flux_midpoints(NUM_WAVELENGTH_BINS) = &
      [320e12_kind_phys, 280e12_kind_phys, 240e12_kind_phys, 200e12_kind_phys]
    real(kind_phys)          :: extraterrestrial_flux_midpoints(NUM_WAVELENGTH_BINS)
    type(grid_t),    pointer :: wavelength_grid
    type(profile_t), pointer :: profile
    type(error_t)            :: error
    character(len=512)       :: errmsg
    integer                  :: errcode
    integer                  :: i

    wavelength_grid => create_wavelength_grid( host_interfaces, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(wavelength_grid))

    profile => create_extraterrestrial_flux_profile( wavelength_grid, &
                host_interfaces, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated( profile ))
    ASSERT(allocated( host_wavelength_grid_interfaces_ ))
    ASSERT(host_wavelength_grid_interfaces_%size == NUM_WAVELENGTH_GRID_INTERFACES)
    do i = 1, NUM_WAVELENGTH_GRID_INTERFACES
      ASSERT_NEAR(host_wavelength_grid_interfaces_%interfaces(i), host_interfaces(i), ABS_ERROR)
    end do
    ASSERT(allocated( tuvx_wavelength_grid_interfaces_ ))
    ASSERT(tuvx_wavelength_grid_interfaces_%size == NUM_WAVELENGTH_GRID_INTERFACES)
    do i = 1, NUM_WAVELENGTH_GRID_INTERFACES
      ASSERT_NEAR(tuvx_wavelength_grid_interfaces_%interfaces(i), host_interfaces(i), ABS_ERROR)
    end do

    call set_extraterrestrial_flux_values( profile, extraterrestrial_flux, errmsg, errcode )
    ASSERT(errcode == 0)

    call profile%get_midpoint_values( extraterrestrial_flux_midpoints, error )
    ASSERT(error%is_success())
    do i = 1, size(extraterrestrial_flux_midpoints)
      ASSERT_NEAR(extraterrestrial_flux_midpoints(i), expected_extraterrestrial_flux_midpoints(i), ABS_ERROR)
    end do

    deallocate( wavelength_grid )
    deallocate( profile )
    call cleanup_photolysis_wavelength_grid_interfaces()

  end subroutine test_update_extraterrestrial_flux

end program test_tuvx_extraterrestrial_flux