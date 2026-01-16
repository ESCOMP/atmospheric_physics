! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_tuvx_surface_albedo

  use musica_ccpp_tuvx_surface_albedo

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_update_surface_albedo()

contains

  subroutine test_update_surface_albedo()
    use musica_ccpp_tuvx_wavelength_grid, only: create_wavelength_grid
    use musica_util,                      only: error_t
    use musica_tuvx_grid,                 only: grid_t
    use musica_tuvx_profile,              only: profile_t
    use ccpp_kinds,                       only: kind_phys

    integer, parameter       :: NUM_WAVELENGTH_BINS = 4
    real, parameter          :: ABS_ERROR = 1e-6
    real(kind_phys)          :: wavelength_grid_interfaces(NUM_WAVELENGTH_BINS + 1) = &
        [200.0e-9_kind_phys, 210.0e-9_kind_phys, 240.0e-9_kind_phys, 300.0e-9_kind_phys, 400.0e-9_kind_phys]
    real(kind_phys)          :: host_surface_albedo = 0.3_kind_phys
    real(kind_phys)          :: expected_surface_albedo_interfaces(NUM_WAVELENGTH_BINS + 1) = 0.3_kind_phys
    real(kind_phys)          :: surface_albedo_interfaces(NUM_WAVELENGTH_BINS + 1)
    type(grid_t),    pointer :: wavelength_grid
    type(profile_t), pointer :: profile
    type(error_t)            :: error
    character(len=512)       :: errmsg
    integer                  :: errcode
    integer                  :: i

    wavelength_grid => create_wavelength_grid(wavelength_grid_interfaces, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(wavelength_grid))

    profile => create_surface_albedo_profile( wavelength_grid, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(profile))

    call set_surface_albedo_values( profile, host_surface_albedo, errmsg, errcode )
    ASSERT(errcode == 0)

    call profile%get_edge_values( surface_albedo_interfaces, error)
    ASSERT(error%is_success())
    do i = 1, size(surface_albedo_interfaces)
      ASSERT_NEAR(surface_albedo_interfaces(i), expected_surface_albedo_interfaces(i), ABS_ERROR)
    end do

    deallocate( profile )
    deallocate( wavelength_grid )

  end subroutine test_update_surface_albedo

end program test_tuvx_surface_albedo