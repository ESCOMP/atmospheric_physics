! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_tuvx_wavelength_grid

  use musica_ccpp_tuvx_wavelength_grid

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_create_wavelength_grid()

contains

  subroutine test_create_wavelength_grid()
    use musica_util,      only: error_t
    use musica_tuvx_grid, only: grid_t
    use ccpp_kinds,       only: kind_phys

    integer, parameter    :: NUM_WAVELENGTH_GRID_MIDPOINTS = 2
    integer, parameter    :: NUM_WAVELENGTH_GRID_INTERFACES = 3
    real, parameter       :: ABS_ERROR = 1e-5
    real(kind_phys)       :: host_interfaces(NUM_WAVELENGTH_GRID_INTERFACES) = [180.0_kind_phys, 200.0_kind_phys, 240.0_kind_phys] ! nm
    real(kind_phys)       :: expected_midpoints(NUM_WAVELENGTH_GRID_MIDPOINTS) = [190.0_kind_phys, 220.0_kind_phys]
    real(kind_phys)       :: interfaces(NUM_WAVELENGTH_GRID_INTERFACES)
    real(kind_phys)       :: midpoints(NUM_WAVELENGTH_GRID_MIDPOINTS)
    type(grid_t), pointer :: wavelength_grid => null()
    character(len=512)    :: errmsg
    integer               :: errcode
    type(error_t)         :: error
    integer               :: i

    wavelength_grid => create_wavelength_grid(host_interfaces, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(wavelength_grid))

    call wavelength_grid%get_edges(interfaces, error)
    ASSERT(error%is_success())
    do i = 1, NUM_WAVELENGTH_GRID_INTERFACES
      ASSERT_NEAR(interfaces(i), host_interfaces(i), ABS_ERROR)
    end do

    call wavelength_grid%get_midpoints(midpoints, error)
    ASSERT(error%is_success())
    do i = 1, NUM_WAVELENGTH_GRID_MIDPOINTS
      ASSERT_NEAR(midpoints(i), expected_midpoints(i), ABS_ERROR)
    end do

    deallocate(wavelength_grid)

  end subroutine test_create_wavelength_grid

end program test_tuvx_wavelength_grid