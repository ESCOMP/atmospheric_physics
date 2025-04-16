! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_musica_grid

  use musica_ccpp_grid, only: grid_t, GRID_INVALID

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  implicit none

  call test_grid_constructor_interfaces()
  call test_grid_constructor_interfaces_centers()
  call test_grid_constructor_evenly_spaced()

contains

  !> @brief Test the grid_constructor_interfaces function
  subroutine test_grid_constructor_interfaces()
    use ccpp_kinds, only: rk => kind_phys
    type(grid_t) :: grid
    real(rk) :: interfaces(3)
    character(len=512) :: error
    integer :: error_code
    interfaces = [1.0_rk, 2.0_rk, 3.0_rk]

    ! Check for valid interfaces
    grid = grid_t(interfaces, error, error_code)
    ASSERT( error_code == 0 )
    ASSERT( size(grid%interfaces_) == 3 )
    ASSERT_NEAR( grid%interfaces_(1), 1.0_rk, 1.0e-6_rk )
    ASSERT_NEAR( grid%interfaces_(2), 2.0_rk, 1.0e-6_rk )
    ASSERT_NEAR( grid%interfaces_(3), 3.0_rk, 1.0e-6_rk )
    ASSERT( size(grid%centers_) == 2 )
    ASSERT_NEAR( grid%centers_(1), 1.5_rk, 1.0e-6_rk )
    ASSERT_NEAR( grid%centers_(2), 2.5_rk, 1.0e-6_rk )

    ! Check for invalid interfaces
    grid = grid_t([1.0_rk, 3.0_rk, 2.0_rk], error, error_code)
    ASSERT( error_code == GRID_INVALID )
  end subroutine test_grid_constructor_interfaces

  !> @brief Test the grid_constructor_interfaces_centers function
  subroutine test_grid_constructor_interfaces_centers()
    use ccpp_kinds, only: rk => kind_phys
    type(grid_t) :: grid
    real(rk) :: interfaces(3)
    real(rk) :: centers(2)
    character(len=512) :: error
    integer :: error_code
    interfaces = [1.0_rk, 2.0_rk, 3.0_rk]
    centers = [1.6_rk, 2.7_rk]

    ! Check for valid interfaces and centers
    grid = grid_t(interfaces, centers, error, error_code)
    ASSERT( error_code == 0 )
    ASSERT( size(grid%interfaces_) == 3 )
    ASSERT_NEAR( grid%interfaces_(1), 1.0_rk, 1.0e-6_rk )
    ASSERT_NEAR( grid%interfaces_(2), 2.0_rk, 1.0e-6_rk )
    ASSERT_NEAR( grid%interfaces_(3), 3.0_rk, 1.0e-6_rk )
    ASSERT( size(grid%centers_) == 2 )
    ASSERT_NEAR( grid%centers_(1), 1.6_rk, 1.0e-6_rk )
    ASSERT_NEAR( grid%centers_(2), 2.7_rk, 1.0e-6_rk )

    ! Check for bad dimensions
    grid = grid_t(interfaces, [1.5_rk], error, error_code)
    ASSERT( error_code == GRID_INVALID)

    ! Check for invalid interfaces
    grid = grid_t([1.0_rk, 3.0_rk, 2.0_rk], centers, error, error_code)
    ASSERT( error_code == GRID_INVALID )

    ! Check for invalid centers
    grid = grid_t(interfaces, [1.5_rk, 1.9_rk], error, error_code)
    ASSERT( error_code == GRID_INVALID )
  end subroutine test_grid_constructor_interfaces_centers

  !> @brief Test the grid_constructor_evenly_spaced function
  subroutine test_grid_constructor_evenly_spaced()
    use ccpp_kinds, only: rk => kind_phys
    type(grid_t) :: grid
    real(rk) :: start, end
    integer :: number_of_sections
    character(len=512) :: error
    integer :: error_code
    start = 1.0_rk
    end = 3.0_rk
    number_of_sections = 2

    ! Check for valid start, end, and number of sections
    grid = grid_t(start, end, number_of_sections, error, error_code)
    ASSERT( error_code == 0 )
    ASSERT( size(grid%interfaces_) == 3 )
    ASSERT_NEAR( grid%interfaces_(1), 1.0_rk, 1.0e-6_rk )
    ASSERT_NEAR( grid%interfaces_(2), 2.0_rk, 1.0e-6_rk )
    ASSERT_NEAR( grid%interfaces_(3), 3.0_rk, 1.0e-6_rk )
    ASSERT( size(grid%centers_) == 2 )
    ASSERT_NEAR( grid%centers_(1), 1.5_rk, 1.0e-6_rk )
    ASSERT_NEAR( grid%centers_(2), 2.5_rk, 1.0e-6_rk )

    ! Check for invalid start, end, and delta
    grid = grid_t(3.0_rk, 1.0_rk, 1, error, error_code)
    ASSERT( error_code == GRID_INVALID )
  end subroutine test_grid_constructor_evenly_spaced

end program test_musica_grid