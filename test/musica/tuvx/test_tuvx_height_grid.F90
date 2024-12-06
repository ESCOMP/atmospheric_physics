program test_tuvx_height_grid

  use musica_ccpp_tuvx_height_grid

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  real, parameter :: ABS_ERROR = 1e-5

  call test_create_height_grid()
  call test_calculate_height_grid_values()

contains

  subroutine test_create_height_grid()
    use musica_util,      only: error_t
    use musica_tuvx_grid, only: grid_t
    use ccpp_kinds,       only: kind_phys

    integer, parameter    :: NUM_HOST_MIDPOINTS = 2
    integer, parameter    :: NUM_HOST_INTERFACES = 3
    real(kind_phys)       :: host_midpoints(NUM_HOST_MIDPOINTS) = [200.8_kind_phys, 100.6_kind_phys]
    real(kind_phys)       :: host_interfaces(NUM_HOST_INTERFACES) = [250.3_kind_phys, 150.2_kind_phys, 50.1_kind_phys]
    real(kind_phys)       :: expected_midpoints(NUM_HOST_MIDPOINTS+1) = [(100.6 + 50.1) * 0.5, 150.2, (250.3 + 200.8) * 0.5]
    real(kind_phys)       :: expected_interfaces(NUM_HOST_INTERFACES+1) = [50.1_kind_phys, 100.6_kind_phys, 200.8_kind_phys, 250.3_kind_phys]
    real(kind_phys)       :: expected_height_deltas(NUM_HOST_INTERFACES)
    real(kind_phys)       :: midpoints(NUM_HOST_MIDPOINTS+1)
    real(kind_phys)       :: interfaces(NUM_HOST_INTERFACES+1)
    real(kind_phys)       :: height_deltas(NUM_HOST_INTERFACES)
    type(grid_t), pointer :: height_grid => null()
    type(error_t)         :: error
    character(len=512)    :: errmsg
    integer               :: errcode
    integer               :: i

    expected_height_deltas(:) = expected_interfaces(2:size(expected_interfaces)) &
                              - expected_interfaces(1:size(expected_interfaces)-1)

    height_grid => create_height_grid(-1, 0, errmsg, errcode)
    ASSERT(errcode == 1)
    ASSERT(.not. associated(height_grid))

    height_grid => create_height_grid(3, 3, errmsg, errcode)
    ASSERT(errcode == 1)
    ASSERT(.not. associated(height_grid))

    height_grid => create_height_grid(NUM_HOST_MIDPOINTS, NUM_HOST_INTERFACES, &
                                      errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(height_grid))

    call set_height_grid_values(height_grid, host_midpoints, host_interfaces, &
                                height_deltas, errmsg, errcode)
    ASSERT(errcode == 0)

    ASSERT(height_grid%number_of_sections(error) == NUM_HOST_MIDPOINTS + 1)
    ASSERT(error%is_success())

    do i = 1, size(height_deltas)
      ASSERT_NEAR(height_deltas(i), expected_height_deltas(i), ABS_ERROR)
    end do

    call height_grid%get_midpoints(midpoints, error)
    ASSERT(error%is_success())
    do i = 1, size(midpoints)
      ASSERT_NEAR(midpoints(i), expected_midpoints(i), ABS_ERROR)
    end do

    call height_grid%get_edges(interfaces, error)
    ASSERT(error%is_success())
    do i = 1, size(interfaces)
      ASSERT_NEAR(interfaces(i), expected_interfaces(i), ABS_ERROR)
    end do

    deallocate( height_grid )

  end subroutine test_create_height_grid

  subroutine test_calculate_height_grid_values()
    use ccpp_kinds, only: kind_phys

    integer, parameter :: NUM_LAYERS = 2
    real(kind_phys)    :: geopotential_height_wrt_surface_at_midpoint(NUM_LAYERS) = [2000.0_kind_phys, 500.0_kind_phys]
    real(kind_phys)    :: geopotential_height_wrt_surface_at_interface(NUM_LAYERS+1) = [3000.0_kind_phys, 1000.0_kind_phys, 0.0_kind_phys]
    real(kind_phys)    :: surface_geopotential = 100.0_kind_phys
    real(kind_phys)    :: reciprocal_of_gravitational_acceleration = 10.0_kind_phys
    real(kind_phys)    :: expected_height_midpoints(NUM_LAYERS) = [3.0_kind_phys, 1.5_kind_phys]
    real(kind_phys)    :: expected_height_interfaces(NUM_LAYERS+1) = [4.0_kind_phys, 2.0_kind_phys, 1.0_kind_phys]
    real(kind_phys)    :: height_midpoints(NUM_LAYERS)
    real(kind_phys)    :: height_interfaces(NUM_LAYERS+1)
    integer            :: i

    call calculate_heights(geopotential_height_wrt_surface_at_midpoint, &
                           geopotential_height_wrt_surface_at_interface, &
                           surface_geopotential, reciprocal_of_gravitational_acceleration, &
                           height_midpoints, height_interfaces)

    do i = 1, size(height_midpoints)
      ASSERT_NEAR(height_midpoints(i), expected_height_midpoints(i), ABS_ERROR)
    end do

    do i = 1, size(height_interfaces)
      ASSERT_NEAR(height_interfaces(i), expected_height_interfaces(i), ABS_ERROR)
    end do

  end subroutine test_calculate_height_grid_values

end program test_tuvx_height_grid