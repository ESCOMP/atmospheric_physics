program test_tuvx_height_grid

  use musica_ccpp_tuvx_height_grid

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_create_height_grid()
  call test_calculate_height_grid_values()

contains

  subroutine test_create_height_grid()
    use musica_util,      only: error_t
    use musica_tuvx_grid, only: grid_t
    use ccpp_kinds,       only: kind_phys

    integer, parameter      :: NUM_HOST_MIDPOINTS = 2
    integer, parameter      :: NUM_HOST_INTERFACES = 3
    real(kind_phys), target :: host_midpoints(NUM_HOST_MIDPOINTS) = [200.8_kind_phys, 100.6_kind_phys]
    real(kind_phys), target :: host_interfaces(NUM_HOST_INTERFACES) = [250.3_kind_phys, 150.2_kind_phys, 50.1_kind_phys]
    type(grid_t), pointer   :: height_grid
    character(len=512)      :: errmsg
    integer                 :: errcode
    real(kind_phys)         :: abs_error = 1e-5
    integer                 :: i

    ! local variables
    real(kind_phys), dimension(NUM_HOST_MIDPOINTS+1)  :: midpoints
    real(kind_phys), dimension(NUM_HOST_INTERFACES+1) :: interfaces
    type(error_t)                                     :: error

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
                                errmsg, errcode)
    ASSERT(errcode == 0)

    ASSERT(height_grid%number_of_sections(error) == NUM_HOST_MIDPOINTS + 1)
    ASSERT(error%is_success())

    call height_grid%get_midpoints(midpoints, error)
    ASSERT(error%is_success())

    call height_grid%get_edges(interfaces, error)
    ASSERT(error%is_success())

    ASSERT_NEAR(midpoints(1), (100.6 + 50.1) * 0.5, abs_error)
    ASSERT_NEAR(midpoints(2), 150.2, abs_error)
    ASSERT_NEAR(midpoints(3), (250.3 + 200.8) * 0.5, abs_error)
    ASSERT_NEAR(interfaces(1), 50.1, abs_error)
    ASSERT_NEAR(interfaces(2), 100.6, abs_error)
    ASSERT_NEAR(interfaces(3), 200.8, abs_error)
    ASSERT_NEAR(interfaces(4), 250.3, abs_error)

    deallocate( height_grid )

  end subroutine test_create_height_grid

  subroutine test_calculate_height_grid_values()

    use ccpp_kinds, only: kind_phys

    integer, parameter                       :: NUM_LAYERS = 2
    real(kind_phys), dimension(NUM_LAYERS)   :: geopotential_height_wrt_surface_at_midpoint  ! m
    real(kind_phys), dimension(NUM_LAYERS+1) :: geopotential_height_wrt_surface_at_interface ! m
    real(kind_phys)                          :: surface_geopotential ! m2 s-2
    real(kind_phys)                          :: reciprocal_of_gravitational_acceleration ! s2 m-1
    real(kind_phys), dimension(NUM_LAYERS)   :: height_midpoints  ! km
    real(kind_phys), dimension(NUM_LAYERS+1) :: height_interfaces ! km

    geopotential_height_wrt_surface_at_midpoint(:) = (/ 2000.0_kind_phys, 500.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(:) = (/ 3000.0_kind_phys, 1000.0_kind_phys, 0.0_kind_phys /)
    surface_geopotential = 100.0_kind_phys
    reciprocal_of_gravitational_acceleration = 10.0_kind_phys

    call calculate_heights(geopotential_height_wrt_surface_at_midpoint, &
                           geopotential_height_wrt_surface_at_interface, &
                           surface_geopotential, reciprocal_of_gravitational_acceleration, &
                           height_midpoints, height_interfaces)

    ASSERT_NEAR(height_midpoints(1), 3.0, 1e-5)
    ASSERT_NEAR(height_midpoints(2), 1.5, 1e-5)
    ASSERT_NEAR(height_interfaces(1), 4.0, 1e-5)
    ASSERT_NEAR(height_interfaces(2), 2.0, 1e-5)
    ASSERT_NEAR(height_interfaces(3), 1.0, 1e-5)

  end subroutine test_calculate_height_grid_values

end program test_tuvx_height_grid
