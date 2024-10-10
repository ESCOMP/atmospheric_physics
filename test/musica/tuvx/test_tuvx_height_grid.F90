program test_tuvx_height_grid

  use iso_c_binding
  use musica_ccpp_tuvx_height_grid
  use ccpp_kinds, only: kind_phys

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_create_height_grid()

contains

  subroutine test_create_height_grid()

    use musica_util,      only: error_t
    use musica_tuvx_grid, only: grid_t
    use iso_c_binding,    only: c_double
    use ccpp_kinds,       only: kind_phys

    integer, parameter      :: NUM_HOST_MIDPOINTS = 2
    integer, parameter      :: NUM_HOST_EDGES = 3
    real(kind_phys), target :: host_midpoints(NUM_HOST_MIDPOINTS) = [200.8_kind_phys, 100.6_kind_phys]
    real(kind_phys), target :: host_edges(NUM_HOST_EDGES) = [250.3_kind_phys, 150.2_kind_phys, 50.1_kind_phys]
    type(grid_t), pointer   :: height_grid
    character(len=512)      :: errmsg
    integer                 :: errcode
    real(kind_phys)         :: abs_error = 1e-5
    integer                 :: i

    ! local variables
    real(kind_phys), dimension(NUM_HOST_MIDPOINTS+1) :: midpoints
    real(kind_phys), dimension(NUM_HOST_EDGES+1) :: edges
    type(error_t) :: error

    height_grid => create_height_grid(-1, 0, errmsg, errcode)
    ASSERT(errcode == 1)
    ASSERT(.not. associated(height_grid))

    height_grid => create_height_grid(3, 3, errmsg, errcode)
    ASSERT(errcode == 1)
    ASSERT(.not. associated(height_grid))

    height_grid => create_height_grid(NUM_HOST_MIDPOINTS, NUM_HOST_EDGES, &
                                      errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(height_grid))

    call set_height_grid_values(height_grid, host_midpoints, host_edges, &
                                errmsg, errcode)
    ASSERT(errcode == 0)

    ASSERT(height_grid%number_of_sections(error) == NUM_HOST_MIDPOINTS + 1)
    ASSERT(error%is_success())

    call height_grid%get_midpoints(midpoints, error)
    ASSERT(error%is_success())

    call height_grid%get_edges(edges, error)
    ASSERT(error%is_success())

    ASSERT_NEAR(midpoints(1), (100.6 + 50.1) * 0.5, abs_error)
    ASSERT_NEAR(midpoints(2), 150.2, abs_error)
    ASSERT_NEAR(midpoints(3), (250.3 + 200.8) * 0.5, abs_error)
    ASSERT_NEAR(edges(1), 50.1, abs_error)
    ASSERT_NEAR(edges(2), 100.6, abs_error)
    ASSERT_NEAR(edges(3), 200.8, abs_error)
    ASSERT_NEAR(edges(4), 250.3, abs_error)

    deallocate( height_grid )

  end subroutine test_create_height_grid

end program test_tuvx_height_grid
