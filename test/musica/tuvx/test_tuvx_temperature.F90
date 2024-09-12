program test_tuvx_temperature

  use musica_ccpp_tuvx_temperature
  use musica_ccpp_tuvx_height_grid, only: create_height_grid

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif
  
  call test_update_temperature()
  
contains

  subroutine test_update_temperature()

    use musica_util,         only: error_t
    use musica_tuvx_grid,    only: grid_t
    use musica_tuvx_profile, only: profile_t
    use ccpp_kinds,          only: kind_phys

    integer, parameter       :: NUM_HOST_MIDPOINTS = 5
    integer, parameter       :: NUM_HOST_EDGES = 3
    real(kind_phys), target  :: host_temperature_mid(NUM_HOST_MIDPOINTS)
    real(kind_phys), target  :: host_surface_temperature = 300.3_kind_phys
    
    type(grid_t), pointer    :: height_grid
    type(profile_t), pointer :: profile
    character(len=512)       :: errmsg
    integer                  :: errcode
    real(kind_phys)          :: abs_error = 1e-5
    integer                  :: i

    ! local variables
    real(kind_phys), dimension(NUM_HOST_MIDPOINTS+2) :: temperature_edges
    type(error_t) :: error

    host_temperature_mid = (/ 800.8_kind_phys, 700.7_kind_phys, 600.6_kind_phys, 500.5_kind_phys, 400.4_kind_phys /)

    height_grid => create_height_grid(NUM_HOST_MIDPOINTS, NUM_HOST_EDGES, &
                                      errmsg, errcode)
    profile => create_temperature_profile( height_grid, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(profile))

    call set_temperature_values( profile, host_temperature_mid, &
                                 host_surface_temperature, errmsg, errcode )
    ASSERT(error%is_success())

    call profile%get_edge_values( temperature_edges, error)
    ASSERT(error%is_success())

    ! ASSERT_NEAR(temperature_edges(1), (100.6 + 50.1) * 0.5, abs_error)
    ! ASSERT_NEAR(temperature_edges(2), 150.2, abs_error)
    ! ASSERT_NEAR(temperature_edges(3), (250.3 + 200.8) * 0.5, abs_error)
    ! ASSERT_NEAR(temperature_edges(4), (250.3 + 200.8) * 0.5, abs_error)
    ! ASSERT_NEAR(temperature_edges(5), (250.3 + 200.8) * 0.5, abs_error)
    ! ASSERT_NEAR(temperature_edges(6), (250.3 + 200.8) * 0.5, abs_error)
    ! ASSERT_NEAR(temperature_edges(7), (250.3 + 200.8) * 0.5, abs_error)
  
    write(*,*) "temperature_edges(1)", temperature_edges(1)
    write(*,*) "temperature_edges(2)", temperature_edges(2)
    write(*,*) "temperature_edges(3)", temperature_edges(3)
    write(*,*) "temperature_edges(4)", temperature_edges(4)
    write(*,*) "temperature_edges(5)", temperature_edges(5)
    write(*,*) "temperature_edges(6)", temperature_edges(6)
    write(*,*) "temperature_edges(7)", temperature_edges(7)

    deallocate( height_grid )
    deallocate( profile )
  
  end subroutine test_update_temperature

end program test_tuvx_temperature