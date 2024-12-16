program test_tuvx_temperature

  use musica_ccpp_tuvx_temperature

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_update_temperature()

contains

  subroutine test_update_temperature()
    use musica_ccpp_tuvx_height_grid, only: create_height_grid
    use musica_util,                  only: error_t
    use musica_tuvx_grid,             only: grid_t
    use musica_tuvx_profile,          only: profile_t
    use ccpp_kinds,                   only: kind_phys

    integer, parameter       :: NUM_HOST_MIDPOINTS = 5
    integer, parameter       :: NUM_HOST_INTERFACES = 6
    real, parameter          :: ABS_ERROR = 1e-4
    real(kind_phys)          :: host_midpoint_temperature(NUM_HOST_MIDPOINTS) = &
                                [800.8_kind_phys, 700.7_kind_phys, 600.6_kind_phys, 500.5_kind_phys, 400.4_kind_phys]
    real(kind_phys)          :: host_surface_temperature = 300.3_kind_phys
    real(kind_phys)          :: expected_temperature_interfaces(NUM_HOST_MIDPOINTS+2) = &
       [300.3_kind_phys, 400.4_kind_phys, 500.5_kind_phys, 600.6_kind_phys, 700.7_kind_phys, 800.8_kind_phys, 800.8_kind_phys]
    real(kind_phys)          :: temperature_interfaces(NUM_HOST_MIDPOINTS+2)
    type(grid_t),    pointer :: height_grid
    type(profile_t), pointer :: profile
    character(len=512)       :: errmsg
    integer                  :: errcode
    type(error_t)            :: error
    integer                  :: i

    height_grid => create_height_grid(NUM_HOST_MIDPOINTS, NUM_HOST_INTERFACES, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(height_grid))

    profile => create_temperature_profile( height_grid, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(profile))

    call set_temperature_values( profile, host_midpoint_temperature, &
                                 host_surface_temperature, errmsg, errcode )
    ASSERT(errcode == 0)

    call profile%get_edge_values( temperature_interfaces, error)
    ASSERT(error%is_success())
    do i = 1, size(temperature_interfaces)
      ASSERT_NEAR(temperature_interfaces(i), expected_temperature_interfaces(i), ABS_ERROR)
    end do

    deallocate( profile )
    deallocate( height_grid )

  end subroutine test_update_temperature

end program test_tuvx_temperature