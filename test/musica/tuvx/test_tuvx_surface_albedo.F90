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

    integer, parameter       :: NUM_WAVELENGTH_BIN = 4
    type(grid_t),    pointer :: wavelength_grid
    type(profile_t), pointer :: profile
    real(kind_phys), target  :: host_surface_albedo = 500.5_kind_phys
    real(kind_phys)          :: wavelength_grid_interfaces(NUM_WAVELENGTH_BIN + 1) = &
        [200.0e-9_kind_phys, 210.0e-9_kind_phys, 240.0e-9_kind_phys, 300.0e-9_kind_phys, 400.0e-9_kind_phys]
    real(kind_phys)          :: surface_albedos(NUM_WAVELENGTH_BIN + 1)
    character(len=512)       :: errmsg
    integer                  :: errcode
    type(error_t)            :: error
    real(kind_phys)          :: abs_error = 1e-4
    integer                  :: i

    wavelength_grid => create_wavelength_grid(wavelength_grid_interfaces, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(associated(wavelength_grid))

    profile => create_surface_albedo_profile( wavelength_grid, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(profile))

    call set_surface_albedo_values( profile, host_surface_albedo, errmsg, errcode )
    ASSERT(errcode == 0)

    call profile%get_edge_values( surface_albedos, error)
    ASSERT(error%is_success())
    ASSERT_NEAR(surface_albedos(1), 500.5, abs_error)
    ASSERT_NEAR(surface_albedos(2), 500.5, abs_error)
    ASSERT_NEAR(surface_albedos(3), 500.5, abs_error)
    ASSERT_NEAR(surface_albedos(4), 500.5, abs_error)
    ASSERT_NEAR(surface_albedos(5), 500.5, abs_error)

    deallocate( profile )
    deallocate( wavelength_grid )

  end subroutine test_update_surface_albedo

end program test_tuvx_surface_albedo