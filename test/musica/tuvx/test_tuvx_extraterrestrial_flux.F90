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
    integer, parameter       :: NUM_PHOTOLYSIS_WAVELENGTH_GRID_SECTIONS = 8
    real, parameter          :: ABS_ERROR = 1e-6
    real, parameter          :: MAGNITUDE_REDUCER = 1e-15 ! Its purpose is to adjust the magnitude of the values to reduce absolute errors
    real(kind_phys)          :: wavelength_grid_interfaces(NUM_WAVELENGTH_BINS + 1) = &
        [200.0e-9_kind_phys, 220.0e-9_kind_phys, 240.0e-9_kind_phys, 260.0e-9_kind_phys, 280.0e-9_kind_phys] ! m
    real(kind_phys)          :: photolysis_wavelength_grid_interfaces(NUM_PHOTOLYSIS_WAVELENGTH_GRID_SECTIONS + 1) = &
        [200.0_kind_phys, 210.0_kind_phys, 220.0_kind_phys, 230.0_kind_phys, 240.0_kind_phys, &
         250.0_kind_phys, 260.0_kind_phys, 270.0_kind_phys, 280.0_kind_phys]                                 ! nm
    real(kind_phys)          :: extraterrestrial_flux(NUM_PHOTOLYSIS_WAVELENGTH_GRID_SECTIONS) = &
        [1.5e13_kind_phys, 1.5e13_kind_phys, 1.4e13_kind_phys, 1.4e13_kind_phys, &
         1.3e13_kind_phys, 1.2e13_kind_phys, 1.1e13_kind_phys, 1.0e13_kind_phys]
    real(kind_phys)          :: expected_extraterrestrial_flux_midpoints(NUM_WAVELENGTH_BINS) = &
        [3.0e14_kind_phys, 2.8e14_kind_phys, 2.5e14_kind_phys, 2.1e14_kind_phys]
    real(kind_phys)          :: extraterrestrial_flux_midpoints(NUM_WAVELENGTH_BINS)
    type(grid_t),    pointer :: wavelength_grid
    type(profile_t), pointer :: profile
    type(error_t)            :: error
    character(len=512)       :: errmsg
    integer                  :: errcode
    integer                  :: i

    wavelength_grid => create_wavelength_grid (wavelength_grid_interfaces, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(wavelength_grid))

    profile => create_extraterrestrial_flux_profile( wavelength_grid, wavelength_grid_interfaces, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(profile))

    call set_extraterrestrial_flux_values( profile, photolysis_wavelength_grid_interfaces, &
                                           extraterrestrial_flux, errmsg, errcode )
    ASSERT(errcode == 0)

    call profile%get_midpoint_values( extraterrestrial_flux_midpoints, error )
    ASSERT(error%is_success())

    do i = 1, size(extraterrestrial_flux_midpoints)
      ASSERT_NEAR(extraterrestrial_flux_midpoints(i) * MAGNITUDE_REDUCER, expected_extraterrestrial_flux_midpoints(i) * MAGNITUDE_REDUCER, ABS_ERROR)
    end do

    call deallocate_photolysis_wavelength_grid_interfaces()
    deallocate( wavelength_grid )
    deallocate( profile )

  end subroutine test_update_extraterrestrial_flux

end program test_tuvx_extraterrestrial_flux