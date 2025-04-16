! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program mock_host

  use musica_ccpp_aerosol_model, only: aerosol_model_t
  use musica_ccpp_aerosol_state, only: aerosol_state_t

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  implicit none

  class(aerosol_model_t), pointer :: aerosol_model => null()
  class(aerosol_state_t), pointer :: aerosol_state => null()
  integer, parameter :: NUMBER_OF_COLUMNS = 1
  integer, parameter :: NUMBER_OF_LEVELS = 2
  integer :: i_time_step

  call initialize_aerosol_model(aerosol_model, aerosol_state)
  do i_time_step = 1, 10
    ! Do something with the aerosol model and state
    call calculate_radiation(aerosol_model, aerosol_state)
    ! Maybe do something else with the aerosol model and state
  end do
  call finalize(aerosol_model, aerosol_state)

contains

  subroutine initialize_aerosol_model(aerosol_model, aerosol_state)
    use musica_ccpp_stub_aerosol_model, only: stub_aerosol_model_t, &
                                              stub_aerosol_model_parameters_t

    class(aerosol_model_t), pointer, intent(inout) :: aerosol_model
    class(aerosol_state_t), pointer, intent(inout) :: aerosol_state
    
    type(stub_aerosol_model_parameters_t) :: parameters
    character(len=512) :: error_message
    integer :: error_code

    ! initialize the aerosol model
    parameters = stub_aerosol_model_parameters_t()
    aerosol_model => stub_aerosol_model_t(parameters, error_message, error_code)
    ASSERT( error_code == 0 )

    ! create an aerosol state for the host model grid
    aerosol_state => aerosol_model%create_state(NUMBER_OF_COLUMNS, &
        NUMBER_OF_LEVELS, error_message, error_code)
    ASSERT( error_code == 0 )
  end subroutine initialize_aerosol_model

  subroutine calculate_radiation(aerosol_model, aerosol_state)
    use ccpp_kinds, only: rk => kind_phys
    use musica_ccpp_grid, only: grid_t
    use musica_ccpp_aerosol_state, only: aerosol_state_t

    class(aerosol_model_t), pointer, intent(inout) :: aerosol_model
    class(aerosol_state_t), pointer, intent(inout) :: aerosol_state
    
    integer, parameter :: SW_WAVELENGTH_BINS = 3
    integer, parameter :: LW_WAVELENGTH_BINS = 5
    type(grid_t) :: short_wavelengths, long_wavelengths
    real(rk) :: sw_ext(1,2,3), sw_abs(1,2,3), sw_sca(1,2,3), sw_asym(1,2,3)
    real(rk) :: lw_ext(1,2,5)
    character(len=512) :: error_message
    integer :: error_code

    ! define the shortwave and longwave wavelengths
    short_wavelengths = grid_t([0.5_rk, 1.0_rk, 1.5_rk, 2.0_rk])
    long_wavelengths = grid_t([5.0_rk, 10.0_rk, 15.0_rk, 20.0_rk, 25.0_rk, 30.0_rk])

    ! compute the shortwave optical properties
    call aerosol_model%optical_properties(aerosol_state, short_wavelengths, &
        error_message, error_code, &
        extinction = sw_ext, &
        absorption = sw_abs, &
        scattering = sw_sca, &
        asymmetry_factor = sw_asym)
    ASSERT( error_code == 0 )
    ! Do something with the shortwave optical properties

    ! compute the longwave optical properties
    call aerosol_model%optical_properties(aerosol_state, long_wavelengths, &
        error_message, error_code, &
        extinction = lw_ext)
    ASSERT( error_code == 0 )
    ! Do something with the longwave optical properties
  end subroutine calculate_radiation

  subroutine finalize(aerosol_model, aerosol_state)

    class(aerosol_model_t), pointer, intent(inout) :: aerosol_model
    class(aerosol_state_t), pointer, intent(inout) :: aerosol_state

    ! Clean up
    if (associated(aerosol_state)) then
      deallocate(aerosol_state)
      aerosol_state => null()
    end if
    if (associated(aerosol_model)) then
      deallocate(aerosol_model)
      aerosol_model => null()
    end if
  end subroutine finalize

end program mock_host