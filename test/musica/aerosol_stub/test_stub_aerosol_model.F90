! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_stub_aerosol_model

  use musica_ccpp_stub_aerosol_model, only: stub_aerosol_model_t, &
                                            stub_aerosol_model_parameters_t, &
                                            STUB_AEROSOL_INVALID_DIMENSION

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  implicit none

  call test_stub_aerosol_model_create_state()
  call test_stub_aerosol_model_optical_properties()

contains

  !> @brief Test the stub_aerosol_model_create_state function
  subroutine test_stub_aerosol_model_create_state()
    use musica_ccpp_aerosol_state, only: aerosol_state_t
    use musica_ccpp_stub_aerosol_state, only: stub_aerosol_state_t
    type(stub_aerosol_model_t), pointer :: model
    type(stub_aerosol_model_parameters_t) :: parameters
    class(aerosol_state_t), pointer :: state
    character(len=512) :: error_message
    integer :: error_code
    parameters = stub_aerosol_model_parameters_t()
    model => stub_aerosol_model_t(parameters, error_message, error_code)
    ASSERT( error_code == 0 )
    state => model%create_state(1, 2, error_message, error_code)
    ASSERT( error_code == 0 )
    select type(state)
    type is (stub_aerosol_state_t)
      ASSERT( state%number_of_columns() == 1 )
      ASSERT( state%number_of_levels() == 2 )
    class default
      ASSERT(.false.)
    end select
    deallocate(model)
    deallocate(state)
  end subroutine test_stub_aerosol_model_create_state

  !> @brief Test the stub_aerosol_model_optical_properties function
  subroutine test_stub_aerosol_model_optical_properties()
    use ccpp_kinds, only: rk => kind_phys
    use musica_ccpp_grid, only: grid_t
    use musica_ccpp_aerosol_state, only: aerosol_state_t
    type(stub_aerosol_model_t), pointer :: model
    type(stub_aerosol_model_parameters_t) :: parameters
    class(aerosol_state_t), pointer :: state
    type(grid_t) :: wavelengths
    real(rk) :: optical_properties(1,2,2)
    real(rk) :: optical_properties_wrong_dims(1,2,5)
    character(len=512) :: error_message
    integer :: error_code
    parameters = stub_aerosol_model_parameters_t()
    model => stub_aerosol_model_t(parameters, error_message, error_code)
    ASSERT( error_code == 0 )
    state => model%create_state(1, 2, error_message, error_code)
    ASSERT( error_code == 0 )
    wavelengths = grid_t([0.5_rk, 1.0_rk, 1.5_rk])
    optical_properties = -999.0_rk
    call model%optical_properties(state, wavelengths, error_message, &
        error_code, extinction = optical_properties)
    ASSERT( error_code == 0 )
    ASSERT( all( optical_properties == 0.0_rk ) )
    call model%optical_properties(state, wavelengths, error_message, &
        error_code, extinction = optical_properties_wrong_dims)
    ASSERT( error_code == STUB_AEROSOL_INVALID_DIMENSION )
    optical_properties = -999.0_rk
    call model%optical_properties(state, wavelengths, error_message, &
        error_code, absorption = optical_properties)
    ASSERT( error_code == 0 )
    ASSERT( all( optical_properties == 0.0_rk ) )
    call model%optical_properties(state, wavelengths, error_message, &
        error_code, absorption = optical_properties_wrong_dims)
    ASSERT( error_code == STUB_AEROSOL_INVALID_DIMENSION )
    optical_properties = -999.0_rk
    call model%optical_properties(state, wavelengths, error_message, &
        error_code, scattering = optical_properties)
    ASSERT( error_code == 0 )
    ASSERT( all( optical_properties == 0.0_rk ) )
    call model%optical_properties(state, wavelengths, error_message, &
        error_code, scattering = optical_properties_wrong_dims)
    ASSERT( error_code == STUB_AEROSOL_INVALID_DIMENSION )
    optical_properties = -999.0_rk
    call model%optical_properties(state, wavelengths, error_message, &
        error_code, asymmetry_factor = optical_properties)
    ASSERT( error_code == 0 )
    ASSERT( all( optical_properties == 0.0_rk ) )
    call model%optical_properties(state, wavelengths, error_message, &
        error_code, asymmetry_factor = optical_properties_wrong_dims)
    ASSERT( error_code == STUB_AEROSOL_INVALID_DIMENSION )
    deallocate(model)
    deallocate(state)
  end subroutine test_stub_aerosol_model_optical_properties

end program test_stub_aerosol_model