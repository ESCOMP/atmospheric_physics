! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_stub_aerosol_model

  use musica_ccpp_aerosol_model, only: aerosol_model_t

  implicit none
  private

  public :: stub_aerosol_model_t, stub_aerosol_model_parameters_t, &
            STUB_AEROSOL_INVALID_DIMENSION, STUB_AEROSOL_INVALID_STATE_TYPE

  !> @brief stub_aerosol_model_parameters_t defines the parameters for the
  !! stub aerosol model. (This model assumes no aerosols are present in
  !! the atmosphere, and therefore has no parameters.)
  type :: stub_aerosol_model_parameters_t
  end type stub_aerosol_model_parameters_t

  !> @brief stub_aerosol_model_t defines the configuration of a simplified
  !! aerosol package, which assumes no aerosols are present in the
  !! atmosphere.
  type, extends(aerosol_model_t) :: stub_aerosol_model_t
  contains
    procedure :: create_state => stub_aerosol_model_create_state
    procedure :: optical_properties => stub_aerosol_model_optical_properties
  end type stub_aerosol_model_t

  interface stub_aerosol_model_t
    module procedure stub_aerosol_model_constructor
  end interface stub_aerosol_model_t

  integer, parameter :: STUB_AEROSOL_INVALID_DIMENSION = 1
  integer, parameter :: STUB_AEROSOL_INVALID_STATE_TYPE = 2

contains

  !> @brief Constructor for stub_aerosol_model_t
  !! @param parameters The parameters for the stub aerosol model.
  !! @param error_message The error message if an error occurs.
  !! @param error_code The error code if an error occurs.
  !! @return The stub aerosol model instance.
  function stub_aerosol_model_constructor(parameters, error_message, &
      error_code) result(model)
    type(stub_aerosol_model_t),             pointer     :: model
    class(stub_aerosol_model_parameters_t), intent(in)  :: parameters
    character(len=512),                     intent(out) :: error_message
    integer,                                intent(out) :: error_code
    error_message = ''
    error_code = 0
    allocate(model, stat=error_code, errmsg=error_message)
    if (error_code == 0) then
      error_message = ''
    end if
  end function stub_aerosol_model_constructor

  !> @brief Create a new aerosol state for the stub aerosol model.
  !! @param this The stub aerosol model instance.
  !! @param number_of_columns The number of columns in the model grid.
  !! @param number_of_levels The number of levels in the model grid.
  !! @param error_message The error message if an error occurs.
  !! @param error_code The error code if an error occurs.
  !! @return The aerosol state instance.
  function stub_aerosol_model_create_state(this, number_of_columns, &
      number_of_levels, error_message, error_code) result(aerosol_state)
    use musica_ccpp_aerosol_state,      only: aerosol_state_t
    use musica_ccpp_stub_aerosol_state, only: stub_aerosol_state_t
    class(stub_aerosol_model_t), intent(in)  :: this
    integer,                     intent(in)  :: number_of_columns
    integer,                     intent(in)  :: number_of_levels
    character(len=512),          intent(out) :: error_message
    integer,                     intent(out) :: error_code
    class(aerosol_state_t),      pointer     :: aerosol_state
    error_message = ''
    error_code = 0
    ! Create a new aerosol state for the stub aerosol model
    aerosol_state => stub_aerosol_state_t(number_of_columns, number_of_levels, &
        error_message, error_code)
  end function stub_aerosol_model_create_state

  !> @brief Compute the optical properties of the aerosol for the stub aerosol model.
  !! @param this The stub aerosol model instance.
  !! @param state The aerosol state instance.
  !! @param wavelengths The wavelengths at which to compute the optical properties.
  !! @param error_message The error message if an error occurs.
  !! @param error_code The error code if an error occurs.
  !! @param extinction Parameterized specific extinction (m^2/kg) [column, level, wavelength].
  !! @param absorption Parameterized specific absorption (m^2/kg) [column, level, wavelength].
  !! @param scattering Single scattering albedo (unitless) [column, level, wavelength].
  !! @param asymmetry_factor Asymmetry factor (unitless) [column, level, wavelength].
  subroutine stub_aerosol_model_optical_properties(this, state, wavelengths, &
      error_message, error_code, extinction, absorption, scattering, &
      asymmetry_factor)
    use ccpp_kinds,                     only: rk => kind_phys
    use musica_ccpp_aerosol_state,      only: aerosol_state_t
    use musica_ccpp_grid,               only: grid_t
    use musica_ccpp_stub_aerosol_state, only: stub_aerosol_state_t
    class(stub_aerosol_model_t), intent(in)  :: this
    class(aerosol_state_t),      intent(in)  :: state
    class(grid_t),               intent(in)  :: wavelengths
    character(len=512),          intent(out) :: error_message
    integer,                     intent(out) :: error_code
    real(rk), optional,          intent(out) :: extinction(:,:,:)
    real(rk), optional,          intent(out) :: absorption(:,:,:)
    real(rk), optional,          intent(out) :: scattering(:,:,:)
    real(rk), optional,          intent(out) :: asymmetry_factor(:,:,:)
    error_message = ''
    error_code = 0
    select type(state)
    class is (stub_aerosol_state_t)
      ! Compute the optical properties of the aerosol
      ! (This model assumes no aerosols are present in the atmosphere,
      ! so the optical properties are set to zero.)
      if (present(extinction)) then
        if (size(extinction, 1) /= state%number_of_columns() .or. &
            size(extinction, 2) /= state%number_of_levels() .or. &
            size(extinction, 3) /= wavelengths%number_of_sections()) then
          error_message = 'Invalid dimensions for extinction'
          error_code = STUB_AEROSOL_INVALID_DIMENSION
          return
        end if
        extinction = 0.0_rk
      end if
      if (present(absorption)) then
        if (size(absorption, 1) /= state%number_of_columns() .or. &
            size(absorption, 2) /= state%number_of_levels() .or. &
            size(absorption, 3) /= wavelengths%number_of_sections()) then
          error_message = 'Invalid dimensions for absorption'
          error_code = STUB_AEROSOL_INVALID_DIMENSION
          return
        end if
        absorption = 0.0_rk
      end if
      if (present(scattering)) then
        if (size(scattering, 1) /= state%number_of_columns() .or. &
            size(scattering, 2) /= state%number_of_levels() .or. &
            size(scattering, 3) /= wavelengths%number_of_sections()) then
          error_message = 'Invalid dimensions for scattering'
          error_code = STUB_AEROSOL_INVALID_DIMENSION
          return
        end if
        scattering = 0.0_rk
      end if
      if (present(asymmetry_factor)) then
        if (size(asymmetry_factor, 1) /= state%number_of_columns() .or. &
            size(asymmetry_factor, 2) /= state%number_of_levels() .or. &
            size(asymmetry_factor, 3) /= wavelengths%number_of_sections()) then
          error_message = 'Invalid dimensions for asymmetry factor'
          error_code = STUB_AEROSOL_INVALID_DIMENSION
          return
        end if
        asymmetry_factor = 0.0_rk
      end if
    class default
      error_message = 'Invalid aerosol state type'
      error_code = STUB_AEROSOL_INVALID_STATE_TYPE
    end select
  end subroutine stub_aerosol_model_optical_properties

end module musica_ccpp_stub_aerosol_model