! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_aerosol_model

  implicit none
  private

  public :: aerosol_model_t

  !> Defines the configuration of any aerosol package (using
  !! any aerosol representation) based on user specification. These values are
  !! set during initialization and do not vary during the simulation.
  !!
  !! Each aerosol package (e.g., MAM, CARMA, etc) must extend the abstract
  !! aerosol_model_t class to define the details of their configuration. Any
  !! package must implement each of the deferred procedures of the abstract
  !! aerosol_model_t class, may include additional private data members and
  !! type-bound procedures, and may override functions of the abstract class.
  !!
  !! Please see the musica_ccpp_stub_aerosol_model module for an example of how the
  !! aerosol_model_t class can be extended for a specific aerosol package.
  type, abstract :: aerosol_model_t
  contains
    procedure(aerosol_model_create_state), deferred :: create_state
    procedure(aerosol_model_optical_properties), deferred :: optical_properties
  end type aerosol_model_t

  abstract interface

    !> Returns a new instance of the aerosol state for the aerosol model.
    !! The aerosol state is used to store the time-and-space varying aerosol
    !! properties for the aerosol model.
    !! @param this The aerosol model instance.
    !! @param number_of_columns The number of columns in the model grid.
    !! @param number_of_levels The number of levels in the model grid.
    !! @param error_message The error message if an error occurs.
    !! @param error_code The error code if an error occurs.
    !! @return The aerosol state instance.
    function aerosol_model_create_state(this, number_of_columns, number_of_levels, &
        error_message, error_code) result(aerosol_state)
      use musica_ccpp_aerosol_state, only: aerosol_state_t
      import :: aerosol_model_t
      class(aerosol_model_t), intent(in)  :: this
      integer,                intent(in)  :: number_of_columns
      integer,                intent(in)  :: number_of_levels
      character(len=512),     intent(out) :: error_message
      integer,                intent(out) :: error_code
      class(aerosol_state_t), pointer     :: aerosol_state
    end function aerosol_model_create_state

    !> Computes the optical properties of the aerosol for the given state and
    !! wavelengths.
    !! @param this The aerosol model instance.
    !! @param state The aerosol state instance.
    !! @param wavelengths The wavelengths at which to compute the optical properties.
    !! @param error_message The error message if an error occurs.
    !! @param error_code The error code if an error occurs.
    !! @param extinction Parameterized specific extinction (m^2/kg) [column, level, wavelength].
    !! @param absorption Parameterized specific absorption (m^2/kg) [column, level, wavelength].
    !! @param scattering Single scattering albedo (unitless) [column, level, wavelength].
    !! @param asymmetry_factor Asymmetry factor (unitless) [column, level, wavelength].
    subroutine aerosol_model_optical_properties(this, state, wavelengths, &
        error_message, error_code, extinction, absorption, scattering, &
        asymmetry_factor)
      use ccpp_kinds, only: rk => kind_phys
      use musica_ccpp_aerosol_state, only: aerosol_state_t
      use musica_ccpp_grid, only: grid_t
      import :: aerosol_model_t
      class(aerosol_model_t),      intent(in)  :: this
      class(aerosol_state_t),      intent(in)  :: state
      class(grid_t),               intent(in)  :: wavelengths
      character(len=512),          intent(out) :: error_message
      integer,                     intent(out) :: error_code
      real(rk), optional,          intent(out) :: extinction(:,:,:)
      real(rk), optional,          intent(out) :: absorption(:,:,:)
      real(rk), optional,          intent(out) :: scattering(:,:,:)
      real(rk), optional,          intent(out) :: asymmetry_factor(:,:,:)
    end subroutine aerosol_model_optical_properties

  end interface

end module musica_ccpp_aerosol_model