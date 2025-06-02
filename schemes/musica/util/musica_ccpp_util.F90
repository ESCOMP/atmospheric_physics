! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_util

  use ccpp_kinds, only: kind_phys

  implicit none

  private
  public :: has_error_occurred, set_constants

  real(kind_phys), parameter, public :: PI = 3.14159265358979323846_kind_phys
  real(kind_phys), parameter, public :: DEGREE_TO_RADIAN = PI / 180.0_kind_phys
  real(kind_phys), public, protected :: MOLAR_MASS_DRY_AIR = -HUGE(1.0_kind_phys) ! kg mol-1

  !> Conversion factor for wavelength interfaces from meters (CAM-SIMA) to nanometers (TUV-x)
  real(kind_phys), parameter, public :: m_to_nm = 1.0e9_kind_phys

contains

  !> @brief Set constants used for MUSICA components
  !> @param[in] molar_mass_dry_air_in The molar mass of dry air (kg mol-1).
  subroutine set_constants(molar_mass_dry_air_in)

    real(kind_phys), intent(in) :: molar_mass_dry_air_in

    MOLAR_MASS_DRY_AIR = molar_mass_dry_air_in
  end subroutine set_constants

  !> @brief Evaluate a MUSICA error for failure and convert to CCPP error data
  !> @param[in] error The error code to evaluate and convert.
  !> @param[out] error_message The CCPP error message.
  !> @param[out] error_code The CCPP error code.
  !> @return True for an error, false for success.
  logical function has_error_occurred(error, error_message, error_code)
    use musica_util, only: error_t

    type(error_t),      intent(in)  :: error
    character(len=512), intent(out) :: error_message
    integer,            intent(out) :: error_code

    character(len=30) :: error_code_str

    if ( error%is_success( ) ) then
      error_code = 0
      error_message = ''
      has_error_occurred = .false.
      return
    end if
    error_code = error%code( )
    write(error_code_str, '(I30)') error%code( )
    error_message = '[MUSICA Error]: ' // error%category( ) // '[' // &
                    trim( adjustl( error_code_str ) ) // ']: ' // error%message( )
    has_error_occurred = .true.

  end function has_error_occurred

end module musica_ccpp_util
