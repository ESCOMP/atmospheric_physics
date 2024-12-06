! Copyright (C) 2024 National Center for Atmospheric Research,
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_util

  use ccpp_kinds, only: kind_phys

  implicit none

  private
  public :: has_error_occurred, calculate_solar_zenith_angle_and_earth_sun_distance

  real(kind_phys), parameter, public :: PI = 3.14159265358979323846_kind_phys
  real(kind_phys), parameter, public :: DEGREE_TO_RADIAN = PI / 180.0_kind_phys

contains

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

  !> Calculate the solar zenith angle and Earth-Sun distance
  !> @param[in] calendar_day Calendar day, including fraction
  !> @param[in] latitude Latitude in radians
  !> @param[in] longitude Longitude in radians
  !> @param[in] earth_eccentricity Earth's eccentricity factor (unitless)
  !> @param[in] earth_obliquity Earth's obliquity in radians
  !> @param[in] perihelion_longitude Earth's mean perihelion longitude at the vernal equinox (radians)
  !> @param[in] moving_vernal_equinox_longitude Earth's moving vernal equinox longitude of perihelion plus pi (radians)
  !> @param[out] solar_zenith_angle Solar zenith angle in radians
  !> @param[out] earth_sun_distance Earth-Sun distance in AU
  !> @param[out] errmsg Error message
  !> @param[out] errcode Error code
  subroutine calculate_solar_zenith_angle_and_earth_sun_distance(calendar_day, &
      latitude, longitude, earth_eccentricity, earth_obliquity, perihelion_longitude, &
      moving_vernal_equinox_longitude, solar_zenith_angle, earth_sun_distance, &
      errmsg, errcode)
    use shr_orb_mod, only: shr_orb_decl, shr_orb_cosz
    use musica_util, only: error_t

    real(kind_phys),    intent(in)  :: calendar_day
    real(kind_phys),    intent(in)  :: latitude(:)
    real(kind_phys),    intent(in)  :: longitude(:)
    real(kind_phys),    intent(in)  :: earth_eccentricity
    real(kind_phys),    intent(in)  :: earth_obliquity
    real(kind_phys),    intent(in)  :: perihelion_longitude
    real(kind_phys),    intent(in)  :: moving_vernal_equinox_longitude
    real(kind_phys),    intent(out) :: solar_zenith_angle(:)
    real(kind_phys),    intent(out) :: earth_sun_distance
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    real(kind_phys) :: delta
    integer :: i_sza

    errcode = 0
    errmsg = ''

    ! Calculate earth/orbit parameters
    call shr_orb_decl(calendar_day, earth_eccentricity, earth_obliquity, &
                      perihelion_longitude, moving_vernal_equinox_longitude, &
                      delta, earth_sun_distance)

    ! Calculate solar zenith angle
    do i_sza = 1, size(solar_zenith_angle)
      solar_zenith_angle(i_sza) = acos(shr_orb_cosz(calendar_day, latitude(i_sza), longitude(i_sza), delta))
    end do

  end subroutine calculate_solar_zenith_angle_and_earth_sun_distance

end module musica_ccpp_util
