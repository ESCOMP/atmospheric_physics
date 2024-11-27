! Copyright (C) 2024 National Science Foundation - National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_musica_ccpp_util

  use ccpp_kinds, only: kind_phys
  use musica_ccpp_util

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_calculate_solar_zenith_angle_and_earth_sun_distance()

contains

  subroutine test_calculate_solar_zenith_angle_and_earth_sun_distance()
    use shr_orb_mod, only: shr_orb_decl
    use musica_util, only: error_t
    implicit none

    integer, parameter :: NUM_COLUMNS = 2
    real(kind_phys), dimension(NUM_COLUMNS) :: latitude
    real(kind_phys), dimension(NUM_COLUMNS) :: longitude
    real(kind_phys), dimension(NUM_COLUMNS) :: solar_zenith_angle
    real(kind_phys) :: calendar_day
    real(kind_phys) :: earth_eccentricity
    real(kind_phys) :: earth_obliquity
    real(kind_phys) :: perihelion_longitude
    real(kind_phys) :: moving_vernal_equinox_longitude
    real(kind_phys) :: earth_sun_distance
    character(len=512) :: errmsg
    integer :: errcode

    ! Greenwich, UK and Wellington, NZ (more or less)
    latitude = (/ 51.5_kind_phys * DEGREE_TO_RADIAN, -41.3_kind_phys * DEGREE_TO_RADIAN /)
    longitude = (/ 0.0_kind_phys * DEGREE_TO_RADIAN, 174.8_kind_phys * DEGREE_TO_RADIAN /)
    calendar_day = 183.5_kind_phys ! noon GMT on July 1
    earth_eccentricity = 0.0167_kind_phys
    earth_obliquity = 23.5_kind_phys * DEGREE_TO_RADIAN
    perihelion_longitude = 102.9_kind_phys * DEGREE_TO_RADIAN
    moving_vernal_equinox_longitude = 182.7_kind_phys * DEGREE_TO_RADIAN ! couldn't find a good value for this so I used what gave the correct SZA for Greenwich

    call calculate_solar_zenith_angle_and_earth_sun_distance(calendar_day, &
        latitude, longitude, earth_eccentricity, earth_obliquity, perihelion_longitude, &
        moving_vernal_equinox_longitude, solar_zenith_angle, earth_sun_distance, &
        errmsg, errcode)

    ! Check Earth-Sun distance is reasonable
    ASSERT_NEAR( earth_sun_distance, 1.0, 0.05 )

    ! Check solar zenith angles are reasonable (using approximate values from https://gml.noaa.gov/grad/solcalc/azel.html)
    ASSERT_NEAR( solar_zenith_angle(1), 0.8792, 5.0 * DEGREE_TO_RADIAN ) ! noon GMT in Greenwich (light)
    ASSERT( solar_zenith_angle(2) > 120.0 * DEGREE_TO_RADIAN )           ! noon GMT in Wellington (dark)

  end subroutine test_calculate_solar_zenith_angle_and_earth_sun_distance

end program test_musica_ccpp_util