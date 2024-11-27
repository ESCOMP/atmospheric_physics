! Copyright (C) 2024 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module shr_orb_mod

  use ccpp_kinds, only: kind_phys

  implicit none

  private
  public :: shr_orb_decl, shr_orb_cosz

  ! This module contains the routines for computing the solar zenith angle and
  ! Earth-Sun distance from https://github.com/ESCOMP/CESM_share
  !
  ! This is a temporary module that will be replaced when a solution for
  ! providing information calculated from CESM shared code is implemented for
  ! CAM-SIMA.
  !
  ! Code is included in the form present in ESCOMP/CESM_share.

  ! Dependencies from other parts of the shared code
  integer, parameter           :: SHR_KIND_R8  = kind_phys
  real(SHR_KIND_R8), parameter :: pi           = 3.14159265358979323846_SHR_KIND_R8  ! pi
  
contains

  !=======================================================================

  SUBROUTINE shr_orb_decl(calday ,eccen ,mvelpp ,lambm0 ,obliqr ,delta ,eccf)

    !-------------------------------------------------------------------------------
    !
    ! Compute earth/orbit parameters using formula suggested by
    ! Duane Thresher.
    !
    !---------------------------Code history----------------------------------------
    !
    ! Original version:  Erik Kluzek
    ! Date:              Oct/1997
    !
    !-------------------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    real   (SHR_KIND_R8),intent(in)  :: calday ! Calendar day, including fraction
    real   (SHR_KIND_R8),intent(in)  :: eccen  ! Eccentricity
    real   (SHR_KIND_R8),intent(in)  :: obliqr ! Earths obliquity in radians
    real   (SHR_KIND_R8),intent(in)  :: lambm0 ! Mean long of perihelion at the
    ! vernal equinox (radians)
    real   (SHR_KIND_R8),intent(in)  :: mvelpp ! moving vernal equinox longitude
    ! of perihelion plus pi (radians)
    real   (SHR_KIND_R8),intent(out) :: delta  ! Solar declination angle in rad
    real   (SHR_KIND_R8),intent(out) :: eccf   ! Earth-sun distance factor (ie. (1/r)**2)

    !---------------------------Local variables-----------------------------
    real   (SHR_KIND_R8),parameter :: dayspy = 365.0_SHR_KIND_R8  ! days per year
    real   (SHR_KIND_R8),parameter :: ve     = 80.5_SHR_KIND_R8   ! Calday of vernal equinox
    ! assumes Jan 1 = calday 1

    real   (SHR_KIND_R8) ::   lambm  ! Lambda m, mean long of perihelion (rad)
    real   (SHR_KIND_R8) ::   lmm    ! Intermediate argument involving lambm
    real   (SHR_KIND_R8) ::   lamb   ! Lambda, the earths long of perihelion
    real   (SHR_KIND_R8) ::   invrho ! Inverse normalized sun/earth distance
    real   (SHR_KIND_R8) ::   sinl   ! Sine of lmm

    ! Compute eccentricity factor and solar declination using
    ! day value where a round day (such as 213.0) refers to 0z at
    ! Greenwich longitude.
    !
    ! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
    ! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
    ! 35:2362-2367.
    !
    ! To get the earths true longitude (position in orbit; lambda in Berger
    ! 1978) which is necessary to find the eccentricity factor and declination,
    ! must first calculate the mean longitude (lambda m in Berger 1978) at
    ! the present day.  This is done by adding to lambm0 (the mean longitude
    ! at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
    ! an increment (delta lambda m in Berger 1978) that is the number of
    ! days past or before (a negative increment) the vernal equinox divided by
    ! the days in a model year times the 2*pi radians in a complete orbit.

    lambm = lambm0 + (calday - ve)*2._SHR_KIND_R8*pi/dayspy
    lmm   = lambm  - mvelpp

    ! The earths true longitude, in radians, is then found from
    ! the formula in Berger 1978:

    sinl  = sin(lmm)
    lamb  = lambm  + eccen*(2._SHR_KIND_R8*sinl + eccen*(1.25_SHR_KIND_R8*sin(2._SHR_KIND_R8*lmm)  &
         &     + eccen*((13.0_SHR_KIND_R8/12.0_SHR_KIND_R8)*sin(3._SHR_KIND_R8*lmm) - 0.25_SHR_KIND_R8*sinl)))

    ! Using the obliquity, eccentricity, moving vernal equinox longitude of
    ! perihelion (plus), and earths true longitude, the declination (delta)
    ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
    ! rho will be used), and thus the eccentricity factor (eccf), can be
    ! calculated from formulas given in Berger 1978.

    invrho = (1._SHR_KIND_R8 + eccen*cos(lamb - mvelpp)) / (1._SHR_KIND_R8 - eccen*eccen)

    ! Set solar declination and eccentricity factor

    delta  = asin(sin(obliqr)*sin(lamb))
    eccf   = invrho*invrho

    return

  END SUBROUTINE shr_orb_decl

  !=======================================================================

  real(SHR_KIND_R8) pure FUNCTION shr_orb_cosz(jday,lat,lon,declin)

    !----------------------------------------------------------------------------
    !
    ! FUNCTION to return the cosine of the solar zenith angle.
    ! Assumes 365.0 days/year.
    !
    !--------------- Code History -----------------------------------------------
    !
    ! Original Author: Brian Kauffman
    ! Date:            Jan/98
    ! History:         adapted from statement FUNCTION in share/orb_cosz.h
    !
    !----------------------------------------------------------------------------

    real   (SHR_KIND_R8),intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
    real   (SHR_KIND_R8),intent(in) :: lat    ! Centered latitude (radians)
    real   (SHR_KIND_R8),intent(in) :: lon    ! Centered longitude (radians)
    real   (SHR_KIND_R8),intent(in) :: declin ! Solar declination (radians)

    !----------------------------------------------------------------------------

    ! perform the calculation of shr_orb_cosz
    shr_orb_cosz = sin(lat)*sin(declin) - cos(lat)*cos(declin) * &
                   cos((jday-floor(jday))*2.0_SHR_KIND_R8*pi + lon)

  END FUNCTION shr_orb_cosz

  !=======================================================================

end module shr_orb_mod