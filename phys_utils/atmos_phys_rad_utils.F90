module atmos_phys_rad_utils
    ! Radiation utility functions

    implicit none
    private

    public :: is_visible
    public :: stratospheric_trace_gas_profile

contains

   pure logical function is_visible(wavenumber)
     ! Returns true if the provided wavenumber is above the visible threshold
     use ccpp_kinds, only: kind_phys

     ! Wavenumber is in the visible if it is above the visible threshold
     ! wavenumber, and in the infrared if it is below the threshold
     ! This function doesn't distinquish between visible and UV.

     ! wavenumber in inverse cm (cm^-1)
     real(kind_phys), intent(in) :: wavenumber

     ! Set threshold between visible and infrared to 0.7 micron, or 14286 cm^-1
     real(kind_phys), parameter :: visible_wavenumber_threshold = 14286._kind_phys  ! cm^-1

     if (wavenumber > visible_wavenumber_threshold) then
        is_visible = .true.
     else
        is_visible = .false.
     end if

   end function is_visible

   pure subroutine stratospheric_trace_gas_profile(trop_mmr, lat, degrees_per_radian, &
        exponent_coefs, pmid, q)
      ! Zonal-mean vertical profile of a prescribed radiatively active trace gas
      ! (J. Kiehl's parameterization from CAM ghg_data::trcmix).
      ! The mixing ratio is uniform (trop_mmr) at and below a climatological tropopause,
      !    ptrop = 250 hPa - 150 hPa * cos(lat)**2,
      ! and decays above it as trop_mmr * (pmid/ptrop)**decay_exponent.
      ! The decay exponent is piecewise linear in dlat = abs(latitude in degrees),
      ! steepening the stratospheric decrease toward the poles:
      !    decay_exponent = exponent_coefs(1) + exponent_coefs(2) * dlat        (dlat <= 45)
      !    decay_exponent = exponent_coefs(3) + exponent_coefs(4) * (dlat - 45) (dlat > 45)

      !REMOVECAM:
      ! Unlike CAM trcmix, the degree conversion used here is exact
      ! (CAM uses a truncated 57.2958 deg/rad for 180/pi) and
      ! the cosine exponent is an integer (CAM uses cos(lat)**2.0),
      ! so profiles differ from CAM's in the last digits.
      ! When CAM is revised to use this subroutine it will cause answer changes.
      !REMOVECAM_END
      use ccpp_kinds, only: kind_phys

      real(kind_phys), intent(in)  :: trop_mmr           ! tropospheric mixing ratio of the gas [kg kg-1]
      real(kind_phys), intent(in)  :: lat                ! latitude [rad]
      real(kind_phys), intent(in)  :: degrees_per_radian ! 180/pi [deg rad-1]
      real(kind_phys), intent(in)  :: exponent_coefs(4)  ! decay exponent coefficients [1]
      real(kind_phys), intent(in)  :: pmid(:)            ! air pressure at layer midpoints [Pa]
      real(kind_phys), intent(out) :: q(:)               ! gas mixing ratio profile [kg kg-1]

      ! latitude of the switch between the two decay exponent branches [deg]
      real(kind_phys), parameter :: exponent_branch_lat = 45._kind_phys

      real(kind_phys) :: abs_lat_degrees ! abs(latitude) [deg]
      real(kind_phys) :: ptrop           ! climatological tropopause pressure [Pa]
      real(kind_phys) :: decay_exponent  ! stratospheric decay exponent [1]
      integer         :: k

      abs_lat_degrees = abs(lat * degrees_per_radian)
      if (abs_lat_degrees <= exponent_branch_lat) then
         decay_exponent = exponent_coefs(1) + exponent_coefs(2) * abs_lat_degrees
      else
         decay_exponent = exponent_coefs(3) + exponent_coefs(4) * (abs_lat_degrees - exponent_branch_lat)
      end if

      ptrop = 250.0e2_kind_phys - 150.0e2_kind_phys * cos(lat)**2

      do k = 1, size(pmid)
         if (pmid(k) >= ptrop) then
            q(k) = trop_mmr
         else
            q(k) = trop_mmr * (pmid(k) / ptrop)**decay_exponent
         end if
      end do

   end subroutine stratospheric_trace_gas_profile

end module atmos_phys_rad_utils
