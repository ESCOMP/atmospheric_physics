module atmos_phys_rad_utils
    ! Radiation utility functions

    implicit none
    private

    public :: is_visible

contains

   pure logical function is_visible(wavenumber)
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

end module atmos_phys_rad_utils
