! Various utilities used in CAM-SIMA chemistry.
module ccpp_chem_utils

  implicit none
  private

  public :: chem_constituent_qmin

contains

  ! Returns the minimum mixing ratio for a given constituent
  ! Used to set appropriate minimum value for various chemical species at register phase.
  function chem_constituent_qmin(constituent_name) result(qmin)
    use ccpp_kinds,   only: kind_phys

    use string_utils, only: to_lower

    character(len=*), intent(in) :: constituent_name  ! Name of the chemical constituent
    real(kind_phys)              :: qmin              ! Minimum mixing ratio

    character(len=len(constituent_name)) :: name_lower

    ! Convert to lowercase for case-insensitive comparison
    name_lower = to_lower(constituent_name) ! impure

    ! Default minimum mixing ratio for chemistry species.
    qmin = 1.e-36_kind_phys

    if (index(name_lower, 'num_a') == 1) then
      ! Aerosol number density.
      qmin = 1.e-5_kind_phys
    else if (trim(name_lower) == 'o3') then
      qmin = 1.e-12_kind_phys
    else if (trim(name_lower) == 'ch4') then
      qmin = 1.e-12_kind_phys
    else if (trim(name_lower) == 'n2o') then
      qmin = 1.e-15_kind_phys
    else if (trim(name_lower) == 'cfc11' .or. trim(name_lower) == 'cfc12') then
      qmin = 1.e-20_kind_phys
    end if

  end function chem_constituent_qmin

end module ccpp_chem_utils
