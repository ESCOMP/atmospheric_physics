! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_species
  ! Mock module for testing
  use ccpp_kinds,   only: kind_phys
  
  implicit none
  private
  public :: micm_indices_constituent_props, micm_molar_mass_array

  integer, parameter :: micm_indices_constituent_props(4) = (/ 1, 2, 3, 4 /)
  real(kind_phys), parameter :: micm_molar_mass_array(4) = &
      (/ 200._kind_phys, 100._kind_phys, 150._kind_phys, 250._kind_phys /)
end module musica_ccpp_species