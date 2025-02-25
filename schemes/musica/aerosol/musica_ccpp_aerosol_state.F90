! Copyright (C) 2025 National Center for Atmospheric Research,
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_aerosol_state

  implicit none
  private

  public :: aerosol_state_t

  !> Defines the state of an aerosol system according to
  !! the aerosol representation of a specific aerosol package.
  type, abstract :: aerosol_state_t
  end type aerosol_state_t

end module musica_ccpp_aerosol_state