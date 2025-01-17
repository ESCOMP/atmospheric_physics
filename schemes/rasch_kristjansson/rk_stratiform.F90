! Copyright (C) 2025 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Rasch and Kristjansson prognostic cloud microphysics and CAM4 macrophysics
! CCPP-ized: Haipeng Lin, January 2025
module rk_stratiform
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: rk_stratiform_init
  public :: rk_stratiform_sedimentation_run
  public :: rk_stratiform_cloud_fractions_run
  public :: rk_stratiform_microphysics_run
  public :: rk_stratiform_prognostic_cloud_water_tendencies_run
  public :: rk_stratiform_microphysics_tendencies_run
  public :: rk_stratiform_cloud_optical_properties_run

  !

  ! temp: convect_shallow_use_shfrc() is not available so set it to
  ! false for now. it is used for UW and UNICON shallow convection schemes
  ! but is unavailable in the pbuf anyway...
  logical :: use_shfrc = .false.



contains

  ! Initialize rk_stratiform
  subroutine rk_stratiform_init(&
    errmsg, errflg)
    ! If qcwat, tcwat, lcwat are not initialized, eventually init them

  end subroutine rk_stratiform_init




end module rk_stratiform
