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
  public :: rk_stratiform_timestep_init
  public :: rk_stratiform_sedimentation_run
  ! public :: rk_stratiform_cloud_fractions_run
  ! public :: rk_stratiform_microphysics_run
  ! public :: rk_stratiform_prognostic_cloud_water_tendencies_run
  ! public :: rk_stratiform_microphysics_tendencies_run
  ! public :: rk_stratiform_cloud_optical_properties_run

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
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    errmsg = ''
    errflg = 0

  end subroutine rk_stratiform_init

  subroutine rk_stratiform_timestep_init(&
    errmsg, errflg)
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    errmsg = ''
    errflg = 0


  end subroutine rk_stratiform_timestep_init


  subroutine rk_stratiform_sedimentation_run(&
    sfliq, snow_sed, &
    prec_sed, &
    prec_str, snow_str, &
    errmsg, errflg)

    real(kind_phys),    intent(in)    :: sfliq(:)     ! stratiform_rain_surface_flux_due_to_sedimentation [kg m-2 s-1]
    real(kind_phys),    intent(in)    :: snow_sed(:)     ! lwe_cloud_ice_sedimentation_rate_at_surface_due_to_microphysics [m s-1]
    real(kind_phys),    intent(out)   :: prec_sed(:)     ! stratiform_cloud_water_surface_flux_due_to_sedimentation [m s-1]
    real(kind_phys),    intent(out)   :: prec_str(:)     ! stratiform_rain_and_snow_surface_flux [m s-1]
    real(kind_phys),    intent(out)   :: snow_str(:)     ! lwe_snow_and_cloud_ice_precipitation_rate_at_surface_due_to_microphysics [m s-1]
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    errmsg = ''
    errflg = 0

    ! Convert rain flux to precip units from mass units
    ! and create cloud water surface flux (rain + snow)
    prec_sed(:ncol) = sfliq(:ncol)/1000._kind_phys + snow_sed(:ncol)

    ! Start accumulation of precipitation and snow flux [m s-1]
    prec_str(:ncol) = 0._kind_phys + prec_sed(:ncol)
    snow_str(:ncol) = 0._kind_phys + snow_sed(:ncol)


  end subroutine rk_stratiform_sedimentation_run




end module rk_stratiform
