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
  public :: rk_stratiform_detrain_convective_condensate_run
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


  subroutine rk_stratiform_sedimentation_run( &
    ncol, &
    sfliq, snow_sed, &
    prec_sed, &
    prec_str, snow_str, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    real(kind_phys),    intent(in)    :: sfliq(:)     ! stratiform_rain_surface_flux_due_to_sedimentation [kg m-2 s-1]
    real(kind_phys),    intent(in)    :: snow_sed(:)     ! sfice = lwe_cloud_ice_sedimentation_rate_at_surface_due_to_microphysics [m s-1]

    ! Output arguments
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

  subroutine rk_stratiform_detrain_convective_condensate_run( &
    ncol, &
    dlf, &
    rliq, &
    prec_str, &
    tend_cldliq, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    real(kind_phys),    intent(in)    :: dlf(:,:)       ! detrainment_of_water_due_to_all_convection [kg kg-1 s-1]
    real(kind_phys),    intent(in)    :: rliq(:,:)      ! vertically_integrated_liquid_water_tendency_due_to_all_convection_to_be_applied_later_in_time_loop [kg kg-1 s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: prec_str(:)     ! stratiform_rain_and_snow_surface_flux [m s-1]

    ! Output arguments
    real(kind_phys),    intent(out)   :: tend_cldliq(:,:) ! tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1 s-1]
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    errmsg = ''
    errflg = 0

    ! Apply detrainment tendency to cloud liquid water
    tend_cldliq(:ncol,:) = dlf(:ncol,:)

    ! Accumulate precipitation and snow after
    ! reserved liquid (vertical integral) has now been used
    ! (snow contribution is zero)
    prec_str(:ncol) = prec_str(:ncol) - rliq(:ncol)

  end subroutine rk_stratiform_detrain_convective_condensate_run

  ! Determine tendencies from prognostic cloud water
  subroutine rk_stratiform_prognostic_cloud_water_tendencies_run( &
    ncol, pver, &
    dtime, &
    latvap, latice, &
    cme, fice, &
    evapheat, prfzheat, meltheat, &
    repartht, &
    evapprec, &
    ice2pr, liq2pr, &
    tend_s, &
    tend_q, &
    tend_cldice, &
    tend_cldliq, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: latvap
    real(kind_phys),    intent(in)    :: latice
    real(kind_phys),    intent(in)    :: cme(:,:)          ! net_condensation_rate_due_to_microphysics [s-1]
    real(kind_phys),    intent(in)    :: fice(:,:)         ! mass_fraction_of_ice_content_within_stratiform_cloud [fraction]
    real(kind_phys),    intent(in)    :: evapheat(:,:)     !
    real(kind_phys),    intent(in)    :: prfzheat(:,:)     !
    real(kind_phys),    intent(in)    :: meltheat(:,:)     !
    real(kind_phys),    intent(in)    :: repartht(:,:)     ! (from microphysical tend)
    real(kind_phys),    intent(in)    :: evapprec(:,:)     !
    real(kind_phys),    intent(in)    :: ice2pr(:,:)       !
    real(kind_phys),    intent(in)    :: liq2pr(:,:)       !

    ! Output arguments
    real(kind_phys),    intent(out)   :: cmeheat(:,:)      ! ... [J kg-1 s-1]
    real(kind_phys),    intent(out)   :: cmeice(:,:)       ! ... [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: cmeliq(:,:)       ! ... [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: tend_s(:,:)       ! tendency_of_dry_air_enthalpy_at_constant_pressure [J kg-1 s-1]
    real(kind_phys),    intent(out)   :: tend_q(:,:)       ! tendency_of_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: tend_cldice(:,:)  ! tendency_of_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: tend_cldliq(:,:)  ! tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1 s-1]
    character(len=512), intent(out)   :: errmsg            ! error message
    integer,            intent(out)   :: errflg            ! error flag

    errmsg = ''
    errflg = 0

    do k = 1, pver
      do i = 1, ncol
        ! Heating from cond-evap within the cloud [J kg-1 s-1]
        cmeheat(i,k)     =   cme(i,k)*(latvap + latice*fice(i,k))

        ! Rate of cond-evap of ice within the cloud [kg kg-1 s-1]
        cmeice(i,k)      =   cme(i,k)*fice(i,k)

        ! Rate of cond-evap of liq within the cloud [kg kg-1 s-1]
        cmeliq(i,k)      =   cme(i,k)*(1._kind_phys-fice(i,k))

        tend_s(i,k)      =   cmeheat(i,k) + &
                                      evapheat(i,k) + prfzheat(i,k) + meltheat(i,k) + repartht(i,k)
        tend_q(i,k)      = - cme(i,k) + evapprec(i,k)
        tend_cldice(i,k) =   cmeice(i,k) - ice2pr(i,k)
        tend_cldliq(i,k) =   cmeliq(i,k) - liq2pr(i,k)
      end do
   end do



  end subroutine rk_stratiform_prognostic_cloud_water_tendencies_run




end module rk_stratiform
