! Copyright (C) 2025 University Corporation for Atmospheric Research
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
  !   note: cloud_fraction_perturbation_run calls the compute_cloud_fraction
  !         scheme run phase for perturbation with a modified rhpert_flag = .true.
  !
  ! refer to test SDF suite_rasch_kristjansson.xml for total order of operations,
  ! as the full RK-stratiform requires other schemes not included in this module.
  public :: rk_stratiform_check_qtlcwat_run
  ! -- cloud_particle_sedimentation --
  public :: rk_stratiform_sedimentation_run
  public :: rk_stratiform_detrain_convective_condensate_run
  public :: rk_stratiform_cloud_fraction_perturbation_run         ! see note.
  public :: rk_stratiform_external_forcings_run
  public :: rk_stratiform_condensate_repartioning_run
  ! -- prognostic_cloud_water --
  public :: rk_stratiform_prognostic_cloud_water_tendencies_run
  public :: rk_stratiform_cloud_optical_properties_run
  public :: rk_stratiform_save_qtlcwat_run

  !

contains

!> \section arg_table_rk_stratiform_check_qtlcwat_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_check_qtlcwat_run.html
  subroutine rk_stratiform_check_qtlcwat_run( &
    ncol, pver, &
    t, q_wv, cldice, cldliq, &
    qcwat, tcwat, lcwat, &  ! from end of last microphysics/macrophysics call.
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: t(:,:)         ! air_temperature [K]
    real(kind_phys),    intent(in)    :: q_wv(:, :)     ! adv: water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: qcwat(:,:)     ! [kg kg-1]
    real(kind_phys),    intent(inout) :: tcwat(:,:)     ! [K]
    real(kind_phys),    intent(inout) :: lcwat(:,:)     ! [kg kg-1]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    errmsg = ''
    errflg = 0

    ! Check that qcwat and tcwat were initialized - if not then initialize
    ! this is made as a separate "run" scheme so it does not have to be used in current CAM

    ! lcwat is initialized from initial conditions of cldice, cldliq
    ! TODO: check if this is always done (appears to be from physpkg.F90) or should be read from snapshot.
    if(qcwat(1,1) < 0._kind_phys .or. tcwat(1,1) < 0._kind_phys) then
      lcwat(:ncol,:) = cldice(:ncol,:) + cldliq(:ncol,:)
    endif

    ! The registry will set default negative values if not read from snapshot.
    if(qcwat(1,1) < 0._kind_phys) then
      qcwat(:ncol,:) = q_wv(:ncol,:)
    endif

    if(tcwat(1,1) < 0._kind_phys) then
      tcwat(:ncol,:) = t(:ncol,:)
    endif

  end subroutine rk_stratiform_check_qtlcwat_run

!> \section arg_table_rk_stratiform_sedimentation_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_sedimentation_run.html
  subroutine rk_stratiform_sedimentation_run( &
    ncol, &
    sfliq, snow_sed, &
    prec_sed, &
    prec_str, snow_str, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    real(kind_phys),    intent(in)    :: sfliq(:)        ! stratiform_rain_flux_at_surface_due_to_sedimentation [kg m-2 s-1]
    real(kind_phys),    intent(in)    :: snow_sed(:)     ! sfice = lwe_cloud_ice_sedimentation_rate_at_surface_due_to_microphysics [m s-1]

    ! Output arguments
    real(kind_phys),    intent(out)   :: prec_sed(:)     ! stratiform_cloud_water_surface_flux_due_to_sedimentation [m s-1]
    real(kind_phys),    intent(out)   :: prec_str(:)     ! lwe_large_scale_precipitation_rate_at_surface [m s-1]
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

!> \section arg_table_rk_stratiform_detrain_convective_condensate_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_detrain_convective_condensate_run.html
  subroutine rk_stratiform_detrain_convective_condensate_run( &
    ncol, &
    dlf, &
    rliq, &
    prec_str, &
    tend_cldliq, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    real(kind_phys),    intent(in)    :: dlf(:,:)       ! detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_all_convection [kg kg-1 s-1]
    real(kind_phys),    intent(in)    :: rliq(:)        ! vertically_integrated_cloud_liquid_water_tendency_due_to_all_convection_to_be_applied_later_in_time_loop [m s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: prec_str(:)     ! lwe_large_scale_precipitation_rate_at_surface [m s-1]

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

  ! Call perturbed cloud fraction and compute perturbation threshold criteria
  ! necessary for prognostic_cloud_water scheme
!> \section arg_table_rk_stratiform_cloud_fraction_perturbation_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_cloud_fraction_perturbation_run.html
  subroutine rk_stratiform_cloud_fraction_perturbation_run( &
    ncol, pver, &
    cappa, gravit, rair, tmelt, &
    top_lev_cloudphys, &
    pmid, ps, temp, sst, &
    q_wv, cldice, &
    phis, &
    shallowcu, deepcu, concld, & ! inputs from convective_cloud_cover
    landfrac, ocnfrac, snowh, &
    cloud, relhum, rhu00, & ! inputs from unperturbed compute_cloud_fraction
    rhdfda, & ! output for prognostic_cloud_water
    errmsg, errflg)

    ! Dependency: compute_cloud_fraction CCPPized scheme run phase.
    ! this scheme is called with an altered rhpert_flag = .true.
    ! then the outputs are combined with the "regular" output of the compute_cloud_fraction
    ! CCPP scheme to get the perturbed quantities for prognostic_cloud_water.
    use compute_cloud_fraction, only: compute_cloud_fraction_run

    ! Input arguments
    integer,         intent(in) :: ncol
    integer,         intent(in) :: pver
    real(kind_phys), intent(in) :: cappa
    real(kind_phys), intent(in) :: gravit
    real(kind_phys), intent(in) :: rair
    real(kind_phys), intent(in) :: tmelt
    integer,         intent(in) :: top_lev_cloudphys ! vertical_layer_index_of_cloud_fraction_top [index]
    real(kind_phys), intent(in) :: pmid(:, :)        ! air_pressure [Pa]
    real(kind_phys), intent(in) :: ps(:)             ! surface_air_pressure [Pa]
    real(kind_phys), intent(in) :: temp(:, :)        ! air_temperature [K]
    real(kind_phys), intent(in) :: sst(:)            ! sea_surface_temperature [K]
    real(kind_phys), intent(in) :: q_wv(:, :)        ! adv: water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys), intent(in) :: cldice(:, :)      ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys), intent(in) :: phis(:)           ! surface_geopotential [m2 s-2]
    real(kind_phys), intent(in) :: shallowcu(:, :)   ! shallow convective cloud fraction
    real(kind_phys), intent(in) :: deepcu(:, :)      ! deep convective cloud fraction
    real(kind_phys), intent(in) :: concld(:, :)      ! convective_cloud_area_fraction [fraction]
    real(kind_phys), intent(in) :: landfrac(:)       ! land_area_fraction [fraction]
    real(kind_phys), intent(in) :: ocnfrac(:)        ! ocean_area_fraction [fraction]
    real(kind_phys), intent(in) :: snowh(:)          ! lwe_surface_snow_depth_over_land [m]

    real(kind_phys), intent(in) :: cloud(:, :)       ! cloud_area_fraction [fraction]
    real(kind_phys), intent(in) :: relhum(:, :)      ! RH for prognostic cldwat [percent]

    ! Input/output arguments
    ! Note: CAM4 intentionally mutates the top level of rhu00 so this is inout here.
    real(kind_phys),  intent(inout) :: rhu00(:, :)      ! RH threshold for cloud [fraction]

    ! Output arguments
    real(kind_phys),    intent(out) :: rhdfda(:, :)     ! derivative of RH w.r.t. cloud fraction for prognostic cloud water [percent]
    character(len=512), intent(out) :: errmsg           ! error message
    integer,            intent(out) :: errflg           ! error flag

    ! Local variables
    integer :: i, k

    ! Local variables (outputs from perturbed compute_cloud_fraction)
    real(kind_phys)  :: cloud2(ncol, pver)

    ! Dummy outputs (unused)
    real(kind_phys)  :: rhcloud2(ncol, pver)
    real(kind_phys)  :: cldst2(ncol, pver)
    real(kind_phys)  :: rhu002(ncol, pver)
    real(kind_phys)  :: icecldf2(ncol, pver)
    real(kind_phys)  :: liqcldf2(ncol, pver)
    real(kind_phys)  :: relhum2(ncol, pver)

    errmsg = ''
    errflg = 0

    ! Call perturbed version of compute_cloud_fraction
    ! WARN: This is NOT CCPP-compliant!
    call compute_cloud_fraction_run( &
      ncol              = ncol,                 &
      pver              = pver,                 &
      cappa             = cappa,                &
      gravit            = gravit,               &
      rair              = rair,                 &
      tmelt             = tmelt,                &
      top_lev_cloudphys = top_lev_cloudphys,    & ! CAM4 macrophysics - top lev is 1
      pmid              = pmid(:ncol,:),        &
      ps                = ps(:ncol),            &
      temp              = temp(:ncol,:),        &
      sst               = sst(:ncol),           &
      q                 = q_wv(:ncol,:),        &
      cldice            = cldice(:ncol,:),      &
      phis              = phis(:ncol),          &
      shallowcu         = shallowcu(:ncol,:),   &
      deepcu            = deepcu(:ncol,:),      &
      concld            = concld(:ncol,:),      &
      landfrac          = landfrac(:ncol),      &
      ocnfrac           = ocnfrac(:ncol),       &
      snowh             = snowh(:ncol),         &
      rhpert_flag       = .true.,               & ! ** apply perturbation here **
      cloud             = cloud2(:ncol, :),     &
      rhcloud           = rhcloud2(:ncol, :),   &
      cldst             = cldst2(:ncol,:),      &
      rhu00             = rhu002(:ncol,:),      &
      icecldf           = icecldf2(:ncol,:),    &
      liqcldf           = liqcldf2(:ncol,:),    &
      relhum            = relhum2(:ncol,:),     &
      errmsg            = errmsg,               &
      errflg            = errflg)

    ! Compute rhdfda (derivative of RH w.r.t. cloud fraction)
    ! for use in the prognostic_cloud_water scheme.
    rhu00(:ncol,1) = 2.0_kind_phys   ! arbitrary number larger than 1 (100%)
    do k = 1, pver
      do i = 1, ncol
         if( relhum(i,k) < rhu00(i,k) ) then
            rhdfda(i,k) = 0.0_kind_phys
         elseif( relhum(i,k) >= 1.0_kind_phys ) then
            rhdfda(i,k) = 0.0_kind_phys
         else
            ! Under certain circumstances, rh+ cause cld not to changed
            ! when at an upper limit, or w/ strong subsidence
            if( ( cloud2(i,k) - cloud(i,k) ) < 1.e-4_kind_phys ) then
               rhdfda(i,k) = 0.01_kind_phys*relhum(i,k)*1.e+4_kind_phys
            else
               rhdfda(i,k) = 0.01_kind_phys*relhum(i,k)/(cloud2(i,k)-cloud(i,k))
            endif
         endif
      enddo
    enddo

  end subroutine rk_stratiform_cloud_fraction_perturbation_run


  ! Compute non-micro and non-macrophysical external forcings
  ! for computing of net condensation rate.
  ! Note: advective forcing of condensate is aggregated into liquid phase.
!> \section arg_table_rk_stratiform_external_forcings_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_external_forcings_run.html
  subroutine rk_stratiform_external_forcings_run( &
    ncol, pver, &
    dtime, &
    t, &
    q_wv, cldice, cldliq, &
    qcwat, tcwat, lcwat, &  ! from end of last physics timestep.
    qtend, ttend, ltend, &  ! output for prognostic_cloud_water
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: dtime
    real(kind_phys),    intent(in)    :: t(:,:)         ! air_temperature [K]
    real(kind_phys),    intent(in)    :: q_wv(:,:)      ! adv: water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]

    real(kind_phys),    intent(in)    :: qcwat(:,:)     ! [kg kg-1]
    real(kind_phys),    intent(in)    :: tcwat(:,:)     ! [K]
    real(kind_phys),    intent(in)    :: lcwat(:,:)     ! [kg kg-1]

    ! Output arguments (for prognostic_cloud_water)
    real(kind_phys),    intent(out)   :: qtend(:,:)     ! not due to micro/macrophysics [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: ttend(:,:)     ! not due to micro/macrophysics [K s-1]
    real(kind_phys),    intent(out)   :: ltend(:,:)     ! not due to micro/macrophysics [kg kg-1 s-1]
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local arguments
    real(kind_phys)  :: totcw(ncol, pver)

    errmsg = ''
    errflg = 0

    totcw(:ncol,:)     = cldice(:ncol,:) + cldliq(:ncol,:)

    qtend(:ncol,:pver) = 1.0_kind_phys / dtime * (q_wv  (:ncol,:pver) - qcwat(:ncol,:pver))
    ttend(:ncol,:pver) = 1.0_kind_phys / dtime * (t     (:ncol,:pver) - tcwat(:ncol,:pver))
    ltend(:ncol,:pver) = 1.0_kind_phys / dtime * (totcw (:ncol,:pver) - lcwat(:ncol,:pver))

  end subroutine rk_stratiform_external_forcings_run

  ! Repartitioning of stratiform condensate,
  ! and compute repartition heating from change in cloud ice
!> \section arg_table_rk_stratiform_condensate_repartioning_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_condensate_repartioning_run.html
  subroutine rk_stratiform_condensate_repartioning_run( &
    ncol, pver, &
    dtime, &
    latice, &
    cldice, cldliq, &
    fice, &   ! from cloud_fraction_fice
    repartht, &
    tend_cldice, &
    tend_cldliq, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: dtime
    real(kind_phys),    intent(in)    :: latice
    real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: fice(:,:)      ! mass_fraction_of_ice_content_within_stratiform_cloud [fraction]

    ! Input/output arguments

    ! Output arguments
    real(kind_phys),    intent(out)   :: repartht(:,:)     ! [J kg-1 s-1]
    real(kind_phys),    intent(out)   :: tend_cldice(:,:)  ! tendency_of_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: tend_cldliq(:,:)  ! tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1 s-1]
    character(len=512), intent(out)   :: errmsg            ! error message
    integer,            intent(out)   :: errflg            ! error flag

    ! Local arguments
    real(kind_phys)  :: totcw(ncol, pver)

    errmsg = ''
    errflg = 0

    totcw(:ncol,:) = cldice(:ncol,:) + cldliq(:ncol,:)
    tend_cldice(:ncol,:) = 1.0_kind_phys / dtime * ( totcw(:ncol,:)*fice(:ncol,:)                 - cldice(:ncol,:) )
    tend_cldliq(:ncol,:) = 1.0_kind_phys / dtime * ( totcw(:ncol,:)*(1.0_kind_phys-fice(:ncol,:)) - cldliq(:ncol,:) )

    repartht(:ncol,:pver) = latice * tend_cldice(:ncol,:pver)

  end subroutine rk_stratiform_condensate_repartioning_run

  ! Determine tendencies from prognostic cloud water
!> \section arg_table_rk_stratiform_prognostic_cloud_water_tendencies_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_prognostic_cloud_water_tendencies_run.html
  subroutine rk_stratiform_prognostic_cloud_water_tendencies_run( &
    ncol, pver, &
    dtime, &
    latvap, latice, &
    qme, fice, &
    evapheat, prfzheat, meltheat, &
    repartht, &
    evapprec, &
    ice2pr, liq2pr, &
    prec_pcw, snow_pcw, &
    prec_str, snow_str, &
    cmeheat, cmeice, cmeliq, &
    tend_s, &
    tend_q, &
    tend_cldice, &
    tend_cldliq, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: dtime
    real(kind_phys),    intent(in)    :: latvap
    real(kind_phys),    intent(in)    :: latice
    real(kind_phys),    intent(in)    :: qme(:,:)          ! net_condensation_rate_due_to_microphysics [s-1]
    real(kind_phys),    intent(in)    :: fice(:,:)         ! mass_fraction_of_ice_content_within_stratiform_cloud [fraction]
    real(kind_phys),    intent(in)    :: evapheat(:,:)     !
    real(kind_phys),    intent(in)    :: prfzheat(:,:)     !
    real(kind_phys),    intent(in)    :: meltheat(:,:)     !
    real(kind_phys),    intent(in)    :: repartht(:,:)     ! (from microphysical tend)
    real(kind_phys),    intent(in)    :: evapprec(:,:)     !
    real(kind_phys),    intent(in)    :: ice2pr(:,:)       !
    real(kind_phys),    intent(in)    :: liq2pr(:,:)       !
    real(kind_phys),    intent(in)    :: prec_pcw(:)       ! lwe_stratiform_precipitation_rate_at_surface [m s-1]
    real(kind_phys),    intent(in)    :: snow_pcw(:)       ! lwe_snow_precipitation_rate_at_surface_due_to_microphysics [m s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: prec_str(:)    ! lwe_large_scale_precipitation_rate_at_surface [m s-1]
    real(kind_phys),    intent(inout) :: snow_str(:)    ! lwe_snow_and_cloud_ice_precipitation_rate_at_surface_due_to_microphysics [m s-1]

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

    integer :: i, k
    errmsg = ''
    errflg = 0

    do k = 1, pver
      do i = 1, ncol
        ! Heating from cond-evap within the cloud [J kg-1 s-1]
        cmeheat(i,k)     =   qme(i,k)*(latvap + latice*fice(i,k))

        ! Rate of cond-evap of ice within the cloud [kg kg-1 s-1]
        cmeice(i,k)      =   qme(i,k)*fice(i,k)

        ! Rate of cond-evap of liq within the cloud [kg kg-1 s-1]
        cmeliq(i,k)      =   qme(i,k)*(1._kind_phys-fice(i,k))

        ! Tendencies from after prognostic_cloud_water...
        tend_s(i,k)      =   cmeheat(i,k) + &
                                      evapheat(i,k) + prfzheat(i,k) + meltheat(i,k) + repartht(i,k)
        tend_q(i,k)      = - qme(i,k) + evapprec(i,k)
        tend_cldice(i,k) =   cmeice(i,k) - ice2pr(i,k)
        tend_cldliq(i,k) =   cmeliq(i,k) - liq2pr(i,k)
      end do
   end do

   prec_str(:ncol) = prec_str(:ncol) + prec_pcw(:ncol)
   snow_str(:ncol) = snow_str(:ncol) + snow_pcw(:ncol)

  end subroutine rk_stratiform_prognostic_cloud_water_tendencies_run

  ! Save Q, T, cloud water at end of stratiform microphysics for use in next timestep
  ! to determine non-microphysical/macrophysical tendencies
!> \section arg_table_rk_stratiform_save_qtlcwat_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_save_qtlcwat_run.html
  subroutine rk_stratiform_save_qtlcwat_run( &
    ncol, pver, &
    t, &
    q_wv, cldice, cldliq, &
    qcwat, tcwat, lcwat, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver

    real(kind_phys),    intent(in)    :: t(:,:)         ! air_temperature [K]
    real(kind_phys),    intent(in)    :: q_wv(:, :)     ! adv: water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]

    ! Output arguments
    real(kind_phys),    intent(out)   :: qcwat(:,:)     ! [kg kg-1]
    real(kind_phys),    intent(out)   :: tcwat(:,:)     ! [K]
    real(kind_phys),    intent(out)   :: lcwat(:,:)     ! [kg kg-1]
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    integer :: k

    errmsg = ''
    errflg = 0

    do k = 1, pver
      qcwat(:ncol,k) = q_wv(:ncol,k)
      tcwat(:ncol,k) = t(:ncol,k)
      lcwat(:ncol,k) = cldice(:ncol,k) + cldliq(:ncol,k)
    enddo

  end subroutine rk_stratiform_save_qtlcwat_run

  ! Compute and save cloud water and ice particle sizes for radiation
!> \section arg_table_rk_stratiform_cloud_optical_properties_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_cloud_optical_properties_run.html
  subroutine rk_stratiform_cloud_optical_properties_run( &
    ncol, pver, &
    tmelt, &
    landfrac, icefrac, snowh, landm, &
    t, ps, pmid, &
    rel, rei, &
    errmsg, errflg)

    ! Dependency: to_be_ccppized
    use cloud_optical_properties, only: cldefr

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver

    real(kind_phys),    intent(in)    :: tmelt
    real(kind_phys),    intent(in)    :: landfrac(:)    ! land_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: icefrac(:)     ! sea_ice_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: snowh(:)       ! lwe_surface_snow_depth_over_land [m]
    real(kind_phys),    intent(in)    :: landm(:)       ! smoothed_land_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: t(:,:)         ! air_temperature [K]
    real(kind_phys),    intent(in)    :: ps(:)          ! surface_air_pressure [Pa]
    real(kind_phys),    intent(in)    :: pmid(:,:)      ! air_pressure [Pa]

    ! Output arguments
    real(kind_phys),    intent(out)   :: rel(:,:)       ! effective_radius_of_stratiform_cloud_liquid_water_particle [um]
    real(kind_phys),    intent(out)   :: rei(:,:)       ! effective_radius_of_stratiform_cloud_ice_particle [um]
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    errmsg = ''
    errflg = 0

    call cldefr( &
      ncol = ncol, &
      pver = pver, &
      tmelt = tmelt, &
      landfrac = landfrac(:ncol), &
      icefrac = icefrac(:ncol), &
      snowh = snowh(:ncol), &
      landm = landm(:ncol), &
      t = t(:ncol,:), &
      ps = ps(:ncol), &
      pmid = pmid(:ncol,:), & ! below output:
      rel = rel(:ncol,:), &
      rei = rei(:ncol,:))

  end subroutine rk_stratiform_cloud_optical_properties_run

end module rk_stratiform
