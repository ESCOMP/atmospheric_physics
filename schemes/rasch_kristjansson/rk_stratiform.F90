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
  !   note: cloud_fraction_perturbation_run is NOT a fully CCPP-compliant scheme
  !         because it calls the compute_cloud_fraction scheme for perturbation.
  !         this is currently not possible to elegantly handle in the framework.
  public :: rk_stratiform_init
  public :: rk_stratiform_timestep_init
  public :: rk_stratiform_sedimentation_run
  public :: rk_stratiform_detrain_convective_condensate_run
  public :: rk_stratiform_cloud_fraction_perturbation_run         ! see note.
  ! public :: rk_stratiform_microphysics_run
  public :: rk_stratiform_prognostic_cloud_water_tendencies_run
  ! public :: rk_stratiform_ice_and_liquid_water_content_run
  ! public :: rk_stratiform_cloud_optical_properties_run
  ! public :: rk_stratiform_save_qtlcwat_run

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
    real(kind_phys),    intent(in)    :: rliq(:)        ! vertically_integrated_liquid_water_tendency_due_to_all_convection_to_be_applied_later_in_time_loop [kg kg-1 s-1]

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

  ! Call perturbed cloud fraction and compute perturbation threshold criteria
  ! necessary for prognostic_cloud_water scheme
  subroutine rk_stratiform_cloud_fraction_perturbation_run( &
    ncol, pver, &
    cappa, gravit, rair, tmelt, &
    top_lev_cloudphys, &
    pmid, ps, temp, sst, &
    q_wv, cldice, &
    phis, &
    shallowcu, deepcu, concld, & ! inputs from convective_cloud_cover
    landfrac, ocnfrac, snowh, &
    cloud, relhum, & ! inputs from unperturbed compute_cloud_fraction
    rhdfda, & ! output for prognostic_cloud_water
    errmsg, errflg)

    ! WARN: This is NOT CCPP-compliant!
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
    real(kind_phys), intent(in) :: q_wv(:, :)        ! water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys), intent(in) :: cldice(:, :)      ! cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
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
    real(kind_phys),  intent(inout) :: rhu00(:, :)      ! RH threshold for cloud

    ! Output arguments
    real(kind_phys),    intent(out) :: rhdfda(:, :)     ! derivative of RH w.r.t. cloud fraction for prognostic cloud water [percent]
    character(len=512), intent(out) :: errmsg           ! error message
    integer,            intent(out) :: errflg           ! error flag

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
      temp              = t(:ncol,:),           &
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
    rhu00(:ncol,1) = 2.0_kind_phys
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

  ! Determine tendencies from prognostic cloud water
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
    real(kind_phys),    intent(inout) :: prec_pcw(:)       ! lwe_stratiform_precipitation_rate_at_surface [m s-1]
    real(kind_phys),    intent(inout) :: snow_pcw(:)       ! lwe_snow_precipitation_rate_at_surface_due_to_microphysics [m s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: prec_str(:)    ! stratiform_rain_and_snow_surface_flux [m s-1]
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

  ! might be better suited for diagnostics ...
  ! subroutine rk_stratiform_ice_and_liquid_water_content_run( &
  !   ncol, pver, &
  !   pmid, &
  !   t, &
  !   cldice, cldliq, &
  !   rhcloud, &
  !   iwc, lwc, &
  !   icimr, icwmr, &
  !   errmsg, errflg)

  !   ! Input arguments
  !   integer,            intent(in)    :: ncol
  !   integer,            intent(in)    :: pver
  !   real(kind_phys),    intent(in)    :: pmid(:,:)      ! air_pressure [Pa]
  !   real(kind_phys),    intent(in)    :: t(:,:)         ! air_temperature [K]
  !   real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
  !   real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
  !   real(kind_phys),    intent(in)    :: rhcloud(:,:)   ! ? [fraction]

  !   ! Output arguments
  !   real(kind_phys),    intent(out)   :: iwc(:)         ! stratiform_cloud_ice_water_content [kg m-3]
  !   real(kind_phys),    intent(out)   :: lwc(:)         ! stratiform_cloud_liquid_water_content [kg m-3]
  !   real(kind_phys),    intent(out)   :: icimr(:)       ! stratiform_cloud_ice_water_mixing_ratio [kg kg-1]
  !   real(kind_phys),    intent(out)   :: icwmr(:)       ! stratiform_cloud_liquid_water_mixing_ratio [kg kg-1]
  !   character(len=512), intent(out)   :: errmsg         ! error message
  !   integer,            intent(out)   :: errflg         ! error flag

  !   integer :: i, k
  !   errmsg = ''
  !   errflg = 0

  !   do k = 1, pver
  !     do i = 1, ncol
  !        iwc(i,k)   = cldice(i,k)*pmid(i,k)/(287.15_kind_phys*t(i,k))
  !        lwc(i,k)   = cldliq(i,k)*pmid(i,k) / &
  !                     (287.15_kind_phys*t(i,k))
  !        icimr(i,k) = cldice(i,k) / max(0.01_kind_phys,rhcloud(i,k))
  !        icwmr(i,k) = cldliq(i,k) / max(0.01_kind_phys,rhcloud(i,k))
  !     end do
  !  end do

  ! end subroutine rk_stratiform_ice_and_liquid_water_content_run


end module rk_stratiform
