! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Compute cloud fractions using RH threshold and other methods.
! CCPP-ized: Haipeng Lin, February 2025
module compute_cloud_fraction
  use ccpp_kinds, only: kind_phys
  private
  save

  public :: compute_cloud_fraction_init
  public :: compute_cloud_fraction_timestep_init
  public :: compute_cloud_fraction_run

  ! Namelist variables
  logical         :: cldfrc_freeze_dry ! flag for Vavrus correction
  logical         :: cldfrc_ice        ! flag to compute ice cloud fraction
  integer         :: iceopt            ! option for ice cloud closure
                                       ! 1=wang & sassen 2=schiller (iciwc)
                                       ! 3=wood & field, 4=Wilson (based on smith)
  real(kind_phys) :: rhminl            ! minimum rh for low stable clouds
  real(kind_phys) :: rhminl_adj_land   ! rhminl adjustment for snowfree land
  real(kind_phys) :: rhminh            ! minimum rh for high stable clouds
  real(kind_phys) :: premit            ! top pressure bound for mid level cloud
  real(kind_phys) :: premib            ! bottom pressure bound for mid level cloud
  real(kind_phys) :: icecrit           ! critical RH for ice clouds in Wilson and Ballard closure (smaller = more ice clouds)

  ! flags
  logical         :: inversion_cld_off ! turn off stratification-based cloud fraction (inversion_cld)

  integer         :: k700              ! model level nearest to 700 mb [index]

  ! constants
  real(kind_phys), parameter :: pnot = 1.e5_kind_phys    ! reference pressure
  real(kind_phys), parameter :: lapse = 6.5e-3_kind_phys ! U.S. Standard Atmosphere lapse rate
  real(kind_phys), parameter :: pretop = 1.0e2_kind_phys ! pressure bounding high cloud

contains
  ! Initialize cloud_fraction from namelist parameters
!> \section arg_table_compute_cloud_fraction_init Argument Table
!! \htmlinclude compute_cloud_fraction_init.html
  subroutine compute_cloud_fraction_init( &
    amIRoot, iulog, &
    pver, &
    pref_mid, &
    inversion_cld_off_in, &
    cldfrc_freeze_dry_in, cldfrc_ice_in, iceopt_in, &
    rhminl_in, rhminl_adj_land_in, &
    rhminh_in, &
    premit_in, premib_in, &
    icecrit_in, &
    errmsg, errflg)

    ! Input arguments
    logical,         intent(in)  :: amIRoot
    integer,         intent(in)  :: iulog
    integer,         intent(in)  :: pver
    real(kind_phys), intent(in)  :: pref_mid(:)          ! reference_pressure [Pa]
    logical,         intent(in)  :: inversion_cld_off_in ! flag for turn off inversion_cld?
    logical,         intent(in)  :: cldfrc_freeze_dry_in ! flag for Vavrus correction
    logical,         intent(in)  :: cldfrc_ice_in        ! flag to compute ice cloud fraction
    integer,         intent(in)  :: iceopt_in            ! option for ice cloud closure
    real(kind_phys), intent(in)  :: rhminl_in            ! minimum rh for low stable clouds
    real(kind_phys), intent(in)  :: rhminl_adj_land_in   ! rhminl adjustment for snowfree land
    real(kind_phys), intent(in)  :: rhminh_in            ! minimum rh for high stable clouds
    real(kind_phys), intent(in)  :: premit_in            ! top pressure bound for mid level cloud
    real(kind_phys), intent(in)  :: premib_in            ! bottom pressure bound for mid level cloud
    real(kind_phys), intent(in)  :: icecrit_in           ! critical RH for ice clouds in Wilson and Ballard closure

    ! Output arguments
    character(len=512), intent(out) :: errmsg            ! error message
    integer,            intent(out) :: errflg            ! error flag

    errflg = 0
    errmsg = ''

    ! Set module variables from input arguments
    cldfrc_freeze_dry = cldfrc_freeze_dry_in
    cldfrc_ice        = cldfrc_ice_in
    iceopt            = iceopt_in
    rhminl            = rhminl_in
    rhminl_adj_land   = rhminl_adj_land_in
    rhminh            = rhminh_in
    premit            = premit_in
    premib            = premib_in
    icecrit           = icecrit_in

    ! Find vertical level nearest 700 mb.
    k700 = minloc(abs(pref_mid(:) - 7.e4_kind_phys), 1)

    if(amIRoot) then
      write(iulog,*) 'compute_cloud_fraction: model level nearest 700 mb is',k700,'which is',pref_mid(k700),'pascals'
    endif

  end subroutine compute_cloud_fraction_init

  ! Timestep initialization for cloud fraction
  ! Specify the perturbation to RH to none by default.
  ! If a scheme has to perturb RH this quantity should be set by an interstitial
  ! before calling cloud_fraction a second time.
!> \section arg_table_compute_cloud_fraction_timestep_init Argument Table
!! \htmlinclude compute_cloud_fraction_timestep_init.html
  subroutine compute_cloud_fraction_timestep_init(rhpert_flag, errmsg, errflg)
    logical,            intent(out) :: rhpert_flag
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    rhpert_flag = .false.
  end subroutine compute_cloud_fraction_timestep_init

  ! Compute cloud fraction.
  !
  ! Calculates cloud fraction using a relative humidity threshold
  ! The threshold depends upon pressure, and upon the presence or absence
  ! of convection as defined by a reasonably large vertical mass flux
  ! entering that layer from below.
  !
  ! Original authors: Many, including Jim McCaa.

  ! previous call routine
  ! subroutine cldfrc(lchnk   ,ncol    , pbuf,  &
  !      pmid    ,temp    ,q       ,omga    , phis, &
  !      shfrc   ,use_shfrc, &
  !      cloud   ,rhcloud, clc     ,pdel    , &
  !      cmfmc   ,cmfmc2  ,landfrac,snowh   ,concld  ,cldst   , &
  !      ts      ,sst     ,ps      ,zdu     ,ocnfrac ,&
  !      rhu00   ,cldice  ,icecldf ,liqcldf ,relhum  ,dindex )
!> \section arg_table_compute_cloud_fraction_run Argument Table
!! \htmlinclude compute_cloud_fraction_run.html
  subroutine compute_cloud_fraction_run( &
    ncol, pver, &
    cappa, gravit, rair, tmelt, &
    top_lev_cloudphys, &
    pmid, ps, temp, sst, &
    q, cldice, &
    phis, &
    shallowcu, deepcu, concld, & ! convective cloud cover
    landfrac, ocnfrac, snowh, &
    rhpert_flag, & ! todo: decide what to do with this
    cloud, rhcloud, &
    cldst, &
    rhu00, icecldf, liqcldf, relhum, &
    errmsg, errflg)

    ! To-be-ccppized dependencies
    use wv_saturation, only: qsat, qsat_water, svp_ice_vect

    ! Arguments
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
    real(kind_phys), intent(in) :: q(:, :)           ! water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys), intent(in) :: cldice(:, :)      ! cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys), intent(in) :: phis(:)           ! surface_geopotential [m2 s-2]
    real(kind_phys), intent(in) :: shallowcu(:, :)   ! shallow convective cloud fraction
    real(kind_phys), intent(in) :: deepcu(:, :)      ! deep convective cloud fraction
    real(kind_phys), intent(in) :: concld(:, :)      ! convective_cloud_area_fraction [fraction]
    real(kind_phys), intent(in) :: landfrac(:)       ! land_area_fraction [fraction]
    real(kind_phys), intent(in) :: ocnfrac(:)        ! ocean_area_fraction [fraction]
    real(kind_phys), intent(in) :: snowh(:)          ! lwe_surface_snow_depth_over_land [m]
    logical,         intent(in) :: rhpert_flag       ! 0 or 1 to perturb rh

    ! Output arguments
    real(kind_phys), intent(out) :: cloud(:, :)      ! cloud_area_fraction [fraction]
    real(kind_phys), intent(out) :: rhcloud(:, :)    ! cloud fraction
    real(kind_phys), intent(out) :: cldst(:, :)      ! stratiform_cloud_area_fraction [fraction]
    real(kind_phys), intent(out) :: rhu00(:, :)      ! RH threshold for cloud
    real(kind_phys), intent(out) :: relhum(:, :)     ! RH for prognostic cldwat [percent]
    real(kind_phys), intent(out) :: icecldf(:, :)    ! ice cloud fraction
    real(kind_phys), intent(out) :: liqcldf(:, :)    ! liquid cloud fraction (combined into cloud)
    character(len=512), intent(out) :: errmsg        ! error message
    integer,            intent(out) :: errflg        ! error flag

    ! Local variables:
    real(kind_phys) :: cld                         ! intermediate scratch variable (low cld)
    real(kind_phys) :: dthdpmn(ncol)               ! most stable lapse rate below 750 mb
    real(kind_phys) :: dthdp                       ! lapse rate (intermediate variable)
    real(kind_phys) :: es(ncol, pver)              ! saturation vapor pressure
    real(kind_phys) :: qs(ncol, pver)              ! saturation specific humidity
    real(kind_phys) :: rhwght                      ! weighting function for rhlim transition
    real(kind_phys) :: rh(ncol, pver)              ! relative humidity
    real(kind_phys) :: rhdif                       ! intermediate scratch variable
    real(kind_phys) :: strat                       ! intermediate scratch variable
    real(kind_phys) :: theta(ncol, pver)           ! potential temperature
    real(kind_phys) :: thetas(ncol)                ! ocean surface potential temperature
    real(kind_phys) :: rhlim                       ! local rel. humidity threshold estimate
    real(kind_phys) :: coef1                       ! coefficient to convert mass flux to mb/d
    real(kind_phys) :: clrsky(ncol)                ! temporary used in random overlap calc
    real(kind_phys) :: rpdeli(ncol, pver - 1)      ! 1./(pmid(k+1)-pmid(k))
    real(kind_phys) :: rhpert                      ! the specified perturbation to rh

    logical  :: cldbnd(ncol)           ! region below high cloud boundary

    integer  :: i, ierror, k           ! column, level indices
    integer  :: kp1, ifld
    integer  :: kdthdp(ncol)
    integer  :: numkcld                ! number of levels in which to allow clouds

    !  In Cloud Ice Content variables
    real(kind_phys) :: a, b, c, as, bs, cs    ! fit parameters
    real(kind_phys) :: Kc                     ! constant for ice cloud calc (wood & field)
    real(kind_phys) :: ttmp                   ! limited temperature
    real(kind_phys) :: icicval                ! empirical iwc value
    real(kind_phys) :: rho                    ! local air density
    real(kind_phys) :: esl(ncol, pver)        ! liq sat vapor pressure
    real(kind_phys) :: esi(ncol, pver)        ! ice sat vapor pressure
    real(kind_phys) :: ncf, phi               ! Wilson and Ballard parameters

    ! Statement functions
    logical land
    land(i) = nint(landfrac(i)) == 1

    errmsg = ''
    errflg = 0

    !==================================================================================
    ! PHILOSOPHY OF PRESENT IMPLEMENTATION
    ! Modification to philosophy for ice supersaturation
    ! philosophy below is based on RH water only. This is 'liquid condensation'
    ! or liquid cloud (even though it will freeze immediately to ice)
    ! The idea is that the RH limits for condensation are strict only for
    ! water saturation
    !
    ! Ice clouds are formed by explicit parameterization of ice nucleation.
    ! Closure for ice cloud fraction is done on available cloud ice, such that
    ! the in-cloud ice content matches an empirical fit
    ! thus, icecldf = min(cldice/icicval,1) where icicval = f(temp,cldice,numice)
    ! for a first cut, icicval=f(temp) only.
    ! Combined cloud fraction is maximum overlap  cloud=max(1,max(icecldf,liqcldf))
    ! No dA/dt term for ice?
    !
    ! There are three co-existing cloud types: convective, inversion related low-level
    ! stratocumulus, and layered cloud (based on relative humidity).  Layered and
    ! stratocumulus clouds do not compete with convective cloud for which one creates
    ! the most cloud.  They contribute collectively to the total grid-box average cloud
    ! amount.  This is reflected in the way in which the total cloud amount is evaluated
    ! (a sum as opposed to a logical "or" operation)
    !
    !==================================================================================
    ! set defaults for rhu00
    rhu00(:, :) = 2.0_kind_phys
    ! define rh perturbation in order to estimate rhdfda
    rhpert = 0.01_kind_phys

    !set Wang and Sassen IWC paramters
    a = 26.87_kind_phys
    b = 0.569_kind_phys
    c = 0.002892_kind_phys
    !set schiller parameters
    as = -68.4202_kind_phys
    bs = 0.983917_kind_phys
    cs = 2.81795_kind_phys
    !set wood and field paramters...
    Kc = 75._kind_phys

    ! Evaluate potential temperature and relative humidity
    ! If not computing ice cloud fraction then hybrid RH, if MG then water RH
    if (cldfrc_ice) then
      do k = top_lev_cloudphys, pver
        call qsat_water(temp(1:ncol, k), pmid(1:ncol, k), esl(1:ncol, k), qs(1:ncol, k), ncol)
        call svp_ice_vect(temp(1:ncol, k), esi(1:ncol, k), ncol)
      end do
    else
      do k = top_lev_cloudphys, pver
        call qsat(temp(1:ncol, k), pmid(1:ncol, k), es(1:ncol, k), qs(1:ncol, k), ncol)
      end do
    end if

    cloud = 0._kind_phys
    icecldf = 0._kind_phys
    liqcldf = 0._kind_phys
    rhcloud = 0._kind_phys
    cldst = 0._kind_phys

    do k = top_lev_cloudphys, pver
      theta(:ncol, k) = temp(:ncol, k)*(pnot/pmid(:ncol, k))**cappa

      do i = 1, ncol
        if(.not. rhpert_flag) then
          ! default: no RH perturbation
          rh(i, k) = q(i, k)/qs(i, k)
        else
          rh(i, k) = q(i, k)/qs(i, k)*(1.0_kind_phys + rhpert)
        endif
        ! record relhum, rh itself will later be modified related with concld
        relhum(i, k) = rh(i, k)
      end do
    end do

    ! Initialize other temporary variables
    ierror = 0
    do i = 1, ncol
      ! Adjust thetas(i) in the presence of non-zero ocean heights.
      ! This reduces the temperature for positive heights according to a standard lapse rate.
      if (ocnfrac(i) .gt. 0.01_kind_phys) thetas(i) = &
        (sst(i) - lapse*phis(i)/gravit)*(pnot/ps(i))**cappa
      if (ocnfrac(i) .gt. 0.01_kind_phys .and. sst(i) .lt. 260._kind_phys) ierror = i
    end do
    coef1 = gravit*864.0_kind_phys    ! conversion to millibars/day

    if (ierror > 0) then
      write (iulog, *) 'COLDSST: encountered in cldfrc:', ierror, ocnfrac(ierror), sst(ierror)
    end if

    do k = top_lev_cloudphys, pver - 1
      rpdeli(:ncol, k) = 1._kind_phys/(pmid(:ncol, k + 1) - pmid(:ncol, k))
    end do

    ! shallow and deep convective cloudiness are calculated in the
    ! convective_cloud_cover CCPP scheme.
#ifndef PERGRO
    do k = top_lev_cloudphys, pver
      do i = 1, ncol
        rh(i, k) = (rh(i, k) - concld(i, k))/(1.0_kind_phys - concld(i, k))
      end do
    end do
#endif
    !==================================================================================
    !
    !          ****** Compute layer cloudiness ******
    !
    !====================================================================
    ! Begin the evaluation of layered cloud amount based on (modified) RH
    !====================================================================
    !
    numkcld = pver
    do k = top_lev_cloudphys + 1, numkcld
      kp1 = min(k + 1, pver)
      do i = 1, ncol

        ! This is now designed to apply FOR LIQUID CLOUDS (condensation > RH water)

        cldbnd(i) = pmid(i, k) .ge. pretop

        if (pmid(i, k) .ge. premib) then
          !==============================================================
          ! This is the low cloud (below premib) block
          !==============================================================
          ! enhance low cloud activation over land with no snow cover
          if (land(i) .and. (snowh(i) <= 0.000001_kind_phys)) then
            rhlim = rhminl - rhminl_adj_land
          else
            rhlim = rhminl
          end if

          rhdif = (rh(i, k) - rhlim)/(1.0_kind_phys - rhlim)
          rhcloud(i, k) = min(0.999_kind_phys, (max(rhdif, 0.0_kind_phys))**2)

          ! SJV: decrease cloud amount if very low water vapor content
          ! (thus very cold): "freeze dry"
          if (cldfrc_freeze_dry) then
            rhcloud(i, k) = rhcloud(i, k)*max(0.15_kind_phys, min(1.0_kind_phys, q(i, k)/0.0030_kind_phys))
          end if

        else if (pmid(i, k) .lt. premit) then
          !==============================================================
          ! This is the high cloud (above premit) block
          !==============================================================
          !
          rhlim = rhminh
          !
          rhdif = (rh(i, k) - rhlim)/(1.0_kind_phys - rhlim)
          rhcloud(i, k) = min(0.999_kind_phys, (max(rhdif, 0.0_kind_phys))**2)
        else
          !==============================================================
          ! This is the middle cloud block
          !==============================================================
          !
          !       linear rh threshold transition between thresholds for low & high cloud
          !
          rhwght = (premib - (max(pmid(i, k), premit)))/(premib - premit)

          if (land(i) .and. (snowh(i) <= 0.000001_kind_phys)) then
            rhlim = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_kind_phys - rhwght)
          else
            rhlim = rhminh*rhwght + rhminl*(1.0_kind_phys - rhwght)
          end if
          rhdif = (rh(i, k) - rhlim)/(1.0_kind_phys - rhlim)
          rhcloud(i, k) = min(0.999_kind_phys, (max(rhdif, 0.0_kind_phys))**2)
        end if
        !==================================================================================
        ! WE NEED TO DOCUMENT THE PURPOSE OF THIS TYPE OF CODE (ASSOCIATED WITH 2ND CALL)
        !==================================================================================
        !      !
        !      ! save rhlim to rhu00, it handles well by itself for low/high cloud
        !      !
        rhu00(i, k) = rhlim
        !==================================================================================

        if (cldfrc_ice) then

          ! Evaluate ice cloud fraction based on in-cloud ice content

          !--------ICE CLOUD OPTION 1--------Wang & Sassen 2002
          !         Evaluate desired in-cloud water content
          !               icicval = f(temp,cldice,numice)
          !         Start with a function of temperature.
          !         Wang & Sassen 2002 (JAS), based on ARM site MMCR (midlat cirrus)
          !           parameterization valid for 203-253K
          !           icival > 0 for t>195K
          if (iceopt .lt. 3) then
            if (iceopt .eq. 1) then
              ttmp = max(195._kind_phys, min(temp(i, k), 253._kind_phys)) - 273.16_kind_phys
              icicval = a + b*ttmp + c*ttmp**2._kind_phys
              !convert units
              rho = pmid(i, k)/(rair*temp(i, k))
              icicval = icicval*1.e-6_kind_phys/rho
            else
              !--------ICE CLOUD OPTION 2--------Schiller 2008 (JGR)
              !          Use a curve based on FISH measurements in
              !          tropics, mid-lats and arctic. Curve is for 180-250K (raise to 273K?)
              !          use median all flights

              ttmp = max(190._kind_phys, min(temp(i, k), 273.16_kind_phys))
              icicval = 10._kind_phys**(as*bs**ttmp + cs)
              !convert units from ppmv to kg/kg
              icicval = icicval*1.e-6_kind_phys*18._kind_phys/28.97_kind_phys
            end if
            !set icecldfraction  for OPTION 1 or OPTION2
            icecldf(i, k) = max(0._kind_phys, min(cldice(i, k)/icicval, 1._kind_phys))

          else if (iceopt .eq. 3) then

            !--------ICE CLOUD OPTION 3--------Wood & Field 2000 (JAS)
            ! eq 6: cloud fraction = 1 - exp (-K * qc/qsati)

            icecldf(i, k) = 1._kind_phys - exp(-Kc*cldice(i, k)/(qs(i, k)*(esi(i, k)/esl(i, k))))
            icecldf(i, k) = max(0._kind_phys, min(icecldf(i, k), 1._kind_phys))
          else
            !--------ICE CLOUD OPTION 4--------Wilson and ballard 1999
            ! inversion of smith....
            !       ncf = cldice / ((1-RHcrit)*qs)
            ! then a function of ncf....
            ncf = cldice(i, k)/((1._kind_phys - icecrit)*qs(i, k))
            if (ncf .le. 0._kind_phys) then
              icecldf(i, k) = 0._kind_phys
            else if (ncf .gt. 0._kind_phys .and. ncf .le. 1._kind_phys/6._kind_phys) then
              icecldf(i, k) = 0.5_kind_phys*(6._kind_phys*ncf)**(2._kind_phys/3._kind_phys)
            else if (ncf .gt. 1._kind_phys/6._kind_phys .and. ncf .lt. 1._kind_phys) then
              phi = (acos(3._kind_phys*(1._kind_phys - ncf)/2._kind_phys**(3._kind_phys/2._kind_phys)) + 4._kind_phys*3.1415927_kind_phys)/3._kind_phys
              icecldf(i, k) = (1._kind_phys - 4._kind_phys*cos(phi)*cos(phi))
            else
              icecldf(i, k) = 1._kind_phys
            end if
            icecldf(i, k) = max(0._kind_phys, min(icecldf(i, k), 1._kind_phys))
          end if
          !TEST: if ice present, icecldf=1.
          !          if (cldice(i,k).ge.1.e-8_kind_phys) then
          !             icecldf(i,k) = 0.99_kind_phys
          !          endif

               !!          if ((cldice(i,k) .gt. icicval) .or. ((cldice(i,k) .gt. 0._kind_phys) .and. (icecldf(i,k) .eq. 0._kind_phys))) then
          !          if (cldice(i,k) .gt. 1.e-8_kind_phys) then
          !             write(iulog,*) 'i,k,pmid,rho,t,cldice,icicval,icecldf,rhcloud: ', &
          !                i,k,pmid(i,k),rho,temp(i,k),cldice(i,k),icicval,icecldf(i,k),rhcloud(i,k)
          !          endif

          !         Combine ice and liquid cloud fraction assuming maximum overlap.
          ! Combined cloud fraction is maximum overlap
          !          cloud(i,k)=min(1._kind_phys,max(icecldf(i,k),rhcloud(i,k)))

          liqcldf(i, k) = (1._kind_phys - icecldf(i, k))*rhcloud(i, k)
          cloud(i, k) = liqcldf(i, k) + icecldf(i, k)
        else
          ! For RK microphysics
          cloud(i, k) = rhcloud(i, k)
        end if
      end do
    end do

    !
    ! Add in the marine strat
    ! MARINE STRATUS SHOULD BE A SPECIAL CASE OF LAYERED CLOUD
    ! CLOUD CURRENTLY CONTAINS LAYERED CLOUD DETERMINED BY RH CRITERIA
    ! TAKE THE MAXIMUM OF THE DIAGNOSED LAYERED CLOUD OR STRATOCUMULUS
    !
    !===================================================================================
    !
    !  SOME OBSERVATIONS ABOUT THE FOLLOWING SECTION OF CODE (missed in earlier look)
    !  K700 IS SET AS A CONSTANT BASED ON HYBRID COORDINATE: IT DOES NOT DEPEND ON
    !  LOCAL PRESSURE; THERE IS NO PRESSURE RAMP => LOOKS LEVEL DEPENDENT AND
    !  DISCONTINUOUS IN SPACE (I.E., STRATUS WILL END SUDDENLY WITH NO TRANSITION)
    !
    !  IT APPEARS THAT STRAT IS EVALUATED ACCORDING TO KLEIN AND HARTMANN; HOWEVER,
    !  THE ACTUAL STRATUS AMOUNT (CLDST) APPEARS TO DEPEND DIRECTLY ON THE RH BELOW
    !  THE STRONGEST PART OF THE LOW LEVEL INVERSION.
    !PJR answers: 1) the rh limitation is a physical/mathematical limitation
    !             cant have more cloud than there is RH
    !             allowed the cloud to exist two layers below the inversion
    !             because the numerics frequently make 50% relative humidity
    !             in level below the inversion which would allow no cloud
    !             2) since  the cloud is only allowed over ocean, it should
    !             be very insensitive to surface pressure (except due to
    !             spectral ringing, which also causes so many other problems
    !             I didnt worry about it.
    !
    !==================================================================================
    if (.not. inversion_cld_off) then
      !
      ! Find most stable level below 750 mb for evaluating stratus regimes
      !
      do i = 1, ncol
        ! Nothing triggers unless a stability greater than this minimum threshold is found
        dthdpmn(i) = -0.125_kind_phys
        kdthdp(i) = 0
      end do

      do k = top_lev_cloudphys + 1, pver
        do i = 1, ncol
          if (pmid(i, k) >= premib .and. ocnfrac(i) .gt. 0.01_kind_phys) then
            ! I think this is done so that dtheta/dp is in units of dg/mb (JJH)
            dthdp = 100.0_kind_phys*(theta(i, k) - theta(i, k - 1))*rpdeli(i, k - 1)
            if (dthdp < dthdpmn(i)) then
              dthdpmn(i) = dthdp
              kdthdp(i) = k     ! index of interface of max inversion
            end if
          end if
        end do
      end do

      ! Also check between the bottom layer and the surface
      ! Only perform this check if the criteria were not met above
      do i = 1, ncol
        if (kdthdp(i) .eq. 0 .and. ocnfrac(i) .gt. 0.01_kind_phys) then
          dthdp = 100.0_kind_phys*(thetas(i) - theta(i, pver))/(ps(i) - pmid(i, pver))
          if (dthdp < dthdpmn(i)) then
            dthdpmn(i) = dthdp
            kdthdp(i) = pver     ! index of interface of max inversion
          end if
        end if
      end do

      do i = 1, ncol
        if (kdthdp(i) /= 0) then
          k = kdthdp(i)
          kp1 = min(k + 1, pver)
          ! Note: strat will be zero unless ocnfrac > 0.01
          strat = min(1._kind_phys, max(0._kind_phys, ocnfrac(i)*((theta(i, k700) - thetas(i))*.057_kind_phys - .5573_kind_phys)))
          !
          ! assign the stratus to the layer just below max inversion
          ! the relative humidity changes so rapidly across the inversion
          ! that it is not safe to just look immediately below the inversion
          ! so limit the stratus cloud by rh in both layers below the inversion
          !
          cldst(i, k) = min(strat, max(rh(i, k), rh(i, kp1)))
        end if
      end do
    end if  ! .not.inversion_cld_off

    do k = top_lev_cloudphys, pver
      do i = 1, ncol
        ! which is greater; standard layered cloud amount or stratocumulus diagnosis
        cloud(i, k) = max(rhcloud(i, k), cldst(i, k))

        ! add in the contributions of convective cloud (determined separately and accounted
        ! for by modifications to the large-scale relative humidity.
        cloud(i, k) = min(cloud(i, k) + concld(i, k), 1.0_kind_phys)
      end do
    end do

  end subroutine compute_cloud_fraction_run
end module compute_cloud_fraction
