! Copyright (C) 2024-2025 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Hack shallow convective scheme.
! The main subroutine was formerly named "cmfmca", and its initialization "mfinti".
!
! Original Author: J. Hack
! CCPPized: Haipeng Lin, October 2024
module hack_convect_shallow
  use ccpp_kinds,           only: kind_phys
  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: hack_convect_shallow_init
  public :: hack_convect_shallow_run

  ! namelist variables for tuning of hack shallow convective scheme
  real(kind_phys) :: cmftau                     ! characteristic adjustment time scale [s]
  real(kind_phys) :: c0                         ! rain water autoconversion coefficient [m-1]

  ! host-model physical constants and shorthands
  real(kind_phys) :: cp                         ! specific heat of dry air [J K-1 kg-1]
  real(kind_phys) :: rgas                       ! gas constant for dry air [J K-1 kg-1]
  real(kind_phys) :: grav                       ! gravitational constant [m s-2]
  real(kind_phys) :: hlat                       ! latent heat of vaporization [J kg-1]
  real(kind_phys) :: rhoh2o                     ! density of liquid water at STP [kg m-3]

  real(kind_phys) :: rcp                        ! reciprocal of cp
  real(kind_phys) :: rgrav                      ! reciprocal of grav
  real(kind_phys) :: rhlat                      ! reciprocal of hlat

  integer         :: limcnv                     ! top vertical interface level limit for convection [index]
                                                ! derived from reference pressures to below 40 mb

  ! internal parameters
  real(kind_phys) :: betamn = 0.10_kind_phys    ! minimum overshoot parameter [???]
  real(kind_phys) :: dzmin  = 0.0_kind_phys     ! minimum convective depth for precipitation [m]
  logical         :: rlxclm = .true.            ! control for relaxing column versus cloud triplet (default: true)
                                                ! true: relaxation timescale should be applied to column as opposed to triplets individually
  real(kind_phys) :: ssfac  = 1.001_kind_phys   ! detrained air supersaturation bound [???]

  ! internal parameters for tolerance
  real(kind_phys) :: tiny   = 1.0e-36_kind_phys ! arbitrary small num in scalar transport estimates
  real(kind_phys) :: eps    = 1.0e-13_kind_phys ! machine dependent convergence criteria
  real(kind_phys) :: tpmax  = 1.50_kind_phys    ! maximum acceptable T perturbation [K]
  real(kind_phys) :: shpmax = 1.50e-3_kind_phys ! maximum acceptable Q perturbation [g g-1]

  ! diagnostic only
  logical         :: debug_verbose = .false.    ! control for debug messages


contains
  ! Initialization of moist convective mass procedure including namelist read.
!> \section arg_table_hack_convect_shallow_init Argument Table
!! \htmlinclude hack_convect_shallow_init.html
  subroutine hack_convect_shallow_init( &
    pver, &
    amIRoot, iulog, &
    cmftau_in, c0_in, &
    rair, cpair, gravit, latvap, rhoh2o_in, &
    pref_edge, &
    use_shfrc, shfrc, &
    top_lev, &
    errmsg, errflg)

    integer,            intent(in)  :: pver         ! number of vertical levels
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog        ! log output unit
    real(kind_phys),    intent(in)  :: cmftau_in    ! characteristic adjustment time scale [s]
    real(kind_phys),    intent(in)  :: c0_in        ! rain water autoconversion coefficient [m-1]
    real(kind_phys),    intent(in)  :: rair         ! gas constant for dry air [J K-1 kg-1]
    real(kind_phys),    intent(in)  :: cpair        ! specific heat of dry air [J K-1 kg-1]
    real(kind_phys),    intent(in)  :: gravit       ! gravitational constant [m s-2]
    real(kind_phys),    intent(in)  :: latvap       ! latent heat of vaporization [J kg-1]
    real(kind_phys),    intent(in)  :: rhoh2o_in    ! density of liquid water [kg m-3]
    real(kind_phys),    intent(in)  :: pref_edge(:) ! reference pressures at interface [Pa]

    logical,            intent(out) :: use_shfrc    ! this shallow scheme provides convective cloud fractions? [flag]
    real(kind_phys),    intent(out) :: shfrc(:,:)   ! (dummy) shallow convective cloud fractions calculated in-scheme [fraction]

    integer,            intent(out) :: top_lev      ! top level for cloud fraction [index]

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! local variables
    integer :: k

    errmsg = ''
    errflg = 0

    ! namelist variables
    cmftau  = cmftau_in
    c0      = c0_in

    if(amIRoot) then
      write(iulog,*) 'tuning parameters hack_convect_shallow: cmftau',cmftau
      write(iulog,*) 'tuning parameters hack_convect_shallow: c0',c0
    endif

    ! host model physical constants
    cp      = cpair
    rcp     = 1.0_kind_phys/cp
    hlat    = latvap
    rhlat   = 1.0_kind_phys/hlat
    grav    = gravit
    rgrav   = 1.0_kind_phys/gravit
    rgas    = rair
    rhoh2o  = rhoh2o_in

    ! determine limit of shallow convection: regions below 40 mb
    ! logic ported from convect_shallow_init with note that this calculation is repeated in the deep
    ! convection interface.
    if(pref_edge(1) >= 4.e3_kind_phys) then
      limcnv = 1
    else
      limcnv = pver + 1
      do k = 1, pver
        if(pref_edge(k) < 4.e3_kind_phys .and. pref_edge(k+1) >= 4.e3_kind_phys) then
          limcnv = k
        endif
      enddo
    endif

    if(amIRoot) then
      write(iulog,*) "hack_convect_shallow_init: convection will be capped at interface ", limcnv, &
                     "which is ", pref_edge(limcnv), " pascals"
    endif

    ! flags for whether this shallow convection scheme
    ! calculates and provides convective cloud fractions
    ! to convective cloud cover scheme.
    !
    ! the Hack scheme does not provide this.
    ! a dummy shfrc is provided and is never used.
    use_shfrc = .false.
    shfrc(:,:) = 0._kind_phys

    ! for Hack shallow convection (CAM4 physics), do not limit cloud fraction
    ! (extend all the way to model top)
    top_lev = 1
  end subroutine hack_convect_shallow_init

  ! Moist convective mass flux procedure.
  !
  ! If stratification is unstable to nonentraining parcel ascent,
  ! complete an adjustment making successive use of a simple cloud model
  ! consisting of three layers (sometimes referred to as a triplet)
  !
  ! Code generalized to allow specification of parcel ("updraft")
  ! properties, as well as convective transport of an arbitrary
  ! number of passive constituents (see q array).  The code
  ! is written so the water vapor field is passed independently
  ! in the calling list from the block of other transported
  ! constituents, even though as currently designed, it is the
  ! first component in the constituents field.
  !
  ! Reports tendencies in cmfdt and dq instead of updating profiles.
  !
  ! Original author: J. Hack, BAB
!> \section arg_table_hack_convect_shallow_run Argument Table
!! \htmlinclude hack_convect_shallow_run.html
  subroutine hack_convect_shallow_run( &
    ncol, pver, pcnst, &
    iulog, &
    const_props, &
    ztodt, &
    pmid, pmiddry, &
    pdel, pdeldry, rpdel, rpdeldry, &
    zm, &
    qpert_in, &
    phis, &
    pblh, &
    t, &
    q, & ! ... below are output arguments:
    dq, &
    qc_sh, &
    cmfdt, &
    cmfmc_sh, &
    cmfdqr, &
    cmfsl, &
    cmflq, &
    precc, &
    cnt_sh, &
    cnb_sh, &
    icwmr, &
    rliq_sh, &
    scheme_name, &
    flx_cnd, &
    errmsg, errflg &
  )
    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! dependency to get constituent index
    use ccpp_const_utils,          only: ccpp_const_get_idx

    ! to_be_ccppized
    use wv_saturation,             only: qsat

    ! Input arguments
    integer,         intent(in)     :: ncol               ! number of atmospheric columns
    integer,         intent(in)     :: pver               ! number of vertical levels
    integer,         intent(in)     :: pcnst              ! number of ccpp constituents
    integer,         intent(in)     :: iulog              ! log output unit
    type(ccpp_constituent_prop_ptr_t), &
                     intent(in)     :: const_props(:)     ! ccpp constituent properties pointer
    real(kind_phys), intent(in)     :: ztodt              ! physics timestep [s]

    real(kind_phys), intent(in)     :: pmid(:,:)          ! midpoint pressures [Pa]
    real(kind_phys), intent(in)     :: pmiddry(:,:)       ! dry pressure at midpoints [Pa]
    real(kind_phys), intent(in)     :: pdel(:,:)          ! layer thickness (delta-p) [Pa]
    real(kind_phys), intent(in)     :: pdeldry(:,:)       ! dry layer thickness [Pa]
    real(kind_phys), intent(in)     :: rpdel(:,:)         ! 1.0 / pdel
    real(kind_phys), intent(in)     :: rpdeldry(:,:)      ! 1.0 / pdeldry

    real(kind_phys), intent(in)     :: zm(:,:)            ! geopotential height at midpoints [m]
    real(kind_phys), intent(in)     :: qpert_in(:)        ! PBL perturbation specific humidity (convective humidity excess) [kg kg-1]
    real(kind_phys), intent(in)     :: phis(:)            ! surface geopotential [m2 s-2]
    real(kind_phys), intent(in)     :: pblh(:)            ! PBL height [m]
    real(kind_phys), intent(in)     :: t(:,:)             ! temperature [K]
    real(kind_phys), intent(in)     :: q(:,:,:)           ! constituents [kg kg-1]

    ! Output arguments
    real(kind_phys), intent(out)    :: dq(:,:,:)          ! constituent tendencies [kg kg-1 s-1]
    real(kind_phys), intent(out)    :: qc_sh(:,:)         ! dq/dt due to export of cloud water / shallow reserved cloud condensate [kg kg-1 s-1]
    real(kind_phys), intent(out)    :: cmfdt(:,:)         ! heating rate (to ptend%s) [J kg-1 s-1]
    real(kind_phys), intent(out)    :: cmfmc_sh(:,:)      ! convective updraft mass flux, shallow [kg s-1 m-2]
    real(kind_phys), intent(out)    :: cmfdqr(:,:)        ! q tendency due to shallow convective rainout [kg kg-1 s-1]
    real(kind_phys), intent(out)    :: cmfsl(:,:)         ! moist shallow convection liquid water static energy flux [W m-2]
    real(kind_phys), intent(out)    :: cmflq(:,:)         ! moist shallow convection total water flux [W m-2]
    real(kind_phys), intent(out)    :: precc(:)           ! shallow convective precipitation rate [m s-1]
    real(kind_phys), intent(out)    :: cnt_sh(:)          ! top level of shallow convective activity [index]
    real(kind_phys), intent(out)    :: cnb_sh(:)          ! bottom level of shallow convective activity [index]
    real(kind_phys), intent(out)    :: icwmr(:,:)         ! shallow convection in-cloud water mixing ratio [kg kg-1]
    real(kind_phys), intent(out)    :: rliq_sh(:)         ! vertically-integrated shallow reserved cloud condensate [m s-1]

    character(len=64),  intent(out) :: scheme_name        ! scheme name
    real(kind_phys), intent(out)    :: flx_cnd(:)         ! net_liquid_and_lwe_ice_fluxes_through_top_and_bottom_of_atmosphere_column [m s-1] for check_energy_chng

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: tpert(ncol)        ! PBL perturbation temperature (convective temperature excess) [K]

    character(len=256) :: const_standard_name ! temp: constituent standard name
    logical            :: const_is_dry        ! temp: constituent is dry flag
    integer            :: const_wv_idx        ! temp: index of water vapor

    real(kind_phys) :: pm(ncol,pver)       ! pressure [Pa]
    real(kind_phys) :: pd(ncol,pver)       ! delta-p [Pa]
    real(kind_phys) :: rpd(ncol,pver)      ! 1./pdel [Pa-1]
    real(kind_phys) :: cmfdq(ncol,pver)    ! dq(wv)/dt due to moist convection (later copied to dq(:,:,const_wv_idx)) [kg kg-1 s-1]
    real(kind_phys) :: gam(ncol,pver)      ! 1/cp (d(qsat)/dT) change in saturation mixing ratio with temp
    real(kind_phys) :: sb(ncol,pver)       ! dry static energy (s bar) [J kg-1]
    real(kind_phys) :: hb(ncol,pver)       ! moist static energy (h bar) [J kg-1]
    real(kind_phys) :: shbs(ncol,pver)     ! sat. specific humidity (sh bar star)
    real(kind_phys) :: hbs(ncol,pver)      ! sat. moist static energy (h bar star)
    real(kind_phys) :: shbh(ncol,pver+1)   ! specific humidity on interfaces
    real(kind_phys) :: sbh(ncol,pver+1)    ! s bar on interfaces
    real(kind_phys) :: hbh(ncol,pver+1)    ! h bar on interfaces
    real(kind_phys) :: cmrh(ncol,pver+1)   ! interface constituent mixing ratio
    real(kind_phys) :: prec(ncol)          ! instantaneous total precipitation
    real(kind_phys) :: dzcld(ncol)         ! depth of convective layer (m)
    real(kind_phys) :: beta(ncol)          ! overshoot parameter (fraction)
    real(kind_phys) :: betamx(ncol)        ! local maximum on overshoot
    real(kind_phys) :: eta(ncol)           ! convective mass flux (kg/m^2 s)
    real(kind_phys) :: etagdt(ncol)        ! eta*grav*dt
    real(kind_phys) :: cldwtr(ncol)        ! cloud water (mass)
    real(kind_phys) :: rnwtr(ncol)         ! rain water  (mass)
    real(kind_phys) :: totcond(ncol)       ! total condensate; mix of precip and cloud water (mass)
    real(kind_phys) :: sc  (ncol)          ! dry static energy   ("in-cloud")
    real(kind_phys) :: shc (ncol)          ! specific humidity   ("in-cloud")
    real(kind_phys) :: hc  (ncol)          ! moist static energy ("in-cloud")
    real(kind_phys) :: cmrc(ncol)          ! constituent mix rat ("in-cloud")
    real(kind_phys) :: dq1(ncol)           ! shb  convective change (lower lvl)
    real(kind_phys) :: dq2(ncol)           ! shb  convective change (mid level)
    real(kind_phys) :: dq3(ncol)           ! shb  convective change (upper lvl)
    real(kind_phys) :: ds1(ncol)           ! sb   convective change (lower lvl)
    real(kind_phys) :: ds2(ncol)           ! sb   convective change (mid level)
    real(kind_phys) :: ds3(ncol)           ! sb   convective change (upper lvl)
    real(kind_phys) :: dcmr1(ncol)         ! q convective change (lower lvl)
    real(kind_phys) :: dcmr2(ncol)         ! q convective change (mid level)
    real(kind_phys) :: dcmr3(ncol)         ! q convective change (upper lvl)
    real(kind_phys) :: estemp(ncol,pver)   ! saturation vapor pressure (scratch)
    real(kind_phys) :: vtemp1(2*ncol)      ! intermediate scratch vector
    real(kind_phys) :: vtemp2(2*ncol)      ! intermediate scratch vector
    real(kind_phys) :: vtemp3(2*ncol)      ! intermediate scratch vector
    real(kind_phys) :: vtemp4(2*ncol)      ! intermediate scratch vector
    real(kind_phys) :: vtemp5(2*ncol)      ! intermediate scratch vector
    integer         :: indx1(ncol)         ! longitude indices for condition true
    logical         :: etagt0              ! true if eta > 0.0
    real(kind_phys) :: cats                ! modified characteristic adj. time
    real(kind_phys) :: rtdt                ! 1./ztodt
    real(kind_phys) :: qprime              ! modified specific humidity pert.
    real(kind_phys) :: tprime              ! modified thermal perturbation
    real(kind_phys) :: pblhgt              ! bounded pbl height (max[pblh,1m])
    real(kind_phys) :: fac1                ! intermediate scratch variable
    real(kind_phys) :: shprme              ! intermediate specific humidity pert.
    real(kind_phys) :: qsattp              ! sat mix rat for thermally pert PBL parcels
    real(kind_phys) :: dz                  ! local layer depth
    real(kind_phys) :: temp1               ! intermediate scratch variable
    real(kind_phys) :: b1                  ! bouyancy measure in detrainment lvl
    real(kind_phys) :: b2                  ! bouyancy measure in condensation lvl
    real(kind_phys) :: temp2               ! intermediate scratch variable
    real(kind_phys) :: temp3               ! intermediate scratch variable
    real(kind_phys) :: g                   ! bounded vertical gradient of hb
    real(kind_phys) :: tmass               ! total mass available for convective exch
    real(kind_phys) :: denom               ! intermediate scratch variable
    real(kind_phys) :: qtest1              ! used in negative q test (middle lvl)
    real(kind_phys) :: qtest2              ! used in negative q test (lower lvl)
    real(kind_phys) :: fslkp               ! flux lw static energy (bot interface)
    real(kind_phys) :: fslkm               ! flux lw static energy (top interface)
    real(kind_phys) :: fqlkp               ! flux total water (bottom interface)
    real(kind_phys) :: fqlkm               ! flux total water (top interface)
    real(kind_phys) :: botflx              ! bottom constituent mixing ratio flux
    real(kind_phys) :: topflx              ! top constituent mixing ratio flux
    real(kind_phys) :: efac1               ! ratio q to convectively induced chg (btm lvl)
    real(kind_phys) :: efac2               ! ratio q to convectively induced chg (mid lvl)
    real(kind_phys) :: efac3               ! ratio q to convectively induced chg (top lvl)
    real(kind_phys) :: tb(ncol,pver)       ! working storage for temp (t bar)
    real(kind_phys) :: shb(ncol,pver)      ! working storage for spec hum (sh bar)
    real(kind_phys) :: adjfac              ! adjustment factor (relaxation related)
    real(kind_phys) :: rktp
    real(kind_phys) :: rk
    integer         :: i,k                 ! longitude, level indices
    integer         :: ii                  ! index on "gathered" vectors
    integer         :: len1                ! vector length of "gathered" vectors
    integer         :: m                   ! constituent index
    integer         :: ktp                 ! tmp indx used to track top of convective layer

    ! debug use quantities
    real(kind_phys) :: rh                  ! relative humidity
    real(kind_phys) :: es                  ! sat vapor pressure
    real(kind_phys) :: hsum1               ! moist static energy integral
    real(kind_phys) :: qsum1               ! total water integral
    real(kind_phys) :: hsum2               ! final moist static energy integral
    real(kind_phys) :: qsum2               ! final total water integral
    real(kind_phys) :: fac                 ! intermediate scratch variable
    integer         :: n                   ! vertical index     (diagnostics)
    integer         :: kp                  ! vertical index     (diagnostics)
    integer         :: kpp                 ! index offset, kp+1 (diagnostics)
    integer         :: kpm1                ! index offset, kp-1 (diagnostics)

    errmsg = ''
    errflg = 0

    scheme_name = 'hack_convect_shallow'

    !---------------------------------------------------
    ! Initialize output tendencies
    !---------------------------------------------------
    cmfdt   (:ncol,:)     = 0._kind_phys
    cmfdq   (:ncol,:)     = 0._kind_phys
    cmfmc_sh(:ncol,:)     = 0._kind_phys
    cmfdqr  (:ncol,:)     = 0._kind_phys
    cmfsl   (:ncol,:)     = 0._kind_phys
    cmflq   (:ncol,:)     = 0._kind_phys
    qc_sh   (:ncol,:)     = 0._kind_phys
    rliq_sh (:ncol)       = 0._kind_phys

    ! Check constituents list and locate water vapor index
    ! (not assumed to be 1)
    call ccpp_const_get_idx(const_props, &
         'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', &
         const_wv_idx, errmsg, errflg)

    !---------------------------------------------------
    ! copy q to dq for passive tracer transport.
    ! this is NOT an initialization. the dq at this point
    ! is not physical (used as temporary here) only at the end
    ! dq is updated to be an actual tendency.
    !---------------------------------------------------
    if(pcnst > 1) then
      ! set dq for passive tracer transport from q as temporary...
      dq(:ncol,:,:) = q(:ncol,:,:)

      ! except for water vapor
      dq(:ncol,:,const_wv_idx) = 0._kind_phys
    endif

    !---------------------------------------------------
    ! Quantity preparations from convect_shallow.F90.
    !---------------------------------------------------

    ! convect_shallow.F90 is not linked to pbuf tpert and always sets to zero.
    ! "This field probably should reference the pbuf tpert field but it doesnt"
    tpert(:ncol) = 0.0_kind_phys

    !---------------------------------------------------
    ! Preparation of working arrays
    !---------------------------------------------------
    ! Ensure that characteristic adjustment time scale (cmftau) assumed
    ! in estimate of eta isn't smaller than model time scale (ztodt)
    ! The time over which the convection is assumed to act (the adjustment
    ! time scale) can be applied with each application of the three-level
    ! cloud model, or applied to the column tendencies after a "hard"
    ! adjustment (i.e., on a 2-delta t time scale) is evaluated
    if (rlxclm) then
       cats   = ztodt             ! relaxation applied to column
       adjfac = ztodt/(max(ztodt,cmftau))
    else
       cats   = max(ztodt,cmftau) ! relaxation applied to triplet
       adjfac = 1.0_kind_phys
    endif
    rtdt = 1.0_kind_phys/ztodt

    ! Move temperature and moisture into working storage
    do k=limcnv,pver
       do i=1,ncol
          tb (i,k) = t(i,k)
          shb(i,k) = q(i,k,const_wv_idx)
       end do
    end do
    do k=1,pver
       do i=1,ncol
          icwmr(i,k) = 0._kind_phys
       end do
    end do

    ! Compute sb,hb,shbs,hbs
    do k = limcnv,pver
       call qsat(tb(1:ncol,k), pmid(1:ncol,k), &
            estemp(1:ncol,k), shbs(1:ncol,k), ncol, &
            gam=gam(1:ncol,k))
    end do

    do k=limcnv,pver
       do i=1,ncol
          sb (i,k) = cp*tb(i,k) + zm(i,k)*grav + phis(i)
          hb (i,k) = sb(i,k) + hlat*shb(i,k)
          hbs(i,k) = sb(i,k) + hlat*shbs(i,k)
       end do
    end do

    ! Compute sbh, shbh
    do k=limcnv+1,pver
       do i=1,ncol
          sbh (i,k) = 0.5_kind_phys*(sb(i,k-1) + sb(i,k))
          shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
          hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
       end do
    end do

    ! Specify properties at top of model (not used, but filling anyway)
    do i=1,ncol
       sbh (i,limcnv) = sb(i,limcnv)
       shbh(i,limcnv) = shb(i,limcnv)
       hbh (i,limcnv) = hb(i,limcnv)
    end do

    ! Zero vertically independent control, tendency & diagnostic arrays
    do i=1,ncol
       prec(i)  = 0.0_kind_phys
       dzcld(i) = 0.0_kind_phys
       cnb_sh(i)= 0.0_kind_phys
       cnt_sh(i)= real(pver+1,kind_phys)
    end do

    if(debug_verbose) then
      ! DEBUG DIAGNOSTICS - Output initial thermodynamic profile
      do i=1,ncol
        if(i == 1) then
          ! Approximate vertical integral of moist static energy
          ! and total precipitable water
          hsum1 = 0.0_kind_phys
          qsum1 = 0.0_kind_phys
          do k=limcnv,pver
            hsum1 = hsum1 + pdel(i,k)*rgrav*hb(i,k)
            qsum1 = qsum1 + pdel(i,k)*rgrav*shb(i,k)
          end do

          write(iulog,8010)
          fac = grav*864._kind_phys
          do k=limcnv,pver
            rh = shb(i,k)/shbs(i,k)
            write(iulog,8020) shbh(i,k),sbh(i,k),hbh(i,k),fac*cmfmc_sh(i,k),cmfsl(i,k), cmflq(i,k)
            write(iulog,8040) tb(i,k),shb(i,k),rh,sb(i,k),hb(i,k),hbs(i,k),ztodt*cmfdt(i,k), &
                          ztodt*cmfdq(i,k),ztodt*cmfdqr(i,k)
          end do
          write(iulog, 8000) prec(i)
        end if
      end do
    endif

    !---------------------------------------------------
    ! Begin moist convective mass flux adjustment procedure.
    ! Formalism ensures that negative cloud liquid water can never occur
    !---------------------------------------------------
    kloop: do k = pver-1,limcnv+1,-1
      do i = 1, ncol
        etagdt(i) = 0.0_kind_phys
        eta   (i) = 0.0_kind_phys
        beta  (i) = 0.0_kind_phys
        ds1   (i) = 0.0_kind_phys
        ds2   (i) = 0.0_kind_phys
        ds3   (i) = 0.0_kind_phys
        dq1   (i) = 0.0_kind_phys
        dq2   (i) = 0.0_kind_phys
        dq3   (i) = 0.0_kind_phys
        ! Specification of "cloud base" conditions
        qprime    = 0.0_kind_phys
        tprime    = 0.0_kind_phys

        ! Assign tprime within the PBL to be proportional to the quantity
        ! tpert (which will be bounded by tpmax), passed to this routine by
        ! the PBL routine.  Don't allow perturbation to produce a dry
        ! adiabatically unstable parcel.  Assign qprime within the PBL to be
        ! an appropriately modified value of the quantity qpert (which will be
        ! bounded by shpmax) passed to this routine by the PBL routine.  The
        ! quantity qprime should be less than the local saturation value
        ! (qsattp=qsat[t+tprime,p]).  In both cases, tpert and qpert are
        ! linearly reduced toward zero as the PBL top is approached.
        pblhgt = max(pblh(i),1.0_kind_phys)
        if ( (zm(i,k+1) <= pblhgt) .and. dzcld(i) == 0.0_kind_phys ) then
           fac1   = max(0.0_kind_phys,1.0_kind_phys-zm(i,k+1)/pblhgt)
           tprime = min(tpert(i),tpmax)*fac1
           qsattp = shbs(i,k+1) + cp*rhlat*gam(i,k+1)*tprime
           shprme = min(min(qpert_in(i),shpmax)*fac1,max(qsattp-shb(i,k+1),0.0_kind_phys))
           qprime = max(qprime,shprme)
        else
           tprime = 0.0_kind_phys
           qprime = 0.0_kind_phys
        end if

        ! Specify "updraft" (in-cloud) thermodynamic properties
        sc (i)    = sb (i,k+1) + cp*tprime
        shc(i)    = shb(i,k+1) + qprime
        hc (i)    = sc (i    ) + hlat*shc(i)
        vtemp4(i) = hc(i) - hbs(i,k)
        dz        = pdel(i,k)*rgas*tb(i,k)*rgrav/pmid(i,k)
        if (vtemp4(i) > 0.0_kind_phys) then
           dzcld(i) = dzcld(i) + dz
        else
           dzcld(i) = 0.0_kind_phys
        end if
      enddo

      if(debug_verbose) then
        ! DEBUG DIAGNOSTICS - output thermodynamic perturbation information
        do i=1,ncol
          if(i == 1) then
            write(iulog,8090) k+1,sc(i),shc(i),hc(i)
          end if
        enddo
      endif


      ! Check on moist convective instability
      ! Build index vector of points where instability exists
      len1 = 0
      do i=1,ncol
         if (vtemp4(i) > 0.0_kind_phys) then
            len1 = len1 + 1
            indx1(len1) = i
         end if
      end do

      if (len1 <= 0) cycle kloop

      ! Current level just below top level => no overshoot
      if (k <= limcnv+1) then
         do ii=1,len1
            i = indx1(ii)
            temp1     = vtemp4(i)/(1.0_kind_phys + gam(i,k))
            cldwtr(i) = max(0.0_kind_phys,(sb(i,k) - sc(i) + temp1))
            beta(i)   = 0.0_kind_phys
            vtemp3(i) = (1.0_kind_phys + gam(i,k))*(sc(i) - sbh(i,k))
         end do
      else
        ! First guess at overshoot parameter using crude buoyancy closure
        ! 10% overshoot assumed as a minimum and 1-c0*dz maximum to start
        ! If pre-existing supersaturation in detrainment layer, beta=0
        ! cldwtr is temporarily equal to hlat*l (l=> liquid water)
        do ii=1,len1
          i = indx1(ii)
          temp1     = vtemp4(i)/(1.0_kind_phys + gam(i,k))
          cldwtr(i) = max(0.0_kind_phys,(sb(i,k)-sc(i)+temp1))
          betamx(i) = 1.0_kind_phys - c0*max(0.0_kind_phys,(dzcld(i)-dzmin))
          b1        = (hc(i) - hbs(i,k-1))*pdel(i,k-1)
          b2        = (hc(i) - hbs(i,k  ))*pdel(i,k  )
          beta(i)   = max(betamn,min(betamx(i), 1.0_kind_phys + b1/b2))
          if (hbs(i,k-1) <= hb(i,k-1)) beta(i) = 0.0_kind_phys

          ! Bound maximum beta to ensure physically realistic solutions
          !
          ! First check constrains beta so that eta remains positive
          ! (assuming that eta is already positive for beta equal zero)
          vtemp1(i) = -(hbh(i,k+1) - hc(i))*pdel(i,k)*rpdel(i,k+1)+ &
                      (1.0_kind_phys + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i))
          vtemp2(i) = (1.0_kind_phys + gam(i,k))*(sc(i) - sbh(i,k))
          vtemp3(i) = vtemp2(i)
          if ((beta(i)*vtemp2(i) - vtemp1(i)) > 0._kind_phys) then
            betamx(i) = 0.99_kind_phys*(vtemp1(i)/vtemp2(i))
            beta(i)   = max(0.0_kind_phys,min(betamx(i),beta(i)))
          end if
        end do

        ! Second check involves supersaturation of "detrainment layer"
        ! small amount of supersaturation acceptable (by ssfac factor)
        do ii=1,len1
          i = indx1(ii)
          if (hb(i,k-1) < hbs(i,k-1)) then
            vtemp1(i) = vtemp1(i)*rpdel(i,k)
            temp2 = gam(i,k-1)*(sbh(i,k) - sc(i) + cldwtr(i)) -  &
                    hbh(i,k) + hc(i) - sc(i) + sbh(i,k)
            temp3 = vtemp3(i)*rpdel(i,k)
            vtemp2(i) = (ztodt/cats)*(hc(i) - hbs(i,k))*temp2/ &
                        (pdel(i,k-1)*(hbs(i,k-1) - hb(i,k-1))) + temp3
            if ((beta(i)*vtemp2(i) - vtemp1(i)) > 0._kind_phys) then
              betamx(i) = ssfac*(vtemp1(i)/vtemp2(i))
              beta(i)   = max(0.0_kind_phys,min(betamx(i),beta(i)))
            end if
          else
             beta(i) = 0.0_kind_phys
          end if
        end do

        ! Third check to avoid introducing 2 delta x thermodynamic
        ! noise in the vertical ... constrain adjusted h (or theta e)
        ! so that the adjustment doesn't contribute to "kinks" in h
        do ii=1,len1
           i = indx1(ii)
           g = min(0.0_kind_phys,hb(i,k) - hb(i,k-1))
           temp1 = (hb(i,k) - hb(i,k-1) - g)*(cats/ztodt)/(hc(i) - hbs(i,k))
           vtemp1(i) = temp1*vtemp1(i) + (hc(i) - hbh(i,k+1))*rpdel(i,k)
           vtemp2(i) = temp1*vtemp3(i)*rpdel(i,k) + (hc(i) - hbh(i,k) - cldwtr(i))* &
                       (rpdel(i,k) + rpdel(i,k+1))
           if ((beta(i)*vtemp2(i) - vtemp1(i)) > 0._kind_phys) then
              if (vtemp2(i) /= 0.0_kind_phys) then
                betamx(i) = vtemp1(i)/vtemp2(i)
              else
                betamx(i) = 0.0_kind_phys
              end if
              beta(i) = max(0.0_kind_phys,min(betamx(i),beta(i)))
           end if
        end do
      end if ! (k <= limcnv+1) Current level just below top level => no overshoot


      ! Calculate mass flux required for stabilization.
      !
      ! Ensure that the convective mass flux, eta, is positive by
      ! setting negative values of eta to zero..
      ! Ensure that estimated mass flux cannot move more than the
      ! minimum of total mass contained in either layer k or layer k+1.
      ! Also test for other pathological cases that result in non-
      ! physical states and adjust eta accordingly.
      do ii=1,len1
        i = indx1(ii)
        beta(i) = max(0.0_kind_phys,beta(i))
        temp1 = hc(i) - hbs(i,k)
        temp2 = ((1.0_kind_phys + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i)) - &
                  beta(i)*vtemp3(i))*rpdel(i,k) - (hbh(i,k+1) - hc(i))*rpdel(i,k+1)
        eta(i) = temp1/(temp2*grav*cats)
        tmass = min(pdel(i,k),pdel(i,k+1))*rgrav
        if (eta(i) > tmass*rtdt .or. eta(i) <= 0.0_kind_phys) eta(i) = 0.0_kind_phys

        ! Check on negative q in top layer (bound beta)
        if (shc(i)-shbh(i,k) < 0.0_kind_phys .and. beta(i)*eta(i) /= 0.0_kind_phys) then
           denom = eta(i)*grav*ztodt*(shc(i) - shbh(i,k))*rpdel(i,k-1)
           beta(i) = max(0.0_kind_phys,min(-0.999_kind_phys*shb(i,k-1)/denom,beta(i)))
        end if

        ! Check on negative q in middle layer (zero eta)
        qtest1 = shb(i,k) + eta(i)*grav*ztodt*((shc(i) - shbh(i,k+1)) - &
                 (1.0_kind_phys - beta(i))*cldwtr(i)*rhlat - beta(i)*(shc(i) - shbh(i,k)))* &
           rpdel(i,k)
        if (qtest1 <= 0.0_kind_phys) eta(i) = 0.0_kind_phys

        ! Check on negative q in lower layer (bound eta)
        fac1 = -(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
        qtest2 = shb(i,k+1) - eta(i)*grav*ztodt*fac1
        if (qtest2 < 0.0_kind_phys) then
           eta(i) = 0.99_kind_phys*shb(i,k+1)/(grav*ztodt*fac1)
        end if
        etagdt(i) = eta(i)*grav*ztodt
      end do

      if(debug_verbose) then
        do i=1,ncol
          if (i == 1) then
            write(iulog,8080) beta(i), eta(i)
          endif
        enddo
      endif

      ! Calculate cloud water, rain water, and thermodynamic changes
      do ii=1,len1
        i = indx1(ii)
        icwmr(i,k) = cldwtr(i)*rhlat
        cldwtr(i) = etagdt(i)*cldwtr(i)*rhlat*rgrav

        ! JJH changes to facilitate export of cloud liquid water --------------------------------
        totcond(i) = (1.0_kind_phys - beta(i))*cldwtr(i)
        rnwtr(i) = min(totcond(i),c0*(dzcld(i)-dzmin)*cldwtr(i))
        ds1(i) = etagdt(i)*(sbh(i,k+1) - sc(i))*rpdel(i,k+1)
        dq1(i) = etagdt(i)*(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
        ds2(i) = (etagdt(i)*(sc(i) - sbh(i,k+1)) +  &
                 hlat*grav*cldwtr(i) - beta(i)*etagdt(i)*(sc(i) - sbh(i,k)))*rpdel(i,k)

        ! JJH change for export of cloud liquid water; must use total condensate
        ! since rainwater no longer represents total condensate
        dq2(i) = (etagdt(i)*(shc(i) - shbh(i,k+1)) - grav*totcond(i) - beta(i)* &
                 etagdt(i)*(shc(i) - shbh(i,k)))*rpdel(i,k)
        ds3(i) = beta(i)*(etagdt(i)*(sc(i) - sbh(i,k)) - hlat*grav*cldwtr(i))* &
                 rpdel(i,k-1)
        dq3(i) = beta(i)*etagdt(i)*(shc(i) - shbh(i,k))*rpdel(i,k-1)

        ! Isolate convective fluxes for later diagnostics
        fslkp = eta(i)*(sc(i) - sbh(i,k+1))
        fslkm = beta(i)*(eta(i)*(sc(i) - sbh(i,k)) - hlat*cldwtr(i)*rtdt)
        fqlkp = eta(i)*(shc(i) - shbh(i,k+1))
        fqlkm = beta(i)*eta(i)*(shc(i) - shbh(i,k))

        ! Update thermodynamic profile (update sb, hb, & hbs later)
        tb (i,k+1) = tb(i,k+1)  + ds1(i)*rcp
        tb (i,k  ) = tb(i,k  )  + ds2(i)*rcp
        tb (i,k-1) = tb(i,k-1)  + ds3(i)*rcp
        shb(i,k+1) = shb(i,k+1) + dq1(i)
        shb(i,k  ) = shb(i,k  ) + dq2(i)
        shb(i,k-1) = shb(i,k-1) + dq3(i)

        ! ** Update diagnostic information for final budget **
        ! Tracking precipitation, temperature & specific humidity tendencies,
        ! rainout term, convective mass flux, convective liquid
        ! water static energy flux, and convective total water flux
        ! The variable afac makes the necessary adjustment to the
        ! diagnostic fluxes to account for adjustment time scale based on
        ! how relaxation time scale is to be applied (column vs. triplet)
        prec(i)    = prec(i) + (rnwtr(i)/rhoh2o)*adjfac

        ! The following variables have units of "units"/second
        cmfdt (i,k+1) = cmfdt (i,k+1) + ds1(i)*rtdt*adjfac
        cmfdt (i,k  ) = cmfdt (i,k  ) + ds2(i)*rtdt*adjfac
        cmfdt (i,k-1) = cmfdt (i,k-1) + ds3(i)*rtdt*adjfac
        cmfdq (i,k+1) = cmfdq (i,k+1) + dq1(i)*rtdt*adjfac
        cmfdq (i,k  ) = cmfdq (i,k  ) + dq2(i)*rtdt*adjfac
        cmfdq (i,k-1) = cmfdq (i,k-1) + dq3(i)*rtdt*adjfac

        ! JJH changes to export cloud liquid water --------------------------------
        qc_sh   (i,k  ) = (grav*(totcond(i)-rnwtr(i))*rpdel(i,k))*rtdt*adjfac
        cmfdqr  (i,k  ) = cmfdqr(i,k  ) + (grav*rnwtr(i)*rpdel(i,k))*rtdt*adjfac
        cmfmc_sh(i,k+1) = cmfmc_sh(i,k+1) + eta(i)*adjfac
        cmfmc_sh(i,k  ) = cmfmc_sh(i,k  ) + beta(i)*eta(i)*adjfac

        ! The following variables have units of w/m**2
        cmfsl (i,k+1) = cmfsl (i,k+1) + fslkp*adjfac
        cmfsl (i,k  ) = cmfsl (i,k  ) + fslkm*adjfac
        cmflq (i,k+1) = cmflq (i,k+1) + hlat*fqlkp*adjfac
        cmflq (i,k  ) = cmflq (i,k  ) + hlat*fqlkm*adjfac
      enddo

      ! Next, convectively modify passive constituents
      ! For now, when applying relaxation time scale to thermal fields after
      ! entire column has undergone convective overturning, constituents will
      ! be mixed using a "relaxed" value of the mass flux determined above
      ! Although this will be inconsistant with the treatment of the thermal
      ! fields, it's computationally much cheaper, no more-or-less justifiable,
      ! and consistent with how the history tape mass fluxes would be used in
      ! an off-line mode (i.e., using an off-line transport model)
      const_modify_loop: do m = 1, pcnst
        ! Water vapor needs to be skipped in the loop.
        if (m == const_wv_idx) then
          cycle const_modify_loop
        endif

        ! assign pd, rpd, pm temporary properties based on constituent dry/moist mixing ratio
        call const_props(m)%is_dry(const_is_dry, errflg, errmsg)
        if(const_is_dry) then
          pd (:ncol,:) = pdeldry (:ncol,:)
          rpd(:ncol,:) = rpdeldry(:ncol,:)
          pm (:ncol,:) = pmiddry (:ncol,:)
        else
          pd (:ncol,:) = pdel    (:ncol,:)
          rpd(:ncol,:) = rpdel   (:ncol,:)
          pm (:ncol,:) = pmid    (:ncol,:)
        endif

        pcl1loop: do ii=1,len1
          i = indx1(ii)

          ! If any of the reported values of the constituent is negative in
          ! the three adjacent levels, nothing will be done to the profile
          if ((dq(i,k+1,m) < 0.0_kind_phys) .or. (dq(i,k,m) < 0.0_kind_phys) .or. (dq(i,k-1,m) < 0.0_kind_phys)) cycle pcl1loop

          ! Specify constituent interface values (linear interpolation)
          cmrh(i,k  ) = 0.5_kind_phys*(dq(i,k-1,m) + dq(i,k  ,m))
          cmrh(i,k+1) = 0.5_kind_phys*(dq(i,k  ,m) + dq(i,k+1,m))

          ! Specify perturbation properties of constituents in PBL
          pblhgt = max(pblh(i),1.0_kind_phys)
          if ( (zm(i,k+1) <= pblhgt) .and. dzcld(i) == 0.0_kind_phys ) then
              fac1 = max(0.0_kind_phys,1.0_kind_phys-zm(i,k+1)/pblhgt)
              ! cmrc(i) = dq(i,k+1,m) + qpert(i,m)*fac1
              ! hplin - qpert for m>1 is always zero
              cmrc(i) = dq(i,k+1,m)
          else
             cmrc(i) = dq(i,k+1,m)
          end if

          ! Determine fluxes, flux divergence => changes due to convection
          ! Logic must be included to avoid producing negative values. A bit
          ! messy since there are no a priori assumptions about profiles.
          ! Tendency is modified (reduced) when pending disaster detected.

          botflx   = etagdt(i)*(cmrc(i) - cmrh(i,k+1))*adjfac
          topflx   = beta(i)*etagdt(i)*(cmrc(i)-cmrh(i,k))*adjfac
          dcmr1(i) = -botflx*rpd(i,k+1)
          efac1    = 1.0_kind_phys
          efac2    = 1.0_kind_phys
          efac3    = 1.0_kind_phys

          if (dq(i,k+1,m)+dcmr1(i) < 0.0_kind_phys) then
             if ( abs(dcmr1(i)) > 1.e-300_kind_phys ) then
                efac1 = max(tiny,abs(dq(i,k+1,m)/dcmr1(i)) - eps)
             else
                efac1 = tiny
             endif
          end if

          if (efac1 == tiny .or. efac1 > 1.0_kind_phys) efac1 = 0.0_kind_phys
          dcmr1(i) = -efac1*botflx*rpd(i,k+1)
          dcmr2(i) = (efac1*botflx - topflx)*rpd(i,k)

          if (dq(i,k,m)+dcmr2(i) < 0.0_kind_phys) then
             if ( abs(dcmr2(i)) > 1.e-300_kind_phys ) then
                efac2 = max(tiny,abs(dq(i,k  ,m)/dcmr2(i)) - eps)
             else
                efac2 = tiny
             endif
          end if

          if (efac2 == tiny .or. efac2 > 1.0_kind_phys) efac2 = 0.0_kind_phys
          dcmr2(i) = (efac1*botflx - efac2*topflx)*rpd(i,k)
          dcmr3(i) = efac2*topflx*rpd(i,k-1)

          if ( (dq(i,k-1,m)+dcmr3(i) < 0.0_kind_phys ) ) then
             if  ( abs(dcmr3(i)) > 1.e-300_kind_phys ) then
                efac3 = max(tiny,abs(dq(i,k-1,m)/dcmr3(i)) - eps)
             else
                efac3 = tiny
             endif
          end if

          if (efac3 == tiny .or. efac3 > 1.0_kind_phys) efac3 = 0.0_kind_phys
          efac3    = min(efac2,efac3)
          dcmr2(i) = (efac1*botflx - efac3*topflx)*rpd(i,k)
          dcmr3(i) = efac3*topflx*rpd(i,k-1)

          dq(i,k+1,m) = dq(i,k+1,m) + dcmr1(i)
          dq(i,k  ,m) = dq(i,k  ,m) + dcmr2(i)
          dq(i,k-1,m) = dq(i,k-1,m) + dcmr3(i)
        end do pcl1loop
      end do const_modify_loop
      ! Constituent modifications complete

      ! This if restructured from a goto
      if (k /= limcnv+1) then
        ! Complete update of thermodynamic structure at integer levels
        ! gather/scatter points that need new values of shbs and gamma
        do ii=1,len1
           i = indx1(ii)
           vtemp1(ii     ) = tb(i,k)
           vtemp1(ii+len1) = tb(i,k-1)
           vtemp2(ii     ) = pmid(i,k)
           vtemp2(ii+len1) = pmid(i,k-1)
        end do
        call qsat(vtemp1(1:2*len1), vtemp2(1:2*len1), &
                  vtemp5(1:2*len1), vtemp3(1:2*len1), 2*len1, gam=vtemp4(1:2*len1))
        do ii=1,len1
           i = indx1(ii)
           shbs(i,k  ) = vtemp3(ii     )
           shbs(i,k-1) = vtemp3(ii+len1)
           gam(i,k  ) = vtemp4(ii     )
           gam(i,k-1) = vtemp4(ii+len1)
           sb (i,k  ) = sb(i,k  ) + ds2(i)
           sb (i,k-1) = sb(i,k-1) + ds3(i)
           hb (i,k  ) = sb(i,k  ) + hlat*shb(i,k  )
           hb (i,k-1) = sb(i,k-1) + hlat*shb(i,k-1)
           hbs(i,k  ) = sb(i,k  ) + hlat*shbs(i,k  )
           hbs(i,k-1) = sb(i,k-1) + hlat*shbs(i,k-1)
        end do

        ! Update thermodynamic information at half (i.e., interface) levels
        do ii=1,len1
           i = indx1(ii)
           sbh (i,k) = 0.5_kind_phys*(sb(i,k) + sb(i,k-1))
           shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
           hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
           sbh (i,k-1) = 0.5_kind_phys*(sb(i,k-1) + sb(i,k-2))
           shbh(i,k-1) = qhalf(shb(i,k-2),shb(i,k-1),shbs(i,k-2),shbs(i,k-1))
           hbh (i,k-1) = sbh(i,k-1) + hlat*shbh(i,k-1)
        end do
      end if ! k /= limcnv+1

      ! Ensure that dzcld is reset if convective mass flux zero
      ! specify the current vertical extent of the convective activity
      ! top of convective layer determined by size of overshoot param.
      do i=1,ncol
        etagt0 = eta(i).gt.0.0_kind_phys
        if ( .not. etagt0) dzcld(i) = 0.0_kind_phys
        if (etagt0 .and. beta(i) > betamn) then
           ktp = k-1
        else
           ktp = k
        end if
        if (etagt0) then
           rk = k
           rktp = ktp
           cnt_sh(i) = min(cnt_sh(i),rktp)
           cnb_sh(i) = max(cnb_sh(i),rk)
        end if
      end do
    end do kloop

    !---------------------------------------------------
    ! apply final thermodynamic tendencies
    !---------------------------------------------------
    ! Set output q tendencies...
    ! ...for water vapor
    dq(:ncol,:,const_wv_idx) = cmfdq(:ncol,:)

    ! ...for other tracers from passive tracer transport
    do m = 1, pcnst
      if (m .ne. const_wv_idx) then
        dq(:ncol,:,m) = (dq(:ncol,:,m) - q(:ncol,:,m))/ztodt
      endif
    enddo

    ! Kludge to prevent cnb_sh-cnt_sh from being zero (in the event
    ! someone decides that they want to divide by this quantity)
    do i=1,ncol
       if (cnb_sh(i) /= 0.0_kind_phys .and. cnb_sh(i) == cnt_sh(i)) then
          cnt_sh(i) = cnt_sh(i) - 1.0_kind_phys
       end if
    end do

    do i=1,ncol
       precc(i) = prec(i)*rtdt
    end do

    ! Compute reserved liquid (not yet in cldliq) for energy integrals.
    ! Treat rliq_sh as flux out bottom, to be added back later.
    do k = 1, pver
       do i = 1, ncol
          rliq_sh(i) = rliq_sh(i) + qc_sh(i,k)*pdel(i,k)/grav
       end do
    end do

    ! rliq_sh is converted to precipitation units [m s-1]
    rliq_sh(:ncol) = rliq_sh(:ncol) / 1000._kind_phys

    ! Prepare boundary fluxes for check_energy [m s-1]
    flx_cnd(:ncol) = precc(:ncol) + rliq_sh(:ncol)

    if(debug_verbose) then
      ! DEBUG DIAGNOSTICS - show final result
      do i=1,ncol
        if (i == 1) then
          fac = grav*864._kind_phys
          write(iulog, 8010)
          do k=limcnv,pver
            rh = shb(i,k)/shbs(i,k)
            write(iulog, 8020) shbh(i,k),sbh(i,k),hbh(i,k),fac*cmfmc_sh(i,k), &
                               cmfsl(i,k), cmflq(i,k)
            write(iulog, 8040) tb(i,k),shb(i,k),rh   ,sb(i,k),hb(i,k), &
                               hbs(i,k), ztodt*cmfdt(i,k),ztodt*cmfdq(i,k), &
                               ztodt*cmfdqr(i,k)
          end do
          write(iulog, 8000) prec(i)

          ! approximate vertical integral of moist static energy and
          ! total preciptable water after adjustment and output changes
          hsum2 = 0.0_kind_phys
          qsum2 = 0.0_kind_phys
          do k=limcnv,pver
            hsum2 = hsum2 + pdel(i,k)*rgrav*hb(i,k)
            qsum2 = qsum2 + pdel(i,k)*rgrav*shb(i,k)
          end do
          write(iulog,8070) hsum1, hsum2, abs(hsum2-hsum1)/hsum2, &
                            qsum1, qsum2, abs(qsum2-qsum1)/qsum2
        end if
      enddo
    endif

    ! Diagnostic use format strings
8000              format(///,10x,'PREC = ',3pf12.6,/)
8010              format('1**        TB      SHB      RH       SB', &
                        '       HB      HBS      CAH      CAM       PRECC ', &
                        '     ETA      FSL       FLQ     **', /)
8020              format(' ----- ',     9x,3p,f7.3,2x,2p,     9x,-3p,f7.3,2x, &
                        f7.3, 37x, 0p,2x,f8.2,  0p,2x,f8.2,2x,f8.2, ' ----- ')
8030              format(' ----- ',  0p,82x,f8.2,  0p,2x,f8.2,2x,f8.2, &
                         ' ----- ')
8040              format(' - - - ',f7.3,2x,3p,f7.3,2x,2p,f7.3,2x,-3p,f7.3,2x, &
                        f7.3, 2x,f8.3,2x,0p,f7.3,3p,2x,f7.3,2x,f7.3,30x, &
                         ' - - - ')
8050              format(' ----- ',110x,' ----- ')
8060              format('1 K =>',  i4,/, &
                           '           TB      SHB      RH       SB', &
                           '       HB      HBS      CAH      CAM       PREC ', &
                           '     ETA      FSL       FLQ', /)
8070              format(' VERTICALLY INTEGRATED MOIST STATIC ENERGY BEFORE, AFTER', &
                        ' AND PERCENTAGE DIFFERENCE => ',1p,2e15.7,2x,2p,f7.3,/, &
                        ' VERTICALLY INTEGRATED MOISTURE            BEFORE, AFTER', &
                        ' AND PERCENTAGE DIFFERENCE => ',1p,2e15.7,2x,2p,f7.3,/)
8080              format(' BETA, ETA => ', 1p,2e12.3)
8090              format (' k+1, sc, shc, hc => ', 1x, i2, 1p, 3e12.4)
  end subroutine hack_convect_shallow_run

  ! qhalf computes the specific humidity at interface levels between two model layers (interpolate moisture)
  pure function qhalf(sh1,sh2,shbs1,shbs2) result(qh)
    real(kind_phys), intent(in) :: sh1    ! humidity of layer 1 [kg kg-1]
    real(kind_phys), intent(in) :: sh2    ! humidity of layer 2 [kg kg-1]
    real(kind_phys), intent(in) :: shbs1  ! saturation specific humidity of layer 1 [kg kg-1]
    real(kind_phys), intent(in) :: shbs2  ! saturation specific humidity of layer 2 [kg kg-1]
    real(kind_phys) :: qh
    qh = min(max(sh1,sh2),(shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))
  end function qhalf
end module hack_convect_shallow
