! Hack shallow convective scheme.
!
! Original Author: J. Hack
! CCPPized: Haipeng Lin, October 2024
module cmfmca
  use ccpp_kinds,           only : kind_phys
  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: cmfmca_init
  public :: cmfmca_run

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

  integer         :: limcnv                     ! top interface level limit for convection [index]
                                                ! derived from reference pressures to below 40 mb

  ! internal parameters
  real(kind_phys) :: betamn                     ! minimum overshoot parameter [???]
  real(kind_phys) :: dzmin                      ! minimum convective depth for precipitation [m]
  logical         :: rlxclm                     ! control for relaxing column versus cloud triplet
  real(kind_phys) :: ssfac = 1.001_kind_phys    ! detrained air supersaturation bound [???]

  ! internal parameters for tolerance
  real(kind_phys) :: tiny = 1.0e-36_kind_phys   ! arbitrary small num in scalar transport estimates
  real(kind_phys) :: eps  = 1.0e-13_kind_phys   ! machine dependent convergence criteria
  real(kind_phys) :: tpmax = 1.50_kind_phys     ! maximum acceptable T perturbation [K]
  real(kind_phys) :: shpmax = 1.50e-3_kind_phys ! maximum acceptable Q perturbation [g/g]

  ! diagnostic only


contains
  ! Initialization of moist convective mass procedure including namelist read.
!> \section arg_table_cmfmca_init Argument Table
!! \htmlinclude cmfmca_init.html
  subroutine cmfmca_init( &
    pver, &
    cmftau_in, c0_in, &
    rair, cpair, gravit, latvap, rhoh2o_in, &
    pref_edge, &
    errmsg, errflg &
  )
    integer,            intent(in)  :: pver         ! number of vertical levels
    real(kind_phys),    intent(in)  :: cmftau_in    ! characteristic adjustment time scale [s]
    real(kind_phys),    intent(in)  :: c0_in        ! rain water autoconversion coefficient [m-1]
    real(kind_phys),    intent(in)  :: rair         ! gas constant for dry air [J K-1 kg-1]
    real(kind_phys),    intent(in)  :: cpair        ! specific heat of dry air [J K-1 kg-1]
    real(kind_phys),    intent(in)  :: gravit       ! gravitational constant [m s-2]
    real(kind_phys),    intent(in)  :: latvap       ! latent heat of vaporization [J kg-1]
    real(kind_phys),    intent(in)  :: rhoh2o_in    ! density of liquid water [kg m-3]
    real(kind_phys),    intent(in)  :: pref_edge(:) ! reference pressures at interface [Pa]

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! local variables
    integer :: k

    ! namelist variables
    cmftau  = cmftau_in
    c0      = c0_in

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
      do k = 1, pver ! was plev
        if(pref_edge(k) < 4.e3_kind_phys .and. pref_edge(k+1) >= 4.e3_kind_phys) then
          limcnv = k
        endif
      enddo
      limcnv = pver + 1 ! was plevp
    endif

    ! in-module parameters (hardcoded)
    dzmin = 0.0_kind_phys
    betamn = 0.10_kind_phys

    ! specify that relaxation timescale should be applied to column as opposed to triplets individually
    rlxclm = .true.
  end subroutine cmfmca_init


!> \section arg_table_cmfmca_run Argument Table
!! \htmlinclude cmfmca_run.html
  subroutine cmfmca_run( &
    ncol, pver, pcnst, &
    nstep, &
    ztodt, &
    pmid, pmiddry, &
    pdel, pdeldry, rpdel, rpdeldry, &
    zm, &
    tpert, qpert, &
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
    rliq_sh &
  )
    ! to_be_ccppized
    use wv_saturation,   only: qsat

    ! FIXME hplin: dependencies on CAM constituents properties ...?
    ! cnst_get_type_byind --> dry/water ...
    use constituents,    only: cnst_get_type_byind

    ! Input arguments
    integer,         intent(in)    :: ncol               ! number of atmospheric columns
    integer,         intent(in)    :: pver               ! number of vertical levels
    integer,         intent(in)    :: pcnst              ! number of ccpp constituents
    integer,         intent(in)    :: nstep              ! current time step index
    real(kind_phys), intent(in)    :: ztodt              ! 2 delta t (model time increment) [s]

    real(kind_phys), intent(in)    :: pmid(:,:)          ! midpoint pressures [Pa]
    real(kind_phys), intent(in)    :: pmiddry(:,:)       ! dry pressure at midpoints [Pa]
    real(kind_phys), intent(in)    :: pdel(:,:)          ! layer thickness (delta-p) [Pa]
    real(kind_phys), intent(in)    :: pdeldry(:,:)       ! dry layer thickness [Pa]
    real(kind_phys), intent(in)    :: rpdel(:,:)         ! 1.0 / pdel
    real(kind_phys), intent(in)    :: rpdeldry(:,:)      ! 1.0 / pdeldry

    real(kind_phys), intent(in)    :: zm(:,:)            ! geopotential height at midpoints [m]
    real(kind_phys), intent(in)    :: qpert(:,:)         ! PBL perturbation specific humidity (convective humidity excess) [kg/kg]
    real(kind_phys), intent(in)    :: phis(:)            ! surface geopotential [m2 s-2]
    real(kind_phys), intent(in)    :: pblh(:)            ! PBL height [m]
    real(kind_phys), intent(in)    :: t(:,:)             ! temperature [K]
    real(kind_phys), intent(in)    :: q(:,:,:)           ! specific humidity and constituents [kg/kg?]

    ! Output arguments
    real(kind_phys), intent(out)   :: dq(:,:,:)          ! constituent tendencies [kg kg-1 s-1]
    real(kind_phys), intent(out)   :: qc_sh(:,:)         ! dq/dt due to export of cloud water / shallow reserved cloud condensate [kg kg-1 s-1]
    real(kind_phys), intent(out)   :: cmfdt(:,:)         ! heating rate (to ptend%s) [J kg-1 s-1]
    real(kind_phys), intent(out)   :: cmfmc_sh(:,:)      ! convective updraft mass flux, shallow [kg s-1 m-2]
    real(kind_phys), intent(out)   :: cmfdqr(:,:)        ! q tendency due to shallow convective rainout [kg kg-1 s-1]
    real(kind_phys), intent(out)   :: cmfsl(:,:)         ! moist shallow convection liquid water static energy flux [W m-2]
    real(kind_phys), intent(out)   :: cmflq(:,:)         ! moist shallow convection total water flux [W m-2]
    real(kind_phys), intent(out)   :: precc(:)           ! shallow convective precipitation rate [m s-1]
    real(kind_phys), intent(out)   :: cnt_sh(:)          ! top level of shallow convective activity [index]
    real(kind_phys), intent(out)   :: cnb_sh(:)          ! bottom level of shallow convective activity [index]
    real(kind_phys), intent(out)   :: icwmr(:,:)         ! shallow convection in-cloud water mixing ratio [kg kg-1]
    real(kind_phys), intent(out)   :: rliq_sh(:)         ! vertically-integrated shallow reserved cloud condensate [m s-1]

    ! Local variables
    real(kind_phys)                :: tpert(:)           ! PBL perturbation temperature (convective temperature excess) [K]

    real(kind_phys) :: pm(ncol,pver)       ! pressure
    real(kind_phys) :: pd(ncol,pver)       ! delta-p
    real(kind_phys) :: rpd(ncol,pver)      ! 1./pdel
    real(kind_phys) :: cmfdq(ncol,pver)    ! dq/dt due to moist convection (later copied to dq(:,:,1))
    real(kind_phys) :: gam(ncol,pver)      ! 1/cp (d(qsat)/dT)
    real(kind_phys) :: sb(ncol,pver)       ! dry static energy (s bar)
    real(kind_phys) :: hb(ncol,pver)       ! moist static energy (h bar)
    real(kind_phys) :: shbs(ncol,pver)     ! sat. specific humidity (sh bar star)
    real(kind_phys) :: hbs(ncol,pver)      ! sat. moist static energy (h bar star)
    real(kind_phys) :: shbh(ncol,pverp)    ! specific humidity on interfaces
    real(kind_phys) :: sbh(ncol,pverp)     ! s bar on interfaces
    real(kind_phys) :: hbh(ncol,pverp)     ! h bar on interfaces
    real(kind_phys) :: cmrh(ncol,pverp)    ! interface constituent mixing ratio
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

    !---------------------------------------------------
    ! Initialize output tendencies
    !---------------------------------------------------
    ! q is copied to dq. dq is used for passive tracer transport
    if(pcnst > 1) then
      dq    (:ncol,:,2:)  = q(:ncol,:,2:)
    endif
    cmfdt   (:ncol,:)     = 0._kind_phys
    cmfdq   (:ncol,:)     = 0._kind_phys
    cmfmc_sh(:ncol,:)     = 0._kind_phys
    cmfdqr  (:ncol,:)     = 0._kind_phys
    cmfsl   (:ncol,:)     = 0._kind_phys
    cmflq   (:ncol,:)     = 0._kind_phys
    qc_sh   (:ncol,:)     = 0._kind_phys
    rliq_sh (:ncol)       = 0._kind_phys

    !---------------------------------------------------
    ! Quantity preparations from convect_shallow.F90.
    !---------------------------------------------------

    ! convect_shallow.F90 is not linked to pbuf tpert and always sets to zero.
    ! "This field probably should reference the pbuf tpert field but it doesnt"
    tpert(:ncol) = 0.0_kind_phys

    ! Reset PBL perturbation in constituents other than water to zero. (appears to be destructive op.)
    qpert(:ncol,2:pcnst) = 0.0_kind_phys

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
          shb(i,k) = q(i,k,1)
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
       cnt_sh(i)= real(pver+1,r8)
    end do

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
           shprme = min(min(qpert(i,1),shpmax)*fac1,max(qsattp-shb(i,k+1),0.0_kind_phys))
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
      const_modify_loop: do m=2,pcnst    ! note: indexing assumes water is first field
        ! FIXME hplin: use CCPP constituents object here
        if (cnst_get_type_byind(m) .eq. 'dry') then
           pd(:ncol,:) = pdeldry(:ncol,:)
           rpd(:ncol,:) = rpdeldry(:ncol,:)
           pm(:ncol,:) = pmiddry(:ncol,:)
        else
           pd(:ncol,:) = pdel(:ncol,:)
           rpd(:ncol,:) = rpdel(:ncol,:)
           pm(:ncol,:) = pmid(:ncol,:)
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
             cmrc(i) = dq(i,k+1,m) + qpert(i,m)*fac1
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
      if (k == limcnv+1) then
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
      else
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
    end do kloop

    !---------------------------------------------------
    ! apply final thermodynamic tendencies
    !---------------------------------------------------
    ! Set output q tendencies
    dq(:ncol,:,1 ) = cmfdq(:ncol,:)
    dq(:ncol,:,2:) = (dq(:ncol,:,2:) - q(:ncol,:,2:))/ztodt

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
    rliq_sh(:ncol) = rliq_sh(:ncol) / 1000._kind_phys
  end subroutine cmfmca_run

  pure function qhalf(sh1,sh2,shbs1,shbs2) result(qh)
    real(kind_phys), intent(in) :: sh1
    real(kind_phys), intent(in) :: sh2
    real(kind_phys), intent(in) :: shbs1
    real(kind_phys), intent(in) :: shbs2
    real(kind_phys) :: qh
    qh = min(max(sh1,sh2),(shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))
  end function qhalf
end module cmfmca
