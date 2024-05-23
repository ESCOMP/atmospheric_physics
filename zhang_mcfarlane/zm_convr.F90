module zm_convr

  use ccpp_kinds, only:  kind_phys

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zm_convr_init                 ! ZM schemea
  public zm_convr_run                  ! ZM schemea

   real(kind_phys) rl         ! wg latent heat of vaporization.
   real(kind_phys) cpres      ! specific heat at constant pressure in j/kg-degk.
   real(kind_phys) :: capelmt ! namelist configurable:
                       ! threshold value for cape for deep convection.
   real(kind_phys) :: ke           ! Tunable evaporation efficiency set from namelist input zmconv_ke
   real(kind_phys) :: ke_lnd
   real(kind_phys) :: c0_lnd       ! set from namelist input zmconv_c0_lnd
   real(kind_phys) :: c0_ocn       ! set from namelist input zmconv_c0_ocn
   integer  :: num_cin      ! set from namelist input zmconv_num_cin
                            ! The number of negative buoyancy regions that are allowed
                            ! before the convection top and CAPE calculations are completed.
   logical  :: zm_org
   real(kind_phys) tau   ! convective time scale
   real(kind_phys) :: tfreez
   real(kind_phys) :: eps1
   real(kind_phys) :: momcu
   real(kind_phys) :: momcd

   logical :: no_deep_pbl ! default = .false.
                          ! no_deep_pbl = .true. eliminates deep convection entirely within PBL


   real(kind_phys) :: rgrav       ! reciprocal of grav
   real(kind_phys) :: rgas        ! gas constant for dry air
   real(kind_phys) :: grav        ! = gravit
   real(kind_phys) :: cp          ! = cpres = cpair

   integer  limcnv       ! top interface level limit for convection

   logical :: lparcel_pbl     ! Switch to turn on mixing of parcel MSE air, and picking launch level to be the top of the PBL.


   real(kind_phys) :: tiedke_add      ! namelist configurable
   real(kind_phys) :: dmpdz_param     ! namelist configurable

contains



!===============================================================================
!> \section arg_table_zm_convr_init Argument Table
!! \htmlinclude zm_convr_init.html
!!
subroutine zm_convr_init(cpair, epsilo, gravit, latvap, tmelt, rair, &
                    limcnv_in, zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_ke_lnd, &
                    zmconv_momcu, zmconv_momcd, zmconv_num_cin, zmconv_org, &
                    no_deep_pbl_in, zmconv_tiedke_add, &
                    zmconv_capelmt, zmconv_dmpdz, zmconv_parcel_pbl, zmconv_tau, &
                    masterproc, iulog, errmsg, errflg)

   real(kind_phys), intent(in)   :: cpair           ! specific heat of dry air (J K-1 kg-1)
   real(kind_phys), intent(in)   :: epsilo          ! ratio of h2o to dry air molecular weights
   real(kind_phys), intent(in)   :: gravit          ! gravitational acceleration (m s-2)
   real(kind_phys), intent(in)   :: latvap          ! Latent heat of vaporization (J kg-1)
   real(kind_phys), intent(in)   :: tmelt           ! Freezing point of water (K)
   real(kind_phys), intent(in)   :: rair            ! Dry air gas constant     (J K-1 kg-1)
   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection
   integer, intent(in)           :: zmconv_num_cin  ! Number negative buoyancy regions that are allowed
                                                    ! before the convection top and CAPE calculations are completed.
   real(kind_phys),intent(in)           :: zmconv_c0_lnd
   real(kind_phys),intent(in)           :: zmconv_c0_ocn
   real(kind_phys),intent(in)           :: zmconv_ke
   real(kind_phys),intent(in)           :: zmconv_ke_lnd
   real(kind_phys),intent(in)           :: zmconv_momcu
   real(kind_phys),intent(in)           :: zmconv_momcd
   logical, intent(in)           :: zmconv_org
   logical, intent(in)           :: no_deep_pbl_in  ! no_deep_pbl = .true. eliminates ZM convection entirely within PBL
   real(kind_phys),intent(in)           :: zmconv_tiedke_add
   real(kind_phys),intent(in)           :: zmconv_capelmt
   real(kind_phys),intent(in)           :: zmconv_dmpdz
   logical, intent(in)           :: zmconv_parcel_pbl ! Should the parcel properties include PBL mixing?
   real(kind_phys),intent(in)           :: zmconv_tau
   logical, intent(in)                  :: masterproc
   integer, intent(in)                  :: iulog
   character(len=512), intent(out)      :: errmsg
   integer, intent(out)                 :: errflg

   errmsg =''
   errflg = 0

   ! Initialization of ZM constants
   limcnv = limcnv_in
   tfreez = tmelt
   eps1   = epsilo
   rl     = latvap
   cpres  = cpair
   rgrav  = 1.0_kind_phys/gravit
   rgas   = rair
   grav   = gravit
   cp     = cpres

   c0_lnd  = zmconv_c0_lnd
   c0_ocn  = zmconv_c0_ocn
   num_cin = zmconv_num_cin
   ke      = zmconv_ke
   ke_lnd  = zmconv_ke_lnd
   zm_org  = zmconv_org
   momcu   = zmconv_momcu
   momcd   = zmconv_momcd

   tiedke_add = zmconv_tiedke_add
   capelmt = zmconv_capelmt
   dmpdz_param = zmconv_dmpdz
   no_deep_pbl = no_deep_pbl_in
   lparcel_pbl = zmconv_parcel_pbl

   tau = zmconv_tau

   if ( masterproc ) then
      write(iulog,*) 'tuning parameters zm_convr_init: tau',tau
      write(iulog,*) 'tuning parameters zm_convr_init: c0_lnd',c0_lnd, ', c0_ocn', c0_ocn
      write(iulog,*) 'tuning parameters zm_convr_init: num_cin', num_cin
      write(iulog,*) 'tuning parameters zm_convr_init: ke',ke
      write(iulog,*) 'tuning parameters zm_convr_init: no_deep_pbl',no_deep_pbl
      write(iulog,*) 'tuning parameters zm_convr_init: zm_capelmt', capelmt
      write(iulog,*) 'tuning parameters zm_convr_init: zm_dmpdz', dmpdz_param
      write(iulog,*) 'tuning parameters zm_convr_init: zm_tiedke_add', tiedke_add
      write(iulog,*) 'tuning parameters zm_convr_init: zm_parcel_pbl', lparcel_pbl
   endif

   if (masterproc) write(iulog,*)'**** ZM: DILUTE Buoyancy Calculation ****'

end subroutine zm_convr_init


!===============================================================================
!> \section arg_table_zm_convr_run Argument Table
!! \htmlinclude zm_convr_run.html
!!
subroutine zm_convr_run(     ncol    ,pver    , &
                    pverp,   gravit  ,latice  ,cpwv    ,cpliq   , rh2o, &
                    t       ,qh      ,prec    , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    ,dpp     , &
                    delt    ,mcon    ,cme     ,cape    , &
                    tpert   ,dlf     ,zdu     ,rprd    , &
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    ql      ,rliq    ,landfrac,                   &
                    org     ,orgt    ,org2d   ,  &
                    rice   ,errmsg  ,errflg)
!-----------------------------------------------------------------------
!
! Purpose:
! Main driver for zhang-mcfarlane convection scheme
!
! Method:
! performs deep convective adjustment based on mass-flux closure
! algorithm.
!
! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
!
! This is contributed code not fully standardized by the CAM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version
! of the CAM where it will include a more straightforward formulation
! and will make use of the standard CAM nomenclature
!
!-----------------------------------------------------------------------
!
! ************************ index of variables **********************
!
!  wg * alpha    array of vertical differencing used (=1. for upstream).
!  w  * cape     convective available potential energy.
!  wg * capeg    gathered convective available potential energy.
!  c  * capelmt  threshold value for cape for deep convection.
!  ic  * cpres    specific heat at constant pressure in j/kg-degk.
!  i  * dpp
!  ic  * delt     length of model time-step in seconds.
!  wg * dp       layer thickness in mbs (between upper/lower interface).
!  wg * dqdt     mixing ratio tendency at gathered points.
!  wg * dsdt     dry static energy ("temp") tendency at gathered points.
!  wg * dudt     u-wind tendency at gathered points.
!  wg * dvdt     v-wind tendency at gathered points.
!  wg * dsubcld  layer thickness in mbs between lcl and maxi.
!  ic  * grav     acceleration due to gravity in m/sec2.
!  wg * du       detrainment in updraft. specified in mid-layer
!  wg * ed       entrainment in downdraft.
!  wg * eu       entrainment in updraft.
!  wg * hmn      moist static energy.
!  wg * hsat     saturated moist static energy.
!  w  * ideep    holds position of gathered points vs longitude index.
!  ic  * pver     number of model levels.
!  wg * j0       detrainment initiation level index.
!  wg * jd       downdraft   initiation level index.
!  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
!  wg * jt       top  level index of deep cumulus convection.
!  w  * lcl      base level index of deep cumulus convection.
!  wg * lclg     gathered values of lcl.
!  w  * lel      index of highest theoretical convective plume.
!  wg * lelg     gathered values of lel.
!  w  * lon      index of onset level for deep convection.
!  w  * maxi     index of level with largest moist static energy.
!  wg * maxg     gathered values of maxi.
!  wg * mb       cloud base mass flux.
!  wg * mc       net upward (scaled by mb) cloud mass flux.
!  wg * md       downward cloud mass flux (positive up).
!  wg * mu       upward   cloud mass flux (positive up). specified
!                at interface
!  ic  * msg      number of missing moisture levels at the top of model.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  i  * pblt     row of pbl top indices.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  wg * qhat     grid slice of upper interface mixing ratio.
!  wg * ql       grid slice of cloud liquid water.
!  wg * qs       grid slice of saturation mixing ratio.
!  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
!  wg * qstpg    grid slice of gathered values of qstp.
!  wg * qu       grid slice of mixing ratio in updraft.
!  ic  * rgas     dry air gas constant.
!  wg * rl       latent heat of vaporization.
!  w  * s        grid slice of scaled dry static energy (t+gz/cp).
!  wg * sd       grid slice of dry static energy in downdraft.
!  wg * sg       grid slice of gathered values of s.
!  wg * shat     grid slice of upper interface dry static energy.
!  wg * su       grid slice of dry static energy in updraft.
!  i/o * t
!  wg * tg       grid slice of gathered values of t.
!  w  * tl       row of parcel temperature at lcl.
!  wg * tlg      grid slice of gathered values of tl.
!  w  * tp       grid slice of parcel temperatures.
!  wg * tpg      grid slice of gathered values of tp.
!  i/o * u        grid slice of u-wind (real).
!  wg * ug       grid slice of gathered values of u.
!  i/o * utg      grid slice of u-wind tendency (real).
!  i/o * v        grid slice of v-wind (real).
!  w  * va       work array re-used by called subroutines.
!  wg * vg       grid slice of gathered values of v.
!  i/o * vtg      grid slice of v-wind tendency (real).
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!  wg * zfg      grid slice of gathered values of zf.
!  wg * zg       grid slice of gathered values of z.
!
!-----------------------------------------------------------------------
!
! multi-level i/o fields:
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.
!
! input arguments
!
   integer, intent(in) :: ncol                    ! number of atmospheric columns
   integer, intent(in) :: pver, pverp

   real(kind_phys), intent(in) :: gravit          ! gravitational acceleration (m s-2)
   real(kind_phys), intent(in) :: latice          ! Latent heat of fusion (J kg-1)
   real(kind_phys), intent(in) :: cpwv            ! specific heat of water vapor (J K-1 kg-1)
   real(kind_phys), intent(in) :: cpliq           ! specific heat of fresh h2o (J K-1 kg-1)
   real(kind_phys), intent(in) :: rh2o            ! Water vapor gas constant (J K-1 kg-1)

   real(kind_phys), intent(in) :: t(:,:)          ! grid slice of temperature at mid-layer.           (ncol,pver)
   real(kind_phys), intent(in) :: qh(:,:)         ! grid slice of specific humidity.                  (ncol,pver)
   real(kind_phys), intent(in) :: pap(:,:)        !                                                   (ncol,pver)
   real(kind_phys), intent(in) :: paph(:,:)     !                                                     (ncol,pver+1)
   real(kind_phys), intent(in) :: dpp(:,:)        ! local sigma half-level thickness (i.e. dshj).     (ncol,pver)
   real(kind_phys), intent(in) :: zm(:,:)       !                                                     (ncol,pver)
   real(kind_phys), intent(in) :: geos(:)       !                                                     (ncol)
   real(kind_phys), intent(in) :: zi(:,:)      !                                                      (ncol,pver+1)
   real(kind_phys), intent(in) :: pblh(:)     !                                                       (ncol)
   real(kind_phys), intent(in) :: tpert(:)      !                                                     (ncol)
   real(kind_phys), intent(in) :: landfrac(:) ! RBN Landfrac                                          (ncol)

! output arguments
!
   real(kind_phys), intent(out) :: qtnd(:,:)           ! specific humidity tendency (kg/kg/s)             (ncol,pver)
   real(kind_phys), intent(out) :: heat(:,:)           ! heating rate (dry static energy tendency, W/kg)  (ncol,pver)
   real(kind_phys), intent(out) :: mcon(:,:)  !   (ncol,pverp)
   real(kind_phys), intent(out) :: dlf(:,:)    ! scattrd version of the detraining cld h2o tend (ncol,pver)
   real(kind_phys), intent(out) :: cme(:,:)    !                                                          (ncol,pver)
   real(kind_phys), intent(out) :: cape(:)        ! w  convective available potential energy.             (ncol)
   real(kind_phys), intent(out) :: zdu(:,:)    ! (ncol,pver)
   real(kind_phys), intent(out) :: rprd(:,:)     ! rain production rate (ncol,pver)

! move these vars from local storage to output so that convective
! transports can be done in outside of conv_cam.
   real(kind_phys), intent(out) :: mu(:,:)  !                                                                 (ncol,pver)
   real(kind_phys), intent(out) :: eu(:,:)  !                                                                 (ncol,pver)
   real(kind_phys), intent(out) :: du(:,:)  !                                                                 (ncol,pver)
   real(kind_phys), intent(out) :: md(:,:)  !                                                                 (ncol,pver)
   real(kind_phys), intent(out) :: ed(:,:)  !                                                                 (ncol,pver)
   real(kind_phys), intent(out) :: dp(:,:)       ! wg layer thickness in mbs (between upper/lower interface). (ncol,pver)
   real(kind_phys), intent(out) :: dsubcld(:)       ! wg layer thickness in mbs between lcl and maxi.         (ncol)
   real(kind_phys), intent(out) :: prec(:)  !                                                                 (ncol)
   real(kind_phys), intent(out) :: rliq(:) ! reserved liquid (not yet in cldliq) for energy integrals         (ncol)
   real(kind_phys), intent(out) :: rice(:) ! reserved ice (not yet in cldce) for energy integrals             (ncol)

   integer,  intent(out) :: ideep(:)  ! column indices of gathered points                              (ncol)

   character(len=512), intent(out)      :: errmsg
   integer, intent(out)                 :: errflg

   real(kind_phys), intent(in)  :: org(:,:)     ! Only used if zm_org is true  ! in
   real(kind_phys), intent(out) :: orgt(:,:)   ! Only used if zm_org is true   ! out
   real(kind_phys), intent(out) :: org2d(:,:)  ! Only used if zm_org is true   ! out

   ! Local variables

   real(kind_phys) zs(ncol)
   real(kind_phys) dlg(ncol,pver)    ! gathrd version of the detraining cld h2o tend
   real(kind_phys) cug(ncol,pver)    ! gathered condensation rate

   real(kind_phys) evpg(ncol,pver)   ! gathered evap rate of rain in downdraft
   real(kind_phys) orgavg(ncol)
   real(kind_phys) dptot(ncol)

   real(kind_phys) mumax(ncol)

!CACNOTE - Figure out real intent for jt and maxg
   integer, intent(inout) :: jt(ncol)                          ! wg top  level index of deep cumulus convection.
   integer, intent(inout) :: maxg(ncol)                        ! wg gathered values of maxi.

   integer lengath
!     diagnostic field used by chem/wetdep codes

!CACNOTE - Figure out real intent for ql
   real(kind_phys),intent(inout):: ql(ncol,pver)                    ! wg grid slice of cloud liquid water.
!
   real(kind_phys) pblt(ncol)           ! i row of pbl top indices.




!
!-----------------------------------------------------------------------
!
! general work fields (local variables):
!
   real(kind_phys) q(ncol,pver)              ! w  grid slice of mixing ratio.
   real(kind_phys) p(ncol,pver)              ! w  grid slice of ambient mid-layer pressure in mbs.
   real(kind_phys) z(ncol,pver)              ! w  grid slice of ambient mid-layer height in metres.
   real(kind_phys) s(ncol,pver)              ! w  grid slice of scaled dry static energy (t+gz/cp).
   real(kind_phys) tp(ncol,pver)             ! w  grid slice of parcel temperatures.
   real(kind_phys) zf(ncol,pver+1)           ! w  grid slice of ambient interface height in metres.
   real(kind_phys) pf(ncol,pver+1)           ! w  grid slice of ambient interface pressure in mbs.
   real(kind_phys) qstp(ncol,pver)           ! w  grid slice of parcel temp. saturation mixing ratio.

   real(kind_phys) tl(ncol)                  ! w  row of parcel temperature at lcl.

   integer lcl(ncol)                  ! w  base level index of deep cumulus convection.
   integer lel(ncol)                  ! w  index of highest theoretical convective plume.
   integer lon(ncol)                  ! w  index of onset level for deep convection.
   integer maxi(ncol)                 ! w  index of level with largest moist static energy.

   real(kind_phys) precip
!
! gathered work fields:
!
   real(kind_phys) qg(ncol,pver)             ! wg grid slice of gathered values of q.
   real(kind_phys) tg(ncol,pver)             ! w  grid slice of temperature at interface.
   real(kind_phys) pg(ncol,pver)             ! wg grid slice of gathered values of p.
   real(kind_phys) zg(ncol,pver)             ! wg grid slice of gathered values of z.
   real(kind_phys) sg(ncol,pver)             ! wg grid slice of gathered values of s.
   real(kind_phys) tpg(ncol,pver)            ! wg grid slice of gathered values of tp.
   real(kind_phys) zfg(ncol,pver+1)          ! wg grid slice of gathered values of zf.
   real(kind_phys) qstpg(ncol,pver)          ! wg grid slice of gathered values of qstp.
   real(kind_phys) ug(ncol,pver)             ! wg grid slice of gathered values of u.
   real(kind_phys) vg(ncol,pver)             ! wg grid slice of gathered values of v.
   real(kind_phys) cmeg(ncol,pver)

   real(kind_phys) rprdg(ncol,pver)           ! wg gathered rain production rate
   real(kind_phys) capeg(ncol)               ! wg gathered convective available potential energy.
   real(kind_phys) tlg(ncol)                 ! wg grid slice of gathered values of tl.
   real(kind_phys) landfracg(ncol)            ! wg grid slice of landfrac

   integer lclg(ncol)       ! wg gathered values of lcl.
   integer lelg(ncol)
!
! work fields arising from gathered calculations.
!
   real(kind_phys) dqdt(ncol,pver)           ! wg mixing ratio tendency at gathered points.
   real(kind_phys) dsdt(ncol,pver)           ! wg dry static energy ("temp") tendency at gathered points.
   real(kind_phys) sd(ncol,pver)             ! wg grid slice of dry static energy in downdraft.
   real(kind_phys) qd(ncol,pver)             ! wg grid slice of mixing ratio in downdraft.
   real(kind_phys) mc(ncol,pver)             ! wg net upward (scaled by mb) cloud mass flux.
   real(kind_phys) qhat(ncol,pver)           ! wg grid slice of upper interface mixing ratio.
   real(kind_phys) qu(ncol,pver)             ! wg grid slice of mixing ratio in updraft.
   real(kind_phys) su(ncol,pver)             ! wg grid slice of dry static energy in updraft.
   real(kind_phys) qs(ncol,pver)             ! wg grid slice of saturation mixing ratio.
   real(kind_phys) shat(ncol,pver)           ! wg grid slice of upper interface dry static energy.
   real(kind_phys) hmn(ncol,pver)            ! wg moist static energy.
   real(kind_phys) hsat(ncol,pver)           ! wg saturated moist static energy.
   real(kind_phys) qlg(ncol,pver)
   real(kind_phys) dudt(ncol,pver)           ! wg u-wind tendency at gathered points.
   real(kind_phys) dvdt(ncol,pver)           ! wg v-wind tendency at gathered points.

   real(kind_phys) qldeg(ncol,pver)        ! cloud liquid water mixing ratio for detrainment (kg/kg)
   real(kind_phys) mb(ncol)                ! wg cloud base mass flux.

   integer jlcl(ncol)
   integer j0(ncol)                 ! wg detrainment initiation level index.
   integer jd(ncol)                 ! wg downdraft initiation level index.

   real(kind_phys),intent(in):: delt                     ! length of model time-step in seconds.

   integer i
   integer ii
   integer k, kk, l, m

   integer msg                      !  ic number of missing moisture levels at the top of model.
   real(kind_phys) qdifr
   real(kind_phys) sdifr

   real(kind_phys), parameter :: dcon  = 25.e-6_kind_phys
   real(kind_phys), parameter :: mucon = 5.3_kind_phys
   real(kind_phys) negadq
   logical doliq


!
!--------------------------Data statements------------------------------

   errmsg = ''
   errflg = 0
!
! Set internal variable "msg" (convection limit) to "limcnv-1"
!
   msg = limcnv - 1
!
! initialize necessary arrays.
! zero out variables not used in cam
!

   if (zm_org) then
      orgt(:,:) = 0._kind_phys
   end if

   qtnd(:,:) = 0._kind_phys
   heat(:,:) = 0._kind_phys
   mcon(:,:) = 0._kind_phys
   rliq(:ncol)   = 0._kind_phys
   rice(:ncol)   = 0._kind_phys

!
! initialize convective tendencies
!
   prec(:ncol) = 0._kind_phys
   do k = 1,pver
      do i = 1,ncol
         dqdt(i,k)  = 0._kind_phys
         dsdt(i,k)  = 0._kind_phys
         dudt(i,k)  = 0._kind_phys
         dvdt(i,k)  = 0._kind_phys
         cme(i,k)   = 0._kind_phys
         rprd(i,k)  = 0._kind_phys
         zdu(i,k)   = 0._kind_phys
         ql(i,k)    = 0._kind_phys
         qlg(i,k)   = 0._kind_phys
         dlf(i,k)   = 0._kind_phys
         dlg(i,k)   = 0._kind_phys
         qldeg(i,k) = 0._kind_phys

      end do
   end do

   do i = 1,ncol
      pblt(i) = pver
      dsubcld(i) = 0._kind_phys

   end do

  if (zm_org) then
! compute vertical average here
      orgavg(:) = 0._kind_phys
      dptot(:) = 0._kind_phys

      do k = 1, pver
        do i = 1,ncol
          if (org(i,k) .gt. 0) then
            orgavg(i) = orgavg(i)+dpp(i,k)*org(i,k)
            dptot(i) = dptot(i)+dpp(i,k)
          endif
        enddo
      enddo

      do i = 1,ncol
        if (dptot(i) .gt. 0) then
          orgavg(i) = orgavg(i)/dptot(i)
        endif
      enddo

      do k = 1, pver
        do i = 1, ncol
           org2d(i,k) = orgavg(i)
        enddo
      enddo

   endif

!
! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.
!
   do i = 1,ncol
      zs(i) = geos(i)*rgrav
      pf(i,pver+1) = paph(i,pver+1)*0.01_kind_phys
      zf(i,pver+1) = zi(i,pver+1) + zs(i)
   end do
   do k = 1,pver
      do i = 1,ncol
         p(i,k) = pap(i,k)*0.01_kind_phys
         pf(i,k) = paph(i,k)*0.01_kind_phys
         z(i,k) = zm(i,k) + zs(i)
         zf(i,k) = zi(i,k) + zs(i)
      end do
   end do
!
   do k = pver - 1,msg + 1,-1
      do i = 1,ncol
         if (abs(z(i,k)-zs(i)-pblh(i)) < (zf(i,k)-zf(i,k+1))*0.5_kind_phys) pblt(i) = k
      end do
   end do
!
! store incoming specific humidity field for subsequent calculation
! of precipitation (through change in storage).
! define dry static energy (normalized by cp).
!
   do k = 1,pver
      do i = 1,ncol
         q(i,k) = qh(i,k)
         s(i,k) = t(i,k) + (grav/cpres)*z(i,k)
         tp(i,k)=0.0_kind_phys
         shat(i,k) = s(i,k)
         qhat(i,k) = q(i,k)
      end do
   end do

   do i = 1,ncol
      capeg(i) = 0._kind_phys
      lclg(i) = 1
      lelg(i) = pver
      maxg(i) = 1
      tlg(i) = 400._kind_phys
      dsubcld(i) = 0._kind_phys
   end do


   !  Evaluate Tparcel, qs(Tparcel), buoyancy and CAPE,
   !     lcl, lel, parcel launch level at index maxi()=hmax

   call buoyan_dilute(  ncol   ,pver     , &
               cpliq   ,latice  ,cpwv    ,rh2o    ,&
               q       ,t       ,p       ,z       ,pf       , &
               tp      ,qstp    ,tl      ,rl      ,cape     , &
               pblt    ,lcl     ,lel     ,lon     ,maxi     , &
               rgas    ,grav    ,cpres   ,msg     , &
               zi      ,zs      ,tpert   , org2d  , landfrac,&
               errmsg  ,errflg)

!
! determine whether grid points will undergo some deep convection
! (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
! (require cape.gt. 0 and lel<lcl as minimum conditions).
!
   lengath = 0
   ideep   = 0
   do i=1,ncol
      if (cape(i) > capelmt) then
         lengath = lengath + 1
         ideep(lengath) = i
      end if
   end do

   if (lengath.eq.0) return
!
! obtain gathered arrays necessary for ensuing calculations.
!
   do k = 1,pver
      do i = 1,lengath
         dp(i,k) = 0.01_kind_phys*dpp(ideep(i),k)
         qg(i,k) = q(ideep(i),k)
         tg(i,k) = t(ideep(i),k)
         pg(i,k) = p(ideep(i),k)
         zg(i,k) = z(ideep(i),k)
         sg(i,k) = s(ideep(i),k)
         tpg(i,k) = tp(ideep(i),k)
         zfg(i,k) = zf(ideep(i),k)
         qstpg(i,k) = qstp(ideep(i),k)
         ug(i,k) = 0._kind_phys
         vg(i,k) = 0._kind_phys
      end do
   end do

!
   do i = 1,lengath
      zfg(i,pver+1) = zf(ideep(i),pver+1)
   end do
   do i = 1,lengath
      capeg(i) = cape(ideep(i))
      lclg(i) = lcl(ideep(i))
      lelg(i) = lel(ideep(i))
      maxg(i) = maxi(ideep(i))
      tlg(i) = tl(ideep(i))
      landfracg(i) = landfrac(ideep(i))
   end do
!
! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.
!
   do k = msg + 1,pver
      do i = 1,lengath
         if (k >= maxg(i)) then
            dsubcld(i) = dsubcld(i) + dp(i,k)
         end if
      end do
   end do
!
! define array of factors (alpha) which defines interfacial
! values, as well as interfacial values for (q,s) used in
! subsequent routines.
!
   do k = msg + 2,pver
      do i = 1,lengath
         sdifr = 0._kind_phys
         qdifr = 0._kind_phys
         if (sg(i,k) > 0._kind_phys .or. sg(i,k-1) > 0._kind_phys) &
            sdifr = abs((sg(i,k)-sg(i,k-1))/max(sg(i,k-1),sg(i,k)))
         if (qg(i,k) > 0._kind_phys .or. qg(i,k-1) > 0._kind_phys) &
            qdifr = abs((qg(i,k)-qg(i,k-1))/max(qg(i,k-1),qg(i,k)))
         if (sdifr > 1.E-6_kind_phys) then
            shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/(sg(i,k-1)-sg(i,k))
         else
            shat(i,k) = 0.5_kind_phys* (sg(i,k)+sg(i,k-1))
         end if
         if (qdifr > 1.E-6_kind_phys) then
            qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/(qg(i,k-1)-qg(i,k))
         else
            qhat(i,k) = 0.5_kind_phys* (qg(i,k)+qg(i,k-1))
         end if
      end do
   end do
!
! obtain cloud properties.
!

   call cldprp(ncol   ,pver    ,pverp   ,cpliq  , &
               latice  ,cpwv    ,rh2o    ,&
               qg      ,tg      ,ug      ,vg      ,pg      , &
               zg      ,sg      ,mu      ,eu      ,du      , &
               md      ,ed      ,sd      ,qd      ,mc      , &
               qu      ,su      ,zfg     ,qs      ,hmn     , &
               hsat    ,shat    ,qlg     , &
               cmeg    ,maxg    ,lelg    ,jt      ,jlcl    , &
               maxg    ,j0      ,jd      ,rl      ,lengath , &
               rgas    ,grav    ,cpres   ,msg     , &
               evpg    ,cug     ,rprdg   ,limcnv  ,landfracg , &
               qldeg    ,qhat    )


!
! convert detrainment from units of "1/m" to "1/mb".
!

   do k = msg + 1,pver
      do i = 1,lengath
         du   (i,k) = du   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         eu   (i,k) = eu   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         ed   (i,k) = ed   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cug  (i,k) = cug  (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cmeg (i,k) = cmeg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         rprdg(i,k) = rprdg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         evpg (i,k) = evpg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
      end do
   end do

   call closure(ncol   ,pver    , &
                qg      ,tg      ,pg      ,zg      ,sg      , &
                tpg     ,qs      ,qu      ,su      ,mc      , &
                du      ,mu      ,md      ,qd      ,sd      , &
                qhat    ,shat    ,dp      ,qstpg   ,zfg     , &
                qlg     ,dsubcld ,mb      ,capeg   ,tlg     , &
                lclg    ,lelg    ,jt      ,maxg    ,1       , &
                lengath ,rgas    ,grav    ,cpres   ,rl      , &
                msg     ,capelmt    )
!
! limit cloud base mass flux to theoretical upper bound.
!
   do i=1,lengath
      mumax(i) = 0
   end do
   do k=msg + 2,pver
      do i=1,lengath
        mumax(i) = max(mumax(i), mu(i,k)/dp(i,k))
      end do
   end do

   do i=1,lengath
      if (mumax(i) > 0._kind_phys) then
         mb(i) = min(mb(i),0.5_kind_phys/(delt*mumax(i)))
      else
         mb(i) = 0._kind_phys
      endif
   end do
   ! If no_deep_pbl = .true., don't allow convection entirely
   ! within PBL (suggestion of Bjorn Stevens, 8-2000)

   if (no_deep_pbl) then
      do i=1,lengath
         if (zm(ideep(i),jt(i)) < pblh(ideep(i))) mb(i) = 0
      end do
   end if

   do k=msg+1,pver
      do i=1,lengath
         mu   (i,k)  = mu   (i,k)*mb(i)
         md   (i,k)  = md   (i,k)*mb(i)
         mc   (i,k)  = mc   (i,k)*mb(i)
         du   (i,k)  = du   (i,k)*mb(i)
         eu   (i,k)  = eu   (i,k)*mb(i)
         ed   (i,k)  = ed   (i,k)*mb(i)
         cmeg (i,k)  = cmeg (i,k)*mb(i)
         rprdg(i,k)  = rprdg(i,k)*mb(i)
         cug  (i,k)  = cug  (i,k)*mb(i)
         evpg (i,k)  = evpg (i,k)*mb(i)

      end do
   end do
!
! compute temperature and moisture changes due to convection.
!
   call q1q2_pjr(ncol   ,pver    ,latice  , &
                 dqdt    ,dsdt    ,qg      ,qs      ,qu      , &
                 su      ,du      ,qhat    ,shat    ,dp      , &
                 mu      ,md      ,sd      ,qd      ,qldeg   , &
                 dsubcld ,jt      ,maxg    ,1       ,lengath , &
                 cpres   ,rl      ,msg     ,          &
                 dlg     ,evpg    ,cug)

!
! gather back temperature and mixing ratio.
!

   do k = msg + 1,pver
      do i = 1,lengath
!
! q is updated to compute net precip.
!
         q(ideep(i),k) = qh(ideep(i),k) + 2._kind_phys*delt*dqdt(i,k)
         qtnd(ideep(i),k) = dqdt (i,k)
         cme (ideep(i),k) = cmeg (i,k)
         rprd(ideep(i),k) = rprdg(i,k)
         zdu (ideep(i),k) = du   (i,k)
         mcon(ideep(i),k) = mc   (i,k)
         heat(ideep(i),k) = dsdt (i,k)*cpres
         dlf (ideep(i),k) = dlg  (i,k)
         ql  (ideep(i),k) = qlg  (i,k)
      end do
   end do

! Compute precip by integrating change in water vapor minus detrained cloud water
   do k = pver,msg + 1,-1
      do i = 1,ncol
          prec(i) = prec(i) - dpp(i,k)* (q(i,k)-qh(i,k)) - dpp(i,k)*(dlf(i,k))*2._kind_phys*delt
      end do
   end do

! obtain final precipitation rate in m/s.
   do i = 1,ncol
      prec(i) = rgrav*max(prec(i),0._kind_phys)/ (2._kind_phys*delt)/1000._kind_phys
   end do

! Compute reserved liquid (not yet in cldliq) for energy integrals.
! Treat rliq as flux out bottom, to be added back later.
   do k = 1, pver
      do i = 1, ncol
          rliq(i) = rliq(i) + (dlf(i,k))*dpp(i,k)/gravit
          rice(i) = rice(i)/gravit
      end do
   end do
   rliq(:ncol) = rliq(:ncol) /1000._kind_phys
   rice(:ncol) = rice(:ncol) /1000._kind_phys

   return
end subroutine zm_convr_run

subroutine zm_convr_finalize
end subroutine zm_convr_finalize

!=========================================================================================

subroutine buoyan_dilute(  ncol   ,pver    , &
                  cpliq   ,latice  ,cpwv    ,rh2o    ,&
                  q       ,t       ,p       ,z       ,pf      , &
                  tp      ,qstp    ,tl      ,rl      ,cape    , &
                  pblt    ,lcl     ,lel     ,lon     ,mx      , &
                  rd      ,grav    ,cp      ,msg     , &
                  zi      ,zs      ,tpert    ,org    , landfrac,&
                  errmsg  ,errflg)
!-----------------------------------------------------------------------
!
! Purpose:
! Calculates CAPE the lifting condensation level and the convective top
! where buoyancy is first -ve.
!
! Method: Calculates the parcel temperature based on a simple constant
! entraining plume model. CAPE is integrated from buoyancy.
! 09/09/04 - Simplest approach using an assumed entrainment rate for
!            testing (dmpdp).
! 08/04/05 - Swap to convert dmpdz to dmpdp
!
! SCAM Logical Switches - DILUTE:RBN - Now Disabled
! ---------------------
! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
!
! References:
! Raymond and Blythe (1992) JAS
!
! Author:
! Richard Neale - September 2004
!
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! input arguments
!
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: pver
   real(kind_phys), intent(in) :: cpliq
   real(kind_phys), intent(in) :: latice
   real(kind_phys), intent(in) :: cpwv
   real(kind_phys), intent(in) :: rh2o

   real(kind_phys), intent(in) :: q(ncol,pver)        ! spec. humidity
   real(kind_phys), intent(in) :: t(ncol,pver)        ! temperature
   real(kind_phys), intent(in) :: p(ncol,pver)        ! pressure
   real(kind_phys), intent(in) :: z(ncol,pver)        ! height
   real(kind_phys), intent(in) :: pf(ncol,pver+1)     ! pressure at interfaces
   real(kind_phys), intent(in) :: pblt(ncol)          ! index of pbl depth
   real(kind_phys), intent(in) :: tpert(ncol)         ! perturbation temperature by pbl processes

! Use z interface/surface relative values for PBL parcel calculations.
   real(kind_phys), intent(in) :: zi(ncol,pver+1)
   real(kind_phys), intent(in) :: zs(ncol)

!
! output arguments
!

   real(kind_phys), intent(out) :: tp(ncol,pver)       ! parcel temperature
   real(kind_phys), intent(out) :: qstp(ncol,pver)     ! saturation mixing ratio of parcel (only above lcl, just q below).
   real(kind_phys), intent(out) :: tl(ncol)            ! parcel temperature at lcl
   real(kind_phys), intent(out) :: cape(ncol)          ! convective aval. pot. energy.
   integer lcl(ncol)        !
   integer lel(ncol)        !
   integer lon(ncol)        ! level of onset of deep convection
   integer mx(ncol)         ! level of max moist static energy

   real(kind_phys)  :: org(:,:)      ! organization parameter
   real(kind_phys), intent(in) :: landfrac(ncol)
   character(len=512), intent(out)      :: errmsg
   integer, intent(out)                 :: errflg

!
!--------------------------Local Variables------------------------------
!
   real(kind_phys) capeten(ncol,5)     ! provisional value of cape
   real(kind_phys) tv(ncol,pver)       !
   real(kind_phys) tpv(ncol,pver)      !
   real(kind_phys) buoy(ncol,pver)

   real(kind_phys) a1(ncol)
   real(kind_phys) a2(ncol)
   real(kind_phys) estp(ncol)
   real(kind_phys) pl(ncol)
   real(kind_phys) plexp(ncol)
   real(kind_phys) hmax(ncol)
   real(kind_phys) hmn(ncol)
   real(kind_phys) y(ncol)

   logical plge600(ncol)
   integer knt(ncol)
   integer lelten(ncol,5)




! Parcel property variables

  real(kind_phys)           :: hmn_lev(ncol,pver)  ! Vertical profile of moist static energy for each column
  real(kind_phys)           :: dp_lev(ncol,pver)   ! Level dpressure between interfaces
  real(kind_phys)           :: hmn_zdp(ncol,pver)  ! Integrals of hmn_lev*dp_lev at each level
  real(kind_phys)           :: q_zdp(ncol,pver)    ! Integrals of q*dp_lev at each level
  real(kind_phys)           :: dp_zfrac             ! Fraction of vertical grid box below mixing top (usually pblt)
  real(kind_phys)           :: parcel_dz(ncol)     ! Depth of parcel mixing (usually parcel_hscale*parcel_dz)
  real(kind_phys)           :: parcel_ztop(ncol)   ! Height of parcel mixing (usually parcel_ztop+zm(nlev))
  real(kind_phys)           :: parcel_dp(ncol)     ! Pressure integral over parcel mixing depth (usually pblt)
  real(kind_phys)           :: parcel_hdp(ncol)    ! Pressure*MSE integral over parcel mixing depth (usually pblt)
  real(kind_phys)           :: parcel_qdp(ncol)    ! Pressure*q integral over parcel mixing depth (usually pblt)
  real(kind_phys)           :: pbl_dz(ncol)        ! Previously diagnosed PBL height
  real(kind_phys)           :: hpar(ncol)          ! Initial MSE of the parcel
  real(kind_phys)           :: qpar(ncol)          ! Initial humidity of the parcel
  real(kind_phys)           :: ql(ncol)          ! Initial parcel humidity (for ientropy routine)
  integer            :: ipar ! Index for top of parcel mixing/launch level.




   real(kind_phys) cp
   real(kind_phys) e
   real(kind_phys) grav

   integer i
   integer k
   integer msg
   integer n

   real(kind_phys) rd
   real(kind_phys) rl


! Scaling of PBL height to give parcel mixing length for lparcel_pbl=True

   real(kind_phys), parameter :: parcel_hscale  = 0.5_kind_phys


!
!-----------------------------------------------------------------------
!
   do n = 1,5
      do i = 1,ncol
         lelten(i,n) = pver
         capeten(i,n) = 0._kind_phys
      end do
   end do
!
   do i = 1,ncol
      lon(i) = pver
      knt(i) = 0
      lel(i) = pver
      mx(i) = lon(i)
      cape(i) = 0._kind_phys
      hmax(i) = 0._kind_phys
      pbl_dz(i) = z(i,nint(pblt(i)))-zs(i) ! mid-point z (zm) reference to PBL depth
      parcel_dz(i) = max(zi(i,pver),parcel_hscale*pbl_dz(i)) ! PBL mixing depth [parcel_hscale*Boundary, but no thinner than zi(i,pver)]
      parcel_ztop(i) = parcel_dz(i)+zs(i) ! PBL mixing height ztop this is wrt zs=0
      parcel_hdp(i) = 0._kind_phys
      parcel_dp(i) = 0._kind_phys
      parcel_qdp(i) = 0._kind_phys
      hpar(i) = 0._kind_phys
      qpar(i) = 0._kind_phys
   end do

   tp(:ncol,:) = t(:ncol,:)
   qstp(:ncol,:) = q(:ncol,:)
   hmn_lev(:ncol,:) = 0._kind_phys



!!! Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
   tv(:ncol,:) = t(:ncol,:) *(1._kind_phys+1.608_kind_phys*q(:ncol,:))/ (1._kind_phys+q(:ncol,:))
   tpv(:ncol,:) = tv(:ncol,:)
   buoy(:ncol,:) = 0._kind_phys


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mix the parcel over a certain dp or dz and take the launch level as the top level
! of this mixing region and the parcel properties as this mixed value
! Should be well mixed by other processes in the very near PBL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if (lparcel_pbl) then

! Vertical profile of MSE and pressure weighted of the same.
   hmn_lev(:ncol,1:pver) = cp*t(:ncol,1:pver) + grav*z(:ncol,1:pver) + rl*q(:ncol,1:pver)
   dp_lev(:ncol,1:pver) = pf(:ncol,2:pver+1)-pf(:ncol,1:pver)
   hmn_zdp(:ncol,1:pver) = hmn_lev(:ncol,1:pver)*dp_lev(:ncol,1:pver)
   q_zdp(:ncol,1:pver) = q(:ncol,1:pver)*dp_lev(:ncol,1:pver)


! Mix profile over vertical length scale of 0.5*PBLH.

   do i = 1,ncol ! Loop columns
      do k = pver,msg + 1,-1

         if (zi(i,k+1)<= parcel_dz(i)) then ! Has to be relative to near-surface layer center elevation
            ipar = k

            if (k == pver) then ! Always at least the full depth of lowest model layer.
               dp_zfrac = 1._kind_phys
            else
               ! Fraction of grid cell depth (mostly 1, except when parcel_ztop is in between levels.
               dp_zfrac =  min(1._kind_phys,(parcel_dz(i)-zi(i,k+1))/(zi(i,k)-zi(i,k+1)))
            end if

            parcel_hdp(i) = parcel_hdp(i)+hmn_zdp(i,k)*dp_zfrac ! Sum parcel profile up to a certain level.
            parcel_qdp(i) = parcel_qdp(i)+q_zdp(i,k)*dp_zfrac ! Sum parcel profile up to a certain level.
            parcel_dp(i)  = parcel_dp(i)+dp_lev(i,k)*dp_zfrac ! SUM dp's for weighting of parcel_hdp

         end if
      end do
      hpar(i) = parcel_hdp(i)/parcel_dp(i)
      qpar(i) = parcel_qdp(i)/parcel_dp(i)
      mx(i) = ipar
   end do

else ! Default method finding level of MSE maximum (nlev sensitive though)
    !
    ! set "launching" level(mx) to be at maximum moist static energy.
    ! search for this level stops at planetary boundary layer top.
    !
    do k = pver,msg + 1,-1
       do i = 1,ncol
          hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
          if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
             hmax(i) = hmn(i)
             mx(i) = k
          end if
       end do
    end do

end if ! Default method of determining parcel launch properties.





! LCL dilute calculation - initialize to mx(i)
! Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
! Original code actually sets LCL as level above wher condensate forms.
! Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.

if (lparcel_pbl) then

! For parcel dilute need to invert hpar and qpar.
! Now need to supply ql(i) as it is mixed parcel version, just q(i,max(i)) in default

   do i = 1,ncol             ! Initialise LCL variables.
      lcl(i) = mx(i)
      tl(i) = (hpar(i)-rl*qpar(i)-grav*parcel_ztop(i))/cp
      ql(i) = qpar(i)
      pl(i) = p(i,mx(i))
   end do

else

   do i = 1,ncol
      lcl(i) = mx(i)
      tl(i) = t(i,mx(i))
      ql(i) = q(i,mx(i))
      pl(i) = p(i,mx(i))
   end do

end if ! Mixed parcel properties



!
! main buoyancy calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
!!!   RBN 9/9/04   !!!

   call parcel_dilute(ncol, pver, cpliq, cpwv, rh2o, latice, msg, mx, p, t, q, &
   tpert, tp, tpv, qstp, pl, tl, ql, lcl, &
   org, landfrac, errmsg, errflg)


! If lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,ncol
      plge600(i) = pl(i).ge.600._kind_phys ! Just change to always allow buoy calculation.
   end do

!
! Main buoyancy calculation.
!
   do k = pver,msg + 1,-1
      do i=1,ncol
         if (k <= mx(i) .and. plge600(i)) then   ! Define buoy from launch level to cloud top.
            tv(i,k) = t(i,k)* (1._kind_phys+1.608_kind_phys*q(i,k))/ (1._kind_phys+q(i,k))
            buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add  ! +0.5K or not?
         else
            qstp(i,k) = q(i,k)
            tp(i,k)   = t(i,k)
            tpv(i,k)  = tv(i,k)
         endif
      end do
   end do



!-------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
   do k = msg + 2,pver
      do i = 1,ncol
         if (k < lcl(i) .and. plge600(i)) then
            if (buoy(i,k+1) > 0._kind_phys .and. buoy(i,k) <= 0._kind_phys) then
               knt(i) = min(num_cin,knt(i) + 1)
               lelten(i,knt(i)) = k
            end if
         end if
      end do
   end do
!
! calculate convective available potential energy (cape).
!
   do n = 1,num_cin
      do k = msg + 1,pver
         do i = 1,ncol
            if (plge600(i) .and. k <= mx(i) .and. k > lelten(i,n)) then
               capeten(i,n) = capeten(i,n) + rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
         end do
      end do
   end do
!
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
   do n = 1,num_cin
      do i = 1,ncol
         if (capeten(i,n) > cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
         end if
      end do
   end do
!
! put lower bound on cape for diagnostic purposes.
!
   do i = 1,ncol
      cape(i) = max(cape(i), 0._kind_phys)
   end do
!
   return
end subroutine buoyan_dilute

subroutine parcel_dilute (ncol, pver, cpliq, cpwv, rh2o, latice, msg, klaunch, p, t, q, &
  tpert, tp, tpv, qstp, pl, tl, ql, lcl, &
  org, landfrac,errmsg,errflg)

! Routine  to determine
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.

!--------------------
implicit none
!--------------------

integer, intent(in) :: ncol
integer, intent(in) :: pver
real(kind_phys), intent(in) :: cpliq
real(kind_phys), intent(in) :: cpwv
real(kind_phys), intent(in) :: rh2o
real(kind_phys), intent(in) :: latice
integer, intent(in) :: msg

integer, intent(in), dimension(ncol) :: klaunch(ncol)

real(kind_phys), intent(in), dimension(ncol,pver) :: p
real(kind_phys), intent(in), dimension(ncol,pver) :: t
real(kind_phys), intent(in), dimension(ncol,pver) :: q
real(kind_phys), intent(in), dimension(ncol) :: tpert ! PBL temperature perturbation.

real(kind_phys), intent(inout), dimension(ncol,pver) :: tp    ! Parcel temp.
real(kind_phys), intent(inout), dimension(ncol,pver) :: qstp  ! Parcel water vapour (sat value above lcl).
real(kind_phys), intent(inout), dimension(ncol) :: tl         ! Actual temp of LCL.
real(kind_phys), intent(inout), dimension(ncol) :: ql ! Actual humidity of LCL
real(kind_phys), intent(inout), dimension(ncol) :: pl          ! Actual pressure of LCL.

integer, intent(inout), dimension(ncol) :: lcl ! Lifting condesation level (first model level with saturation).

real(kind_phys), intent(out), dimension(ncol,pver) :: tpv   ! Define tpv within this routine.

character(len=512), intent(out)      :: errmsg
integer, intent(out)                 :: errflg



real(kind_phys), dimension(:,:) :: org
real(kind_phys), intent(in), dimension(ncol) :: landfrac
!--------------------

! Have to be careful as s is also dry static energy.


! If we are to retain the fact that CAM loops over grid-points in the internal
! loop then we need to dimension sp,atp,mp,xsh2o with ncol.


real(kind_phys) tmix(ncol,pver)        ! Tempertaure of the entraining parcel.
real(kind_phys) qtmix(ncol,pver)       ! Total water of the entraining parcel.
real(kind_phys) qsmix(ncol,pver)       ! Saturated mixing ratio at the tmix.
real(kind_phys) smix(ncol,pver)        ! Entropy of the entraining parcel.
real(kind_phys) xsh2o(ncol,pver)       ! Precipitate lost from parcel.
real(kind_phys) ds_xsh2o(ncol,pver)    ! Entropy change due to loss of condensate.
real(kind_phys) ds_freeze(ncol,pver)   ! Entropy change sue to freezing of precip.
real(kind_phys) dmpdz2d(ncol,pver)     ! variable detrainment rate

real(kind_phys) mp(ncol)    ! Parcel mass flux.
real(kind_phys) qtp(ncol)   ! Parcel total water.
real(kind_phys) sp(ncol)    ! Parcel entropy.

real(kind_phys) sp0(ncol)    ! Parcel launch entropy.
real(kind_phys) qtp0(ncol)   ! Parcel launch total water.
real(kind_phys) mp0(ncol)    ! Parcel launch relative mass flux.

real(kind_phys) lwmax      ! Maximum condesate that can be held in cloud before rainout.
real(kind_phys) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
real(kind_phys) dmpdz      ! Parcel fractional mass entrainment rate (/m)
real(kind_phys) dpdz,dzdp  ! Hydrstatic relation and inverse of.
real(kind_phys) senv       ! Environmental entropy at each grid point.
real(kind_phys) qtenv      ! Environmental total water "   "   ".
real(kind_phys) penv       ! Environmental total pressure "   "   ".
real(kind_phys) tenv       ! Environmental total temperature "   "   ".
real(kind_phys) new_s      ! Hold value for entropy after condensation/freezing adjustments.
real(kind_phys) new_q      ! Hold value for total water after condensation/freezing adjustments.
real(kind_phys) dp         ! Layer thickness (center to center)
real(kind_phys) tfguess    ! First guess for entropy inversion - crucial for efficiency!
real(kind_phys) tscool     ! Super cooled temperature offset (in degC) (eg -35).

real(kind_phys) qxsk, qxskp1        ! LCL excess water (k, k+1)
real(kind_phys) dsdp, dqtdp, dqxsdp ! LCL s, qt, p gradients (k, k+1)
real(kind_phys) slcl,qtlcl,qslcl    ! LCL s, qt, qs values.
real(kind_phys) org2rkm, org2Tpert
real(kind_phys) dmpdz_lnd, dmpdz_mask

integer rcall       ! Number of ientropy call for errors recording
integer nit_lheat     ! Number of iterations for condensation/freezing loop.
integer i,k,ii   ! Loop counters.

!======================================================================
!    SUMMARY
!
!  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
!           and entrains at each level with a specified entrainment rate.
!
! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.
!
!======================================================================
!
! Set some values that may be changed frequently.
!

if (zm_org) then
   org2rkm = 10._kind_phys
   org2Tpert = 0._kind_phys
endif
nit_lheat = 2 ! iterations for ds,dq changes from condensation freezing.
dmpdz=dmpdz_param       ! Entrainment rate. (-ve for /m)
dmpdz_lnd=-1.e-3_kind_phys
lwmax = 1.e-3_kind_phys    ! Need to put formula in for this.
tscool = 0.0_kind_phys   ! Temp at which water loading freezes in the cloud.

qtmix=0._kind_phys
smix=0._kind_phys

qtenv = 0._kind_phys
senv = 0._kind_phys
tenv = 0._kind_phys
penv = 0._kind_phys

qtp0 = 0._kind_phys
sp0  = 0._kind_phys
mp0 = 0._kind_phys

qtp = 0._kind_phys
sp = 0._kind_phys
mp = 0._kind_phys

new_q = 0._kind_phys
new_s = 0._kind_phys

! **** Begin loops ****

do k = pver, msg+1, -1
   do i=1,ncol

! Initialize parcel values at launch level.

      if (k == klaunch(i)) then

         if (lparcel_pbl) then ! Modifcations to parcel properties if lparcel_pbl set.

            qtp0(i) = ql(i)     ! Parcel launch q (PBL mixed value).
            sp0(i)  = entropy(tl(i),pl(i),qtp0(i),cpliq,cpwv,rh2o) ! Parcel launch entropy could be a mixed parcel.

         else

            qtp0(i) = q(i,k)    ! Parcel launch total water (assuming subsaturated)
            sp0(i)  = entropy(t(i,k),p(i,k),qtp0(i),cpliq,cpwv,rh2o) ! Parcel launch entropy.

         end if

         mp0(i)  = 1._kind_phys       ! Parcel launch relative mass (i.e. 1 parcel stays 1 parcel for dmpdp=0, undilute).
         smix(i,k)  = sp0(i)
         qtmix(i,k) = qtp0(i)
         tfguess = t(i,k)
         rcall = 1
         call ientropy (rcall,i,smix(i,k),p(i,k),qtmix(i,k),tmix(i,k),qsmix(i,k),tfguess,cpliq,cpwv,rh2o,errmsg,errflg)
      end if

! Entraining levels

      if (k < klaunch(i)) then

! Set environmental values for this level.

         dp = (p(i,k)-p(i,k+1)) ! In -ve mb as p decreasing with height - difference between center of layers.
         qtenv = 0.5_kind_phys*(q(i,k)+q(i,k+1))         ! Total water of environment.
         tenv  = 0.5_kind_phys*(t(i,k)+t(i,k+1))
         penv  = 0.5_kind_phys*(p(i,k)+p(i,k+1))

         senv  = entropy(tenv,penv,qtenv,cpliq,cpwv,rh2o)  ! Entropy of environment.

! Determine fractional entrainment rate /pa given value /m.

         dpdz = -(penv*grav)/(rgas*tenv) ! in mb/m since  p in mb.
         dzdp = 1._kind_phys/dpdz                  ! in m/mb
         if (zm_org) then
            dmpdz_mask = landfrac(i) * dmpdz_lnd + (1._kind_phys - landfrac(i)) * dmpdz
            dmpdp = (dmpdz_mask/(1._kind_phys+org(i,k)*org2rkm))*dzdp              ! /mb Fractional entrainment
         else
            dmpdp = dmpdz*dzdp
         endif

! Sum entrainment to current level
! entrains q,s out of intervening dp layers, in which linear variation is assumed
! so really it entrains the mean of the 2 stored values.

         sp(i)  = sp(i)  - dmpdp*dp*senv
         qtp(i) = qtp(i) - dmpdp*dp*qtenv
         mp(i)  = mp(i)  - dmpdp*dp

! Entrain s and qt to next level.

         smix(i,k)  = (sp0(i)  +  sp(i)) / (mp0(i) + mp(i))
         qtmix(i,k) = (qtp0(i) + qtp(i)) / (mp0(i) + mp(i))

! Invert entropy from s and q to determine T and saturation-capped q of mixture.
! t(i,k) used as a first guess so that it converges faster.

         tfguess = tmix(i,k+1)
         rcall = 2
         call ientropy(rcall,i,smix(i,k),p(i,k),qtmix(i,k),tmix(i,k),qsmix(i,k),tfguess,cpliq,cpwv,rh2o,errmsg,errflg)

!
! Determine if this is lcl of this column if qsmix <= qtmix.
! FIRST LEVEL where this happens on ascending.

         if (qsmix(i,k) <= qtmix(i,k) .and. qsmix(i,k+1) > qtmix(i,k+1)) then
            lcl(i) = k
            qxsk   = qtmix(i,k) - qsmix(i,k)
            qxskp1 = qtmix(i,k+1) - qsmix(i,k+1)
            dqxsdp = (qxsk - qxskp1)/dp
            pl(i)  = p(i,k+1) - qxskp1/dqxsdp    ! pressure level of actual lcl.
            dsdp   = (smix(i,k)  - smix(i,k+1))/dp
            dqtdp  = (qtmix(i,k) - qtmix(i,k+1))/dp
            slcl   = smix(i,k+1)  +  dsdp* (pl(i)-p(i,k+1))
            qtlcl  = qtmix(i,k+1) +  dqtdp*(pl(i)-p(i,k+1))

            tfguess = tmix(i,k)
            rcall = 3
            call ientropy (rcall,i,slcl,pl(i),qtlcl,tl(i),qslcl,tfguess,cpliq,cpwv,rh2o,errmsg,errflg)

!            write(iulog,*)' '
!            write(iulog,*)' p',p(i,k+1),pl(i),p(i,lcl(i))
!            write(iulog,*)' t',tmix(i,k+1),tl(i),tmix(i,lcl(i))
!            write(iulog,*)' s',smix(i,k+1),slcl,smix(i,lcl(i))
!            write(iulog,*)'qt',qtmix(i,k+1),qtlcl,qtmix(i,lcl(i))
!            write(iulog,*)'qs',qsmix(i,k+1),qslcl,qsmix(i,lcl(i))

         endif
!
      end if !  k < klaunch


   end do ! Levels loop
end do ! Columns loop

!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Could stop now and test with this as it will provide some estimate of buoyancy
!! without the effects of freezing/condensation taken into account for tmix.

!! So we now have a profile of entropy and total water of the entraining parcel
!! Varying with height from the launch level klaunch parcel=environment. To the
!! top allowed level for the existence of convection.

!! Now we have to adjust these values such that the water held in vaopor is < or
!! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
!! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously
!! provides latent heating to the mixed parcel and so this has to be added back
!! to it. But does this also increase qsmix as well? Also freezing processes


xsh2o = 0._kind_phys
ds_xsh2o = 0._kind_phys
ds_freeze = 0._kind_phys

!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate solution twice for accuracy



do k = pver, msg+1, -1
   do i=1,ncol

! Initialize variables at k=klaunch

      if (k == klaunch(i)) then

! Set parcel values at launch level assume no liquid water.

         tp(i,k)    = tmix(i,k)
         qstp(i,k)  = q(i,k)
         if (zm_org) then
            tpv(i,k)   =  (tp(i,k) + (org2Tpert*org(i,k)+tpert(i))) * (1._kind_phys+1.608_kind_phys*qstp(i,k)) / (1._kind_phys+qstp(i,k))
         else
            tpv(i,k)   =  (tp(i,k) + tpert(i)) * (1._kind_phys+1.608_kind_phys*qstp(i,k)) / (1._kind_phys+qstp(i,k))
         endif

      end if

      if (k < klaunch(i)) then

! Initiaite loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.

! Iterate nit_lheat times for s,qt changes.

         do ii=0,nit_lheat-1

! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).

            xsh2o(i,k) = max (0._kind_phys, qtmix(i,k) - qsmix(i,k) - lwmax)

! Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)

            ds_xsh2o(i,k) = ds_xsh2o(i,k+1) - cpliq * log (tmix(i,k)/tfreez) * max(0._kind_phys,(xsh2o(i,k)-xsh2o(i,k+1)))
!
! Entropy of freezing: latice times amount of water involved divided by T.
!

            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) == 0._kind_phys) then ! One off freezing of condensate.
               ds_freeze(i,k) = (latice/tmix(i,k)) * max(0._kind_phys,qtmix(i,k)-qsmix(i,k)-xsh2o(i,k)) ! Gain of LH
            end if

            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) /= 0._kind_phys) then ! Continual freezing of additional condensate.
               ds_freeze(i,k) = ds_freeze(i,k+1)+(latice/tmix(i,k)) * max(0._kind_phys,(qsmix(i,k+1)-qsmix(i,k)))
            end if

! Adjust entropy and accordingly to sum of ds (be careful of signs).

            new_s = smix(i,k) + ds_xsh2o(i,k) + ds_freeze(i,k)

! Adjust liquid water and accordingly to xsh2o.

            new_q = qtmix(i,k) - xsh2o(i,k)

! Invert entropy to get updated Tmix and qsmix of parcel.

            tfguess = tmix(i,k)
            rcall =4
            call ientropy (rcall,i,new_s, p(i,k), new_q, tmix(i,k), qsmix(i,k), tfguess,cpliq,cpwv,rh2o,errmsg,errflg)

         end do  ! Iteration loop for freezing processes.

! tp  - Parcel temp is temp of mixture.
! tpv - Parcel v. temp should be density temp with new_q total water.

         tp(i,k)    = tmix(i,k)

! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)

         if (new_q > qsmix(i,k)) then  ! Super-saturated so condensate present - reduces buoyancy.
            qstp(i,k) = qsmix(i,k)
         else                          ! Just saturated/sub-saturated - no condensate virtual effects.
            qstp(i,k) = new_q
         end if

         if (zm_org) then
            tpv(i,k) = (tp(i,k)+(org2Tpert*org(i,k)+tpert(i)))* (1._kind_phys+1.608_kind_phys*qstp(i,k)) / (1._kind_phys+ new_q)
         else
            tpv(i,k) = (tp(i,k)+tpert(i))* (1._kind_phys+1.608_kind_phys*qstp(i,k)) / (1._kind_phys+ new_q)
         endif

      end if ! k < klaunch

   end do ! Loop for columns

end do  ! Loop for vertical levels.


return
end subroutine parcel_dilute

!-----------------------------------------------------------------------------------------
real(kind_phys) function entropy(TK,p,qtot,cpliq,cpwv,rh2o)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
! from Raymond and Blyth 1992
!
     real(kind_phys), intent(in) :: p,qtot,TK
     real(kind_phys), intent(in) :: cpliq
     real(kind_phys), intent(in) :: cpwv
     real(kind_phys), intent(in) :: rh2o

     real(kind_phys) :: qv,qst,e,est,L
     real(kind_phys), parameter :: pref = 1000._kind_phys

L = rl - (cpliq - cpwv)*(TK-tfreez)         ! T IN CENTIGRADE

call qsat_hPa(TK, p, est, qst)

qv = min(qtot,qst)                         ! Partition qtot into vapor part only.
e = qv*p / (eps1 +qv)

entropy = (cpres + qtot*cpliq)*log( TK/tfreez) - rgas*log( (p-e)/pref ) + &
        L*qv/TK - qv*rh2o*log(qv/qst)

end FUNCTION entropy

!
!-----------------------------------------------------------------------------------------
SUBROUTINE ientropy (rcall,icol,s,p,qt,T,qst,Tfg,cpliq,cpwv,rh2o,errmsg,errflg)
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg).
! Inverts entropy, pressure and total water qt
! for T and saturated vapor mixing ratio
!

! CACNOTE - Remove this when pass in lat/lon or pass out lchnk,icol (Note: lchnk will not exist in CAM-SIMA)
!  use phys_grid, only: get_rlon_p, get_rlat_p

  integer, intent(in) :: icol, rcall
  real(kind_phys), intent(in)  :: s, p, Tfg, qt
  real(kind_phys), intent(in) :: cpliq
  real(kind_phys), intent(in) :: cpwv
  real(kind_phys), intent(in) :: rh2o
  real(kind_phys), intent(out) :: qst, T
  character(len=512), intent(out)      :: errmsg
  integer, intent(out)                 :: errflg

  real(kind_phys) :: est, this_lat,this_lon
  real(kind_phys) :: a,b,c,d,ebr,fa,fb,fc,pbr,qbr,rbr,sbr,tol1,xm,tol
  integer :: i

  logical :: converged

  ! Max number of iteration loops.
  integer, parameter :: LOOPMAX = 100
  real(kind_phys), parameter :: EPS = 3.e-8_kind_phys

  converged = .false.

  ! Invert the entropy equation -- use Brent's method
  ! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.

  T = Tfg                  ! Better first guess based on Tprofile from conv.

  a = Tfg-10    !low bracket
  b = Tfg+10    !high bracket

  fa = entropy(a, p, qt,cpliq,cpwv,rh2o) - s
  fb = entropy(b, p, qt,cpliq,cpwv,rh2o) - s

  c=b
  fc=fb
  tol=0.001_kind_phys

  converge: do i=0, LOOPMAX
     if ((fb > 0.0_kind_phys .and. fc > 0.0_kind_phys) .or. &
          (fb < 0.0_kind_phys .and. fc < 0.0_kind_phys)) then
        c=a
        fc=fa
        d=b-a
        ebr=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if

     tol1=2.0_kind_phys*EPS*abs(b)+0.5_kind_phys*tol
     xm=0.5_kind_phys*(c-b)
     converged = (abs(xm) <= tol1 .or. fb == 0.0_kind_phys)
     if (converged) exit converge

     if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
        sbr=fb/fa
        if (a == c) then
           pbr=2.0_kind_phys*xm*sbr
           qbr=1.0_kind_phys-sbr
        else
           qbr=fa/fc
           rbr=fb/fc
           pbr=sbr*(2.0_kind_phys*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0_kind_phys))
           qbr=(qbr-1.0_kind_phys)*(rbr-1.0_kind_phys)*(sbr-1.0_kind_phys)
        end if
        if (pbr > 0.0_kind_phys) qbr=-qbr
        pbr=abs(pbr)
        if (2.0_kind_phys*pbr  <  min(3.0_kind_phys*xm*qbr-abs(tol1*qbr),abs(ebr*qbr))) then
           ebr=d
           d=pbr/qbr
        else
           d=xm
           ebr=d
        end if
     else
        d=xm
        ebr=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )

     fb = entropy(b, p, qt,cpliq,cpwv,rh2o) - s

  end do converge

  T = b
  call qsat_hPa(T, p, est, qst)

  if (.not. converged) then
!CACNOTE - Revisit this with Jesse
!     this_lat = get_rlat_p(lchnk, icol)*57.296_kind_phys
!     this_lon = get_rlon_p(lchnk, icol)*57.296_kind_phys
!     write(errmsg,100) 'ZM_CONV: IENTROPY. Details: call#,lchnk,icol= ',rcall,lchnk,icol, &
!          ' lat: ',this_lat,' lon: ',this_lon, &
!          ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qt(g/kg) = ', 1000._kind_phys*qt, &
!          ' qst(g/kg) = ', 1000._kind_phys*qst,', s(J/kg) = ',s
     errflg=1
  end if

100 format (A,I1,I4,I4,7(A,F6.2))

end SUBROUTINE ientropy

subroutine cldprp(ncol   ,pver    ,pverp   ,cpliq   , &
                  latice  ,cpwv    ,rh2o    ,&
                  q       ,t       ,u       ,v       ,p       , &
                  z       ,s       ,mu      ,eu      ,du      , &
                  md      ,ed      ,sd      ,qd      ,mc      , &
                  qu      ,su      ,zf      ,qst     ,hmn     , &
                  hsat    ,shat    ,ql      , &
                  cmeg    ,jb      ,lel     ,jt      ,jlcl    , &
                  mx      ,j0      ,jd      ,rl      ,il2g    , &
                  rd      ,grav    ,cp      ,msg     , &
                  evp     ,cu      ,rprd    ,limcnv  ,landfrac, &
                  qcde     ,qhat  )

!-----------------------------------------------------------------------
!CACNOTE - fill in documentation
!
! Purpose:
! <Say what the routine does>
!
! Method:
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             original version cldprop.
!
! Author: See above, modified by P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments disabled
! these are to facilitate future development. We expect to release a
! cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
!
!-----------------------------------------------------------------------

   implicit none

!------------------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol
   integer, intent(in) :: pver
   integer, intent(in) :: pverp

   real(kind_phys), intent(in) :: cpliq
   real(kind_phys), intent(in) :: latice
   real(kind_phys), intent(in) :: cpwv
   real(kind_phys), intent(in) :: rh2o

   real(kind_phys), intent(in) :: q(ncol,pver)         ! spec. humidity of env
   real(kind_phys), intent(in) :: t(ncol,pver)         ! temp of env
   real(kind_phys), intent(in) :: p(ncol,pver)         ! pressure of env
   real(kind_phys), intent(in) :: z(ncol,pver)         ! height of env
   real(kind_phys), intent(in) :: s(ncol,pver)         ! normalized dry static energy of env
   real(kind_phys), intent(in) :: zf(ncol,pverp)       ! height of interfaces
   real(kind_phys), intent(in) :: u(ncol,pver)         ! zonal velocity of env
   real(kind_phys), intent(in) :: v(ncol,pver)         ! merid. velocity of env

   real(kind_phys), intent(in) :: landfrac(ncol) ! RBN Landfrac

   integer, intent(in) :: jb(ncol)              ! updraft base level
   integer, intent(in) :: lel(ncol)             ! updraft launch level
   integer, intent(in) :: mx(ncol)              ! updraft base level (same is jb)
   integer, intent(out) :: jt(ncol)              ! updraft plume top
   integer, intent(out) :: jlcl(ncol)            ! updraft lifting cond level
   integer, intent(out) :: j0(ncol)              ! level where updraft begins detraining
   integer, intent(out) :: jd(ncol)              ! level of downdraft
   integer, intent(in) :: limcnv                 ! convection limiting level
   integer, intent(in) :: il2g                   !CORE GROUP REMOVE
   integer, intent(in) :: msg                    ! missing moisture vals (always 0)
   real(kind_phys), intent(in) :: rl                    ! latent heat of vap
   real(kind_phys), intent(in) :: shat(ncol,pver)      ! interface values of dry stat energy
   real(kind_phys), intent(in) :: qhat(ncol,pver)      ! wg grid slice of upper interface mixing ratio.

!
! output
!

   real(kind_phys), intent(out) :: rprd(ncol,pver)     ! rate of production of precip at that layer
   real(kind_phys), intent(out) :: du(ncol,pver)       ! detrainement rate of updraft
   real(kind_phys), intent(out) :: ed(ncol,pver)       ! entrainment rate of downdraft
   real(kind_phys), intent(out) :: eu(ncol,pver)       ! entrainment rate of updraft
   real(kind_phys), intent(out) :: hmn(ncol,pver)      ! moist stat energy of env
   real(kind_phys), intent(out) :: hsat(ncol,pver)     ! sat moist stat energy of env
   real(kind_phys), intent(out) :: mc(ncol,pver)       ! net mass flux
   real(kind_phys), intent(out) :: md(ncol,pver)       ! downdraft mass flux
   real(kind_phys), intent(out) :: mu(ncol,pver)       ! updraft mass flux
   real(kind_phys), intent(out) :: qd(ncol,pver)       ! spec humidity of downdraft
   real(kind_phys), intent(out) :: ql(ncol,pver)       ! liq water of updraft
   real(kind_phys), intent(out) :: qst(ncol,pver)      ! saturation mixing ratio of env.
   real(kind_phys), intent(out) :: qu(ncol,pver)       ! spec hum of updraft
   real(kind_phys), intent(out) :: sd(ncol,pver)       ! normalized dry stat energy of downdraft
   real(kind_phys), intent(out) :: su(ncol,pver)       ! normalized dry stat energy of updraft
   real(kind_phys), intent(out) :: qcde(ncol,pver)     ! cloud water mixing ratio for detrainment (kg/kg)

   real(kind_phys) rd                   ! gas constant for dry air
   real(kind_phys) grav                 ! gravity
   real(kind_phys) cp                   ! heat capacity of dry air

!
! Local workspace
!
   real(kind_phys) gamma(ncol,pver)
   real(kind_phys) dz(ncol,pver)
   real(kind_phys) iprm(ncol,pver)
   real(kind_phys) hu(ncol,pver)
   real(kind_phys) hd(ncol,pver)
   real(kind_phys) eps(ncol,pver)
   real(kind_phys) f(ncol,pver)
   real(kind_phys) k1(ncol,pver)
   real(kind_phys) i2(ncol,pver)
   real(kind_phys) ihat(ncol,pver)
   real(kind_phys) i3(ncol,pver)
   real(kind_phys) idag(ncol,pver)
   real(kind_phys) i4(ncol,pver)
   real(kind_phys) qsthat(ncol,pver)
   real(kind_phys) hsthat(ncol,pver)
   real(kind_phys) gamhat(ncol,pver)
   real(kind_phys) cu(ncol,pver)
   real(kind_phys) evp(ncol,pver)
   real(kind_phys) cmeg(ncol,pver)
   real(kind_phys) qds(ncol,pver)
   real(kind_phys) c0mask(ncol)

   real(kind_phys) hmin(ncol)
   real(kind_phys) expdif(ncol)
   real(kind_phys) expnum(ncol)
   real(kind_phys) ftemp(ncol)
   real(kind_phys) eps0(ncol)
   real(kind_phys) rmue(ncol)
   real(kind_phys) zuef(ncol)
   real(kind_phys) zdef(ncol)
   real(kind_phys) epsm(ncol)
   real(kind_phys) ratmjb(ncol)
   real(kind_phys) est(ncol)
   real(kind_phys) totpcp(ncol)
   real(kind_phys) totevp(ncol)
   real(kind_phys) alfa(ncol)
   real(kind_phys) ql1
   real(kind_phys) tu
   real(kind_phys) estu
   real(kind_phys) qstu

   real(kind_phys) small
   real(kind_phys) mdt

   real(kind_phys) fice(ncol,pver)        ! ice fraction in precip production
   real(kind_phys) tug(ncol,pver)

   real(kind_phys) tvuo(ncol,pver)        ! updraft virtual T w/o freezing heating
   real(kind_phys) tvu(ncol,pver)         ! updraft virtual T with freezing heating
   real(kind_phys) totfrz(ncol)
   real(kind_phys) frz (ncol,pver)        ! rate of freezing
   integer  jto(ncol)              ! updraft plume old top
   integer  tmplel(ncol)

   integer  iter, itnum
   integer  m

   integer khighest
   integer klowest
   integer kount
   integer i,k

   logical doit(ncol)
   logical done(ncol)
!
!------------------------------------------------------------------------------
!

   do i = 1,il2g
      ftemp(i) = 0._kind_phys
      expnum(i) = 0._kind_phys
      expdif(i) = 0._kind_phys
      c0mask(i)  = c0_ocn * (1._kind_phys-landfrac(i)) +   c0_lnd * landfrac(i)
   end do
!
!jr Change from msg+1 to 1 to prevent blowup
!
   do k = 1,pver
      do i = 1,il2g
         dz(i,k) = zf(i,k) - zf(i,k+1)
      end do
   end do

!
! initialize many output and work variables to zero
!
   do k = 1,pver
      do i = 1,il2g
         k1(i,k) = 0._kind_phys
         i2(i,k) = 0._kind_phys
         i3(i,k) = 0._kind_phys
         i4(i,k) = 0._kind_phys
         mu(i,k) = 0._kind_phys
         f(i,k) = 0._kind_phys
         eps(i,k) = 0._kind_phys
         eu(i,k) = 0._kind_phys
         du(i,k) = 0._kind_phys
         ql(i,k) = 0._kind_phys
         cu(i,k) = 0._kind_phys
         evp(i,k) = 0._kind_phys
         cmeg(i,k) = 0._kind_phys
         qds(i,k) = q(i,k)
         md(i,k) = 0._kind_phys
         ed(i,k) = 0._kind_phys
         sd(i,k) = s(i,k)
         qd(i,k) = q(i,k)
         mc(i,k) = 0._kind_phys
         qu(i,k) = q(i,k)
         su(i,k) = s(i,k)
         call qsat_hPa(t(i,k), p(i,k), est(i), qst(i,k))

         if ( p(i,k)-est(i) <= 0._kind_phys ) then
            qst(i,k) = 1.0_kind_phys
         end if

         gamma(i,k) = qst(i,k)*(1._kind_phys + qst(i,k)/eps1)*eps1*rl/(rd*t(i,k)**2)*rl/cp
         hmn(i,k) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
         hsat(i,k) = cp*t(i,k) + grav*z(i,k) + rl*qst(i,k)
         hu(i,k) = hmn(i,k)
         hd(i,k) = hmn(i,k)
         rprd(i,k) = 0._kind_phys

         fice(i,k) = 0._kind_phys
         tug(i,k)  = 0._kind_phys
         qcde(i,k)   = 0._kind_phys
         tvuo(i,k) = (shat(i,k) - grav/cp*zf(i,k))*(1._kind_phys + 0.608_kind_phys*qhat(i,k))
         tvu(i,k) = tvuo(i,k)
         frz(i,k)  = 0._kind_phys

      end do
   end do

!
!jr Set to zero things which make this routine blow up
!
   do k=1,msg
      do i=1,il2g
         rprd(i,k) = 0._kind_phys
      end do
   end do
!
! interpolate the layer values of qst, hsat and gamma to
! layer interfaces
!
   do k = 1, msg+1
      do i = 1,il2g
         hsthat(i,k) = hsat(i,k)
         qsthat(i,k) = qst(i,k)
         gamhat(i,k) = gamma(i,k)
      end do
   end do
   do i = 1,il2g
      totpcp(i) = 0._kind_phys
      totevp(i) = 0._kind_phys
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (abs(qst(i,k-1)-qst(i,k)) > 1.E-6_kind_phys) then
            qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*qst(i,k)/ (qst(i,k-1)-qst(i,k))
         else
            qsthat(i,k) = qst(i,k)
         end if
         hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
         if (abs(gamma(i,k-1)-gamma(i,k)) > 1.E-6_kind_phys) then
            gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*gamma(i,k-1)*gamma(i,k)/ &
                                (gamma(i,k-1)-gamma(i,k))
         else
            gamhat(i,k) = gamma(i,k)
         end if
      end do
   end do
!
! initialize cloud top to highest plume top.
!jr changed hard-wired 4 to limcnv+1 (not to exceed pver)
!
   jt(:) = pver
   do i = 1,il2g
      jt(i) = max(lel(i),limcnv+1)
      jt(i) = min(jt(i),pver)
      jd(i) = pver
      jlcl(i) = lel(i)
      hmin(i) = 1.E6_kind_phys
   end do
!
! find the level of minimum hsat, where detrainment starts
!

   do k = msg + 1,pver
      do i = 1,il2g
         if (hsat(i,k) <= hmin(i) .and. k >= jt(i) .and. k <= jb(i)) then
            hmin(i) = hsat(i,k)
            j0(i) = k
         end if
      end do
   end do
   do i = 1,il2g
      j0(i) = min(j0(i),jb(i)-2)
      j0(i) = max(j0(i),jt(i)+2)
!
! Fix from Guang Zhang to address out of bounds array reference
!
      j0(i) = min(j0(i),pver)
   end do
!
! Initialize certain arrays inside cloud
!
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= jb(i)) then
            hu(i,k) = hmn(i,mx(i)) + cp*tiedke_add
            su(i,k) = s(i,mx(i)) + tiedke_add
         end if
      end do
   end do
!
! *********************************************************
! compute taylor series for approximate eps(z) below
! *********************************************************
!
   do k = pver - 1,msg + 1,-1
      do i = 1,il2g
         if (k < jb(i) .and. k >= jt(i)) then
            k1(i,k) = k1(i,k+1) + (hmn(i,mx(i))-hmn(i,k))*dz(i,k)
            ihat(i,k) = 0.5_kind_phys* (k1(i,k+1)+k1(i,k))
            i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
            idag(i,k) = 0.5_kind_phys* (i2(i,k+1)+i2(i,k))
            i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
            iprm(i,k) = 0.5_kind_phys* (i3(i,k+1)+i3(i,k))
            i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
         end if
      end do
   end do
!
! re-initialize hmin array for ensuing calculation.
!
   do i = 1,il2g
      hmin(i) = 1.E6_kind_phys
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i) .and. hmn(i,k) <= hmin(i)) then
            hmin(i) = hmn(i,k)
            expdif(i) = hmn(i,mx(i)) - hmin(i)
         end if
      end do
   end do
!
! *********************************************************
! compute approximate eps(z) using above taylor series
! *********************************************************
!
   do k = msg + 2,pver
      do i = 1,il2g
         expnum(i) = 0._kind_phys
         ftemp(i) = 0._kind_phys
         if (k < jt(i) .or. k >= jb(i)) then
            k1(i,k) = 0._kind_phys
            expnum(i) = 0._kind_phys
         else
            expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1)*(zf(i,k)-z(i,k)) + &
                        hsat(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))
         end if
         if ((expdif(i) > 100._kind_phys .and. expnum(i) > 0._kind_phys) .and. &
            k1(i,k) > expnum(i)*dz(i,k)) then
            ftemp(i) = expnum(i)/k1(i,k)
            f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 + &
                     (2._kind_phys*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2* &
                     ftemp(i)**3 + (-5._kind_phys*k1(i,k)*i2(i,k)*i3(i,k)+ &
                     5._kind_phys*i2(i,k)**3+k1(i,k)**2*i4(i,k))/ &
                     k1(i,k)**3*ftemp(i)**4
            f(i,k) = max(f(i,k),0._kind_phys)
            f(i,k) = min(f(i,k),0.0002_kind_phys)
         end if
      end do
   end do
   do i = 1,il2g
      if (j0(i) < jb(i)) then
         if (f(i,j0(i)) < 1.E-6_kind_phys .and. f(i,j0(i)+1) > f(i,j0(i))) j0(i) = j0(i) + 1
      end if
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= j0(i)) then
            f(i,k) = max(f(i,k),f(i,k-1))
         end if
      end do
   end do
   do i = 1,il2g
      eps0(i) = f(i,j0(i))
      eps(i,jb(i)) = eps0(i)
   end do
!
! This is set to match the Rasch and Kristjansson paper
!
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i)) then
            eps(i,k) = f(i,j0(i))
         end if
      end do
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k < j0(i) .and. k >= jt(i)) eps(i,k) = f(i,k)
      end do
   end do

   itnum = 1
   do iter=1, itnum

!
! specify the updraft mass flux mu, entrainment eu, detrainment du
! and moist static energy hu.
! here and below mu, eu,du, md and ed are all normalized by mb
!
      do i = 1,il2g
         if (eps0(i) > 0._kind_phys) then
            mu(i,jb(i)) = 1._kind_phys
            eu(i,jb(i)) = mu(i,jb(i))/dz(i,jb(i))
         end if

         tmplel(i) = jt(i)
      end do
      do k = pver,msg + 1,-1
         do i = 1,il2g
            if (eps0(i) > 0._kind_phys .and. (k >= tmplel(i) .and. k < jb(i))) then
               zuef(i) = zf(i,k) - zf(i,jb(i))
               rmue(i) = (1._kind_phys/eps0(i))* (exp(eps(i,k+1)*zuef(i))-1._kind_phys)/zuef(i)
               mu(i,k) = (1._kind_phys/eps0(i))* (exp(eps(i,k  )*zuef(i))-1._kind_phys)/zuef(i)
               eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
               du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
            end if
         end do
      end do

      khighest = pverp
      klowest = 1
      do i=1,il2g
         khighest = min(khighest,lel(i))
         klowest = max(klowest,jb(i))
      end do
      do k = klowest-1,khighest,-1
         do i = 1,il2g
            if (k <= jb(i)-1 .and. k >= lel(i) .and. eps0(i) > 0._kind_phys) then
               if (mu(i,k) < 0.02_kind_phys) then
                  hu(i,k) = hmn(i,k)
                  mu(i,k) = 0._kind_phys
                  eu(i,k) = 0._kind_phys
                  du(i,k) = mu(i,k+1)/dz(i,k)
               else
                  hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) + &
                            dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k)- du(i,k)*hsat(i,k))
               end if
            end if
         end do
      end do
!
! reset cloud top index beginning from two layers above the
! cloud base (i.e. if cloud is only one layer thick, top is not reset
!
      do i=1,il2g
         doit(i) = .true.
         totfrz(i)= 0._kind_phys
         do k = pver,msg + 1,-1
            totfrz(i)= totfrz(i)+ frz(i,k)*dz(i,k)
         end do
      end do
      do k=klowest-2,khighest-1,-1
         do i=1,il2g
            if (doit(i) .and. k <= jb(i)-2 .and. k >= lel(i)-1) then
               if (hu(i,k) <= hsthat(i,k) .and. hu(i,k+1) > hsthat(i,k+1) &
                  .and. mu(i,k) >= 0.02_kind_phys) then
                  if (hu(i,k)-hsthat(i,k) < -2000._kind_phys) then
                     jt(i) = k + 1
                     doit(i) = .false.
                  else
                     jt(i) = k
                     doit(i) = .false.
                  end if
               else if ( (hu(i,k) > hu(i,jb(i)) .and. totfrz(i)<=0._kind_phys) .or. mu(i,k) < 0.02_kind_phys) then
                  jt(i) = k + 1
                  doit(i) = .false.
               end if
            end if
         end do
      end do

      if (iter == 1)  jto(:) = jt(:)

      do k = pver,msg + 1,-1
         do i = 1,il2g
            if (k >= lel(i) .and. k <= jt(i) .and. eps0(i) > 0._kind_phys) then
               mu(i,k) = 0._kind_phys
               eu(i,k) = 0._kind_phys
               du(i,k) = 0._kind_phys
               hu(i,k) = hmn(i,k)
            end if
            if (k == jt(i) .and. eps0(i) > 0._kind_phys) then
               du(i,k) = mu(i,k+1)/dz(i,k)
               eu(i,k) = 0._kind_phys
               mu(i,k) = 0._kind_phys
            end if
         end do
      end do

      do i = 1,il2g
         done(i) = .false.
      end do
      kount = 0
      do k = pver,msg + 2,-1
         do i = 1,il2g
            if (k == jb(i) .and. eps0(i) > 0._kind_phys) then
               qu(i,k) = q(i,mx(i))
               su(i,k) = (hu(i,k)-rl*qu(i,k))/cp
            end if
            if (( .not. done(i) .and. k > jt(i) .and. k < jb(i)) .and. eps0(i) > 0._kind_phys) then
               su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) + &
                         dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)
               qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)- &
                               du(i,k)*qst(i,k))
               tu = su(i,k) - grav/cp*zf(i,k)
               call qsat_hPa(tu, (p(i,k)+p(i,k-1))/2._kind_phys, estu, qstu)
               if (qu(i,k) >= qstu) then
                  jlcl(i) = k
                  kount = kount + 1
                  done(i) = .true.
               end if
            end if
         end do
         if (kount >= il2g) goto 690
      end do
690   continue
      do k = msg + 2,pver
         do i = 1,il2g
            if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0._kind_phys) then
               su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/(cp* (1._kind_phys+gamhat(i,k)))
               qu(i,k) = qsthat(i,k) + gamhat(i,k)*(hu(i,k)-hsthat(i,k))/ &
                        (rl* (1._kind_phys+gamhat(i,k)))
            end if
         end do
      end do

! compute condensation in updraft
      tmplel(:il2g) = jb(:il2g)

      do k = pver,msg + 2,-1
         do i = 1,il2g
             if (k >= jt(i) .and. k < tmplel(i) .and. eps0(i) > 0._kind_phys) then

               cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/ &
                         dz(i,k)- (eu(i,k)-du(i,k))*s(i,k))/(rl/cp)
               if (k == jt(i)) cu(i,k) = 0._kind_phys
               cu(i,k) = max(0._kind_phys,cu(i,k))
            end if
         end do
      end do


! compute condensed liquid, rain production rate
! accumulate total precipitation (condensation - detrainment of liquid)
! Note ql1 = ql(k) + rprd(k)*dz(k)/mu(k)
! The differencing is somewhat strange (e.g. du(i,k)*ql(i,k+1)) but is
! consistently applied.
!    mu, ql are interface quantities
!    cu, du, eu, rprd are midpoint quantites

         do k = pver,msg + 2,-1
            do i = 1,il2g
               rprd(i,k) = 0._kind_phys
               if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0._kind_phys .and. mu(i,k) >= 0.0_kind_phys) then
                  if (mu(i,k) > 0._kind_phys) then
                     ql1 = 1._kind_phys/mu(i,k)* (mu(i,k+1)*ql(i,k+1)- &
                           dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))
                     ql(i,k) = ql1/ (1._kind_phys+dz(i,k)*c0mask(i))
                  else
                     ql(i,k) = 0._kind_phys
                  end if
                  totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1))
                  rprd(i,k) = c0mask(i)*mu(i,k)*ql(i,k)
                  qcde(i,k) = ql(i,k)


               end if
            end do
         end do
!
   end do   !iter
!
! specify downdraft properties (no downdrafts if jd.ge.jb).
! scale down downward mass flux profile so that net flux
! (up-down) at cloud base in not negative.
!
   do i = 1,il2g
!
! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
!
      alfa(i) = 0.1_kind_phys
      jt(i) = min(jt(i),jb(i)-1)
      jd(i) = max(j0(i),jt(i)+1)
      jd(i) = min(jd(i),jb(i))
      hd(i,jd(i)) = hmn(i,jd(i)-1)
      if (jd(i) < jb(i) .and. eps0(i) > 0._kind_phys) then
         epsm(i) = eps0(i)
         md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
      end if
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k > jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._kind_phys) then
            zdef(i) = zf(i,jd(i)) - zf(i,k)
            md(i,k) = -alfa(i)/ (2._kind_phys*eps0(i))*(exp(2._kind_phys*epsm(i)*zdef(i))-1._kind_phys)/zdef(i)
         end if
      end do
   end do

   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= jb(i)) .and. eps0(i) > 0._kind_phys .and. jd(i) < jb(i)) then
            ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1._kind_phys)
            md(i,k) = md(i,k)*ratmjb(i)
         end if
      end do
   end do

   small = 1.e-20_kind_phys
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= pver) .and. eps0(i) > 0._kind_phys) then
            ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
            mdt = min(md(i,k),-small)
            hd(i,k) = (md(i,k-1)*hd(i,k-1) - dz(i,k-1)*ed(i,k-1)*hmn(i,k-1))/mdt
         end if
      end do
   end do
!
! calculate updraft and downdraft properties.
!
   do k = msg + 2,pver
      do i = 1,il2g
         if ((k >= jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._kind_phys .and. jd(i) < jb(i)) then
            qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/ &
               (rl*(1._kind_phys + gamhat(i,k)))
         end if
      end do
   end do

   do i = 1,il2g
      qd(i,jd(i)) = qds(i,jd(i))
      sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cp
   end do
!
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jd(i) .and. k < jb(i) .and. eps0(i) > 0._kind_phys) then
            qd(i,k+1) = qds(i,k+1)
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._kind_phys)
            mdt = min(md(i,k+1),-small)
            sd(i,k+1) = ((rl/cp*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) + md(i,k)*sd(i,k))/mdt
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   end do
   do i = 1,il2g
      totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) - md(i,jb(i))*qd(i,jb(i))
   end do
!!$   if (.true.) then
   if (.false.) then
      do i = 1,il2g
         k = jb(i)
         if (eps0(i) > 0._kind_phys) then
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._kind_phys)
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   endif

   do i = 1,il2g
      totpcp(i) = max(totpcp(i),0._kind_phys)
      totevp(i) = max(totevp(i),0._kind_phys)
   end do
!
   do k = msg + 2,pver
      do i = 1,il2g
         if (totevp(i) > 0._kind_phys .and. totpcp(i) > 0._kind_phys) then
            md(i,k)  = md (i,k)*min(1._kind_phys, totpcp(i)/(totevp(i)+totpcp(i)))
            ed(i,k)  = ed (i,k)*min(1._kind_phys, totpcp(i)/(totevp(i)+totpcp(i)))
            evp(i,k) = evp(i,k)*min(1._kind_phys, totpcp(i)/(totevp(i)+totpcp(i)))
         else
            md(i,k) = 0._kind_phys
            ed(i,k) = 0._kind_phys
            evp(i,k) = 0._kind_phys
         end if
! cmeg is the cloud water condensed - rain water evaporated
! rprd is the cloud water converted to rain - (rain evaporated)
         cmeg(i,k) = cu(i,k) - evp(i,k)
         rprd(i,k) = rprd(i,k)-evp(i,k)
      end do
   end do

!
   do k = msg + 1,pver
      do i = 1,il2g
         mc(i,k) = mu(i,k) + md(i,k)
      end do
   end do
!
   return
end subroutine cldprp

subroutine closure(ncol   ,pver, &
                   q       ,t       ,p       ,z       ,s       , &
                   tp      ,qs      ,qu      ,su      ,mc      , &
                   du      ,mu      ,md      ,qd      ,sd      , &
                   qhat    ,shat    ,dp      ,qstp    ,zf      , &
                   ql      ,dsubcld ,mb      ,cape    ,tl      , &
                   lcl     ,lel     ,jt      ,mx      ,il1g    , &
                   il2g    ,rd      ,grav    ,cp      ,rl      , &
                   msg     ,capelmt )
!-----------------------------------------------------------------------
!CACNOTE - fill in documentation
!
! Purpose:
! <Say what the routine does>
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: G. Zhang and collaborators. CCM contact:P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
!
!-----------------------------------------------------------------------

!
!-----------------------------Arguments---------------------------------
!
   integer, intent(in) :: ncol
   integer, intent(in) :: pver

   real(kind_phys), intent(inout) :: q(ncol,pver)        ! spec humidity
   real(kind_phys), intent(inout) :: t(ncol,pver)        ! temperature
   real(kind_phys), intent(inout) :: p(ncol,pver)        ! pressure (mb)
   real(kind_phys), intent(inout) :: mb(ncol)            ! cloud base mass flux
   real(kind_phys), intent(in) :: z(ncol,pver)        ! height (m)
   real(kind_phys), intent(in) :: s(ncol,pver)        ! normalized dry static energy
   real(kind_phys), intent(in) :: tp(ncol,pver)       ! parcel temp
   real(kind_phys), intent(in) :: qs(ncol,pver)       ! sat spec humidity
   real(kind_phys), intent(in) :: qu(ncol,pver)       ! updraft spec. humidity
   real(kind_phys), intent(in) :: su(ncol,pver)       ! normalized dry stat energy of updraft
   real(kind_phys), intent(in) :: mc(ncol,pver)       ! net convective mass flux
   real(kind_phys), intent(in) :: du(ncol,pver)       ! detrainment from updraft
   real(kind_phys), intent(in) :: mu(ncol,pver)       ! mass flux of updraft
   real(kind_phys), intent(in) :: md(ncol,pver)       ! mass flux of downdraft
   real(kind_phys), intent(in) :: qd(ncol,pver)       ! spec. humidity of downdraft
   real(kind_phys), intent(in) :: sd(ncol,pver)       ! dry static energy of downdraft
   real(kind_phys), intent(in) :: qhat(ncol,pver)     ! environment spec humidity at interfaces
   real(kind_phys), intent(in) :: shat(ncol,pver)     ! env. normalized dry static energy at intrfcs
   real(kind_phys), intent(in) :: dp(ncol,pver)       ! pressure thickness of layers
   real(kind_phys), intent(in) :: qstp(ncol,pver)     ! spec humidity of parcel
   real(kind_phys), intent(in) :: zf(ncol,pver+1)     ! height of interface levels
   real(kind_phys), intent(in) :: ql(ncol,pver)       ! liquid water mixing ratio

   real(kind_phys), intent(in) :: cape(ncol)          ! available pot. energy of column
   real(kind_phys), intent(in) :: tl(ncol)
   real(kind_phys), intent(in) :: dsubcld(ncol)       ! thickness of subcloud layer

   integer, intent(in) :: lcl(ncol)        ! index of lcl
   integer, intent(in) :: lel(ncol)        ! index of launch leve
   integer, intent(in) :: jt(ncol)         ! top of updraft
   integer, intent(in) :: mx(ncol)         ! base of updraft
!
!--------------------------Local variables------------------------------
!
   real(kind_phys) dtpdt(ncol,pver)
   real(kind_phys) dqsdtp(ncol,pver)
   real(kind_phys) dtmdt(ncol,pver)
   real(kind_phys) dqmdt(ncol,pver)
   real(kind_phys) dboydt(ncol,pver)
   real(kind_phys) thetavp(ncol,pver)
   real(kind_phys) thetavm(ncol,pver)

   real(kind_phys) dtbdt(ncol),dqbdt(ncol),dtldt(ncol)
   real(kind_phys) beta
   real(kind_phys) capelmt
   real(kind_phys) cp
   real(kind_phys) dadt(ncol)
   real(kind_phys) debdt
   real(kind_phys) dltaa
   real(kind_phys) eb
   real(kind_phys) grav

   integer i
   integer il1g
   integer il2g
   integer k, kmin, kmax
   integer msg

   real(kind_phys) rd
   real(kind_phys) rl
! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.
!
   do i = il1g,il2g
      mb(i) = 0._kind_phys
      eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
      dtbdt(i) = (1._kind_phys/dsubcld(i))* (mu(i,mx(i))*(shat(i,mx(i))-su(i,mx(i)))+ &
                  md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
      dqbdt(i) = (1._kind_phys/dsubcld(i))* (mu(i,mx(i))*(qhat(i,mx(i))-qu(i,mx(i)))+ &
                 md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
      debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
      dtldt(i) = -2840._kind_phys* (3.5_kind_phys/t(i,mx(i))*dtbdt(i)-debdt/eb)/ &
                 (3.5_kind_phys*log(t(i,mx(i)))-log(eb)-4.805_kind_phys)**2
   end do
!
!   dtmdt and dqmdt are cumulus heating and drying.
!
   do k = msg + 1,pver
      do i = il1g,il2g
         dtmdt(i,k) = 0._kind_phys
         dqmdt(i,k) = 0._kind_phys
      end do
   end do
!
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k == jt(i)) then
            dtmdt(i,k) = (1._kind_phys/dp(i,k))*(mu(i,k+1)* (su(i,k+1)-shat(i,k+1)- &
                          rl/cp*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-shat(i,k+1)))
            dqmdt(i,k) = (1._kind_phys/dp(i,k))*(mu(i,k+1)* (qu(i,k+1)- &
                         qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*(qd(i,k+1)-qhat(i,k+1)))
         end if
      end do
   end do
!
   beta = 0._kind_phys
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k > jt(i) .and. k < mx(i)) then
            dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+mc(i,k+1)* (s(i,k)-shat(i,k+1)))/ &
                         dp(i,k) - rl/cp*du(i,k)*(beta*ql(i,k)+ (1-beta)*ql(i,k+1))

            dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+cp/rl* (su(i,k+1)-s(i,k)))- &
                          mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*(su(i,k)-s(i,k)))+md(i,k+1)* &
                         (qd(i,k+1)-qhat(i,k+1)+cp/rl*(sd(i,k+1)-s(i,k)))-md(i,k)* &
                         (qd(i,k)-qhat(i,k)+cp/rl*(sd(i,k)-s(i,k))))/dp(i,k) + &
                          du(i,k)* (beta*ql(i,k)+(1-beta)*ql(i,k+1))
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k >= lel(i) .and. k <= lcl(i)) then
            thetavp(i,k) = tp(i,k)* (1000._kind_phys/p(i,k))** (rd/cp)*(1._kind_phys+1.608_kind_phys*qstp(i,k)-q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000._kind_phys/p(i,k))** (rd/cp)*(1._kind_phys+0.608_kind_phys*q(i,k))
            dqsdtp(i,k) = qstp(i,k)* (1._kind_phys+qstp(i,k)/eps1)*eps1*rl/(rd*tp(i,k)**2)
!
! dtpdt is the parcel temperature change due to change of
! subcloud layer properties during convection.
!
            dtpdt(i,k) = tp(i,k)/ (1._kind_phys+rl/cp* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))* &
                        (dtbdt(i)/t(i,mx(i))+rl/cp* (dqbdt(i)/tl(i)-q(i,mx(i))/ &
                         tl(i)**2*dtldt(i)))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1._kind_phys/(1._kind_phys+1.608_kind_phys*qstp(i,k)-q(i,mx(i)))* &
                          (1.608_kind_phys * dqsdtp(i,k) * dtpdt(i,k) -dqbdt(i))) - (dtmdt(i,k)/t(i,k)+0.608_kind_phys/ &
                          (1._kind_phys+0.608_kind_phys*q(i,k))*dqmdt(i,k)))*grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k > lcl(i) .and. k < mx(i)) then
            thetavp(i,k) = tp(i,k)* (1000._kind_phys/p(i,k))** (rd/cp)*(1._kind_phys+0.608_kind_phys*q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000._kind_phys/p(i,k))** (rd/cp)*(1._kind_phys+0.608_kind_phys*q(i,k))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = (dtbdt(i)/t(i,mx(i))+0.608_kind_phys/ (1._kind_phys+0.608_kind_phys*q(i,mx(i)))*dqbdt(i)- &
                          dtmdt(i,k)/t(i,k)-0.608_kind_phys/ (1._kind_phys+0.608_kind_phys*q(i,k))*dqmdt(i,k))* &
                          grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do

!
! buoyant energy change is set to 2/3*excess cape per 3 hours
!
   dadt(il1g:il2g)  = 0._kind_phys
   kmin = minval(lel(il1g:il2g))
   kmax = maxval(mx(il1g:il2g)) - 1
   do k = kmin, kmax
      do i = il1g,il2g
         if ( k >= lel(i) .and. k <= mx(i) - 1) then
            dadt(i) = dadt(i) + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
         endif
      end do
   end do
   do i = il1g,il2g
      dltaa = -1._kind_phys* (cape(i)-capelmt)
      if (dadt(i) /= 0._kind_phys) mb(i) = max(dltaa/tau/dadt(i),0._kind_phys)
   end do
!
   return
end subroutine closure

subroutine q1q2_pjr(ncol   ,pver    ,latice  ,&
                    dqdt    ,dsdt    ,q       ,qs      ,qu      , &
                    su      ,du      ,qhat    ,shat    ,dp      , &
                    mu      ,md      ,sd      ,qd      ,ql      , &
                    dsubcld ,jt      ,mx      ,il1g    ,il2g    , &
                    cp      ,rl      ,msg     ,          &
                    dl      ,evp     ,cu)

   implicit none

!-----------------------------------------------------------------------
! CACNOTE -fill in documentation
!
! Purpose:
! <Say what the routine does>
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: phil rasch dec 19 1995
!
!-----------------------------------------------------------------------


   real(kind_phys), intent(in) :: cp

   integer, intent(in) :: ncol
   integer, intent(in) :: pver
   real(kind_phys), intent(in) :: latice
   integer, intent(in) :: il1g
   integer, intent(in) :: il2g
   integer, intent(in) :: msg

   real(kind_phys), intent(in) :: q(ncol,pver)
   real(kind_phys), intent(in) :: qs(ncol,pver)
   real(kind_phys), intent(in) :: qu(ncol,pver)
   real(kind_phys), intent(in) :: su(ncol,pver)
   real(kind_phys), intent(in) :: du(ncol,pver)
   real(kind_phys), intent(in) :: qhat(ncol,pver)
   real(kind_phys), intent(in) :: shat(ncol,pver)
   real(kind_phys), intent(in) :: dp(ncol,pver)
   real(kind_phys), intent(in) :: mu(ncol,pver)
   real(kind_phys), intent(in) :: md(ncol,pver)
   real(kind_phys), intent(in) :: sd(ncol,pver)
   real(kind_phys), intent(in) :: qd(ncol,pver)
   real(kind_phys), intent(in) :: ql(ncol,pver)
   real(kind_phys), intent(in) :: evp(ncol,pver)
   real(kind_phys), intent(in) :: cu(ncol,pver)
   real(kind_phys), intent(in) :: dsubcld(ncol)

   real(kind_phys),intent(out) :: dqdt(ncol,pver),dsdt(ncol,pver)
   real(kind_phys),intent(out) :: dl(ncol,pver)

   integer kbm
   integer ktm
   integer jt(ncol)
   integer mx(ncol)
!
! work fields:
!
   integer i
   integer k

   real(kind_phys) emc
   real(kind_phys) rl
!-------------------------------------------------------------------
   do k = msg + 1,pver
      do i = il1g,il2g
         dsdt(i,k) = 0._kind_phys
         dqdt(i,k) = 0._kind_phys
         dl(i,k) = 0._kind_phys
      end do
   end do

!
! find the highest level top and bottom levels of convection
!
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

   do k = ktm,pver-1
      do i = il1g,il2g
         emc = -cu (i,k)               &         ! condensation in updraft
               +evp(i,k)                         ! evaporating rain in downdraft

         dsdt(i,k) = -rl/cp*emc &
                     + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1)) &
                        -mu(i,k)*   (su(i,k)-shat(i,k)) &
                        +md(i,k+1)* (sd(i,k+1)-shat(i,k+1)) &
                        -md(i,k)*   (sd(i,k)-shat(i,k)) &
                       )/dp(i,k)

         dqdt(i,k) = emc + &
                    (+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)) &
                     -mu(i,k)*   (qu(i,k)-qhat(i,k)) &
                     +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1)) &
                     -md(i,k)*   (qd(i,k)-qhat(i,k)) &
                    )/dp(i,k)

         dl(i,k) = du(i,k)*ql(i,k+1)

      end do
   end do

!
   do k = kbm,pver
      do i = il1g,il2g
         if (k == mx(i)) then
            dsdt(i,k) = (1._kind_phys/dsubcld(i))* &
                        (-mu(i,k)* (su(i,k)-shat(i,k)) &
                         -md(i,k)* (sd(i,k)-shat(i,k)) &
                        )
            dqdt(i,k) = (1._kind_phys/dsubcld(i))* &
                        (-mu(i,k)*(qu(i,k)-qhat(i,k)) &
                         -md(i,k)*(qd(i,k)-qhat(i,k)) &
                        )
         else if (k > mx(i)) then
            dsdt(i,k) = dsdt(i,k-1)
            dqdt(i,k) = dqdt(i,k-1)
         end if
      end do
   end do
!
   return
end subroutine q1q2_pjr


! Wrapper for qsat_water that does translation between Pa and hPa
! qsat_water uses Pa internally, so get it right, need to pass in Pa.
! Afterward, set es back to hPa.
subroutine qsat_hPa(t, p, es, qm)
!CACNOTE - Need to figure out how to handle this
  use wv_saturation, only: qsat_water

  ! Inputs
  real(kind_phys), intent(in) :: t    ! Temperature (K)
  real(kind_phys), intent(in) :: p    ! Pressure (hPa)
  ! Outputs
  real(kind_phys), intent(out) :: es  ! Saturation vapor pressure (hPa)
  real(kind_phys), intent(out) :: qm  ! Saturation mass mixing ratio
                               ! (vapor mass over dry mass, kg/kg)

  call qsat_water(t, p*100._kind_phys, es, qm)

  es = es*0.01_kind_phys

end subroutine qsat_hPa

end module zm_convr
