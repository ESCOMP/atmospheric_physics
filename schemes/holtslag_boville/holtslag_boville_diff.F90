! Module to compute mixing coefficients associated with turbulence in the 
! planetary boundary layer and elsewhere.  PBL coefficients are based on Holtslag 
! and Boville, 1991.
!
! Original authors:
! Standardized: J. Rosinski, June 1992
! Reviewed:     P. Rasch, B. Boville, August 1992
! Reviewed:     P. Rasch, April 1996
! Reviewed:     B. Boville, April 1996
! Rewritten:    B. Boville, May 2000
! Rewritten:    B. Stevens, August 2000
! Modularized:  J. McCaa, September 2004
! Functional Updates from M. Waxmonsky, February 2025
! CCPP-ized:    Haipeng Lin, May 2025
module holtslag_boville_diff
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: holtslag_boville_diff_init
  public :: hb_pbl_independent_coefficients_run        ! formerly trbintd
  public :: hb_pbl_dependent_coefficients_run          ! formerly pblintd (also used by CLUBB)
  public :: hb_diff_exchange_coefficients_run          ! formerly austausch_pbl and compute_hb_diff logic
  public :: hb_diff_free_atm_exchange_coefficients_run ! formerly compute_hb_free_atm_diff logic
  public :: holtslag_boville_diff_finalize

  ! tuning parameters used by HB and related subroutines
  ! PBL limits
  real(kind_phys), parameter :: pblmaxp = 4.e4_kind_phys            ! PBL max depth [Pa]
  real(kind_phys), parameter :: zkmin   = 0.01_kind_phys            ! Minimum kneutral*f(ri)

  ! PBL parameters
  real(kind_phys), parameter :: onet    = 1._kind_phys/3._kind_phys ! 1/3 power in wind gradient expression
  real(kind_phys), parameter :: betam   = 15.0_kind_phys            ! Constant in wind gradient expression
  real(kind_phys), parameter :: betas   =  5.0_kind_phys            ! Constant in surface layer gradient expression
  real(kind_phys), parameter :: betah   = 15.0_kind_phys            ! Constant in temperature gradient expression
  real(kind_phys), parameter :: fakn    =  7.2_kind_phys            ! Constant in turbulent prandtl number
  real(kind_phys), parameter :: fak     =  8.5_kind_phys            ! Constant in surface temperature excess
  real(kind_phys), parameter :: ricr    =  0.3_kind_phys            ! Critical richardson number
  real(kind_phys), parameter :: sffrac  =  0.1_kind_phys            ! Surface layer fraction of boundary layer
  real(kind_phys), parameter :: binm    = betam*sffrac              ! betam * sffrac
  real(kind_phys), parameter :: binh    = betah*sffrac              ! betah * sffrac

  ! derived parameters from initialization
  real(kind_phys), allocatable :: ml2(:)      ! mixing lengths squared [m^2] at interfaces
  real(kind_phys)              :: ccon        ! fak * sffrac * karman
  integer                      :: npbl        ! maximum # of levels in PBL from surface
  integer                      :: ntop_turb   ! top    level to which turbulent vertical diffusion is applied
  integer                      :: nbot_turb   ! bottom level to which turbulent vertical diffusion is applied

contains

   ! Initialization subroutine to set module parameters, including mixing range
!> \section arg_table_holtslag_boville_diff_init Argument Table
!! \htmlinclude arg_table_holtslag_boville_diff_init.html
  subroutine holtslag_boville_diff_init( &
    amIRoot, iulog, &
    pver, pverp, &
    karman, &
    pref_mid, &
    is_hbr_pbl_scheme, &
    ntop_turb_in, &
    errmsg, errflg)

    ! Input arguments
    logical,            intent(in)    :: amIRoot           ! are we on the MPI root task?
    integer,            intent(in)    :: iulog             ! log output unit
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: pverp
    real(kind_phys),    intent(in)    :: karman            ! von_karman_constant [1]
    real(kind_phys),    intent(in)    :: pref_mid(:)       ! reference_pressure_in_atmosphere_layer [Pa]
    logical,            intent(in)    :: is_hbr_pbl_scheme ! is HBR = true; is HB = false [flag]
    integer,            intent(in)    :: ntop_turb_in      ! top turbulence level [index]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg            ! error message
    integer,            intent(out)   :: errflg            ! error flag

    ! Local variables
    integer :: k

    errmsg = ''
    errflg = 0

    ccon = fak * sffrac * karman

    ! ntop_turb must be 1 or <= nbot_molec (lowest vertical level with molecular diffusion)
    ! In WACCM-X, when 'waccmx_mode' is 'ionosphere' or 'neutral', this should be set to
    ! press_lim_idx(ntop_eddy_pres, top=.true.) where ntop_eddy_pres = 1.e-7_kind_phys [Pa]
    !
    ! FIXME: this is not as clean as I would like it to be
    ! but I want to avoid depending on host model press_lim_idx; so the index is passed
    ! into this initialization subroutine. This also avoids having to query something like
    ! 'waccmx_is' or 'waccmx_mode' that is host model specific inside this scheme
    ! to make it portable, but this means the flow on the host model side has to
    ! specify this top turbulence index. (hplin, 5/8/25)
    ntop_turb = ntop_turb_in

    nbot_turb = pver

    ! allocate and populate mixing lengths squared to grid vertical interfaces
    allocate(ml2(pverp), stat=errflg, errmsg=errmsg)
    if(errflg /= 0) return

    ml2(ntop_turb) = 0._kind_phys
    do k = ntop_turb+1, nbot_turb
      ml2(k) = 30.0_kind_phys**2                 ! HB scheme:  length scale = 30m
      if(is_hbr_pbl_scheme) then
         ml2(k) = 1.0_kind_phys**2               ! HBR scheme: length scale = 1m
      end if
    end do
    ml2(nbot_turb+1) = 0._kind_phys

    ! Limit pbl height to regions below 400 mb
    ! npbl = max # of levels (from bottom) in pbl
    npbl = 0
    do k = nbot_turb,ntop_turb,-1
       if (pref_mid(k) >= pblmaxp) then
          npbl = npbl + 1
       end if
    end do
    npbl = max(npbl, 1)

    if(amIRoot) then
      write(iulog,*) 'Holtslag-Boville PBL: PBL height will be limited to bottom ', npbl, &
                     ' model levels. Top is ', pref_mid(pver+1-npbl), ' pascals'
    end if

  end subroutine holtslag_boville_diff_init

  ! Computation of time-dependent, PBL height-independent variables for HB mixing.
  ! Original author: B. Stevens, rewrite August 2000
!> \section arg_table_hb_pbl_independent_coefficients_run Argument Table
!! \htmlinclude arg_table_hb_pbl_independent_coefficients_run.html
  pure subroutine hb_pbl_independent_coefficients_run( &
    ncol, pver, &
    zvir, rair, cpair, gravit, karman, &
    exner, t, &
    q_wv, &
    z, &
    pmid, &
    u, v, &
    taux, tauy, &
    shflx, q_wv_flx, &
    ! below output
    thv, ustar, &
    khfs, kqfs, kbfs, &
    obklen, s2, ri, &
    errmsg, errflg)
    use atmos_phys_pbl_utils, only: calc_virtual_temperature, calc_friction_velocity, calc_obukhov_length, &
                                    calc_ideal_gas_rrho, &
                                    calc_kinematic_heat_flux, calc_kinematic_water_vapor_flux, &
                                    calc_kinematic_buoyancy_flux

    ! Input arguments
    integer,         intent(in)  :: ncol                     ! # of atmospheric columns
    integer,         intent(in)  :: pver                     ! # of vertical levels
    real(kind_phys), intent(in)  :: zvir
    real(kind_phys), intent(in)  :: rair
    real(kind_phys), intent(in)  :: cpair
    real(kind_phys), intent(in)  :: gravit
    real(kind_phys), intent(in)  :: karman
    real(kind_phys), intent(in)  :: exner(:,:)               ! exner [1]
    real(kind_phys), intent(in)  :: t   (:,:)                ! temperature [K]
    real(kind_phys), intent(in)  :: q_wv(:,:)                ! specific humidity [kg kg-1]
    real(kind_phys), intent(in)  :: z   (:,:)                ! height above surface [m]
    real(kind_phys), intent(in)  :: u   (:,:)                ! zonal velocity [m s-1]
    real(kind_phys), intent(in)  :: v   (:,:)                ! meridional velocity [m s-1]
    real(kind_phys), intent(in)  :: taux (:)                 ! zonal stress [N m-2]
    real(kind_phys), intent(in)  :: tauy (:)                 ! meridional stress [N m-2]
    real(kind_phys), intent(in)  :: shflx(:)                 ! sensible heat flux [W m-2]
    real(kind_phys), intent(in)  :: q_wv_flx(:)              ! upward water vapor flux at surface [kg m-2 s-1]
    real(kind_phys), intent(in)  :: pmid(:,:)                ! midpoint pressures

    ! Output arguments
    real(kind_phys), intent(out) :: thv(:,:)                 ! virtual potential temperature [K]
    real(kind_phys), intent(out) :: ustar(:)                 ! surface friction velocity [m s-1]
    real(kind_phys), intent(out) :: khfs(:)                  ! kinematic surface heat flux [K m s-1]
    real(kind_phys), intent(out) :: kqfs(:)                  ! kinematic surface water vapor flux [kg kg-1 m s-1]
    real(kind_phys), intent(out) :: kbfs(:)                  ! surface kinematic buoyancy flux [m^2 s-3]
    real(kind_phys), intent(out) :: obklen(:)                ! Obukhov length [m]
    real(kind_phys), intent(out) :: s2(:,:)                  ! shear squared [s-2]
    real(kind_phys), intent(out) :: ri(:,:)                  ! richardson number: n2/s2 [1]
    character(len=512), intent(out)   :: errmsg              ! error message
    integer,            intent(out)   :: errflg              ! error flag

    integer :: i, k
    real(kind_phys) :: dvdz2                                 ! velocity shear squared [m^2 s-2]
    real(kind_phys) :: dz                                    ! delta z between midpoints [m]
    real(kind_phys) :: rrho(ncol)                            ! 1 / bottom level density [m^3 kg-1]
    real(kind_phys) :: n2(ncol, pver)                        ! brunt vaisaila frequency [s-2]
    real(kind_phys) :: th(ncol, pver)                        ! potential temperature [K]

    errmsg = ''
    errflg = 0

    ! calculate potential temperature [K]
    th(:ncol,:) = t(:ncol,:) * exner(:ncol,:)

    ! virtual temperature
    thv(:,:) = 0._kind_phys
    thv(:ncol,ntop_turb:) = calc_virtual_temperature(th(:ncol,ntop_turb:), q_wv(:ncol,ntop_turb:), zvir)

    ! Compute ustar, Obukhov length, and kinematic surface fluxes.
    rrho(:ncol)   = calc_ideal_gas_rrho(rair, t(:ncol,pver), pmid(:ncol,pver))
    ustar(:ncol)  = calc_friction_velocity(taux(:ncol),tauy(:ncol), rrho(:ncol))
    khfs(:ncol)   = calc_kinematic_heat_flux(shflx(:ncol), rrho(:ncol), cpair)
    kqfs(:ncol)   = calc_kinematic_water_vapor_flux(q_wv_flx(:ncol), rrho(:ncol))
    kbfs(:ncol)   = calc_kinematic_buoyancy_flux(khfs(:ncol), zvir, th(:ncol,pver), kqfs(:ncol))
    obklen(:ncol) = calc_obukhov_length(thv(:ncol,pver), ustar(:ncol), gravit, karman, kbfs(:ncol))

    ! Compute s^2 (shear squared), n^2 (brunt vaisaila frequency), and ri (Richardson number = n^2/s^2)
    ! using virtual temperature
    ! (formerly trbintd)
    do k = ntop_turb,nbot_turb-1
       do i = 1,ncol
          dvdz2   = (u(i,k)-u(i,k+1))**2 + (v(i,k)-v(i,k+1))**2
          dvdz2   = max(dvdz2,1.e-36_kind_phys)
          dz      = z(i,k) - z(i,k+1)
          s2(i,k) = dvdz2/(dz**2)
          n2(i,k) = gravit*2.0_kind_phys*(thv(i,k) - thv(i,k+1))/((thv(i,k) + thv(i,k+1))*dz)
          ri(i,k) = n2(i,k)/s2(i,k)
       end do
    end do
  end subroutine hb_pbl_independent_coefficients_run

  ! Compute time-dependent, PBL height-dependent variables
  ! (formerly pblintd)
  ! The PBL depth follows:
  !    Holtslag, A.A.M., and B.A. Boville, 1993:
  !    Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate Model.
  !    J. Clim., vol. 6., p. 1825--1842. https://doi.org/10.1175/1520-0442(1993)006<1825:LVNBLD>2.0.CO;2
  !
  ! Updated by Holtslag and Hack to exclude the surface layer from the
  ! definition of the boundary layer Richardson number. Ri is now defined
  ! across the outer layer of the pbl (between the top of the surface
  ! layer and the pbl top) instead of the full pbl (between the surface and
  ! the pbl top). For simiplicity, the surface layer is assumed to be the
  ! region below the first model level (otherwise the boundary layer depth
  ! determination would require iteration).
  !
  ! Original author: B. Stevens, August 2000, extracted from pbldiff
!> \section arg_table_hb_pbl_dependent_coefficients_run Argument Table
!! \htmlinclude arg_table_hb_pbl_dependent_coefficients_run.html
  pure subroutine hb_pbl_dependent_coefficients_run( &
    ncol, pver, pverp, &
    gravit, &
    z, zi, &
    u, v, &
    cldn, &
    ! below inputs from pbl_independent_coefficients
    thv, ustar, kbfs, obklen, &
    ! below output
    pblh, &
    wstar, bge, &
    errmsg, errflg)

    ! Input arguments
    integer,         intent(in)  :: ncol                     ! # of atmospheric columns
    integer,         intent(in)  :: pver                     ! # of vertical levels
    integer,         intent(in)  :: pverp                    ! # of vertical level interfaces
    real(kind_phys), intent(in)  :: gravit
    real(kind_phys), intent(in)  :: z   (:,:)                ! height above surface [m]
    real(kind_phys), intent(in)  :: zi  (:,:)                ! height above surface [m], interfaces
    real(kind_phys), intent(in)  :: u   (:,:)                ! zonal velocity [m s-1]
    real(kind_phys), intent(in)  :: v   (:,:)                ! meridional velocity [m s-1]
    real(kind_phys), intent(in)  :: cldn(:,:)                ! stratiform cloud fraction [fraction]

    ! Input arguments (output from hb_pbl_independent_coefficients)
    real(kind_phys), intent(in)  :: thv (:,:)                ! virtual potential temperature [K]
    real(kind_phys), intent(in)  :: ustar(:)                 ! surface friction velocity [m s-1]
    real(kind_phys), intent(in)  :: kbfs(:)                  ! surface kinematic buoyancy flux [m^2 s-3]
    real(kind_phys), intent(in)  :: obklen(:)                ! Obukhov length [m]

    ! Output arguments
    real(kind_phys), intent(out) :: pblh(:)                  ! boundary-layer height [m]
    real(kind_phys), intent(out) :: wstar(:)                 ! convective scale velocity [m s-1]
    real(kind_phys), intent(out) :: bge(:)                   ! buoyancy gradient enhancement [m s-2]
    character(len=512), intent(out)   :: errmsg              ! error message
    integer,            intent(out)   :: errflg              ! error flag

    ! Local variables
    integer :: i, k
    real(kind_phys) :: phiminv(ncol)   ! inverse phi function for momentum
    real(kind_phys) :: rino(ncol,pver) ! bulk Richardson no. from level to ref lev
    real(kind_phys) :: tlv(ncol)       ! ref. level potential tmp + tmp excess
    real(kind_phys) :: vvk             ! velocity magnitude squared

    !logical  :: unstable(ncol)         ! points with unstable pbl (positive virtual heat flux)
    logical  :: check(ncol)            ! false if Richardson number > critical
    logical  :: ocncldcheck(ncol)      ! true if ocean surface (not implemented) and cloud in lowest layer

    ! Local parameters
    real(kind_phys), parameter :: tiny = 1.e-36_kind_phys   ! lower bound for wind magnitude
    real(kind_phys), parameter :: fac  = 100._kind_phys     ! ustar parameter in height diagnosis

    errmsg = ''
    errflg = 0

    ! Compute Obukhov length virtual temperature flux and various arrays for use later:
    do i=1,ncol
       check(i)     = .true.
       rino(i,pver) = 0.0_kind_phys
       pblh(i)      = z(i,pver)
    end do

    ! PBL height calculation:  Scan upward until the Richardson number between
    ! the first level and the current level exceeds the "critical" value.
    do k=pver-1,pver-npbl+1,-1
       do i=1,ncol
          if (check(i)) then
             vvk = (u(i,k) - u(i,pver))**2 + (v(i,k) - v(i,pver))**2 + fac*ustar(i)**2
             vvk = max(vvk,tiny)
             rino(i,k) = gravit*(thv(i,k) - thv(i,pver))*(z(i,k)-z(i,pver))/(thv(i,pver)*vvk)
             ! Modified for boundary layer height diagnosis: Bert Holtslag, June 1994
             ! >>>>>>>>>  (Use ricr = 0.3 in this formulation)
             if (rino(i,k) >= ricr) then
                pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/(rino(i,k) - rino(i,k+1)) * &
                     (z(i,k) - z(i,k+1))
                check(i) = .false.
             end if
          end if
       end do
    end do

    ! Estimate an effective surface temperature to account for surface fluctuations
    do i=1,ncol
       if (check(i)) pblh(i) = z(i,pverp-npbl)
       !unstable(i) = (kbfs(i) > 0._kind_phys)
       check(i)  = (kbfs(i) > 0._kind_phys)
       if (check(i)) then
          phiminv(i)   = (1._kind_phys - binm*pblh(i)/obklen(i))**onet
          rino(i,pver) = 0.0_kind_phys
          tlv(i)       = thv(i,pver) + kbfs(i)*fak/(ustar(i)*phiminv(i))
       end if
    end do

    ! Improve pblh estimate for unstable conditions using the convective temperature excess:
    do i = 1,ncol
      bge(i) = 1.e-8_kind_phys
    end do

    do k=pver-1,pver-npbl+1,-1
      do i=1,ncol
        if (check(i)) then
          vvk = (u(i,k) - u(i,pver))**2 + (v(i,k) - v(i,pver))**2 + fac*ustar(i)**2
          vvk = max(vvk,tiny)
          rino(i,k) = gravit*(thv(i,k) - tlv(i))*(z(i,k)-z(i,pver))/(thv(i,pver)*vvk)
          if (rino(i,k) >= ricr) then
            pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/(rino(i,k) - rino(i,k+1)) * &
                 (z(i,k) - z(i,k+1))
            bge(i) = 2._kind_phys*gravit/(thv(i,k)+thv(i,k+1))*(thv(i,k)-thv(i,k+1))/(z(i,k)-z(i,k+1))*pblh(i)
            if (bge(i).lt.0._kind_phys) then
              bge(i) = 1.e-8_kind_phys
            endif
            check(i) = .false.
          end if
        end if
      end do
    end do

    ! PBL height must be greater than some minimum mechanical mixing depth
    ! Several investigators have proposed minimum mechanical mixing depth
    ! relationships as a function of the local friction velocity, u*.  We
    ! make use of a linear relationship of the form h = c u* where c=700.
    ! The scaling arguments that give rise to this relationship most often
    ! represent the coefficient c as some constant over the local coriolis
    ! parameter.  Here we make use of the experimental results of Koracin
    ! and Berkowicz (1988) [BLM, Vol 43]: https://doi.org/10.1007/BF00153969
    ! for wich they recommend 0.07/f
    ! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
    ! latitude value for f so that c = 0.07/f = 700. Also, do not allow
    ! PBL to exceed some maximum (npbl) number of allowable points
    do i=1,ncol
       if (check(i)) pblh(i) = z(i,pverp-npbl)
       pblh(i) = max(pblh(i),700.0_kind_phys*ustar(i))
       wstar(i) = (max(0._kind_phys,kbfs(i))*gravit*pblh(i)/thv(i,pver))**onet
    end do

    ! Final requirement on PBL height is that it must be greater than the depth
    ! of the lowest model level over ocean if there is any cloud diagnosed in
    ! the lowest model level.  This is to deal with the inadequacies of the
    ! current "dry" formulation of the boundary layer, where this test is
    ! used to identify circumstances where there is marine stratus in the
    ! lowest level, and to provide a weak ventilation of the layer to avoid
    ! a pathology in the cloud scheme (locking in low-level stratiform cloud)
    ! If over an ocean surface, and any cloud is diagnosed in the
    ! lowest level, set pblh to 50 meters higher than top interface of lowest level
    !
    !  jrm This is being applied everywhere (not just ocean)!
    do i=1,ncol
       ocncldcheck(i) = .false.
       if (cldn(i,pver).ge.0.0_kind_phys) ocncldcheck(i) = .true.
       if (ocncldcheck(i)) pblh(i) = max(pblh(i),zi(i,pver) + 50._kind_phys)
    end do
  end subroutine hb_pbl_dependent_coefficients_run

  ! Atmosphere boundary layer computation.
  !
  ! Nonlocal scheme that determines eddy diffusivities based on a
  ! specified boundary layer height and a turbulent velocity scale;
  ! also, countergradient effects for heat and moisture, and constituents
  ! are included, along with temperature and humidity perturbations which
  ! measure the strength of convective thermals in the lower part of the
  ! atmospheric boundary layer.
  !
  ! For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
  ! Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate Model.
  ! J. Clim., vol. 6., p. 1825--1842. https://doi.org/10.1175/1520-0442(1993)006<1825:LVNBLD>2.0.CO;2
  !
  ! Updated by Holtslag and Hack to exclude the surface layer from the
  ! definition of the boundary layer Richardson number. Ri is now defined
  ! across the outer layer of the pbl (between the top of the surface
  ! layer and the pbl top) instead of the full pbl (between the surface and
  ! the pbl top). For simiplicity, the surface layer is assumed to be the
  ! region below the first model level (otherwise the boundary layer depth
  ! determination would require iteration).
  !
  ! Original authors: B. Boville, B. Stevens, rewrite August 2000
!> \section arg_table_hb_diff_exchange_coefficients_run Argument Table
!! \htmlinclude arg_table_hb_diff_exchange_coefficients_run.html
  pure subroutine hb_diff_exchange_coefficients_run( &
    ncol, pver, pverp, &
    karman, cpair, &
    z, &
    is_hbr_pbl_scheme, &
    ! input from hb_pbl_independent_coefficients
    kqfs, khfs, kbfs, &
    ustar, obklen, &
    s2, ri, &
    ! input from hb_pbl_dependent_coefficients
    pblh, wstar, bge, &
    ! below output
    kvm, kvh, kvq, &
    cgh, cgs, &
    tpert, qpert, &
    tke, &
    errmsg, errflg)

    use atmos_phys_pbl_utils, only: calc_eddy_flux_coefficient

    ! Input arguments
    integer,         intent(in)  :: ncol                     ! # of atmospheric columns
    integer,         intent(in)  :: pver                     ! # of vertical levels
    integer,         intent(in)  :: pverp                    ! # of vertical level interfaces
    real(kind_phys), intent(in)  :: karman
    real(kind_phys), intent(in)  :: cpair
    real(kind_phys), intent(in)  :: z   (:,:)                ! height above surface [m]
    logical,         intent(in)  :: is_hbr_pbl_scheme

    real(kind_phys), intent(in)  :: khfs(:)                  ! kinematic surface heat flux [K m s-1]
    real(kind_phys), intent(in)  :: kqfs(:)                  ! kinematic surface constituent flux [kg kg-1 m s-1]
    real(kind_phys), intent(in)  :: kbfs(:)                  ! surface kinematic buoyancy flux [m^2 s-3]
    real(kind_phys), intent(in)  :: ustar(:)                 ! surface friction velocity [m s-1]
    real(kind_phys), intent(in)  :: obklen(:)                ! Obukhov length [m]
    real(kind_phys), intent(in)  :: s2(:,:)                  ! shear squared [s-2]
    real(kind_phys), intent(in)  :: ri(:,:)                  ! richardson number: n2/s2 [1]

    real(kind_phys), intent(in)  :: pblh(:)                  ! boundary-layer height [m]
    real(kind_phys), intent(in)  :: wstar(:)                 ! convective scale velocity [m s-1]
    real(kind_phys), intent(in)  :: bge(:)                   ! buoyancy gradient enhancement [m s-2]

    ! Output variables
    real(kind_phys), intent(out) :: kvm(:,:)                 ! eddy diffusivity for momentum [m^2 s-1], interfaces
    real(kind_phys), intent(out) :: kvh(:,:)                 ! eddy diffusivity for heat [m^2 s-1], interfaces
    real(kind_phys), intent(out) :: kvq(:,:)                 ! eddy diffusivity for constituents [m^2 s-1], interfaces
    real(kind_phys), intent(out) :: cgh(:,:)                 ! counter-gradient term for heat [J kg-1 m-1], interfaces
    real(kind_phys), intent(out) :: cgs(:,:)                 ! counter-gradient star (cg/flux) [s m-2], interfaces
    real(kind_phys), intent(out) :: tpert(:)                 ! convective temperature excess [K]
    real(kind_phys), intent(out) :: qpert(:)                 ! convective humidity excess [kg kg-1]
    real(kind_phys), intent(out) :: tke(:,:)                 ! (estimated) turbulent kinetic energy [m^2 s-2], interfaces
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    integer :: i, k
    real(kind_phys) :: kvf(ncol,pverp) ! free atmospheric eddy diffusivity [m^2 s-1]
    integer  :: ktopbl(ncol)           ! index of first midpoint inside PBL (diagnostic) [index]

    real(kind_phys) :: phiminv(ncol)   ! inverse phi function for momentum
    real(kind_phys) :: phihinv(ncol)   ! inverse phi function for heat
    real(kind_phys) :: wm(ncol)        ! turbulent velocity scale for momentum
    real(kind_phys) :: zp(ncol)        ! current level height + one level up
    real(kind_phys) :: fak1(ncol)      ! k*ustar*pblh
    real(kind_phys) :: fak2(ncol)      ! k*wm*pblh
    real(kind_phys) :: fak3(ncol)      ! fakn*wstar/wm
    real(kind_phys) :: pblk(ncol)      ! level eddy diffusivity for momentum
    real(kind_phys) :: pr(ncol)        ! Prandtl number for eddy diffusivities
    real(kind_phys) :: zl(ncol)        ! zmzp / Obukhov length
    real(kind_phys) :: zh(ncol)        ! zmzp / pblh
    real(kind_phys) :: zzh(ncol)       ! (1-(zmzp/pblh))**2
    real(kind_phys) :: zmzp            ! level height halfway between zm and zp
    real(kind_phys) :: term            ! intermediate calculation
    real(kind_phys) :: kve             ! diffusivity at entrainment layer in unstable cases [m s-2]

    logical  :: unstable(ncol)         ! points with unstable pbl (positive virtual heat flux) [count]
    logical  :: pblpt(ncol)            ! points within pbl

    errmsg = ''
    errflg = 0

    ! Get atmosphere exchange coefficients
    kvf(:ncol,:) = 0.0_kind_phys
    do k = ntop_turb, nbot_turb-1
       do i = 1, ncol
          kvf(i,k+1) = calc_eddy_flux_coefficient(ml2(k), ri(i, k), s2(i, k))
       end do
    end do

    ! Compute PBL exchange coefficients (formerly austausch_pbl)
    kvm = 0._kind_phys
    kvh = 0._kind_phys
    kve = 0._kind_phys
    cgh = 0._kind_phys
    cgs = 0._kind_phys
    tpert = 0._kind_phys
    qpert = 0._kind_phys
    ktopbl = 0._kind_phys
    tke = 0._kind_phys

    ! Initialize height independent arrays
    do i=1,ncol
      unstable(i) = (kbfs(i) > 0._kind_phys)
      pblk(i) = 0.0_kind_phys
      fak1(i) = ustar(i)*pblh(i)*karman
      if (unstable(i)) then
        phiminv(i) = (1._kind_phys - binm*pblh(i)/obklen(i))**onet
        phihinv(i) = sqrt(1._kind_phys - binh*pblh(i)/obklen(i))
        wm(i)      = ustar(i)*phiminv(i)
        fak2(i)    = wm(i)*pblh(i)*karman
        fak3(i)    = fakn*wstar(i)/wm(i)
        tpert(i)   = max(khfs(i)*fak/wm(i),0._kind_phys)
        qpert(i)   = max(kqfs(i)*fak/wm(i),0._kind_phys)
      else
        tpert(i)   = max(khfs(i)*fak/ustar(i),0._kind_phys)
        qpert(i)   = max(kqfs(i)*fak/ustar(i),0._kind_phys)
      end if
    end do

    ! Initialize output arrays with free atmosphere values
    do k=1,pverp
       do i=1,ncol
          kvm(i,k) = kvf(i,k)
          kvh(i,k) = kvf(i,k)
          cgh(i,k) = 0.0_kind_phys
          cgs(i,k) = 0.0_kind_phys
       end do
    end do

    ! Main level loop to compute the diffusivities and counter-gradient terms. These terms are
    ! only calculated at points determined to be in the interior of the pbl (pblpt(i)==.true.),
    ! and then calculations are directed toward regime: stable vs unstable, surface vs outer
    ! layer.
    do k=pver,pver-npbl+2,-1
      do i=1,ncol
        pblpt(i) = (z(i,k) < pblh(i))
        if (pblpt(i)) then
          ktopbl(i) = k
          zp(i)  = z(i,k-1)
          if (zkmin == 0.0_kind_phys .and. zp(i) > pblh(i)) then
            zp(i) = pblh(i)
          endif
          zmzp    = 0.5_kind_phys*(z(i,k) + zp(i)) ! we think this is an approximation to the interface height (where KVs are calculated)
          zh(i)   = zmzp/pblh(i)
          zl(i)   = zmzp/obklen(i)
          zzh(i)  = zh(i)*max(0._kind_phys,(1._kind_phys - zh(i)))**2
          if (unstable(i)) then
            if (zh(i) < sffrac) then
               term     = (1._kind_phys - betam*zl(i))**onet
               pblk(i)  = fak1(i)*zzh(i)*term
               pr(i)    = term/sqrt(1._kind_phys - betah*zl(i))
            else
               pblk(i)  = fak2(i)*zzh(i)
               pr(i)    = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
               cgs(i,k) = fak3(i)/(pblh(i)*wm(i))
               cgh(i,k) = khfs(i)*cgs(i,k)*cpair
            end if
          else
            if (zl(i) <= 1._kind_phys) then
              pblk(i) = fak1(i)*zzh(i)/(1._kind_phys + betas*zl(i))
            else
              pblk(i) = fak1(i)*zzh(i)/(betas + zl(i))
            end if
            pr(i)    = 1._kind_phys
          end if
          kvm(i,k) = max(pblk(i),kvf(i,k))
          kvh(i,k) = max(pblk(i)/pr(i),kvf(i,k))
        end if
      end do
    end do

    ! Check whether last allowed midpoint is within PBL
    ! HBR scheme only
    if(is_hbr_pbl_scheme) then
      ! apply new diffusivity at entrainment zone
      do i = 1,ncol
        if (bge(i) > 1.e-7_kind_phys) then
          k = ktopbl(i)
          kve = 0.2_kind_phys*(wstar(i)**3+5._kind_phys*ustar(i)**3)/bge(i)
          kvm(i,k) = kve
          kvh(i,k) = kve
        end if
      end do
    end if

    ! Crude estimate of tke (tke=0 above boundary layer)
    do k = max(pverp-npbl,2),pverp
      do i = 1, ncol
        if (z(i,k-1) < pblh(i)) then
          tke(i,k) = (kvm(i,k) / pblh(i))**2
        endif
      end do
    end do

    ! Finalize
    kvq(:ncol,:) = kvh(:ncol,:)
  end subroutine hb_diff_exchange_coefficients_run

  ! A version of hb_diff_exchange_coefficients
  ! that only computes free atmosphere exchanges (no PBL computations)
  !
  ! Original authors: B. Stevens, rewrite August 2000
  !                   Thomas Toniazzo, Peter H. Lauritzen, June 2023
!> \section arg_table_hb_diff_free_atm_exchange_coefficients_run Argument Table
!! \htmlinclude arg_table_hb_diff_free_atm_exchange_coefficients_run.html
  pure subroutine hb_diff_free_atm_exchange_coefficients_run( &
    ncol, pver, pverp, &
    ! input from hb_pbl_independent_coefficients
    s2, ri, &
    ! top of CLUBB (or layer above where HB applies)
    bottom_boundary, &
    ! below output
    kvm, kvh, kvq, &
    cgh, cgs, &
    errmsg, errflg)

    use atmos_phys_pbl_utils, only: calc_free_atm_eddy_flux_coefficient

    ! Input arguments
    integer,         intent(in)  :: ncol         ! # of atmospheric columns
    integer,         intent(in)  :: pver         ! # of vertical levels
    integer,         intent(in)  :: pverp        ! # of vertical level interfaces
    real(kind_phys), intent(in)  :: s2(:,:)      ! shear squared [s-2]
    real(kind_phys), intent(in)  :: ri(:,:)      ! richardson number: n2/s2 [1]

    !REMOVECAM: change to integer once CAM snapshots are no longer used.
    real(kind_phys), intent(in)  :: bottom_boundary(:)  ! level above which HB runs [index]

    ! Output variables
    real(kind_phys), intent(out) :: kvm(:,:)     ! eddy diffusivity for momentum [m^2 s-1], interfaces
    real(kind_phys), intent(out) :: kvh(:,:)     ! eddy diffusivity for heat [m^2 s-1], interfaces
    real(kind_phys), intent(out) :: kvq(:,:)     ! eddy diffusivity for constituents [m^2 s-1], interfaces
    real(kind_phys), intent(out) :: cgh(:,:)     ! counter-gradient term for heat [J kg-1 m-1], interfaces
    real(kind_phys), intent(out) :: cgs(:,:)     ! counter-gradient star (cg/flux) [s m-2], interfaces
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    ! Local variables
    integer :: i, k
    real(kind_phys) :: kvf(ncol,pverp)           ! free atmospheric eddy diffusivity [m^2 s-1]
    integer :: bottom_boundary_int(ncol)

    errmsg = ''
    errflg = 0

    ! Get free atmosphere exchange coefficients
    kvf(:ncol,:) = 0.0_kind_phys
    do k = ntop_turb, nbot_turb-1
       do i = 1, ncol
          kvf(i,k+1) = calc_free_atm_eddy_flux_coefficient(ml2(k), ri(i, k), s2(i, k))
       end do
    end do

    kvq(:ncol,:) = kvf(:ncol,:)
    kvm(:ncol,:) = kvf(:ncol,:)
    kvh(:ncol,:) = kvf(:ncol,:)
    cgh(:ncol,:) = 0._kind_phys
    cgs(:ncol,:) = 0._kind_phys

    ! HB coefficients will be zeroed out in the layers below
    ! the bottom_boundary, which is the topmost layer where CLUBB is active.
    !   1     --- TOA ---
    !   2  .. HB (free atm) ..
    !   3
    !  ...
    !    --- bottom_boundary ---    | levels here
    !  ...     .. CLUBB ..          | have HB coefficients
    !  ...                          | zeroed out as CLUBB is active.
    ! pverp  --- surface ---
    bottom_boundary_int(:ncol) = int(bottom_boundary(:ncol))
    do i = 1, ncol
       do k = bottom_boundary_int(i), pverp
         kvm(i,k) = 0._kind_phys
         kvh(i,k) = 0._kind_phys
         kvq(i,k) = 0._kind_phys
       end do
    end do
  end subroutine hb_diff_free_atm_exchange_coefficients_run

!> \section arg_table_holtslag_boville_diff_finalize Argument Table
!! \htmlinclude arg_table_holtslag_boville_diff_finalize.html
  subroutine holtslag_boville_diff_finalize(errmsg, errflg)
    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    errmsg = ''
    errflg = 0

    if(allocated(ml2)) then
      deallocate(ml2)
    endif
  end subroutine holtslag_boville_diff_finalize

end module holtslag_boville_diff
