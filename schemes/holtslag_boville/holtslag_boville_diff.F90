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
  real(kind_phys)            :: ml2         ! mixing length squared [m2]
  real(kind_phys)            :: ccon        ! fak * sffrac * karman
  integer                    :: npbl        ! maximum # of levels in PBL from surface
  integer                    :: ntop_turb   ! top    level to which turbulent vertical diffusion is applied
  integer                    :: nbot_turb   ! bottom level to which turbulent vertical diffusion is applied

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

    ! mixing lengths squared (valid from [ntop_turb+1, nbot_turb] interfaces)
    ml2 = 30.0_kind_phys**2                 ! HB scheme:  length scale = 30m
    if(is_hbr_pbl_scheme) then
      ml2 = 1.0_kind_phys**2                ! HBR scheme: length scale = 1m
    end if

    ! Limit pbl height to regions below 400 mb
    ! npbl = max # of levels (from bottom) in pbl
    npbl = max(count(pref_mid(ntop_turb:nbot_turb) >= pblmaxp), 1)

    if(amIRoot) then
      write(iulog,*) 'Holtslag-Boville PBL: PBL height will be limited to bottom ', npbl, &
                     ' model levels. Top is ', pref_mid(pver+1-npbl), ' pascals'
      write(iulog,*) 'Holtslag-Boville PBL: top level of turbulence is ', ntop_turb, &
                     ' and bottom is ', nbot_turb
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

    ! Parameters:
    real(kind_phys), parameter :: minimum_velocity_shear_squared = 1.e-36_kind_phys

    ! Local variables:
    integer :: i, k
    real(kind_phys) :: rrho(ncol)                            ! 1 / bottom level density [m^3 kg-1]
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

    ! Compute s^2 (shear squared), n^2 (brunt vaisaila frequency), and ri (Richardson number = n^2/s^2) using virtual temperature
    ! (formerly trbintd)
    do k = ntop_turb,nbot_turb-1
      do i = 1,ncol
        s2(i,k) = calc_shear_squared(u(i,k), u(i,k+1), &
                                     v(i,k), v(i,k+1), &
                                     z(i,k), z(i,k+1), minimum_velocity_shear_squared)
        ri(i,k) = calc_bulk_richardson_number(thv(i,k), thv(i,k+1), &
                                              z(i,k),   z(i,k+1),   &
                                              s2(i,k),  gravit)
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
    real(kind_phys) :: scaled_phiminv  ! inverse phi function for momentum
    real(kind_phys) :: rino(ncol,pver) ! bulk Richardson no. from level to ref lev
    real(kind_phys) :: tlv(ncol)       ! ref. level potential tmp + tmp excess

    logical  :: check(ncol)            ! false if Richardson number > critical

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
             rino(i,k) = calc_richardson_number_at_height(u(i,k), u(i,pver), &
                                                          v(i,k), v(i,pver), &
                                                          z(i,k), z(i,pver), &
                                                          thv(i,k), thv(i,pver), ustar(i), gravit)
             ! Modified for boundary layer height diagnosis: Bert Holtslag, June 1994
             ! >>>>>>>>>  (Use ricr = 0.3 in this formulation)
             if (rino(i,k) >= ricr) then
                pblh(i) = linear_interpolate_height_wrt_richardson(z(i,k:k+1), rino(i,k:k+1))
                check(i) = .false.
             end if
          end if
       end do
    end do

    ! Estimate an effective surface temperature to account for surface fluctuations
    do i=1,ncol
       if (check(i)) pblh(i) = z(i,pverp-npbl)
       check(i)  = (kbfs(i) > 0._kind_phys)
       if (check(i)) then
          scaled_phiminv = comp_unstable_scaled_phiminv(pblh(i), obklen(i))
          rino(i,pver)   = 0.0_kind_phys
          tlv(i)         = compute_appropriate_temperature_at_z(thv(i,pver), kbfs(i), ustar(i), scaled_phiminv)
       end if
    end do

    ! Improve pblh estimate for unstable conditions using the convective temperature excess:
    do i = 1,ncol
      bge(i) = 1.e-8_kind_phys
    end do

    do k=pver-1,pver-npbl+1,-1
      do i=1,ncol
        if (check(i)) then
          rino(i,k) = calc_modified_richardson_number_at_height(u(i,k), u(i,pver), &
                                                                v(i,k), v(i,pver), &
                                                                z(i,k), z(i,pver), &
                                                                thv(i,k), thv(i,pver), tlv(i), ustar(i), gravit)

          if (rino(i,k) >= ricr) then
            pblh(i) = linear_interpolate_height_wrt_richardson(z(i,k:k+1), rino(i,k:k+1))

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
       pblh(i)  = max(pblh(i),700.0_kind_phys*ustar(i))
       wstar(i) = comp_wstar( max(0._kind_phys, kbfs(i)), pblh(i), thv(i,pver), gravit)
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
       if (cldn(i,pver) >= 0.0_kind_phys) then
         pblh(i) = max(pblh(i),zi(i,pver) + 50._kind_phys)
       end if
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

    real(kind_phys), parameter :: minimum_eddy_flux_coefficient = 0.01_kind_phys ! CCM1 2.f.14
    real(kind_phys), parameter :: a                             =  7.2_kind_phys ! Constant in turbulent prandtl number (CCM2, page 76)
    real(kind_phys), parameter :: b                             =  8.5_kind_phys ! Constant in surface temperature excess (CCM2, page 76)

    ! Local variables
    integer :: i, k
    real(kind_phys) :: kvf(ncol,pverp)      ! free atmospheric eddy diffusivity [m^2 s-1]
    integer         :: ktopbl(ncol)         ! index of first midpoint inside PBL (diagnostic) [index]
    real(kind_phys) :: scaled_phiminv(ncol) ! inverse phi function for momentum
    real(kind_phys) :: scaled_phihinv(ncol) ! inverse phi function for heat
    real(kind_phys) :: wm(ncol)             ! turbulent velocity scale for momentum
    real(kind_phys) :: fak1(ncol)           ! k*ustar*pblh
    real(kind_phys) :: fak2(ncol)           ! k*wm*pblh
    real(kind_phys) :: prandtl_factor(ncol) ! a*wstar/wm (Part of CCM2, page 76, equation 4.e.33)
    real(kind_phys) :: pblk                 ! level eddy diffusivity for momentum
    real(kind_phys) :: pr                   ! Prandtl number for eddy diffusivities
    real(kind_phys) :: zh                   ! zmzp / pblh
    real(kind_phys) :: zzh                  ! (1-(zmzp/pblh))**2
    real(kind_phys) :: zmzp                 ! level height halfway between zm and zp
    real(kind_phys) :: phiminv              ! intermediate calculation
    real(kind_phys) :: kve                  ! diffusivity at entrainment layer in unstable cases [m s-2]

    logical         :: unstable(ncol)       ! points with unstable pbl (positive virtual heat flux) [count]

    errmsg = ''
    errflg = 0

    ! Get atmosphere exchange coefficients
    kvf(:ncol,:) = 0.0_kind_phys
    do k = ntop_turb+1, nbot_turb-1
       do i = 1, ncol
          kvf(i,k+1) = max(calc_eddy_flux_coefficient(ml2, ri(i, k), s2(i, k)), &
                           minimum_eddy_flux_coefficient)
       end do
    end do
    kvf(:ncol,ntop_turb+1) = minimum_eddy_flux_coefficient

    ! Compute PBL exchange coefficients (formerly austausch_pbl)
    kvm = 0._kind_phys
    kvh = 0._kind_phys
    kve = 0._kind_phys
    cgh = 0._kind_phys
    cgs = 0._kind_phys
    tpert = 0._kind_phys
    qpert = 0._kind_phys
    tke = 0._kind_phys

    ! Mark columns with positive virtual heat flux:
    unstable = (kbfs > 0._kind_phys)

    ! Initialize height independent arrays
    do i=1,ncol
      fak1(i) = ustar(i)*pblh(i)*karman
      if (unstable(i)) then
        scaled_phiminv(i) = comp_unstable_scaled_phiminv(pblh(i), obklen(i))
        scaled_phihinv(i) = comp_unstable_scaled_phihinv(pblh(i), obklen(i))
        wm(i)             = ustar(i)*scaled_phiminv(i)
        fak2(i)           = wm(i)*pblh(i)*karman
        prandtl_factor(i) = a*wstar(i)/wm(i)
        tpert(i)          = max(khfs(i)*b/wm(i),    0._kind_phys)
        qpert(i)          = max(kqfs(i)*b/wm(i),    0._kind_phys)
      else
        tpert(i)          = max(khfs(i)*b/ustar(i), 0._kind_phys)
        qpert(i)          = max(kqfs(i)*b/ustar(i), 0._kind_phys)
      end if
    end do

    ! Initialize output arrays with free atmosphere values
    do k=1, pverp
       do i=1, ncol
          kvm(i,k) = kvf(i,k)
          kvh(i,k) = kvf(i,k)
       end do
    end do

    ! Main level loop to compute the diffusivities and counter-gradient terms. These terms are
    ! only calculated at points determined to be in the interior of the pbl (z(i,k) < pblh(i))
    ! and then calculations are directed toward regime: stable vs unstable, surface vs outer
    ! layer.
    do k = pver, pver-npbl+2, -1
      do i = 1, ncol
        if (z(i,k) < pblh(i)) then
          zmzp    = 0.5_kind_phys*(z(i,k) + z(i,k-1)) ! we think this is an approximation to the interface height (where KVs are calculated)
          zh      = zmzp/pblh(i)
          zzh     = zh*max(0._kind_phys,(1._kind_phys - zh))**2

          ktopbl(i) = k

          if (unstable(i)) then
            if (zh < sffrac) then
              phiminv = comp_unstable_phiminv(zmzp, obklen(i))
              pblk    = fak1(i)*zzh*phiminv
              pr      = phiminv/comp_unstable_phihinv(zmzp, obklen(i))
            else
              pblk     = fak2(i)*zzh
              pr       = scaled_phiminv(i)/scaled_phihinv(i) + ccon*prandtl_factor(i)/fak
              cgs(i,k) = prandtl_factor(i)/(pblh(i)*wm(i))
              cgh(i,k) = khfs(i)*cgs(i,k)*cpair
            end if
          else
            pblk = (ustar(i)*pblh(i)*karman)*zzh/comp_stable_phih(zmzp, obklen(i))
            pr   = 1._kind_phys
          end if
          kvm(i,k) = max(pblk, kvf(i,k))
          kvh(i,k) = max(pblk/pr, kvf(i,k))
        end if
      end do
    end do

    ! HBR scheme only:
    ! Check whether last allowed midpoint is within PBL
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
    do k = ntop_turb+1, nbot_turb-1
       do i = 1, ncol
          kvf(i,k+1) = calc_free_atm_eddy_flux_coefficient(ml2, ri(i, k), s2(i, k))
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

  end subroutine holtslag_boville_diff_finalize

  ! Utility pure elemental functions used in computation
  pure elemental function calc_bulk_richardson_number(thv1, thv2, z1, z2, s2, g) result(ri)
    real(kind_phys), intent(in) :: thv1
    real(kind_phys), intent(in) :: thv2
    real(kind_phys), intent(in) :: z1
    real(kind_phys), intent(in) :: z2
    real(kind_phys), intent(in) :: s2
    real(kind_phys), intent(in) :: g
    real(kind_phys)             :: ri
    real(kind_phys)             :: n2

    n2 = calc_brunt_vaisaila_frequency(thv1, thv2, z1, z2, g)
    ri = n2/s2
  end function calc_bulk_richardson_number

  pure elemental function calc_brunt_vaisaila_frequency(thv1, thv2, z1, z2, g) result(n2)
    real(kind_phys), intent(in) :: thv1
    real(kind_phys), intent(in) :: thv2
    real(kind_phys), intent(in) :: z1
    real(kind_phys), intent(in) :: z2
    real(kind_phys), intent(in) :: g
    real(kind_phys) :: n2

    n2 = g*2.0_kind_phys*(thv1-thv2)/((thv1+thv2)*(z1-z2))
  end function calc_brunt_vaisaila_frequency

  pure elemental function calc_shear_squared(u1, u2, v1, v2, z1, z2, minimum_velocity_shear_squared) result(s2)
    ! Hack, J., Boville, B., Briegleb, B., Kiehl, J., & Williamson, D. (1993).
    ! Description of the NCAR Community Climate Model (CCM2). University Corporation for Atmospheric Research.
    ! https://doi.org/10.5065/D6QZ27XV (Original work published 1993)
    ! Page 71, Equation 4.e.10
    real(kind_phys), intent(in) :: u1
    real(kind_phys), intent(in) :: u2
    real(kind_phys), intent(in) :: v1
    real(kind_phys), intent(in) :: v2
    real(kind_phys), intent(in) :: z1
    real(kind_phys), intent(in) :: z2
    real(kind_phys), intent(in) :: minimum_velocity_shear_squared
    real(kind_phys) :: s2

    s2 = max((u1-u2)**2 + (v1-v2)**2, minimum_velocity_shear_squared)/((z1-z2)**2)
  end function calc_shear_squared

  pure elemental function comp_wstar(kbfs, pblh, thv, gravity) result(wstar)
    ! Hack, J., Boville, B., Briegleb, B., Kiehl, J., & Williamson, D. (1993).
    ! Description of the NCAR Community Climate Model (CCM2). University Corporation for Atmospheric Research.
    ! https://doi.org/10.5065/D6QZ27XV (Original work published 1993)
    ! Page 75, Equation 4.e.32
    real(kind_phys), intent(in)  :: kbfs
    real(kind_phys), intent(in)  :: pblh
    real(kind_phys), intent(in)  :: thv
    real(kind_phys), intent(in)  :: gravity
    real(kind_phys)              :: wstar

    wstar = (kbfs*gravity*pblh/thv)**(1._kind_phys/3._kind_phys)
  end function comp_wstar

  pure elemental function calc_modified_richardson_number_at_height(uh, us, vh, vs, h, zs, thvh, thvs, tlv, ustar, g) result(rih)
    ! Modified version of base richardson number equation that incorporates reference level potential temperature
    real(kind_phys), intent(in) :: uh, us, vh, vs, h, zs, thvh, thvs, tlv, ustar, g
    real(kind_phys) :: rih
    real(kind_phys), parameter   :: b  = 100._kind_phys     ! Derivation of b comes from page 253
    real(kind_phys), parameter   :: tiny = 1.e-36_kind_phys ! Prevents division by 0

    rih = g * (thvh-tlv)*(h-zs) / &
          (thvs * max((uh-us)**2 + (vh-vs)**2 + b*ustar**2, tiny))
  end function calc_modified_richardson_number_at_height

  pure elemental function calc_richardson_number_at_height(uh, us, vh, vs, h, zs, thvh, thvs, ustar, g) result(rih)
    ! Vogelezang, D.H.P., Holtslag, A.A.M. Evaluation and model impacts of alternative
    ! boundary-layer height formulations. Boundary-Layer Meteorol 81, 245–269 (1996).
    ! https://doi.org/10.1007/BF02430331
    ! Equation 3, page 251
    real(kind_phys), intent(in) :: uh, us, vh, vs, h, zs, thvh, thvs, ustar, g
    real(kind_phys) :: rih
    real(kind_phys), parameter   :: b  = 100._kind_phys     ! Derivation of b comes from page 253
    real(kind_phys), parameter   :: tiny = 1.e-36_kind_phys ! Prevents division by 0

    rih = g * (thvh-thvs)*(h-zs) / &
          (thvs * max((uh-us)**2 + (vh-vs)**2 + b*ustar**2, tiny))
  end function calc_richardson_number_at_height

  pure function linear_interpolate_height_wrt_richardson(z, rino) result(pblh)
    real(kind_phys), intent(in) :: z(2)
    real(kind_phys), intent(in) :: rino(2)
    real(kind_phys)             :: pblh

    pblh = z(2) + (ricr - rino(2))/(rino(1) - rino(2)) * (z(1) - z(2))
  end function linear_interpolate_height_wrt_richardson

  pure elemental function compute_appropriate_temperature_at_z(thv, kbfs, ustar, phiminv) result(tlv)
    ! Holtslag, A.A.M., Van Meijgaard, E. & De Rooy, W.C. A comparison of boundary layer diffusion
    ! schemes in unstable conditions over land. Boundary-Layer Meteorol 76, 69–95 (1995).
    ! https://doi.org/10.1007/BF00710891
    ! Equation 10, page 72
    real(kind_phys), intent(in) :: thv
    real(kind_phys), intent(in) :: kbfs
    real(kind_phys), intent(in) :: ustar
    real(kind_phys), intent(in) :: phiminv
    real(kind_phys)             :: tlv
    real(kind_phys), parameter  :: b = 8.5_kind_phys

    tlv = thv + kbfs * b /(ustar * phiminv)
  end function compute_appropriate_temperature_at_z

  pure elemental function comp_stable_phih(z, L) result(phih)
    ! Hack, J., Boville, B., Briegleb, B., Kiehl, J., & Williamson, D. (1993).
    ! Description of the NCAR Community Climate Model (CCM2). University Corporation for Atmospheric Research.
    ! https://doi.org/10.5065/D6QZ27XV (Original work published 1993)
    ! Page 75, Equation 4.e.23 and 4.e.25
    real(kind_phys), intent(in) :: z ! height above surface [m]
    real(kind_phys), intent(in) :: L ! Obukhov length [m]
    real(kind_phys) :: zL
    real(kind_phys) :: phih          ! Vertical temperature gradient [1]

    zL = z/L
    if (zL <= 1._kind_phys) then
      phih = 1._kind_phys + 5._kind_phys*zL
    else
      phih = 5._kind_phys + zL
    endif
  end function comp_stable_phih

  pure elemental function comp_unstable_phiminv(z, L) result(phiminv)
    ! Hack, J., Boville, B., Briegleb, B., Kiehl, J., & Williamson, D. (1993).
    ! Description of the NCAR Community Climate Model (CCM2). University Corporation for Atmospheric Research.
    ! https://doi.org/10.5065/D6QZ27XV (Original work published 1993)
    ! Equation 4.e.28, page 75
    real(kind_phys), intent(in) :: z       ! height above surface [m]
    real(kind_phys), intent(in) :: L       ! Obukhov length [m]
    real(kind_phys)             :: phiminv ! Wind gradient [1]

    phiminv = (1._kind_phys - 15.0_kind_phys*(z/L))**(1._kind_phys/3._kind_phys)
  end function comp_unstable_phiminv

  pure elemental function comp_unstable_scaled_phiminv(z, L) result(phiminv)
    ! Modified version of:
    ! Hack, J., Boville, B., Briegleb, B., Kiehl, J., & Williamson, D. (1993).
    ! Description of the NCAR Community Climate Model (CCM2). University Corporation for Atmospheric Research.
    ! https://doi.org/10.5065/D6QZ27XV (Original work published 1993)
    ! Equation 4.e.28, page 75
    ! Scales z/L by surface layer fraction of boundary layer.
    ! Uses power of 1/3 instead of -1/3 for optimized usage.
    real(kind_phys), intent(in) :: z       ! height above surface [m]
    real(kind_phys), intent(in) :: L       ! Obukhov length [m]
    real(kind_phys)             :: phiminv ! Wind gradient [1]

    phiminv = (1._kind_phys - (15.0_kind_phys*sffrac)*z/L)**(1._kind_phys/3._kind_phys)
  end function comp_unstable_scaled_phiminv

  pure elemental function comp_unstable_phihinv(z, L) result(phiminv)
    ! Hack, J., Boville, B., Briegleb, B., Kiehl, J., & Williamson, D. (1993).
    ! Description of the NCAR Community Climate Model (CCM2). University Corporation for Atmospheric Research.
    ! https://doi.org/10.5065/D6QZ27XV (Original work published 1993)
    ! Page 75, Equation 4.e.26
    real(kind_phys), intent(in) :: z ! height above surface [m]
    real(kind_phys), intent(in) :: L ! Obukhov length [m]
    real(kind_phys)             :: phiminv

    phiminv = sqrt(1._kind_phys - 15.0_kind_phys*(z/L))
  end function comp_unstable_phihinv

  pure elemental function comp_unstable_scaled_phihinv(z, L) result(scaled_phihinv)
    ! Modified version of:
    ! Hack, J., Boville, B., Briegleb, B., Kiehl, J., & Williamson, D. (1993).
    ! Description of the NCAR Community Climate Model (CCM2). University Corporation for Atmospheric Research.
    ! https://doi.org/10.5065/D6QZ27XV (Original work published 1993)
    ! Page 75, Equation 4.e.26
    ! Scales z/L by surface layer fraction of boundary layer.
    real(kind_phys), intent(in) :: z ! height above surface [m]
    real(kind_phys), intent(in) :: L ! Obukhov length [m]
    real(kind_phys)             :: scaled_phihinv

    scaled_phihinv = sqrt(1._kind_phys - (15.0_kind_phys*sffrac)*z/L)
  end function comp_unstable_scaled_phihinv

end module holtslag_boville_diff
