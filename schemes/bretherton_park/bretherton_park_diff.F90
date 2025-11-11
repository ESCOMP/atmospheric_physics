! University of Washington (UW) moist turbulence scheme (also known as diag_TKE)
! This is the CCPP interface to this scheme.
!
! Original references:
! A new moist turbulence parametrization in the Community Atmosphere Model
!   by Christopher S. Bretherton and Sungsu Park. J. Climate. 2009. 22. 3422-3448
!   https://doi.org/10.1175/2008JCLI2556.1
! The University of Washington shallow convection and moist turbulence schemes
! and their impact on climate simulations with the Community Atmosphere Model
!    by Sungsu Park and Christopher S. Bretherton. J. Climate. 2009. 22. 3449-3469
!    https://doi.org/10.1175/2008JCLI2557.1
!
! Based on eddy_diff_cam originally from Sungsu Park.
module bretherton_park_diff
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! public CCPP-compliant interfaces
  public :: bretherton_park_diff_init
  public :: bretherton_park_diff_run

  ! Tunable parameters:
  ! Number of iterations for solution
  integer,         parameter :: nturb = 5

  real(kind_phys), parameter :: ml2 = 30.0_kind_phys**2 ! mixing lengths squared for computing free air diffusivity
                                                        ! used for interfaces in [ntop_eddy+1, nbot_eddy]

  logical,         parameter :: use_kvf = .false.       ! .true. (.false.) : initialize kvh/kvm = kvf (0)
  real(kind_phys), parameter :: lambda  = 0.5_kind_phys ! Under-relaxation factor (0 < lambda =< 1)

  ! Module storage:
  logical         :: do_diffusion_const_wet(1)
  integer         :: ntop_eddy
  integer         :: nbot_eddy

contains

!> \section arg_table_bretherton_park_diff_init Argument Table
!! \htmlinclude bretherton_park_diff_init.html
  subroutine bretherton_park_diff_init( &
             amIRoot, iulog, &
             pver, pverp, &
             gravit, cpair, rair, zvir, latvap, latice, karman, &
             ntop_eddy_in, &
             pref_mid, &
             ! namelist options:
             eddy_lbulk_max, eddy_leng_max, eddy_max_bot_pressure, eddy_moist_entrain_a2l, &
             ! output:
             ncvmax, &
             errmsg, errflg)

    ! Driver subroutines for UW PBL scheme.
    use eddy_diff, only: init_eddy_diff

    ! Input arguments
    logical,         intent(in)  :: amIRoot
    integer,         intent(in)  :: iulog
    integer,         intent(in)  :: pver
    integer,         intent(in)  :: pverp
    real(kind_phys), intent(in)  :: gravit
    real(kind_phys), intent(in)  :: cpair
    real(kind_phys), intent(in)  :: rair
    real(kind_phys), intent(in)  :: zvir
    real(kind_phys), intent(in)  :: latvap
    real(kind_phys), intent(in)  :: latice
    real(kind_phys), intent(in)  :: karman
    integer,         intent(in)  :: ntop_eddy_in           ! Top interface level to which eddy vertical diffusivity is applied [index]
    real(kind_phys), intent(in)  :: pref_mid(:)            ! reference_pressure [Pa]
    real(kind_phys), intent(in)  :: eddy_lbulk_max         ! Maximum master length [m]
    real(kind_phys), intent(in)  :: eddy_leng_max          ! Maximum dissipation length [m]
    real(kind_phys), intent(in)  :: eddy_max_bot_pressure  ! Bottom pressure level for eddy_leng_max [hPa]
    real(kind_phys), intent(in)  :: eddy_moist_entrain_a2l ! Moist entrainment enhancement parameter [1]

    ! Output arguments
    integer,         intent(out) :: ncvmax                 ! maximum number of convective layers (CLs) [count]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: k
    real(kind_phys) :: leng_max(pver)

    errmsg = ''
    errflg = 0

    ! Only diffuse water vapor
    do_diffusion_const_wet(1) = .true.

    ! from the set_vertical_diffusion_top scheme.
    ntop_eddy = ntop_eddy_in

    ! hardcoded to be surface layer:
    nbot_eddy = pver

    ! Set maximum number of convective layers (CLs) in the UW scheme
    ! this was hardcoded to pver at the time of CCPPization, but could be made a parameter if needed.
    ! this is used as a dimension for many of the CL diagnostic variables in the run phase.
    ncvmax = pver

    do k = 1, pver
      if (pref_mid(k) <= eddy_max_bot_pressure*1.e2_kind_phys) then ! hPa to Pa
        leng_max(k) = eddy_leng_max
      else
        leng_max(k) = 40.e3_kind_phys
      end if
    end do

    if(amIRoot) then
      write(iulog,*) 'eddy diffusivity scheme is UW Moist Turbulence Scheme by Bretherton and Park'
      write(iulog,*) 'Bretherton and Park (UW) PBL: nturb = ', nturb
      write(iulog,*) 'Bretherton and Park (UW) PBL: eddy_leng_max = ', eddy_leng_max,' lbulk_max = ', eddy_lbulk_max
      do k = 1,pver
        write(iulog,*) 'Bretherton and Park (UW) PBL: ',k,pref_mid(k),' leng_max = ',leng_max(k)
      end do
    end if

    call init_eddy_diff(pver, ncvmax, gravit, cpair, rair, zvir, &
         latvap, latice, ntop_eddy, nbot_eddy, karman, &
         eddy_lbulk_max, leng_max, &
         eddy_moist_entrain_a2l, errmsg)

    if(trim(errmsg) /= "") then
      errflg = 1
      return
    end if

  end subroutine bretherton_park_diff_init

  ! Interface to UW PBL scheme to compute eddy diffusivities.
  ! Eddy diffusivities are calculated in a fully implicit way through iterative process.
  ! CL = convective layers; STL = stable turbulent layers.
  !
  ! Original author: Sungsu Park, August 2006, May 2008.
!> \section arg_table_bretherton_park_diff_run Argument Table
!! \htmlinclude bretherton_park_diff_run.html
  subroutine bretherton_park_diff_run( &
             ncol, pver, pverp, pcnst, ncvmax, &
             iulog, &
             dt, &
             const_props, &
             do_iss, am_correction, do_beljaars, &
             is_first_timestep, &
             gravit, cpair, rair, latvap, latice, &
             t, tint, &
             qv, ql, qi, &
             s, &
             p, rhoi, dpidz_sq, &
             cldn, &
             z, zi, &
             pmid, pint, &
             u, v, &
             taux, tauy, &
             shflx, qflx, wstarent, &
             ksrftms, dragblj, &
             qrl, wsedl, &
             ! input/output:
             kvm, kvh, &
             tauresx, tauresy, &
             ! output:
             s2, n2, ri, &
             kvq, rrho, &
             ustar, &
             pblh, pblhp, minpblh, &
             cgh, cgs, tpert, qpert, wpert, &
             tke, tkes, wcap, &
             wsed, turbtype, &
             bprod, sprod, sfi, sfuh, sflh, &
             qlfd, &
             chu, chs, cmu, cms, &
             kbase_o, ktop_o, ncvfin_o, &
             kbase_mg, ktop_mg, ncvfin_mg, &
             kbase_f, ktop_f, ncvfin_f, &
             wet, web, jtbu, jbbu, &
             evhc, jt2slv, n2ht, n2hb, &
             lwp, opt_depth, radinvfrac, radf, &
             wstar, wstar3fact, &
             ebrk, wbrk, lbrk, &
             ricl, ghcl, shcl, smcl, &
             ghi, shi, smi, rii, lengi, &
             errorPBL, &
             errmsg, errflg)

    ! pbl utils dependencies:
    use atmos_phys_pbl_utils, only: calc_eddy_flux_coefficient, calc_ideal_gas_rrho, calc_friction_velocity

    ! to-be-ccppized dependency:
    use coords_1d, only: coords1d
    use wv_saturation, only: qsat

    ! Driver routines for UW PBL scheme.
    use eddy_diff, only: trbintd, caleddy

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! dependency to get constituent index
    use ccpp_const_utils,          only: ccpp_const_get_idx

    ! Driver routines for diffusion solver to be called iteratively within this run phase.
    use diffusion_solver, only: implicit_surface_stress_add_drag_coefficient_run
    use diffusion_stubs,  only: turbulent_mountain_stress_add_drag_coefficient_run
    use diffusion_solver, only: vertical_diffusion_wind_damping_rate_run
    use diffusion_solver, only: vertical_diffusion_diffuse_horizontal_momentum_run
    use diffusion_solver, only: vertical_diffusion_diffuse_dry_static_energy_run
    use diffusion_solver, only: vertical_diffusion_diffuse_tracers_run

    ! Input variables
    integer,         intent(in)    :: ncol
    integer,         intent(in)    :: pver
    integer,         intent(in)    :: pverp
    integer,         intent(in)    :: pcnst
    integer,         intent(in)    :: ncvmax              ! max # of CLs (can set to pver) [count]
    integer,         intent(in)    :: iulog
    real(kind_phys), intent(in)    :: dt                  ! Physics timestep [s]
    type(ccpp_constituent_prop_ptr_t), &
                     intent(in)    :: const_props(:)      ! CCPP constituent properties pointer
    logical,         intent(in)    :: do_iss              ! Use implicit turbulent surface stress computation [flag]
    logical,         intent(in)    :: am_correction       ! Do angular momentum conservation correction [flag]
    logical,         intent(in)    :: do_beljaars         ! Do Beljaars drag in vertical diffusion? [flag]
    logical,         intent(in)    :: is_first_timestep   ! is_first_timestep [flag]
    real(kind_phys), intent(in)    :: gravit
    real(kind_phys), intent(in)    :: cpair
    real(kind_phys), intent(in)    :: rair
    real(kind_phys), intent(in)    :: latvap
    real(kind_phys), intent(in)    :: latice
    real(kind_phys), intent(in)    :: t(:, :)             ! Temperature [K]
    real(kind_phys), intent(in)    :: tint(:, :)          ! Temperature defined on interfaces [K]
    real(kind_phys), intent(in)    :: qv(:, :)            ! Water vapor mixing ratio [kg kg-1]
    real(kind_phys), intent(in)    :: ql(:, :)            ! Liquid water mixing ratio [kg kg-1]
    real(kind_phys), intent(in)    :: qi(:, :)            ! Ice mixing ratio [kg kg-1]
    real(kind_phys), intent(in)    :: s(:, :)             ! Dry static energy [J kg-1]
    type(coords1d),  intent(in)    :: p                   ! Pressure coordinates for solver [Pa]
    real(kind_phys), intent(in)    :: rhoi(:, :)          ! Density at interfaces [kg m-3]
    real(kind_phys), intent(in)    :: dpidz_sq(:, :)      ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4], interfaces
    real(kind_phys), intent(in)    :: cldn(:, :)          ! Stratiform cloud fraction [fraction]
    real(kind_phys), intent(in)    :: z(:, :)             ! Layer mid-point height above surface [m]
    real(kind_phys), intent(in)    :: zi(:, :)            ! Interface height above surface [m]
    real(kind_phys), intent(in)    :: pmid(:, :)          ! Layer mid-point pressure [Pa]
    real(kind_phys), intent(in)    :: pint(:, :)          ! Interface pressure [Pa]
    real(kind_phys), intent(in)    :: u(:, :)             ! Zonal velocity [m s-1]
    real(kind_phys), intent(in)    :: v(:, :)             ! Meridional velocity [m s-1]
    real(kind_phys), intent(in)    :: taux(:)             ! Zonal wind stress at surface [N m-2]
    real(kind_phys), intent(in)    :: tauy(:)             ! Meridional wind stress at surface [N m-2]
    real(kind_phys), intent(in)    :: shflx(:)            ! Sensible heat flux at surface [W m-2]
    real(kind_phys), intent(in)    :: qflx(:, :)          ! Constituent fluxes at surface [kg m-2 s-1]
    logical,         intent(in)    :: wstarent            ! .true. means use the 'wstar' entrainment closure [flag]
    real(kind_phys), intent(in)    :: ksrftms(:)          ! Surface drag coefficient of turbulent mountain stress [kg m-2 s-1]
    real(kind_phys), intent(in)    :: dragblj(:, :)       ! Drag profile from Beljaars SGO form drag [s-1]
    real(kind_phys), intent(in)    :: qrl(:, :)           ! LW radiative cooling rate [K s-1]
    real(kind_phys), intent(in)    :: wsedl(:, :)         ! Sedimentation velocity of stratiform liquid cloud droplet [m s-1]

    ! Input/output variables
    real(kind_phys), intent(inout) :: kvm(:, :)           ! Eddy diffusivity for momentum [m2 s-1]
                                                          ! (in from previous timestep, output to current timestep)
    real(kind_phys), intent(inout) :: kvh(:, :)           ! Eddy diffusivity for heat [m2 s-1]
                                                          ! (in from previous timestep, output to current timestep)
    real(kind_phys), intent(inout) :: tauresx(:)          ! Residual stress to be added in vdiff to correct for turb [N m-2]
    real(kind_phys), intent(inout) :: tauresy(:)          ! Stress mismatch between sfc and atm accumulated in prior timesteps [N m-2]

    ! Output variables, including those for diagnostic output
    real(kind_phys), intent(out)   :: s2(:, :)            ! Shear squared, defined at interfaces except surface [s-2]
    real(kind_phys), intent(out)   :: n2(:, :)            ! Buoyancy frequency, defined at interfaces except surface [s-2]
    real(kind_phys), intent(out)   :: ri(:, :)            ! Richardson number, 'n2/s2', defined at interfaces except surface [s-2]
    real(kind_phys), intent(out)   :: kvq(:, :)           ! Eddy diffusivity for constituents, moisture and tracers [m2 s-1]
    real(kind_phys), intent(out)   :: rrho(:)             ! Reciprocal of density at the lowest layer [m3 kg-1]
    real(kind_phys), intent(out)   :: ustar(:)            ! Surface friction velocity [m s-1]
    real(kind_phys), intent(out)   :: pblh(:)             ! PBL top height [m]
    real(kind_phys), intent(out)   :: pblhp(:)            ! PBL top pressure [Pa]
    real(kind_phys), intent(out)   :: minpblh(:)          ! Minimum PBL height based on surface stress [m]
    real(kind_phys), intent(out)   :: cgh(:, :)           ! Counter-gradient term for heat [J kg-1 m-1]
    real(kind_phys), intent(out)   :: cgs(:, :)           ! Counter-gradient star [cg flux-1]
    real(kind_phys), intent(out)   :: tpert(:)            ! Convective temperature excess [K]
    real(kind_phys), intent(out)   :: qpert(:)            ! Convective humidity excess [kg kg-1]
    real(kind_phys), intent(out)   :: wpert(:)            ! Turbulent velocity excess [m s-1]
    real(kind_phys), intent(out)   :: tke(:, :)           ! Turbulent kinetic energy [m2 s-2]
    real(kind_phys), intent(out)   :: tkes(:)             ! TKE at surface interface [m2 s-2]
    real(kind_phys), intent(out)   :: wcap(:, :)          ! Normalized TKE at all interfaces [m2 s-2]
    real(kind_phys), intent(out)   :: wsed(:, :)          ! Sedimentation velocity at the top of each CL (for sedimentation-entrainment feedback) [m s-1]
    integer,         intent(out)   :: turbtype(:, :)      ! Turbulence type identifier at all interfaces [1]
    real(kind_phys), intent(out)   :: bprod(:, :)         ! Buoyancy production of tke (interfaces) [m2 s-3]
    real(kind_phys), intent(out)   :: sprod(:, :)         ! Shear production [m2 s-3]
    real(kind_phys), intent(out)   :: sfi(:, :)           ! Interfacial layer saturation fraction [fraction]
    real(kind_phys), intent(out)   :: sfuh(:, :)          ! Saturation fraction in upper half-layer [fraction]
    real(kind_phys), intent(out)   :: sflh(:, :)          ! Saturation fraction in lower half-layer [fraction]
    real(kind_phys), intent(out)   :: qlfd(:, :)          ! Liquid water specific humidity for diffusion [kg kg-1]
    ! Buoyancy coefficients : w'b' = ch * w'sl' + cm * w'qt'
    real(kind_phys), intent(out)   :: chu(:, :)           ! Heat buoyancy coef for dry states, interfaces [m s-2 kg J-1]
    real(kind_phys), intent(out)   :: chs(:, :)           ! Heat buoyancy coef for sat states, interfaces [m s-2 kg J-1]
    real(kind_phys), intent(out)   :: cmu(:, :)           ! Moisture buoyancy coef for dry states,  interfaces [m s-2 kg-1 kg]
    real(kind_phys), intent(out)   :: cms(:, :)           ! Moisture buoyancy coef for sat states, interfaces [m s-2 kg-1 kg]
    real(kind_phys), intent(out)   :: kbase_o(:, :)       ! Original external base interface index of CL from 'exacol', ncvmax [index]
    real(kind_phys), intent(out)   :: ktop_o(:, :)        ! Original external top  interface index of CL from 'exacol', ncvmax [index]
    real(kind_phys), intent(out)   :: ncvfin_o(:)         ! Original number of CLs from 'exacol' [count]
    real(kind_phys), intent(out)   :: kbase_mg(:, :)      ! 'kbase' after extending-merging from 'zisocl', ncvmax [index]
    real(kind_phys), intent(out)   :: ktop_mg(:, :)       ! 'ktop' after extending-merging from 'zisocl', ncvmax [index]
    real(kind_phys), intent(out)   :: ncvfin_mg(:)        ! 'ncvfin' after extending-merging from 'zisocl' [count]
    real(kind_phys), intent(out)   :: kbase_f(:, :)       ! Final 'kbase' after extending-merging & including SRCL, ncvmax [index]
    real(kind_phys), intent(out)   :: ktop_f(:, :)        ! Final 'ktop' after extending-merging & including SRCL, ncvmax [index]
    real(kind_phys), intent(out)   :: ncvfin_f(:)         ! Final 'ncvfin' after extending-merging & including SRCL [count]
    real(kind_phys), intent(out)   :: wet(:, :)           ! Entrainment rate at the CL top, ncvmax  [m s-1]
    real(kind_phys), intent(out)   :: web(:, :)           ! Entrainment rate at the CL base, ncvmax [m s-1] (Set to zero if CL is based at surface)
    real(kind_phys), intent(out)   :: jtbu(:, :)          ! Buoyancy jump across the CL top, ncvmax  [m s-2]
    real(kind_phys), intent(out)   :: jbbu(:, :)          ! Buoyancy jump across the CL base, ncvmax [m s-2]
    real(kind_phys), intent(out)   :: evhc(:, :)          ! Evaporative enhancement factor at the CL top, ncvmax [1]
    real(kind_phys), intent(out)   :: jt2slv(:, :)        ! Jump of slv (liquid water virtual static energy) (across two layers)
                                                          ! at CL top used only for evhc (evaporative enhancement factor at CL top), ncvmax [J kg-1]
    real(kind_phys), intent(out)   :: n2ht(:, :)          ! n2 defined at the CL top  interface but using
                                                          ! sfuh(kt)   instead of sfi(kt), ncvmax [s-2]
    real(kind_phys), intent(out)   :: n2hb(:, :)          ! n2 defined at the CL base interface but using
                                                          ! sflh(kb-1) instead of sfi(kb), ncvmax [s-2]
    real(kind_phys), intent(out)   :: lwp(:, :)           ! LWP in the CL top layer, ncvmax [kg m-2]
    real(kind_phys), intent(out)   :: opt_depth(:, :)     ! Optical depth of the CL top layer, ncvmax [1]
    real(kind_phys), intent(out)   :: radinvfrac(:, :)    ! Fraction of radiative cooling confined in the top portion of CL top layer, ncvmax [fraction]
    real(kind_phys), intent(out)   :: radf(:, :)          ! Buoyancy production at the CL top due to LW radiative cooling, ncvmax [m2 s-3]
    real(kind_phys), intent(out)   :: wstar(:, :)         ! Convective velocity in each CL, ncvmax [m s-1]
    real(kind_phys), intent(out)   :: wstar3fact(:, :)    ! Enhancement of 'wstar3' due to entrainment (inverse), ncvmax [1]
    real(kind_phys), intent(out)   :: ebrk(:, :)          ! Net mean TKE of CL including entrainment effect, ncvmax [m2 s-2]
    real(kind_phys), intent(out)   :: wbrk(:, :)          ! Net mean normalized TKE (W) of CL,
                                                          ! 'ebrk/b1' including entrainment effect, ncvmax [m2 s-2]
    real(kind_phys), intent(out)   :: lbrk(:, :)          ! Energetic internal thickness of CL, ncvmax [m]
    real(kind_phys), intent(out)   :: ricl(:, :)          ! CL internal mean Richardson number, ncvmax [1]
    real(kind_phys), intent(out)   :: ghcl(:, :)          ! Half of normalized buoyancy production of CL, ncvmax [1]
    real(kind_phys), intent(out)   :: shcl(:, :)          ! Galperin instability function of heat-moisture of CL, ncvmax [1]
    real(kind_phys), intent(out)   :: smcl(:, :)          ! Galperin instability function of mementum of CL, ncvmax [1]
    real(kind_phys), intent(out)   :: ghi(:, :)           ! Half of normalized buoyancy production at all interfaces [1]
    real(kind_phys), intent(out)   :: shi(:, :)           ! Galperin instability function of heat-moisture at all interfaces [1]
    real(kind_phys), intent(out)   :: smi(:, :)           ! Galperin instability function of heat-moisture at all interfaces [1]
    real(kind_phys), intent(out)   :: rii(:, :)           ! Interfacial Richardson number defined at all interfaces [1]
    real(kind_phys), intent(out)   :: lengi(:, :)         ! Turbulence length scale at all interfaces [m]
    real(kind_phys), intent(out)   :: errorPBL(:)         ! Error function showing whether PBL produced convergent solution or not [m2 s-1]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k, iturb
    integer :: const_wv_idx
    integer :: ipbl(ncol)     ! If 1, PBL is CL, while if 0, PBL is STL.
    integer :: kpblh(ncol)    ! Layer index containing PBL top within or at the base interface (NOT USED)

    character(2048) :: warnstring                ! Warning messages from driver routine

    real(kind_phys) :: went(ncol)                ! Entrainment rate at the PBL top interface [ m/s ] (NOT USED)

    real(kind_phys) :: kvf(ncol, pverp)          ! Free atmospheric eddy diffusivity [ m2/s ]
    real(kind_phys) :: kvm_in(ncol, pverp)       ! Eddy diffusivity for momentum [m2 s-1], previous timestep
    real(kind_phys) :: kvh_in(ncol, pverp)       ! Eddy diffusivity for heat [m2 s-1], previous timestep
    real(kind_phys) :: kvm_ce(ncol, pverp)       ! Eddy diffusivity for momentum [m2 s-1], input into caleddy
    real(kind_phys) :: kvh_ce(ncol, pverp)       ! Eddy diffusivity for heat [m2 s-1], input into caleddy

    real(kind_phys) :: qt(ncol, pver)            ! Total specific humidity [ kg/kg ]
    real(kind_phys) :: sl(ncol, pver)            ! Liquid water static energy [ J/kg ]
    real(kind_phys) :: slv(ncol, pver)           ! Liquid water virtual static energy [ J/kg ]
    real(kind_phys) :: slslope(ncol, pver)       ! Slope of 'sl' in each layer
    real(kind_phys) :: qtslope(ncol, pver)       ! Slope of 'qt' in each layer
    real(kind_phys) :: qvfd(ncol, pver)          ! Specific humidity for diffusion [ kg/kg ]
    real(kind_phys) :: tfd(ncol, pver)           ! Temperature for diffusion [ K ]
    real(kind_phys) :: slfd(ncol, pver)          ! Liquid static energy [ J/kg ]
    real(kind_phys) :: qtfd(ncol, pver, 1)       ! Total specific humidity [ kg/kg ]
    real(kind_phys) :: ufd(ncol, pver)           ! U-wind for diffusion [ m/s ]
    real(kind_phys) :: vfd(ncol, pver)           ! V-wind for diffusion [ m/s ]

    ! Output arguments from iterative calls of diffusion solver
    real(kind_phys) :: ufd_out(ncol, pver)       ! U-wind for diffusion [ m/s ]
    real(kind_phys) :: vfd_out(ncol, pver)       ! V-wind for diffusion [ m/s ]
    real(kind_phys) :: slfd_out(ncol, pver)      ! Liquid static energy [ J/kg ]
    real(kind_phys) :: qtfd_out(ncol, pver, 1)   ! Total specific humidity [ kg/kg ]

    ! Input arguments for CCPPized diffusion solver
    real(kind_phys) :: ksrf(ncol)
    real(kind_phys) :: tau_damp_rate(ncol, pver)
    real(kind_phys) :: ubc_mmr_dummy(ncol, 1)
    logical  :: cnst_fixed_ubc(1)

    real(kind_phys) :: jnk1d(ncol)
    real(kind_phys) :: jnk2d(ncol, pverp)
    real(kind_phys) :: zero(ncol)
    real(kind_phys) :: zero2d(ncol, pverp)
    real(kind_phys) :: es                         ! Saturation vapor pressure
    real(kind_phys) :: qs                         ! Saturation specific humidity
    real(kind_phys) :: ep2, templ, temps

    errmsg = ''
    errflg = 0

    ! Check constituents list and locate water vapor index
    ! (not assumed to be 1)
    call ccpp_const_get_idx(const_props, &
         'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', &
         const_wv_idx, errmsg, errflg)
    if (errflg /= 0) return

    ! Initialize dummy variables to pass into diffusion solver.
    zero(:) = 0._kind_phys
    zero2d(:, :) = 0._kind_phys

    ! Set initial state
    ufd(:ncol, :) = u(:ncol, :)
    vfd(:ncol, :) = v(:ncol, :)
    tfd(:ncol, :) = t(:ncol, :)
    qvfd(:ncol, :) = qv(:ncol, :)
    qlfd(:ncol, :) = ql(:ncol, :)

    ! Save kvm, kvh from previous timestep for use in the driver routines,
    ! as the kvm and kvh inout arguments are overwritten during the iterative
    ! process.
    kvm_in(:, :pverp) = kvm(:, :pverp)
    kvh_in(:, :pverp) = kvh(:, :pverp)

    ! Prepare invariant drag coefficients for diffusion solver
    ! during the iterative process using CCPPized subroutines.
    !
    ! Calculate surface drag rate
    ksrf(:ncol) = 0._kind_phys
    call implicit_surface_stress_add_drag_coefficient_run( &
      ncol   = ncol, &
      pver   = pver, &
      do_iss = do_iss, &
      taux   = taux(:ncol), &
      tauy   = tauy(:ncol), &
      u0     = u(:ncol, :pver), & ! use original state.
      v0     = v(:ncol, :pver), & ! use original state.
      ! below input/output:
      ksrf   = ksrf(:ncol), &
      errmsg = errmsg, &
      errflg = errflg)
    if (errflg /= 0) return

    ! Add TMS surface drag rate
    call turbulent_mountain_stress_add_drag_coefficient_run( &
      ncol = ncol, &
      pver = pver, &
      ksrftms = ksrftms(:ncol), &
      ! below input/output:
      ksrf = ksrf(:ncol), &
      errmsg = errmsg, &
      errflg = errflg)
    if (errflg /= 0) return

    ! Based on the drag coefficients, calculate wind damping rates
    call vertical_diffusion_wind_damping_rate_run( &
      ncol=ncol, &
      pver=pver, &
      gravit=gravit, &
      p=p, & ! Coords1D, pressure coordinates [Pa]
      ksrf=ksrf(:ncol), &
      ! below output:
      tau_damp_rate=tau_damp_rate(:ncol, :pver), &
      errmsg=errmsg, &
      errflg=errflg)
    if (errflg /= 0) return

    ! Iterative loop:
    iturb_loop: do iturb = 1, nturb

      ! Total stress includes 'tms'.
      ! Here, in computing 'tms', we can use either iteratively changed 'ufd,vfd' or the
      ! initially given 'u,v' to the PBL scheme. Note that normal stress, 'taux, tauy'
      ! are not changed by iteration. In order to treat 'tms' in a fully implicit way,
      ! updated wind is used here.

      ! Compute ustar
      rrho(:ncol) = calc_ideal_gas_rrho(rair, tfd(:ncol, pver), pmid(:ncol, pver))
      ustar(:ncol) = calc_friction_velocity(taux(:ncol) - ksrftms(:ncol)*ufd(:ncol, pver), & ! Zonal wind stress
                                            tauy(:ncol) - ksrftms(:ncol)*vfd(:ncol, pver), & ! Meridional wind stress
                                            rrho(:ncol))

      minpblh(:ncol) = 100.0_kind_phys*ustar(:ncol)   ! By construction, 'minpblh' is larger than 1 [m] when 'ustar_min = 0.01'.

      ! Calculate (qt,sl,n2,s2,ri) from a given set of (t,qv,ql,qi,u,v)
      call trbintd( &
        ncol, pver, z, ufd, vfd, tfd, pmid, &
        s2, n2, ri, zi, pint, cldn, qtfd, qvfd, &
        qlfd, qi, sfi, sfuh, sflh, slfd, slv, slslope, &
        qtslope, chs, chu, cms, cmu)

      ! Save initial (i.e., before iterative diffusion) profile of (qt,sl) at each iteration.
      ! Only necessary for (qt,sl) not (u,v) because (qt,sl) are newly calculated variables.
      if (iturb == 1) then
        qt(:ncol, :) = qtfd(:ncol, :, 1)
        sl(:ncol, :) = slfd(:ncol, :)
      end if

      ! Get free atmosphere exchange coefficients. This 'kvf' is not used in UW moist PBL scheme
      if (use_kvf) then
        kvf(:ncol, :) = 0.0_kind_phys
        do k = ntop_eddy, nbot_eddy - 1
          do i = 1, ncol
            ! mixing length squared (ml2) is fixed parameter value from ntop_eddy+1 to nbot_eddy, otherwise 0
            if(k >= ntop_eddy+1 .and. k <= nbot_eddy) then
              kvf(i, k + 1) = calc_eddy_flux_coefficient(ml2, ri(i, k), s2(i, k))
            else
              kvf(i, k + 1) = calc_eddy_flux_coefficient(0.0_kind_phys, ri(i, k), s2(i, k))
            end if
          end do
        end do
      else
        kvf = 0._kind_phys
      end if

      ! caleddy driver subroutine:
      !
      ! Calculate eddy diffusivity (kvh,kvm) and (tke,bprod,sprod) using
      ! a given (kvh,kvm) which are used only for initializing (bprod,sprod)  at
      ! the first part of caleddy. (bprod,sprod) are fully updated at the end of
      ! caleddy after calculating (kvh,kvm)
      !
      ! Depending on model timestep and iteration number, kvh and kvm input to caleddy
      ! as the initial guess of bprod, sprod
      ! differ as necessary for 'wstar-based' entrainment closure.
      if (iturb == 1 .and. is_first_timestep) then
        ! First iteration of first model timestep: Use free tropospheric value or zero.
        kvh_ce = kvf
        kvm_ce = kvf
      else
        ! Further iterations or non-first model timesteps: use value from
        ! previous iteration or previous timestep.
        kvh_ce = kvh
        kvm_ce = kvm
      end if

      call caleddy(ncol, pver, &
                   slfd, qtfd, qlfd, slv, ufd, &
                   vfd, pint, z, zi, &
                   qflx(:,const_wv_idx), shflx, slslope, qtslope, &
                   chu, chs, cmu, cms, sfuh, &
                   sflh, n2, s2, ri, rrho, &
                   pblh, ustar, &
                   kvh_ce, kvm_ce, kvh, kvm, &
                   tpert, qpert, qrl, kvf, tke, &
                   wstarent, bprod, sprod, minpblh, wpert, &
                   tkes, went, turbtype, &
                   kbase_o, ktop_o, ncvfin_o, &
                   kbase_mg, ktop_mg, ncvfin_mg, &
                   kbase_f, ktop_f, ncvfin_f, &
                   wet, web, jtbu, jbbu, &
                   evhc, jt2slv, n2ht, n2hb, &
                   lwp, opt_depth, radinvfrac, radf, &
                   wstar, wstar3fact, &
                   ebrk, wbrk, lbrk, ricl, ghcl, &
                   shcl, smcl, ghi, shi, smi, &
                   rii, lengi, wcap, pblhp, cldn, &
                   ipbl, kpblh, wsedl, wsed, &
                   warnstring, errmsg)

      if (trim(warnstring) /= "") then
        write (iulog, *) "eddy_diff_cam: Messages from caleddy follow."
        write (iulog, *) warnstring
      end if

      if (trim(errmsg) /= "") then
        errflg = 1
        return
      end if

      ! Calculate errorPBL to check whether PBL produced convergent solutions or not.
      if (iturb == nturb) then
        do i = 1, ncol
          errorPBL(i) = 0._kind_phys
          do k = 1, pver
            errorPBL(i) = errorPBL(i) + (kvh(i, k) - kvh_ce(i, k))**2
          end do
          errorPBL(i) = sqrt(errorPBL(i)/pver)
        end do
      end if

      ! Eddy diffusivities which will be used for the initialization of (bprod,
      ! sprod) in 'caleddy' at the next iteration step.
      !
      ! This is updated from the values from the output of caleddy
      ! and from the initial input to caleddy
      if (iturb > 1 .and. iturb < nturb) then
        kvm(:ncol, :) = lambda*kvm(:ncol, :) + (1._kind_phys - lambda)*kvm_ce(:ncol, :)
        kvh(:ncol, :) = lambda*kvh(:ncol, :) + (1._kind_phys - lambda)*kvh_ce(:ncol, :)
      end if

      ! Set nonlocal terms to zero for flux diagnostics, since not used by caleddy.
      cgh(:ncol, :) = 0._kind_phys
      cgs(:ncol, :) = 0._kind_phys

      if (iturb < nturb) then

        ! Each time we diffuse the original state
        slfd(:ncol, :) = sl(:ncol, :)
        qtfd(:ncol, :, 1) = qt(:ncol, :)
        ufd(:ncol, :) = u(:ncol, :)
        vfd(:ncol, :) = v(:ncol, :)

        ! Diffuse initial profile of each time step using a given (kvh_out,kvm_out)
        ! Call the CCPPized subroutines for the diffusion solver
        ! in iteration. This is not specified in the SDF but instead
        ! called internally because it is a iterative process.
        ! Notes:
        ! - there is no molecular diffusion used here.
        ! - in tracers, only water vapor is diffused (ncnst = 1)
        ufd_out(:, :) = 0._kind_phys
        vfd_out(:, :) = 0._kind_phys
        slfd_out(:, :) = 0._kind_phys
        call vertical_diffusion_diffuse_horizontal_momentum_run( &
          ncol=ncol, &
          pver=pver, &
          pverp=pverp, &
          dt=dt, &
          rair=rair, &
          gravit=gravit, &
          do_iss=do_iss, &
          am_correction=am_correction, &
          itaures=.false., &
          t=t(:ncol, :pver), &
          p=p, & ! Coords1D, pressure coordinates [Pa]
          rhoi=rhoi(:ncol, :pverp), &
          taux=taux(:ncol), &
          tauy=tauy(:ncol), &
          tau_damp_rate=tau_damp_rate(:ncol, :pver), & ! tau damp rate from above
          kvm=kvm(:ncol, :pverp), &
          ksrftms=ksrftms(:ncol), &
          dragblj=dragblj(:ncol, :pver), &
          dpidz_sq=dpidz_sq(:ncol, :pverp), & ! moist
          u0=ufd(:ncol, :pver), &
          v0=vfd(:ncol, :pver), &
          dse0=slfd(:ncol, :pver), &
          ! input/output
          ! (not actually updated since itaures = .false. in this internal call.)
          tauresx=tauresx(:ncol), &
          tauresy=tauresy(:ncol), &
          ! below output
          u1=ufd_out(:ncol, :pver), &
          v1=vfd_out(:ncol, :pver), &
          dse1=slfd_out(:ncol, :pver), &
          dtk=jnk2d(:ncol, :), & ! unused dummy.
          tautmsx=jnk1d(:ncol), & ! unused dummy.
          tautmsy=jnk1d(:ncol), & ! unused dummy.
          ! arguments for Beljaars
          do_beljaars=do_beljaars, &
          errmsg=errmsg, &
          errflg=errflg)

        if (errflg /= 0) return

        ! Update u, v, dse with updated iterative values after diffusion solver
        ufd(:, :pver) = ufd_out(:, :pver)
        vfd(:, :pver) = vfd_out(:, :pver)
        slfd(:, :pver) = slfd_out(:, :pver)

        ! Diffuse dry static energy
        call vertical_diffusion_diffuse_dry_static_energy_run( &
          ncol=ncol, &
          pver=pver, &
          dt=dt, &
          gravit=gravit, &
          p=p, & ! Coords1D, pressure coordinates [Pa]
          rhoi=rhoi(:ncol, :pverp), &
          shflx=shflx(:ncol), &
          dse_top=zero(:ncol), & ! = zero
          kvh=kvh(:ncol, :pverp), &
          cgh=cgh(:ncol, :pverp), &
          dpidz_sq=dpidz_sq(:ncol, :pverp), & ! moist
          ! input/output
          dse=slfd(:ncol, :pver), &
          errmsg=errmsg, &
          errflg=errflg)

        if (errflg /= 0) return

        ! Diffuse tracers
        ! Only water vapor is actually diffused here; all constituent indices are subset to just water vapor,
        ! or sized 1, in order to be compatible with the underlying CCPPized subroutine.
        ubc_mmr_dummy(:ncol, :1) = 0._kind_phys
        cnst_fixed_ubc(:1) = .false.

        qtfd_out(:, :, :) = 0._kind_phys
        call vertical_diffusion_diffuse_tracers_run( &
          ncol = ncol, &
          pver = pver, &
          ncnst = 1, & ! only water vapor is diffused here.
          dt = dt, &
          rair = rair, &
          gravit = gravit, &
          do_diffusion_const = do_diffusion_const_wet, & ! moist constituents to diffuse
          p = p, & ! Coords1D, pressure coordinates [Pa]
          t = t(:ncol, :pver), &
          rhoi = rhoi(:ncol, :pverp), &
          cflx = qflx(:ncol, const_wv_idx:const_wv_idx), & ! subset to water vapor only.
          kvh = kvh(:ncol, :pverp), &
          kvq = kvh(:ncol, :pverp), & ! [sic] kvh used for kvq here.
          cgs = cgs(:ncol, :pverp), &
          qmincg = zero(:ncol), &
          dpidz_sq = dpidz_sq(:ncol, :pverp), & ! moist
          ! upper boundary conditions from ubc module
          ubc_mmr = ubc_mmr_dummy(:ncol, :1), &
          cnst_fixed_ubc = cnst_fixed_ubc(:1), & ! = .false.
          q0 = qtfd(:ncol, :pver, :1), &
          q1 = qtfd_out(:ncol, :pver, :1), &
          errmsg = errmsg, &
          errflg = errflg)

        if (errflg /= 0) return

        ! update q with after in iterative process.
        qtfd(:ncol, :pver, :1) = qtfd_out(:ncol, :pver, :1)

        ! TODO (hplin, 5/9/2025): after these are subset to ncol check if we
        ! need to initialize some outs to 0; compute_vdiff did not do this before

        ! Retrieve (tfd,qvfd,qlfd) from (slfd,qtfd) in order to
        ! use 'trbintd' at the next iteration.

        do k = 1, pver
          do i = 1, ncol
            ! ----------------------------------------------------- !
            ! Compute the condensate 'qlfd' in the updated profiles !
            ! ----------------------------------------------------- !
            ! Option.1 : Assume grid-mean condensate is homogeneously diffused by the moist turbulence scheme.
            !            This should be used if 'pseudodiff = .false.' in vertical_diffusion.F90.
            ! Modification : Need to be check whether below is correct in the presence of ice, qi.
            !                I should understand why the variation of ice, qi is neglected during diffusion.
            templ = (slfd(i, k) - gravit*z(i, k))/cpair
            call qsat(templ, pmid(i, k), es, qs)
            ep2 = .622_kind_phys
            temps = templ + (qtfd(i, k, 1) - qs)/(cpair/latvap + latvap*qs/(rair*templ**2))
            call qsat(temps, pmid(i, k), es, qs)
            qlfd(i, k) = max(qtfd(i, k, 1) - qi(i, k) - qs, 0._kind_phys)
            ! Option.2 : Assume condensate is not diffused by the moist turbulence scheme.
            !            This should bs used if 'pseudodiff = .true.'  in vertical_diffusion.F90.
            ! qlfd(i,k) = ql(i,k)
            ! ----------------------------- !
            ! Compute the other 'qvfd, tfd' !
            ! ----------------------------- !
            qvfd(i, k) = max(0._kind_phys, qtfd(i, k, 1) - qi(i, k) - qlfd(i, k))
            tfd(i, k) = (slfd(i, k) + latvap*qlfd(i, k) + (latvap + latice)*qi(i, k) - gravit*z(i, k))/cpair
          end do
        end do
      end if

    end do iturb_loop  ! End of 'iturb' iteration

    kvq(:ncol, :) = kvh(:ncol, :)

  end subroutine bretherton_park_diff_run

end module bretherton_park_diff
