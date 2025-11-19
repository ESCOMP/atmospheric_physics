! Interstitial scheme to prepare inputs for vertical diffusion solver
! from surface coupler fluxes and stress data
module holtslag_boville_diff_interstitials
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! CCPP-compliant public interfaces
  public :: hb_diff_set_vertical_diffusion_top_init
  public :: hb_diff_set_total_surface_stress_run
  public :: hb_diff_prepare_vertical_diffusion_inputs_run
  public :: hb_free_atm_diff_prepare_vertical_diffusion_inputs_run

contains

  ! Interstitial to set top of vertical diffusion.
!> \section arg_table_hb_diff_set_vertical_diffusion_top_init Argument Table
!! \htmlinclude hb_diff_set_vertical_diffusion_top_init.html
  subroutine hb_diff_set_vertical_diffusion_top_init( &
    ntop_eddy, &
    errmsg, errflg)

    ! Host model dependency
    ! See note below. If removal of this dependency is desired in a low-top
    ! model configuration, this call can simply be replaced with ntop_eddy = 1.
    use ref_pres, only: press_lim_idx

    ! Output arguments
    integer,            intent(out) :: ntop_eddy          ! Vertical layer index of vertical diffusion top [index]
    character(len=512), intent(out) :: errmsg             ! Error message
    integer,            intent(out) :: errflg             ! Error flag

    ! Local parameters
    real(kind_phys), parameter :: ntop_eddy_pres = 1.e-7_kind_phys  ! Pressure below which eddy diffusion is not done in WACCM-X [Pa]

    errmsg = ''
    errflg = 0

    ! Set top level for vertical diffusion.
    ! When WACCM-X is active, the top of the model is sufficiently high that
    ! ntop_eddy is not at level 1 (top-of-atmosphere).
    ! In other configurations, the top of the model pressure is larger than
    ! ntop_eddy_pres, thus ntop_eddy computes to 1.
    ntop_eddy = press_lim_idx(ntop_eddy_pres, top=.true.)

  end subroutine hb_diff_set_vertical_diffusion_top_init

  ! Set total surface stresses for input into the HB PBL scheme.
  !
  ! NOTE: Temporarily, TMS and Beljaars are "hard-coupled" into this subroutine because
  ! the current focus is CAM4 and TMS/Beljaars drag have not been CCPPized.
  ! In the future, after TMS/Beljaars is implemented, this subroutine should only initialize
  ! surface stresses from the coupler, then TMS and Beljaars can work on adding the respective
  ! surface stresses to the tautotx/tautoty used by the PBL scheme.
!> \section arg_table_hb_diff_set_total_surface_stress_run Argument Table
!! \htmlinclude hb_diff_set_total_surface_stress_run.html
  subroutine hb_diff_set_total_surface_stress_run( &
    ncol, &
    wsx_from_coupler, wsy_from_coupler, &
    tautmsx, tautmsy, &
    taubljx, taubljy, &
    tautotx, tautoty, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol
    real(kind_phys),    intent(in)  :: wsx_from_coupler(:)      ! Surface eastward wind stress from coupler [Pa]
    real(kind_phys),    intent(in)  :: wsy_from_coupler(:)      ! Surface northward wind stress from coupler [Pa]
    real(kind_phys),    intent(in)  :: tautmsx(:)               ! Eastward turbulent mountain surface stress [Pa]
    real(kind_phys),    intent(in)  :: tautmsy(:)               ! Northward turbulent mountain surface stress [Pa]
    real(kind_phys),    intent(in)  :: taubljx(:)               ! Eastward Beljaars surface stress [Pa]
    real(kind_phys),    intent(in)  :: taubljy(:)               ! Northward Beljaars surface stress [Pa]

    ! Output arguments
    real(kind_phys),    intent(out) :: tautotx(:)               ! Eastward total stress at surface for boundary layer scheme [Pa]
    real(kind_phys),    intent(out) :: tautoty(:)               ! Northward total stress at surface for boundary layer scheme [Pa]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Initialize total surface stresses from coupler
    ! These are used for HB diffusion scheme and later PBL diagnostics but
    ! not for the vertical diffusion solver, which uses surface stresses from the coupler
    ! or just zero (in the case of CLUBB)
    tautotx(:ncol) = wsx_from_coupler(:ncol)
    tautoty(:ncol) = wsy_from_coupler(:ncol)

    ! Add turbulent mountain stress to total surface stress
    tautotx(:ncol) = tautotx(:ncol) + tautmsx(:ncol)
    tautoty(:ncol) = tautoty(:ncol) + tautmsy(:ncol)

    ! Add Beljaars integrated drag to total surface stress
    tautotx(:ncol) = tautotx(:ncol) + taubljx(:ncol)
    tautoty(:ncol) = tautoty(:ncol) + taubljy(:ncol)

  end subroutine hb_diff_set_total_surface_stress_run

  ! Interstitial for full HB (CAM4) to handle inputs from coupler and pass them
  ! to vertical diffusion solver.
!> \section arg_table_hb_diff_prepare_vertical_diffusion_inputs_run Argument Table
!! \htmlinclude hb_diff_prepare_vertical_diffusion_inputs_run.html
  subroutine hb_diff_prepare_vertical_diffusion_inputs_run( &
    ncol, pverp, pcnst, &
    wsx_from_coupler, wsy_from_coupler, &
    shf_from_coupler, &
    cflx_from_coupler, &
    ! below output
    taux, tauy, &
    shflux, &
    cflux, &
    itaures, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol       ! Number of atmospheric columns [count]
    integer,            intent(in)  :: pverp      ! Number of vertical interfaces [count]
    integer,            intent(in)  :: pcnst      ! Number of CCPP constituents [count]
    real(kind_phys),    intent(in)  :: wsx_from_coupler(:)     ! Surface eastward wind stress from coupler [Pa]
    real(kind_phys),    intent(in)  :: wsy_from_coupler(:)     ! Surface northward wind stress from coupler [Pa]
    real(kind_phys),    intent(in)  :: shf_from_coupler(:)     ! Surface upward sensible heat flux from coupler [W m-2]
    real(kind_phys),    intent(in)  :: cflx_from_coupler(:,:)  ! Surface upward constituent fluxes from coupler [kg m-2 s-1]

    ! Output arguments
    real(kind_phys),    intent(out) :: taux(:)                 ! Eastward stress at surface for vertical diffusion [Pa]
    real(kind_phys),    intent(out) :: tauy(:)                 ! Northward stress at surface for vertical diffusion [Pa]
    real(kind_phys),    intent(out) :: shflux(:)               ! Surface upward sensible heat flux for vertical diffusion [W m-2]
    real(kind_phys),    intent(out) :: cflux(:,:)              ! Surface upward constituent fluxes for vertical diffusion [kg m-2 s-1]
    logical,            intent(out) :: itaures                 ! Flag for updating residual stress at surface in vertical diffusion [flag]
    character(len=512), intent(out) :: errmsg                  ! Error message
    integer,            intent(out) :: errflg                  ! Error flag

    errmsg = ''
    errflg = 0

    ! Copy surface fluxes from coupler to vertical diffusion input arrays
    taux(:ncol)    = wsx_from_coupler(:ncol)
    tauy(:ncol)    = wsy_from_coupler(:ncol)
    shflux(:ncol)  = shf_from_coupler(:ncol)
    cflux(:ncol,:pcnst) = cflx_from_coupler(:ncol,:pcnst)

    ! Set flag for updating residual stress to true
    itaures = .true.

  end subroutine hb_diff_prepare_vertical_diffusion_inputs_run

  ! Interstitial for free atmosphere version of HB used above CLUBB which will allow
  ! the diffusion solver to handle non-water vapor surface fluxes (CAM6)
  ! or no surface fluxes (CAM7)
!> \section arg_table_hb_free_atm_diff_prepare_vertical_diffusion_inputs_run Argument Table
!! \htmlinclude hb_free_atm_diff_prepare_vertical_diffusion_inputs_run.html
  subroutine hb_free_atm_diff_prepare_vertical_diffusion_inputs_run( &
    ncol, pverp, pcnst, &
    const_props, &
    apply_nonwv_cflx, &
    cflx_from_coupler, &
    ! below output
    taux, tauy, &
    shflux, &
    cflux, &
    itaures, &
    errmsg, errflg)

    use coords_1d,  only: Coords1D

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! dependency to get constituent index
    use ccpp_const_utils,          only: ccpp_const_get_idx

    ! Input arguments
    integer,            intent(in)  :: ncol       ! Number of atmospheric columns [count]
    integer,            intent(in)  :: pverp      ! Number of vertical interfaces [count]
    integer,            intent(in)  :: pcnst      ! Number of CCPP constituents [count]
    type(ccpp_constituent_prop_ptr_t), &
                        intent(in)  :: const_props(:)           ! CCPP constituent properties pointer
    logical,            intent(in)  :: apply_nonwv_cflx         ! Flag for applying constituent fluxes excluding water vapor [flag]
    real(kind_phys),    intent(in)  :: cflx_from_coupler(:,:)   ! Surface upward constituent fluxes from coupler [kg m-2 s-1]

    ! Output arguments
    real(kind_phys),    intent(out) :: taux(:)                  ! Eastward stress at surface for vertical diffusion [Pa]
    real(kind_phys),    intent(out) :: tauy(:)                  ! Northward stress at surface for vertical diffusion [Pa]
    real(kind_phys),    intent(out) :: shflux(:)                ! Surface upward sensible heat flux for vertical diffusion [W m-2]
    real(kind_phys),    intent(out) :: cflux(:,:)               ! Surface upward constituent fluxes for vertical diffusion [kg m-2 s-1]
    logical,            intent(out) :: itaures                  ! Flag for updating residual stress at surface in vertical diffusion [flag]
    character(len=512), intent(out) :: errmsg                   ! Error message
    integer,            intent(out) :: errflg                   ! Error flag

    ! Local variables
    integer :: const_wv_idx                                     ! Water vapor constituent index

    errmsg = ''
    errflg = 0

    ! Check constituents list and locate water vapor index
    ! (not assumed to be 1)
    call ccpp_const_get_idx(const_props, &
         'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', &
         const_wv_idx, errmsg, errflg)
    if (errflg /= 0) return

    ! CLUBB applies some fluxes itself, but we still want constituent
    ! fluxes applied here (except water vapor).
    taux(:ncol)   = 0._kind_phys
    tauy(:ncol)   = 0._kind_phys
    shflux(:ncol) = 0._kind_phys

    ! If apply_nonwv_cflx is on, then non-water vapor constituent fluxes from the coupler
    ! are applied using vertical diffusion.
    if (apply_nonwv_cflx) then
      ! Copy non-water vapor constituent fluxes from coupler
      cflux(:ncol, :) = cflx_from_coupler(:ncol, :)
      ! But still zero out water vapor flux
      cflux(:ncol, const_wv_idx) = 0._kind_phys
    else
      ! Surface fluxes applied in CLUBB emissions module (CAM7)
      ! so vertical diffusion applies no fluxes.
      cflux(:ncol, :) = 0._kind_phys
    end if

    ! Set flag for updating residual stress to true
    itaures = .true.

  end subroutine hb_free_atm_diff_prepare_vertical_diffusion_inputs_run

end module holtslag_boville_diff_interstitials
