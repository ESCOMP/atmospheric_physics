! Interstitial scheme to prepare inputs for vertical diffusion solver
! common to all PBL schemes.
!
! These include:
!   - the p vertical coordinate (coords1d), including destroying this object at end of timestep.
!   - potential temperature (th)
module vertical_diffusion_interstitials
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! CCPP-compliant public interfaces
  public :: vertical_diffusion_prepare_inputs_run
  public :: vertical_diffusion_prepare_inputs_timestep_final

  public :: compute_kinematic_fluxes_and_obklen_run

contains

  ! Interstitial to prepare moist vertical coordinate for solver,
  ! and compute potential temperature.
!> \section arg_table_vertical_diffusion_prepare_inputs_run Argument Table
!! \htmlinclude vertical_diffusion_prepare_inputs_run.html
  subroutine vertical_diffusion_prepare_inputs_run( &
    ncol, pver, pverp, &
    pint, &
    t, exner, &
    p, &
    th, &
    errmsg, errflg)

    use coords_1d,  only: coords1d

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    integer,            intent(in)  :: pverp
    real(kind_phys),    intent(in)  :: pint(:,:) ! Air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)  :: t(:,:)    ! temperature [K]
    real(kind_phys),    intent(in)  :: exner(:,:)! exner [1]
    type(coords1d),     intent(out) :: p         ! Vertical moist pressure coordinates for vertical diffusion [Pa]
    real(kind_phys),    intent(out) :: th(:,:)   ! Potential temperature [K]
    character(len=512), intent(out) :: errmsg    ! Error message
    integer,            intent(out) :: errflg    ! Error flag

    errmsg = ''
    errflg = 0

    ! Initialize pressure coordinate object for vertical diffusion solver
    p = coords1d(pint(:ncol,:pverp))

    ! Calculate potential temperature [K]
    th(:ncol,:) = t(:ncol,:) * exner(:ncol,:)

  end subroutine vertical_diffusion_prepare_inputs_run

  ! Interstitial to clean up vertical coordinate after use.
!> \section arg_table_vertical_diffusion_prepare_inputs_timestep_final Argument Table
!! \htmlinclude vertical_diffusion_prepare_inputs_timestep_final.html
  subroutine vertical_diffusion_prepare_inputs_timestep_final(p, errmsg, errflg)
    use coords_1d,  only: coords1d

    type(coords1d),     intent(inout) :: p       ! Vertical moist pressure coordinates for vertical diffusion [Pa]
    character(len=512), intent(out)   :: errmsg  ! Error message
    integer,            intent(out)   :: errflg  ! Error flag

    errmsg = ''
    errflg = 0

    call p%finalize()
  end subroutine vertical_diffusion_prepare_inputs_timestep_final

  ! Compute kinematic fluxes (heat, water vapor, buoyancy), and obukhov length
!> \section arg_table_compute_kinematic_fluxes_and_obklen_run Argument Table
!! \htmlinclude compute_kinematic_fluxes_and_obklen_run.html
  subroutine compute_kinematic_fluxes_and_obklen_run( &
    ncol, pver, pcnst, &
    const_props, &
    zvir, cpair, gravit, karman, &
    shf_from_coupler, cflx_from_coupler, &
    q_wv, &
    th, &
    rrho, &
    ustar, &
    ! below output:
    khfs, kqfs, kbfs, &
    obklen, &
    errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! dependency to get constituent index
    use ccpp_const_utils,          only: ccpp_const_get_idx

    use atmos_phys_pbl_utils,      only: calc_virtual_temperature, calc_obukhov_length, &
                                         calc_kinematic_heat_flux, calc_kinematic_water_vapor_flux, &
                                         calc_kinematic_buoyancy_flux

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    integer,            intent(in)  :: pcnst
    type(ccpp_constituent_prop_ptr_t), &
                        intent(in)  :: const_props(:)          ! CCPP constituent properties pointer
    real(kind_phys),    intent(in)  :: zvir
    real(kind_phys),    intent(in)  :: cpair
    real(kind_phys),    intent(in)  :: gravit
    real(kind_phys),    intent(in)  :: karman
    real(kind_phys),    intent(in)  :: shf_from_coupler(:)     ! Surface upward sensible heat flux from coupler [W m-2]
    real(kind_phys),    intent(in)  :: cflx_from_coupler(:,:)  ! Surface upward constituent fluxes from coupler [kg m-2 s-1]
    real(kind_phys),    intent(in)  :: q_wv(:,:)               ! specific humidity [kg kg-1]

    real(kind_phys),    intent(in)  :: th(:,:)                 ! Potential temperature [K]

    ! rrho and ustar from PBL scheme.
    ! for HB, it is computed in hb_pbl_independent_coefficients_run;
    ! for UW, it is iteratively updated in bretherton_park_diff_run and returned as a final result.
    real(kind_phys),    intent(in)  :: rrho(:)                 ! 1 / bottom level density [m3 kg-1]
    real(kind_phys),    intent(in)  :: ustar(:)                ! surface friction velocity [m s-1]

    ! Output arguments:
    real(kind_phys),    intent(out) :: khfs(:)                 ! kinematic surface heat flux [K m s-1]
    real(kind_phys),    intent(out) :: kqfs(:)                 ! kinematic surface water vapor flux [kg kg-1 m s-1]
    real(kind_phys),    intent(out) :: kbfs(:)                 ! surface kinematic buoyancy flux [m2 s-3]
    real(kind_phys),    intent(out) :: obklen(:)               ! Obukhov length [m]

    character(len=512), intent(out) :: errmsg                  ! Error message
    integer,            intent(out) :: errflg                  ! Error flag

    ! Local variables
    integer :: const_wv_idx                                    ! Water vapor constituent index
    real(kind_phys) :: thvs(ncol)                              ! Virtual potential temperature at surface [K]

    ! Check constituents list and locate water vapor index
    ! (not assumed to be 1)
    call ccpp_const_get_idx(const_props, &
         'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', &
         const_wv_idx, errmsg, errflg)
    if (errflg /= 0) return

    thvs  (:ncol) = calc_virtual_temperature(th(:ncol,pver), q_wv(:ncol,pver), zvir)
    khfs  (:ncol) = calc_kinematic_heat_flux(shf_from_coupler(:ncol), rrho(:ncol), cpair)
    kqfs  (:ncol) = calc_kinematic_water_vapor_flux(cflx_from_coupler(:ncol,const_wv_idx), rrho(:ncol))
    kbfs  (:ncol) = calc_kinematic_buoyancy_flux(khfs(:ncol), zvir, th(:ncol,pver), kqfs(:ncol))
    obklen(:ncol) = calc_obukhov_length(thvs(:ncol), ustar(:ncol), gravit, karman, kbfs(:ncol))

  end subroutine compute_kinematic_fluxes_and_obklen_run

end module vertical_diffusion_interstitials
