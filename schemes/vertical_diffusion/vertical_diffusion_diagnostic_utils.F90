! Diagnostic helper functions for vertical diffusion
module vertical_diffusion_diagnostic_utils
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: vertical_diffusion_diagnostic_profiles

contains

  ! Computes profile output before and after vertical diffusion
  ! for PBL analysis
  !
  ! This subroutine was designed to be reusable in current CAM and SIMA
  subroutine vertical_diffusion_diagnostic_profiles( &
             ncol, pver, pverp, &
             ztodt, &
             latvap, latice, zvir, cpair, gravit, rair, &
             pint, pmid, zi, zm, &
             kvh, kvm, cgh, cgs, &
             cam_in_shf, &
             cam_in_cflx_wv, cam_in_cflx_cldliq, cam_in_cflx_cldice, &
             tautotx, tautoty, &
             t0, q0_wv, q0_cldliq, q0_cldice, s0, u0, v0, &
                 q1_wv, q1_cldliq, q1_cldice, s1, u1, v1, &
             tend_s, tend_u, tend_v, tend_q_wv, tend_q_cldliq, tend_q_cldice, &
             qt_pre_PBL, sl_pre_PBL, slv_pre_PBL, ftem_pre_PBL, &
             qt_aft_PBL, sl_aft_PBL, slv_aft_PBL, &
             qv_aft_PBL, ql_aft_PBL, qi_aft_PBL, &
             t_aft_PBL, ftem_aft_PBL, &
             u_aft_PBL, v_aft_PBL, &
             slten, qtten, tten, rhten, &
             slflx, qtflx, uflx, vflx, &
             slflx_cg, qtflx_cg, uflx_cg, vflx_cg)

    ! To-be-ccppized dependency.
    use wv_saturation, only: qsat

    ! Input arguments
    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: pverp
    real(kind_phys),  intent(in)  :: ztodt                   ! Physics timestep [s]
    real(kind_phys),  intent(in)  :: latvap                  ! Latent heat of vaporization [J kg-1]
    real(kind_phys),  intent(in)  :: latice                  ! Latent heat of fusion [J kg-1]
    real(kind_phys),  intent(in)  :: zvir                    ! rh2o/rair - 1 [1]
    real(kind_phys),  intent(in)  :: cpair                   ! Specific heat of dry air at constant pressure [J kg-1 K-1]
    real(kind_phys),  intent(in)  :: gravit                  ! Acceleration due to gravity [m s-2]
    real(kind_phys),  intent(in)  :: rair                    ! Gas constant for dry air [J kg-1 K-1]
    real(kind_phys),  intent(in)  :: pint(:, :)              ! Interface pressure [Pa]
    real(kind_phys),  intent(in)  :: pmid(:, :)              ! Midpoint pressure [Pa]
    real(kind_phys),  intent(in)  :: zi(:, :)                ! Geopotential height at interfaces [m]
    real(kind_phys),  intent(in)  :: zm(:, :)                ! Geopotential height at midpoints [m]

    ! Diffusion coefficients and counter-gradient terms
    real(kind_phys),  intent(in)  :: kvh(:, :)               ! Eddy diffusivity for heat at interfaces [m2 s-1]
    real(kind_phys),  intent(in)  :: kvm(:, :)               ! Eddy diffusivity for momentum at interfaces [m2 s-1]
    real(kind_phys),  intent(in)  :: cgh(:, :)               ! Counter-gradient term for heat [K m s-1]
    real(kind_phys),  intent(in)  :: cgs(:, :)               ! Counter-gradient star [s m-2]

    ! Surface fluxes from coupler
    real(kind_phys),  intent(in)  :: cam_in_shf(:)           ! Surface sensible heat flux [W m-2]

    ! NOTE: These surface fluxes need to be disassembled from cam_in_cflx into wv, cldliq, cldice
    real(kind_phys),  intent(in)  :: cam_in_cflx_wv(:)       ! Surface water vapor flux [kg m-2 s-1]
    real(kind_phys),  intent(in)  :: cam_in_cflx_cldliq(:)   ! Surface cloud liquid flux [kg m-2 s-1]
    real(kind_phys),  intent(in)  :: cam_in_cflx_cldice(:)   ! Surface cloud ice flux [kg m-2 s-1]

    ! Total surface stresses
    ! NOTE: pre-update (not using provisional winds for TMS/Beljaars orographic surface drag)
    real(kind_phys),  intent(in)  :: tautotx(:)              ! Total zonal surface stress [N m-2]
    real(kind_phys),  intent(in)  :: tautoty(:)              ! Total meridional surface stress [N m-2]

    ! Initial state (before vertical diffusion)
    real(kind_phys),  intent(in)  :: t0(:, :)                ! Temperature before diffusion [K]
    real(kind_phys),  intent(in)  :: q0_wv(:, :)             ! Water vapor mixing ratio before diffusion [kg kg-1]
    real(kind_phys),  intent(in)  :: q0_cldliq(:, :)         ! Cloud liquid water mixing ratio before diffusion [kg kg-1]
    real(kind_phys),  intent(in)  :: q0_cldice(:, :)         ! Cloud ice water mixing ratio before diffusion [kg kg-1]
    real(kind_phys),  intent(in)  :: s0(:, :)                ! Dry static energy before diffusion [J kg-1]
    real(kind_phys),  intent(in)  :: u0(:, :)                ! Zonal wind before diffusion [m s-1]
    real(kind_phys),  intent(in)  :: v0(:, :)                ! Meridional wind before diffusion [m s-1]

    ! Final state (after vertical diffusion)
    real(kind_phys),  intent(in)  :: q1_wv(:, :)             ! Water vapor mixing ratio after diffusion [kg kg-1]
    real(kind_phys),  intent(in)  :: q1_cldliq(:, :)         ! Cloud liquid water mixing ratio after diffusion [kg kg-1]
    real(kind_phys),  intent(in)  :: q1_cldice(:, :)         ! Cloud ice water mixing ratio after diffusion [kg kg-1]
    real(kind_phys),  intent(in)  :: s1(:, :)                ! Dry static energy after diffusion [J kg-1]
    real(kind_phys),  intent(in)  :: u1(:, :)                ! Zonal wind after diffusion [m s-1]
    real(kind_phys),  intent(in)  :: v1(:, :)                ! Meridional wind after diffusion [m s-1]

    ! Tendencies attributable to vertical diffusion
    real(kind_phys),  intent(in)  :: tend_s(:, :)            ! Dry static energy tendency [J kg-1 s-1]
    real(kind_phys),  intent(in)  :: tend_u(:, :)            ! Zonal wind tendency [m s-2]
    real(kind_phys),  intent(in)  :: tend_v(:, :)            ! Meridional wind tendency [m s-2]

    ! TODO? hplin -- constituent tendencies were assembled into ccpp_constituent_tendencies
    ! in vertical_diffusion_tendencies_run. The diagnostic CCPP scheme may? have to use
    ! constituent properties to split out the wv, cldliq, cldice tendencies
    real(kind_phys),  intent(in)  :: tend_q_wv(:, :)         ! Water vapor tendency [kg kg-1 s-1]
    real(kind_phys),  intent(in)  :: tend_q_cldliq(:, :)     ! Cloud liquid water tendency [kg kg-1 s-1]
    real(kind_phys),  intent(in)  :: tend_q_cldice(:, :)     ! Cloud ice water tendency [kg kg-1 s-1]

    ! Output diagnostic profiles - before vertical diffusion
    real(kind_phys),  intent(out) :: qt_pre_PBL(:, :)        ! Total water mixing ratio before PBL [kg kg-1]
    real(kind_phys),  intent(out) :: sl_pre_PBL(:, :)        ! Liquid water static energy before PBL [J kg-1]
    real(kind_phys),  intent(out) :: slv_pre_PBL(:, :)       ! Virtual liquid water static energy before PBL [J kg-1]
    real(kind_phys),  intent(out) :: ftem_pre_PBL(:, :)      ! Relative humidity before PBL [percent]

    ! Output diagnostic profiles - after vertical diffusion
    real(kind_phys),  intent(out) :: qt_aft_PBL(:, :)        ! Total water mixing ratio after PBL [kg kg-1]
    real(kind_phys),  intent(out) :: sl_aft_PBL(:, :)        ! Liquid water static energy after PBL [J kg-1]
    real(kind_phys),  intent(out) :: slv_aft_PBL(:, :)       ! Virtual liquid water static energy after PBL [J kg-1]
    real(kind_phys),  intent(out) :: qv_aft_PBL(:, :)        ! Water vapor mixing ratio after PBL [kg kg-1]
    real(kind_phys),  intent(out) :: ql_aft_PBL(:, :)        ! Cloud liquid water mixing ratio after PBL [kg kg-1]
    real(kind_phys),  intent(out) :: qi_aft_PBL(:, :)        ! Cloud ice water mixing ratio after PBL [kg kg-1]
    real(kind_phys),  intent(out) :: t_aft_PBL(:, :)         ! Temperature after PBL [K]
    real(kind_phys),  intent(out) :: ftem_aft_PBL(:, :)      ! Relative humidity after PBL [percent]
    real(kind_phys),  intent(out) :: u_aft_PBL(:, :)         ! Zonal wind after PBL [m s-1]
    real(kind_phys),  intent(out) :: v_aft_PBL(:, :)         ! Meridional wind after PBL [m s-1]

    ! Output tendency diagnostics
    real(kind_phys),  intent(out) :: slten(:, :)             ! Liquid water static energy tendency [J kg-1 s-1]
    real(kind_phys),  intent(out) :: qtten(:, :)             ! Total water mixing ratio tendency [kg kg-1 s-1]
    real(kind_phys),  intent(out) :: tten(:, :)              ! Temperature tendency [K s-1]
    real(kind_phys),  intent(out) :: rhten(:, :)             ! Relative humidity tendency [% s-1]

    ! Output flux diagnostics
    real(kind_phys),  intent(out) :: slflx(:, :)             ! Liquid static energy flux at interfaces [W m-2]
    real(kind_phys),  intent(out) :: qtflx(:, :)             ! Total water flux at interfaces [kg m-2 s-1]
    real(kind_phys),  intent(out) :: uflx(:, :)              ! Zonal momentum flux at interfaces [N m-2]
    real(kind_phys),  intent(out) :: vflx(:, :)              ! Meridional momentum flux at interfaces [N m-2]
    real(kind_phys),  intent(out) :: slflx_cg(:, :)          ! Counter-gradient liquid static energy flux at interfaces [W m-2]
    real(kind_phys),  intent(out) :: qtflx_cg(:, :)          ! Counter-gradient total water flux at interfaces [kg m-2 s-1]
    real(kind_phys),  intent(out) :: uflx_cg(:, :)           ! Counter-gradient zonal momentum flux at interfaces [N m-2]
    real(kind_phys),  intent(out) :: vflx_cg(:, :)           ! Counter-gradient meridional momentum flux at interfaces [N m-2]

    ! Local variables
    integer :: i, k
    real(kind_phys) :: rztodt                                ! Reciprocal of timestep [s-1]
    real(kind_phys) :: rhoair                                ! Air density [kg m-3]
    real(kind_phys) :: ftem_pre(ncol, pver)                  ! Saturation vapor pressure before PBL [Pa]
    real(kind_phys) :: ftem_aft(ncol, pver)                  ! Saturation vapor pressure after PBL [Pa]
    real(kind_phys) :: tem2_pre(ncol, pver)                  ! Saturation specific humidity before PBL [kg kg-1]
    real(kind_phys) :: tem2_aft(ncol, pver)                  ! Saturation specific humidity after PBL [kg kg-1]
    real(kind_phys) :: s_aft_PBL(ncol, pver)                 ! Dry static energy after PBL [J kg-1]

    rztodt = 1.0_kind_phys / ztodt

    ! ====================================================
    ! Compute profiles before vertical diffusion (pre-PBL)
    ! ====================================================

    ! Liquid water static energy before PBL
    sl_pre_PBL(:ncol, :pver) = s0(:ncol, :pver) &
                               - latvap * q0_cldliq(:ncol, :pver) &
                               - (latvap + latice) * q0_cldice(:ncol, :pver)

    ! Total water mixing ratio before PBL
    qt_pre_PBL(:ncol, :pver) = q0_wv(:ncol, :pver) &
                               + q0_cldliq(:ncol, :pver) &
                               + q0_cldice(:ncol, :pver)

    ! Virtual liquid water static energy before PBL
    slv_pre_PBL(:ncol, :pver) = sl_pre_PBL(:ncol, :pver) &
                                * (1.0_kind_phys + zvir * qt_pre_PBL(:ncol, :pver))

    ! Compute saturation vapor pressure and relative humidity before PBL
    do k = 1, pver
       call qsat(t0(:ncol, k), pmid(:ncol, k), tem2_pre(:ncol, k), ftem_pre(:ncol, k), ncol)
    end do
    ftem_pre_PBL(:ncol, :pver) = q0_wv(:ncol, :pver) / ftem_pre(:ncol, :pver) * 100.0_kind_phys

    ! ====================================================
    ! Compute profiles after vertical diffusion (post-PBL)
    ! ====================================================

    ! Apply tendencies to get final state variables
    !
    ! Note: this is indeed computed twice because at the point this utility subroutine runs,
    ! the tendencies have not yet been applied.
    qv_aft_PBL(:ncol, :pver) = q0_wv(:ncol, :pver)     + tend_q_wv(:ncol, :pver)     * ztodt
    ql_aft_PBL(:ncol, :pver) = q0_cldliq(:ncol, :pver) + tend_q_cldliq(:ncol, :pver) * ztodt
    qi_aft_PBL(:ncol, :pver) = q0_cldice(:ncol, :pver) + tend_q_cldice(:ncol, :pver) * ztodt
    u_aft_PBL(:ncol, :pver)  = u0(:ncol, :pver)        + tend_u(:ncol, :pver)        * ztodt
    v_aft_PBL(:ncol, :pver)  = v0(:ncol, :pver)        + tend_v(:ncol, :pver)        * ztodt

    ! Dry static energy after PBL
    s_aft_PBL(:ncol, :pver) = s0(:ncol, :pver) + tend_s(:ncol, :pver) * ztodt

    ! Liquid water static energy after PBL
    sl_aft_PBL(:ncol, :pver) = s1(:ncol, :pver) &
                               - latvap * ql_aft_PBL(:ncol, :pver) & ! or q1_cldliq
                               - (latvap + latice) * qi_aft_PBL(:ncol, :pver) ! or q1_cldice

    ! Total water mixing ratio after PBL
    qt_aft_PBL(:ncol, :pver) = qv_aft_PBL(:ncol, :pver) &
                               + ql_aft_PBL(:ncol, :pver) &
                               + qi_aft_PBL(:ncol, :pver)

    ! Virtual liquid water static energy after PBL
    slv_aft_PBL(:ncol, :pver) = sl_aft_PBL(:ncol, :pver) &
                                * (1.0_kind_phys + zvir * qt_aft_PBL(:ncol, :pver))

    ! Temperature after PBL (derived from dry static energy)
    t_aft_PBL(:ncol, :pver) = (s_aft_PBL(:ncol, :pver) - gravit * zm(:ncol, :pver)) / cpair

    ! Compute saturation vapor pressure and relative humidity after PBL
    do k = 1, pver
       call qsat(t_aft_PBL(:ncol, k), pmid(:ncol, k), tem2_aft(:ncol, k), ftem_aft(:ncol, k), ncol)
    end do
    ftem_aft_PBL(:ncol, :pver) = qv_aft_PBL(:ncol, :pver) / ftem_aft(:ncol, :pver) * 100.0_kind_phys

    ! ====================================================
    ! Compute tendency diagnostics
    ! ====================================================

    ! Liquid water static energy tendency
    slten(:ncol, :pver) = (sl_aft_PBL(:ncol, :pver) - sl_pre_PBL(:ncol, :pver)) * rztodt

    ! Total water mixing ratio tendency
    qtten(:ncol, :pver) = (qt_aft_PBL(:ncol, :pver) - qt_pre_PBL(:ncol, :pver)) * rztodt

    ! Temperature tendency
    tten(:ncol, :pver) = (t_aft_PBL(:ncol, :pver) - t0(:ncol, :pver)) * rztodt

    ! Relative humidity tendency
    rhten(:ncol, :pver) = (ftem_aft_PBL(:ncol, :pver) - ftem_pre_PBL(:ncol, :pver)) * rztodt

    ! ===================================================================
    ! Compute flux diagnostics at interfaces
    ! ===================================================================

    ! Initialize top interface (no flux at TOA)
    slflx(:ncol, 1)    = 0.0_kind_phys
    qtflx(:ncol, 1)    = 0.0_kind_phys
    uflx(:ncol, 1)     = 0.0_kind_phys
    vflx(:ncol, 1)     = 0.0_kind_phys
    slflx_cg(:ncol, 1) = 0.0_kind_phys
    qtflx_cg(:ncol, 1) = 0.0_kind_phys
    uflx_cg(:ncol, 1)  = 0.0_kind_phys
    vflx_cg(:ncol, 1)  = 0.0_kind_phys

    ! Compute fluxes at model interfaces
    do k = 2, pver
       do i = 1, ncol
          ! Air density at interface k [kg m-3]
          rhoair = pint(i, k) / (rair * ((0.5_kind_phys * (slv_aft_PBL(i, k) + slv_aft_PBL(i, k-1)) &
                                        - gravit * zi(i, k)) / cpair))

          ! Liquid static energy flux [W m-2]
          slflx(i, k) = kvh(i, k) * &
                        (cgh(i, k) - rhoair * (sl_aft_PBL(i, k-1) - sl_aft_PBL(i, k)) / (zm(i, k-1) - zm(i, k)))

          ! Total water flux [kg m-2 s-1]
          qtflx(i, k) = kvh(i, k) * &
               ( - rhoair * (qt_aft_PBL(i, k-1) - qt_aft_PBL(i, k)) / (zm(i, k-1) - zm(i, k)) &
                 + rhoair * (cam_in_cflx_wv(i) + cam_in_cflx_cldliq(i) + cam_in_cflx_cldice(i)) * cgs(i, k) )

          ! Zonal momentum flux [N m-2]
          uflx(i, k) = kvm(i, k) * &
                       ( - rhoair * (u_aft_PBL(i, k-1) - u_aft_PBL(i, k)) / (zm(i, k-1) - zm(i, k)) )

          ! Meridional momentum flux [N m-2]
          vflx(i, k) = kvm(i, k) * &
                       ( - rhoair * (v_aft_PBL(i, k-1) - v_aft_PBL(i, k)) / (zm(i, k-1) - zm(i, k)) )

          ! Counter-gradient fluxes
          slflx_cg(i, k) = kvh(i, k) * cgh(i, k)
          qtflx_cg(i, k) = kvh(i, k) * rhoair * (cam_in_cflx_wv(i) + cam_in_cflx_cldliq(i) &
                                                + cam_in_cflx_cldice(i)) * cgs(i, k)
          uflx_cg(i, k)  = 0.0_kind_phys
          vflx_cg(i, k)  = 0.0_kind_phys
       end do
    end do

    ! Set fluxes at bottom interface using surface fluxes
    slflx(:ncol, pverp)    = cam_in_shf(:ncol)
    qtflx(:ncol, pverp)    = cam_in_cflx_wv(:ncol)
    uflx(:ncol, pverp)     = tautotx(:ncol)
    vflx(:ncol, pverp)     = tautoty(:ncol)
    slflx_cg(:ncol, pverp) = 0.0_kind_phys
    qtflx_cg(:ncol, pverp) = 0.0_kind_phys
    uflx_cg(:ncol, pverp)  = 0.0_kind_phys
    vflx_cg(:ncol, pverp)  = 0.0_kind_phys

  end subroutine vertical_diffusion_diagnostic_profiles

end module vertical_diffusion_diagnostic_utils
