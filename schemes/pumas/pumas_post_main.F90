! Post-processing interstitial for PUMAS.
module pumas_post_main
  implicit none
  private

  public :: pumas_post_main_run

contains

  !> \section arg_table_pumas_post_main_run Argument Table
  !! \htmlinclude pumas_post_main_run.html
  subroutine pumas_post_main_run(ncol, nlev, trop_cloud_top_lev, rair, pmid, temp, &
                cldliq, numliq, cldice, numice, snow, numsnow, graupel, numgraup,  &
                strat_cldfrc, effi, dei, pgam, lamc, des, degrau, errmsg, errcode)
    use ccpp_kinds,        only: kind_phys
    use pumas_kinds,       only: pumas_r8=>kind_r8
    use micro_pumas_utils, only: size_dist_param_basic, size_dist_param_liq, &
                                 mg_ice_props, mg_liq_props, avg_diameter,   &
                                 mincld, qsmall, rhog, rhosn, rhows, rhoi

    integer,                 intent(in)  :: ncol
    integer,                 intent(in)  :: nlev
    integer,                 intent(in)  :: trop_cloud_top_lev  ! Index of the top model level for which
                                                                ! tropospheric cloud physics is run (index)
    real(kind_phys),                 intent(in)  :: rair        ! gas constant of dry air (J kg-1 K-1)
    real(kind_phys), dimension(:,:), intent(in)  :: pmid        ! layer midpoint pressure (Pa)
    real(kind_phys), dimension(:,:), intent(in)  :: temp        ! air temperature, before microphysics heating is applied (K)
    real(kind_phys), dimension(:,:), intent(in)  :: cldliq      ! updated cloud liquid water mixing ratio (kg/kg)
    real(kind_phys), dimension(:,:), intent(in)  :: numliq      ! updated cloud droplet number concentration (kg-1)
    real(kind_phys), dimension(:,:), intent(in)  :: cldice      ! updated cloud ice mixing ratio (kg/kg)
    real(kind_phys), dimension(:,:), intent(in)  :: numice      ! updated cloud ice number concentration (kg-1)
    real(kind_phys), dimension(:,:), intent(in)  :: snow        ! updated snow mixing ratio (kg/kg)
    real(kind_phys), dimension(:,:), intent(in)  :: numsnow     ! updated snow number concentration (kg-1)
    real(kind_phys), dimension(:,:), intent(in)  :: graupel     ! updated graupel mixing ratio (kg/kg)
    real(kind_phys), dimension(:,:), intent(in)  :: numgraup    ! updated graupel number concentration (kg-1)
    real(kind_phys), dimension(:,:), intent(in)  :: strat_cldfrc! total stratiform cloud area fraction (= ast)
    real(kind_phys), dimension(:,:), intent(out) :: effi       ! ice effective radius (micron)
    real(kind_phys), dimension(:,:), intent(out) :: dei        ! ice effective diameter for radiation (micron)
    real(kind_phys), dimension(:,:), intent(out) :: pgam       ! droplet size distribution shape parameter for radiation (1)
    real(kind_phys), dimension(:,:), intent(out) :: lamc       ! droplet size distribution slope for radiation (1)
    real(kind_phys), dimension(:,:), intent(out) :: des        ! snow effective diameter for radiation (micron)
    real(kind_phys), dimension(:,:), intent(out) :: degrau     ! graupel effective diameter for radiation (m)
    character(len=*),          intent(out) :: errmsg
    integer,                 intent(out) :: errcode

    integer        :: k
    real(pumas_r8) :: icimrst(ncol, nlev)  ! in-cloud (grid-mean) ice mixing ratio
    real(pumas_r8) :: niic(ncol, nlev)     ! in-cloud (grid-mean) ice number conc
    real(pumas_r8) :: rei(ncol, nlev)      ! ice slope param, then effective radius
    real(pumas_r8) :: icwmrst(ncol, nlev)  ! in-cloud (grid-mean) liquid mixing ratio
    real(pumas_r8) :: ncic(ncol, nlev)     ! in-cloud (grid-mean) droplet number conc
    real(pumas_r8) :: rho(ncol, nlev)      ! air density
    real(pumas_r8) :: mu(ncol, nlev)       ! droplet size distribution shape parameter
    real(pumas_r8) :: lambdac(ncol, nlev)  ! droplet size distribution slope
    real(pumas_r8) :: qs(ncol, nlev)       ! grid-mean snow mixing ratio
    real(pumas_r8) :: ns(ncol, nlev)       ! grid-mean snow number conc
    real(pumas_r8) :: dsout2(ncol, nlev)   ! mean snow particle diameter
    real(pumas_r8) :: desm(ncol, nlev)     ! snow effective diameter (m)
    real(pumas_r8) :: qg(ncol, nlev)       ! grid-mean graupel mixing ratio
    real(pumas_r8) :: ng(ncol, nlev)       ! grid-mean graupel number conc
    real(pumas_r8) :: dgout2(ncol, nlev)   ! mean graupel particle diameter

    errmsg = ' '
    errcode = 0

    ! Ice effective radius is recomputed here from the post-microphysics
    ! grid-mean in-cloud ice, NOT from the raw per-substep value PUMAS returns.
    ! PUMAS uses the total stratiform fraction for ice cloud (icecldf = ast).
    icimrst(:,:) = min(real(cldice(:ncol,:), pumas_r8) / &
                max(mincld, real(strat_cldfrc(:ncol,:), pumas_r8)), 0.005_pumas_r8)
    niic(:,:)   =    real(numice(:ncol,:), pumas_r8) / &
                max(mincld, real(strat_cldfrc(:ncol,:), pumas_r8))

    rei(:,:) = 25._pumas_r8
    do k = trop_cloud_top_lev, nlev
      call size_dist_param_basic(mg_ice_props, icimrst(:,k), niic(:,k), rei(:,k), ncol)
    end do

    where (icimrst(:,:) >= qsmall)
      rei(:,:) = 1.5_pumas_r8 / rei(:,:) * 1.e6_pumas_r8
    elsewhere
      rei(:,:) = 25._pumas_r8
    end where

    effi(:ncol,:) = real(rei(:,:), kind_phys)

    ! Ice effective diameter for radiation follows from the recomputed
    ! effective radius.
    dei(:ncol,:) = real(rei(:,:) * rhoi/rhows * 2._pumas_r8, kind_phys)

    ! Droplet size distribution parameters for radiation are likewise recomputed
    ! from the post-microphysics grid-mean in-cloud liquid.
    ! Air density uses the temperature BEFORE the microphysics heating is applied
    ! The original comment in CAM noted:
    !   "State instead of state_loc to preserve answers for MG1 (and in any
    !    case, it is unlikely to make much difference)."
    ! PUMAS uses the total stratiform fraction for liquid cloud (liqcldf = ast)
    rho(:,:)     = real(pmid(:ncol,:), pumas_r8) / &
                (real(rair, pumas_r8) * real(temp(:ncol,:), pumas_r8))
    icwmrst(:,:) = min(real(cldliq(:ncol,:), pumas_r8) / &
                max(mincld, real(strat_cldfrc(:ncol,:), pumas_r8)), 0.005_pumas_r8)
    ncic(:,:)   =    real(numliq(:ncol,:), pumas_r8) / &
                max(mincld, real(strat_cldfrc(:ncol,:), pumas_r8))

    mu(:,:)      = 0._pumas_r8
    lambdac(:,:) = 0._pumas_r8
    do k = trop_cloud_top_lev, nlev
      call size_dist_param_liq(mg_liq_props, icwmrst(:,k), ncic(:,k), rho(:,k), &
                               mu(:,k), lambdac(:,k), ncol)
    end do

    ! size_dist_param_liq flags no-cloud points with mu = -100 and reset to zero.
    where (icwmrst(:,trop_cloud_top_lev:) < qsmall)
      mu(:,trop_cloud_top_lev:) = 0._pumas_r8
    end where

    pgam(:ncol,:) = real(mu(:,:), kind_phys)
    lamc(:ncol,:) = real(lambdac(:,:), kind_phys)

    ! Snow effective diameter for radiation, from the post-microphysics grid-mean snow
    qs(:,:) = real(snow(:ncol,:), pumas_r8)
    ns(:,:) = real(numsnow(:ncol,:), pumas_r8)

    dsout2(:,:) = 0._pumas_r8
    desm(:,:)   = 0._pumas_r8
    where (qs(:,trop_cloud_top_lev:) >= 1.e-7_pumas_r8)
      dsout2(:,trop_cloud_top_lev:) = avg_diameter( &
           qs(:,trop_cloud_top_lev:), &
           ns(:,trop_cloud_top_lev:) * rho(:,trop_cloud_top_lev:), &
           rho(:,trop_cloud_top_lev:), rhosn)
      desm(:,trop_cloud_top_lev:) = dsout2(:,trop_cloud_top_lev:) * &
           3._pumas_r8 * rhosn/rhows
    end where

    des(:ncol,:) = real(desm(:,:) * 1.e6_pumas_r8, kind_phys)

    ! Graupel effective diameter for radiation, from the post-microphysics grid-mean graupel
    qg(:,:) = real(graupel(:ncol,:), pumas_r8)
    ng(:,:) = real(numgraup(:ncol,:), pumas_r8)

    dgout2(:,:) = 0._pumas_r8
    where (qg(:,trop_cloud_top_lev:) >= 1.e-7_pumas_r8)
      dgout2(:,trop_cloud_top_lev:) = avg_diameter( &
           qg(:,trop_cloud_top_lev:), &
           ng(:,trop_cloud_top_lev:) * rho(:,trop_cloud_top_lev:), &
           rho(:,trop_cloud_top_lev:), rhog)
    end where

    degrau(:ncol,:) = 0._kind_phys
    where (qg(:,trop_cloud_top_lev:) >= 1.e-7_pumas_r8)
      degrau(:ncol,trop_cloud_top_lev:) = real(dgout2(:,trop_cloud_top_lev:) * &
           3._pumas_r8 * rhog/rhows, kind_phys)
    end where

  end subroutine pumas_post_main_run

end module pumas_post_main
