! In-cloud water paths consumed by the RRTMGP cloud optics.
! - in-cloud liquid/ice water paths (iclwp/iciwp), optionally
!   accounting for convective in-cloud water via convective_cloud_water.
! - snow/graupel-adjusted cloud fractions and in-cloud snow/graupel
!   water paths (cldfsnow/icswp, cldfgrau/icgrauwp).
module cloud_water_paths
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: cloud_water_paths_run

  ! Minimum cloud fraction floor for in-cloud conversions [fraction]
  real(kind_phys), parameter :: mincld = 0.0001_kind_phys
  ! Graupel uses a larger cloud fraction floor [fraction]
  ! than the other hydrometeors when converting to an in-cloud graupel water path
  real(kind_phys), parameter :: mincld_grau = 1.e-2_kind_phys

contains

!> \section arg_table_cloud_water_paths_run Argument Table
!! \htmlinclude cloud_water_paths_run.html
  subroutine cloud_water_paths_run( &
    ncol, pver, top_lev, &
    conv_water_in_rad, &
    gravit, rair, &
    pmid, t, pdel, &
    ls_liq, ls_ice, &
    cld, concld, &
    totg_liq, totg_ice, &
    qsout, qgout, &
    iclwp, iciwp, &
    cldfsnow, icswp, &
    cldfgrau, icgrauwp, &
    icimr, icwmr, iwc, lwc, &
    errmsg, errflg)

    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    integer,         intent(in)  :: top_lev           ! troposphere cloud physics top level [index]
    integer,         intent(in)  :: conv_water_in_rad ! control variable, with zero meaning no convective water used [1]

    real(kind_phys), intent(in)  :: gravit            ! gravitational acceleration [m s-2]
    real(kind_phys), intent(in)  :: rair              ! dry air gas constant [J kg-1 K-1]

    real(kind_phys), intent(in)  :: pmid(:, :)        ! (ncol,pver) midpoint pressure [Pa]
    real(kind_phys), intent(in)  :: t(:, :)           ! (ncol,pver) air temperature [K]
    real(kind_phys), intent(in)  :: pdel(:, :)        ! (ncol,pver) pressure thickness [Pa]
    real(kind_phys), intent(in)  :: ls_liq(:, :)      ! (ncol,pver) stratiform cloud liquid (CLDLIQ) [kg kg-1]
    real(kind_phys), intent(in)  :: ls_ice(:, :)      ! (ncol,pver) stratiform cloud ice (CLDICE) [kg kg-1]
    real(kind_phys), intent(in)  :: cld(:, :)         ! (ncol,pver) total cloud fraction [fraction]
    real(kind_phys), intent(in)  :: concld(:, :)      ! (ncol,pver) convective cloud fraction [fraction]
    real(kind_phys), intent(in)  :: totg_liq(:, :)    ! (ncol,pver) grid box total cloud liquid mixing ratio [kg kg-1]
    real(kind_phys), intent(in)  :: totg_ice(:, :)    ! (ncol,pver) grid box total cloud ice mixing ratio [kg kg-1]
    real(kind_phys), intent(in)  :: qsout(:, :)       ! (ncol,pver) grid box snow mixing ratio from microphysics [kg kg-1]
    real(kind_phys), intent(in)  :: qgout(:, :)       ! (ncol,pver) grid box graupel mixing ratio from microphysics [kg kg-1]

    real(kind_phys), intent(out) :: iclwp(:, :)       ! (ncol,pver) in-cloud cloud liquid water path [kg m-2]
    real(kind_phys), intent(out) :: iciwp(:, :)       ! (ncol,pver) in-cloud cloud ice water path [kg m-2]
    real(kind_phys), intent(out) :: cldfsnow(:, :)    ! (ncol,pver) cloud fraction adjusted for snow [fraction]
    real(kind_phys), intent(out) :: icswp(:, :)       ! (ncol,pver) in-cloud snow water path [kg m-2]
    real(kind_phys), intent(out) :: cldfgrau(:, :)    ! (ncol,pver) cloud fraction adjusted for graupel [fraction]
    real(kind_phys), intent(out) :: icgrauwp(:, :)    ! (ncol,pver) in-cloud graupel water path [kg m-2]

    ! Diagnostic-only outputs (CAM history: ICIMR/ICWMR/IWC/LWC)
    real(kind_phys), intent(out) :: icimr(:, :)       ! (ncol,pver) in-cloud ice mixing ratio [kg kg-1]
    real(kind_phys), intent(out) :: icwmr(:, :)       ! (ncol,pver) in-cloud water mixing ratio [kg kg-1]
    real(kind_phys), intent(out) :: iwc(:, :)         ! (ncol,pver) grid box average ice water content [kg m-3]
    real(kind_phys), intent(out) :: lwc(:, :)         ! (ncol,pver) grid box average liquid water content [kg m-3]

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    real(kind_phys) :: allcld_ice(ncol, pver) ! all-cloud (stratiform + convective) ice
    real(kind_phys) :: allcld_liq(ncol, pver) ! all-cloud (stratiform + convective) liquid
    integer         :: i, k

    errmsg = ''
    errflg = 0

    ! ==============================================================
    ! Part A: in-cloud liquid/ice water paths for radiation.
    ! Source: CAM cloud_diagnostics_calc, two-moment branch.
    ! ==============================================================

    ! Adjust in-cloud water values to take account of convective
    ! in-cloud water. It is used to calculate the values of
    ! iclwp and iciwp to pass to the radiation.
    if (conv_water_in_rad /= 0) then
      allcld_ice(:ncol, :) = totg_ice(:ncol, :) ! Grid-avg all cloud liquid
      allcld_liq(:ncol, :) = totg_liq(:ncol, :) ! Grid-avg all cloud ice
    else
      allcld_liq(:ncol, :) = 0._kind_phys
      allcld_ice(:ncol, :) = 0._kind_phys
      allcld_liq(:ncol, top_lev:pver) = ls_liq(:ncol, top_lev:pver)  ! Grid-avg all cloud liquid
      allcld_ice(:ncol, top_lev:pver) = ls_ice(:ncol, top_lev:pver)  ! Grid-avg all cloud ice
    end if

    ! Compute in cloud ice and liquid mixing ratios
    ! Note that 'iclwp, iciwp' are used for radiation computation.
    iciwp = 0._kind_phys
    iclwp = 0._kind_phys
    icimr = 0._kind_phys
    icwmr = 0._kind_phys
    iwc = 0._kind_phys
    lwc = 0._kind_phys

    do k = top_lev, pver
      do i = 1, ncol
        ! Limits for in-cloud mixing ratios consistent with MG microphysics
        ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
        icimr(i, k) = min(allcld_ice(i, k)/max(mincld, cld(i, k)), 0.005_kind_phys)
        icwmr(i, k) = min(allcld_liq(i, k)/max(mincld, cld(i, k)), 0.005_kind_phys)
        iwc(i, k) = allcld_ice(i, k)*pmid(i, k)/(rair*t(i, k))
        lwc(i, k) = allcld_liq(i, k)*pmid(i, k)/(rair*t(i, k))
        ! Calculate total cloud water paths in each layer
        iciwp(i, k) = icimr(i, k)*pdel(i, k)/gravit
        iclwp(i, k) = icwmr(i, k)*pdel(i, k)/gravit
      end do
    end do

    ! ==============================================================
    ! Part B: snow/graupel-adjusted cloud fractions and water paths.
    ! Source: CAM micro_pumas_cam in-cloud diagnostics block.
    ! ==============================================================
    icswp = 0._kind_phys
    cldfsnow = 0._kind_phys
    icgrauwp = 0._kind_phys
    cldfgrau = 0._kind_phys

    do k = top_lev, pver
      do i = 1, ncol
        ! ------------------------------
        ! Adjust cloud fraction for snow
        ! ------------------------------
        cldfsnow(i, k) = cld(i, k)
        ! If cloud and only ice (no convective cloud or ice), then set to 0.
        if ((cldfsnow(i, k) > 1.e-4_kind_phys) .and. &
            (concld(i, k) < 1.e-4_kind_phys) .and. &
            (ls_liq(i, k) < 1.e-10_kind_phys)) then
          cldfsnow(i, k) = 0._kind_phys
        end if
        ! If no cloud and snow, then set to 0.25
        if ((cldfsnow(i, k) <= 1.e-4_kind_phys) .and. (qsout(i, k) > 1.e-6_kind_phys)) then
          cldfsnow(i, k) = 0.25_kind_phys
        end if
        ! Calculate in-cloud snow water path
        icswp(i, k) = qsout(i, k)/max(mincld, cldfsnow(i, k))*pdel(i, k)/gravit

        ! ---------------------------------
        ! Adjust cloud fraction for graupel
        ! ---------------------------------
        cldfgrau(i, k) = cld(i, k)
        ! If cloud and only ice (no convective cloud or ice), then set to 0.
        if ((cldfgrau(i, k) > 1.e-4_kind_phys) .and. &
            (concld(i, k) < 1.e-4_kind_phys) .and. &
            (ls_liq(i, k) < 1.e-10_kind_phys)) then
          cldfgrau(i, k) = 0._kind_phys
        end if
        ! If no cloud and graupel, then set to 0.25
        if ((cldfgrau(i, k) <= 1.e-4_kind_phys) .and. (qgout(i, k) > 1.e-9_kind_phys)) then
          cldfgrau(i, k) = 0.25_kind_phys
        end if

        ! Calculate in-cloud graupel water path
        icgrauwp(i, k) = qgout(i, k)/max(mincld_grau, cldfgrau(i, k))*pdel(i, k)/gravit
      end do
    end do

  end subroutine cloud_water_paths_run

end module cloud_water_paths
