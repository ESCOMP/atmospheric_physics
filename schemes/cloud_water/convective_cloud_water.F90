! Grid-box mean total (stratiform + convective) cloud liquid and ice for radiation:
! blends deep and shallow convective in-cloud condensate with the stratiform condensate.
module convective_cloud_water
  implicit none
  private

  public :: convective_cloud_water_run

  ! conv_water_in_rad options: how convective in-cloud water enters the radiation inputs
  integer, parameter :: conv_water_none = 0      ! stratiform condensate only
  integer, parameter :: conv_water_arith_avg = 1 ! area-weighted arithmetic average
  integer, parameter :: conv_water_emis_avg = 2  ! area-weighted average in emissivity

contains

!> \section arg_table_convective_cloud_water_run Argument Table
!! \htmlinclude convective_cloud_water_run.html
  subroutine convective_cloud_water_run( &
    ncol, pver, &
    conv_water_in_rad, frac_limit, &
    one_mom_clouds, &
    gravit, &
    pdel, ls_liq, ls_ice, &
    sh_icwmr, dp_icwmr, &
    sh_frac, dp_frac, ast, rei, &
    totg_liq, totg_ice, &
    conv_liq, conv_ice, tot_liq, tot_ice, &
    totg_liq_sh, totg_liq_dp, totg_ice_sh, totg_ice_dp, &
    fresh, fredp, frecu, fretot, &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    ! 0==> No; 1==> Yes-Arithmetic average;
    ! 2==> Yes-Average in emissivity.
    integer,         intent(in)  :: conv_water_in_rad
    real(kind_phys), intent(in)  :: frac_limit                  ! minimum cloud fraction [fraction]
    ! one-moment (RK) microphysics: use rei unclamped in the ice absorptivity;
    ! two-moment clamps rei to [13,130] micron
    logical,         intent(in)  :: one_mom_clouds
    real(kind_phys), intent(in)  :: gravit                      ! gravitational acceleration [m s-2]

    real(kind_phys), intent(in)  :: pdel(:, :)                  ! (ncol,pver) moist pressure difference across layer [Pa]
    real(kind_phys), intent(in)  :: ls_liq(:, :)                ! (ncol,pver) large-scale (stratiform) contributions to GBA cloud liq [kg kg-1]
    real(kind_phys), intent(in)  :: ls_ice(:, :)                ! (ncol,pver) large-scale (stratiform) contributions to GBA cloud ice [kg kg-1]

    real(kind_phys), intent(in)  :: sh_icwmr(:, :)              ! (ncol,pver) shallow conv. cloud water [kg kg-1]
    real(kind_phys), intent(in)  :: dp_icwmr(:, :)              ! (ncol,pver) deep conv. cloud water [kg kg-1]
    real(kind_phys), intent(in)  :: sh_frac(:, :)               ! (ncol,pver) shallow convective cloud fraction [fraction]
    real(kind_phys), intent(in)  :: dp_frac(:, :)               ! (ncol,pver) deep convective cloud fraction [fraction]
    real(kind_phys), intent(in)  :: ast(:, :)                   ! (ncol,pver) physical liquid+ice stratus cloud fraction [fraction]
    real(kind_phys), intent(in)  :: rei(:, :)                   ! (ncol,pver) ice crystal effective radius [micron]

    real(kind_phys), intent(out) :: totg_liq(:, :)              ! (ncol,pver) grid box total cloud liquid mixing ratio [kg kg-1]
    real(kind_phys), intent(out) :: totg_ice(:, :)              ! (ncol,pver) grid box total cloud ice mixing ratio [kg kg-1]

    ! Diagnostic-only outputs:
    real(kind_phys), intent(out) :: conv_liq(:, :)              ! (ncol,pver) convective contributions to IC cloud liquid [kg kg-1]
    real(kind_phys), intent(out) :: conv_ice(:, :)              ! (ncol,pver) convective contributions to IC cloud ice [kg kg-1]
    real(kind_phys), intent(out) :: tot_liq(:, :)               ! (ncol,pver) total IC liquid [kg kg-1]
    real(kind_phys), intent(out) :: tot_ice(:, :)               ! (ncol,pver) total IC ice [kg kg-1]
    real(kind_phys), intent(out) :: totg_liq_sh(:, :)           ! (ncol,pver) grid-mean LWP from shallow convective cloud [kg kg-1]
    real(kind_phys), intent(out) :: totg_liq_dp(:, :)           ! (ncol,pver) grid-mean LWP from deep convective cloud [kg kg-1]
    real(kind_phys), intent(out) :: totg_ice_sh(:, :)           ! (ncol,pver) grid-mean IWP from shallow convective cloud [kg kg-1]
    real(kind_phys), intent(out) :: totg_ice_dp(:, :)           ! (ncol,pver) grid-mean IWP from deep convective cloud [kg kg-1]
    real(kind_phys), intent(out) :: fresh(:, :)                 ! (ncol,pver) fractional occurrence of shallow cumulus [1]
    real(kind_phys), intent(out) :: fredp(:, :)                 ! (ncol,pver) fractional occurrence of deep cumulus [1]
    real(kind_phys), intent(out) :: frecu(:, :)                 ! (ncol,pver) fractional occurrence of cumulus [1]
    real(kind_phys), intent(out) :: fretot(:, :)                ! (ncol,pver) fractional occurrence of cloud [1]

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables:
    integer         :: i, k                                          ! Lon, lev indices buff stuff.
    real(kind_phys) :: cu_icwmr                                      ! Convective  water for this grid-box.
    real(kind_phys) :: ls_icwmr                                      ! Large-scale water for this grid-box.
    real(kind_phys) :: tot_icwmr                                     ! Large-scale water for this grid-box.
    real(kind_phys) :: ls_frac                                       ! Large-scale cloud frac for this grid-box.
    real(kind_phys) :: tot0_frac, cu0_frac, dp0_frac, sh0_frac
    real(kind_phys) :: kabs, kabsi, alpha, dp0, sh0
    real(kind_phys) :: ls_ice_mass_frac                              ! Ice mass fraction of the large-scale condensate.

    real(kind_phys), parameter   :: kabsl = 0.090361_kind_phys
    real(kind_phys), parameter   :: ic_limit = 1.e-12_kind_phys

    errmsg = ''
    errflg = 0

    ! The below computations only write output if conv_water_in_rad = 1, 2
    ! For conv_water_in_rad = 0, assign zeros. Otherwise abort.
    if (conv_water_in_rad /= conv_water_arith_avg .and. conv_water_in_rad /= conv_water_emis_avg) then
      totg_liq(:, :) = 0._kind_phys
      totg_ice(:, :) = 0._kind_phys
      conv_liq(:, :) = 0._kind_phys
      conv_ice(:, :) = 0._kind_phys
      tot_liq(:, :) = 0._kind_phys
      tot_ice(:, :) = 0._kind_phys
      totg_liq_sh(:, :) = 0._kind_phys
      totg_liq_dp(:, :) = 0._kind_phys
      totg_ice_sh(:, :) = 0._kind_phys
      totg_ice_dp(:, :) = 0._kind_phys
      fresh(:, :) = 0._kind_phys
      fredp(:, :) = 0._kind_phys
      frecu(:, :) = 0._kind_phys
      fretot(:, :) = 0._kind_phys
      if (conv_water_in_rad /= conv_water_none) then
        errflg = 1
        errmsg = 'convective_cloud_water_run: invalid conv_water_in_rad (must be 0, 1, or 2)'
      end if
      return
    end if

    ! --------------------------------------------------------------- !
    ! Loop through grid-boxes and determine:                          !
    ! 1. Effective mean in-cloud convective ice/liquid (deep+shallow) !
    ! 2. Effective mean in-cloud total ice/liquid (ls+convective)     !
    ! --------------------------------------------------------------- !
    fresh(:, :) = 0._kind_phys
    fredp(:, :) = 0._kind_phys
    frecu(:, :) = 0._kind_phys
    fretot(:, :) = 0._kind_phys

    do k = 1, pver
      do i = 1, ncol

        if (sh_frac(i, k) <= frac_limit .or. sh_icwmr(i, k) <= ic_limit) then
          sh0_frac = 0._kind_phys
        else
          sh0_frac = sh_frac(i, k)
        end if
        if (dp_frac(i, k) <= frac_limit .or. dp_icwmr(i, k) <= ic_limit) then
          dp0_frac = 0._kind_phys
        else
          dp0_frac = dp_frac(i, k)
        end if
        cu0_frac = sh0_frac + dp0_frac

        ! For the moment calculate the emissivity based upon the ls clouds ice fraction
        ls_ice_mass_frac = min(1._kind_phys, max(0._kind_phys, ls_ice(i, k)/(ls_ice(i, k) + ls_liq(i, k) + 1.e-36_kind_phys)))

        if ((cu0_frac < frac_limit) .or. ((sh_icwmr(i, k) + dp_icwmr(i, k)) < ic_limit)) then

          cu0_frac = 0._kind_phys
          cu_icwmr = 0._kind_phys

          ls_frac = ast(i, k)
          if (ls_frac < frac_limit) then
            ls_frac = 0._kind_phys
            ls_icwmr = 0._kind_phys
          else
            ls_icwmr = (ls_liq(i, k) + ls_ice(i, k))/ls_frac ! Convert to IC value.
          end if

          tot0_frac = ls_frac
          tot_icwmr = ls_icwmr

        else

          ! Select radiation constants (effective radii) for emissivity averaging.

          if (one_mom_clouds) then
            kabsi = 0.005_kind_phys + 1._kind_phys/rei(i, k)
          else
            kabsi = 0.005_kind_phys + 1._kind_phys/min(max(13._kind_phys, rei(i, k)), 130._kind_phys)
          end if
          kabs = kabsl*(1._kind_phys - ls_ice_mass_frac) + kabsi*ls_ice_mass_frac
          alpha = -1.66_kind_phys*kabs*pdel(i, k)/gravit*1000.0_kind_phys

          ! Selecting cumulus in-cloud water.

          select case (conv_water_in_rad)
          case (conv_water_arith_avg)
            cu_icwmr = (sh0_frac*sh_icwmr(i, k) + dp0_frac*dp_icwmr(i, k))/max(frac_limit, cu0_frac)
          case (conv_water_emis_avg)
            sh0 = exp(alpha*sh_icwmr(i, k))
            dp0 = exp(alpha*dp_icwmr(i, k))
            cu_icwmr = log((sh0_frac*sh0 + dp0_frac*dp0)/max(frac_limit, cu0_frac))
            cu_icwmr = cu_icwmr/alpha
          end select

          ! Selecting total in-cloud water.
          ! Attribute large-scale/convective area fraction differently from default.

          ls_frac = ast(i, k)
          ls_icwmr = (ls_liq(i, k) + ls_ice(i, k))/max(frac_limit, ls_frac) ! Convert to IC value.
          tot0_frac = (ls_frac + cu0_frac)

          select case (conv_water_in_rad)
          case (conv_water_arith_avg)
            tot_icwmr = (ls_frac*ls_icwmr + cu0_frac*cu_icwmr)/max(frac_limit, tot0_frac)
          case (conv_water_emis_avg)
            tot_icwmr = log((ls_frac*exp(alpha*ls_icwmr) + cu0_frac*exp(alpha*cu_icwmr))/max(frac_limit, tot0_frac))
            tot_icwmr = tot_icwmr/alpha
          end select

        end if

        ! Repartition convective cloud water into liquid and ice phase.
        ! Currently, this partition is made using the ice fraction of stratus condensate.
        ! In future, we should use ice fraction explicitly computed from the convection scheme.

        conv_ice(i, k) = cu_icwmr*ls_ice_mass_frac
        conv_liq(i, k) = cu_icwmr*(1._kind_phys - ls_ice_mass_frac)

        tot_ice(i, k) = tot_icwmr*ls_ice_mass_frac
        tot_liq(i, k) = tot_icwmr*(1._kind_phys - ls_ice_mass_frac)

        totg_ice(i, k) = tot0_frac*tot_icwmr*ls_ice_mass_frac
        totg_liq(i, k) = tot0_frac*tot_icwmr*(1._kind_phys - ls_ice_mass_frac)

        ! Grid-mean convective water
        totg_ice_sh(i, k) = sh0_frac*sh_icwmr(i, k)*ls_ice_mass_frac
        totg_ice_dp(i, k) = dp0_frac*dp_icwmr(i, k)*ls_ice_mass_frac
        totg_liq_sh(i, k) = sh0_frac*sh_icwmr(i, k)*(1._kind_phys - ls_ice_mass_frac)
        totg_liq_dp(i, k) = dp0_frac*dp_icwmr(i, k)*(1._kind_phys - ls_ice_mass_frac)
        if (sh0_frac > frac_limit) then
          fresh(i, k) = 1._kind_phys
        end if
        if (dp0_frac > frac_limit) then
          fredp(i, k) = 1._kind_phys
        end if
        if (cu0_frac > frac_limit) then
          frecu(i, k) = 1._kind_phys
        end if
        if (tot0_frac > frac_limit) then
          fretot(i, k) = 1._kind_phys
        end if

      end do
    end do

  end subroutine convective_cloud_water_run

end module convective_cloud_water
