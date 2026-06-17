! Post-processing interstitial for PUMAS.
module pumas_post_main
  implicit none
  private

  public :: pumas_post_main_run

contains

  !> \section arg_table_pumas_post_main_run Argument Table
  !! \htmlinclude pumas_post_main_run.html
  subroutine pumas_post_main_run(ncol, cldice, numice, strat_cldfrc, effi, errmsg, errcode)
    use ccpp_kinds,        only: kind_phys
    use pumas_kinds,       only: pumas_r8=>kind_r8
    use micro_pumas_utils, only: size_dist_param_basic, mg_ice_props, mincld, qsmall

    integer,                 intent(in)  :: ncol
    real(kind_phys), dimension(:,:), intent(in)  :: cldice      ! updated cloud ice mixing ratio (kg/kg)
    real(kind_phys), dimension(:,:), intent(in)  :: numice      ! updated cloud ice number concentration (kg-1)
    real(kind_phys), dimension(:,:), intent(in)  :: strat_cldfrc! total stratiform cloud area fraction (= ast)
    real(kind_phys), dimension(:,:), intent(out) :: effi       ! ice effective radius (micron)
    character(len=512),          intent(out) :: errmsg
    integer,                 intent(out) :: errcode

    integer        :: k, nlev
    real(pumas_r8) :: icimrst(ncol, size(cldice,2))  ! in-cloud (grid-mean) ice mixing ratio
    real(pumas_r8) :: niic(ncol, size(cldice,2))     ! in-cloud (grid-mean) ice number conc
    real(pumas_r8) :: rei(ncol, size(cldice,2))      ! ice slope param, then effective radius

    errmsg = ' '
    errcode = 0

    nlev = size(cldice, 2)

    ! Ice effective radius is recomputed here from the post-microphysics
    ! grid-mean in-cloud ice, NOT from the raw per-substep value PUMAS returns.
    ! PUMAS uses the total stratiform fraction for ice cloud (icecldf = ast).
    ! Mirrors CAM micro_pumas_cam.F90:2677 (icimrst), :3174 (niic), :3171-3193 (rei).
    icimrst(:,:) = min(real(cldice(:ncol,:), pumas_r8) / &
                max(mincld, real(strat_cldfrc(:ncol,:), pumas_r8)), 0.005_pumas_r8)
    niic(:,:)   =    real(numice(:ncol,:), pumas_r8) / &
                max(mincld, real(strat_cldfrc(:ncol,:), pumas_r8))

    rei(:,:) = 25._pumas_r8
    do k = 1, nlev
      call size_dist_param_basic(mg_ice_props, icimrst(:,k), niic(:,k), rei(:,k), ncol)
    end do

    where (icimrst(:,:) >= qsmall)
      rei(:,:) = 1.5_pumas_r8 / rei(:,:) * 1.e6_pumas_r8
    elsewhere
      rei(:,:) = 25._pumas_r8
    end where

    effi(:ncol,:) = real(rei(:,:), kind_phys)

  end subroutine pumas_post_main_run

end module pumas_post_main
