! Prepare flux variables for the energy checker after the RK stratiform block
! (cloud sedimentation, detrainment of convective condensate, prognostic
! cloud water).
! The reserved convective liquid (rliq) was counted by the convection schemes
! as a flux out of the column (zm_prepare_flux_for_check_energy adds it to
! their precipitation) and re-enters the column here as cloud liquid through
! the detrainment tendency, so it is removed from the water flux the energy
! check sees. The large-scale precipitation rate itself stays physical.
! The resulting net liquid and (lwe) ice fluxes are provided to the
! check_energy_chng CCPPized scheme.
module rk_stratiform_prepare_flux_for_check_energy
  implicit none
  private

  public :: rk_stratiform_prepare_flux_for_check_energy_run

contains

!> \section arg_table_rk_stratiform_prepare_flux_for_check_energy_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_prepare_flux_for_check_energy_run.html
  subroutine rk_stratiform_prepare_flux_for_check_energy_run( &
    ncol, &
    prec_str, rliq, snow_str, &
    scheme_name, &
    flx_cnd, flx_ice, &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    ! Input arguments
    integer,            intent(in)  :: ncol
    real(kind_phys),    intent(in)  :: prec_str(:)            ! Large-scale precipitation rate [m s-1]
    real(kind_phys),    intent(in)  :: rliq(:)                ! Reserved liquid water tendency [m s-1]
    real(kind_phys),    intent(in)  :: snow_str(:)            ! Large-scale snow and cloud ice precipitation rate [m s-1]

    ! Output arguments
    character(len=64),  intent(out) :: scheme_name            ! Scheme name for energy checking
    real(kind_phys),    intent(out) :: flx_cnd(:)             ! Liquid and ice fluxes [m s-1]
    real(kind_phys),    intent(out) :: flx_ice(:)             ! Ice fluxes [m s-1]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Set scheme name for energy checking (CAM's name for this check)
    scheme_name = "cldwat_tend"

    ! Water leaving the column as large-scale precipitation, less the reserved
    ! convective liquid that entered it through detrainment
    flx_cnd(:ncol) = prec_str(:ncol) - rliq(:ncol)

    ! Frozen precipitation represents net ice fluxes through column boundaries
    flx_ice(:ncol) = snow_str(:ncol)

  end subroutine rk_stratiform_prepare_flux_for_check_energy_run

end module rk_stratiform_prepare_flux_for_check_energy
