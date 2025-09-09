! Prepare flux variables for energy checker after ZM deep convection
! Has to run after both zm_convr (provides reserved liquid and reserved ice)
! and zm_conv_evap+interstitial to deep (provides snow_dp and prec_dp) have
! been ran.
! The resulting net liquid (lwe) and ice fluxes are combined and provided to
! the check_energy_chng CCPPized scheme.
module zm_prepare_flux_for_check_energy
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  public :: zm_prepare_flux_for_check_energy_run

contains

!> \section arg_table_zm_prepare_flux_for_check_energy_run Argument Table
!! \htmlinclude zm_prepare_flux_for_check_energy_run.html
  subroutine zm_prepare_flux_for_check_energy_run( &
    ncol, &
    prec_dp, rliq, snow_dp, rice, &
    scheme_name, &
    flx_cnd, flx_ice, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol
    real(kind_phys),    intent(in)  :: prec_dp(:)             ! Deep convection precipitation rate [m s-1]
    real(kind_phys),    intent(in)  :: rliq(:)                ! Reserved liquid water tendency [m s-1]
    real(kind_phys),    intent(in)  :: snow_dp(:)             ! Deep convection frozen precipitation rate [m s-1]
    real(kind_phys),    intent(in)  :: rice(:)                ! Reserved ice tendency [m s-1]

    ! Output arguments
    character(len=16),  intent(out) :: scheme_name            ! Scheme name for energy checking
    real(kind_phys),    intent(out) :: flx_cnd(:)             ! Combined liquid and ice fluxes [m s-1]
    real(kind_phys),    intent(out) :: flx_ice(:)             ! Combined ice fluxes [m s-1]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i

    errmsg = ''
    errflg = 0

    ! Set scheme name for energy checking
    scheme_name = "convect_deep_zm"

    ! Combine precipitation and reserved liquid for total liquid+ice flux
    ! This represents net liquid and ice fluxes through column boundaries
    flx_cnd(:ncol) = prec_dp(:ncol) + rliq(:ncol)

    ! Combine frozen precipitation and reserved ice for total ice flux
    ! This represents net ice fluxes through column boundaries
    flx_ice(:ncol) = snow_dp(:ncol) + rice(:ncol)

  end subroutine zm_prepare_flux_for_check_energy_run

end module zm_prepare_flux_for_check_energy
