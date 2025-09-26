! Save temperature tendency from deep convection
! for use by gravity wave parameterization
!
! This scheme has to be run after all deep convective schemes,
! but before the final tendencies are applied.
module save_ttend_from_convect_deep
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: save_ttend_from_convect_deep_run

contains

!> \section arg_table_save_ttend_from_convect_deep_run Argument Table
!! \htmlinclude save_ttend_from_convect_deep_run.html
  subroutine save_ttend_from_convect_deep_run( &
    ncol, pver, &
    tend_s, cpair, &
    ttend_dp, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol                     ! Number of atmospheric columns [count]
    integer,            intent(in)  :: pver                     ! Number of vertical layers [count]
    real(kind_phys),    intent(in)  :: tend_s(:, :)            ! Enthalpy tendency from deep convection [J kg-1 s-1]
    real(kind_phys),    intent(in)  :: cpair                   ! Specific heat of dry air at constant pressure [J kg-1 K-1]

    ! Output arguments
    real(kind_phys),    intent(out) :: ttend_dp(:, :)          ! Temperature tendency from deep convection [K s-1]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    ! Convert to temperature tendency
    do k = 1, pver
      do i = 1, ncol
        ttend_dp(i, k) = tend_s(i, k) / cpair
      end do
    end do

  end subroutine save_ttend_from_convect_deep_run

end module save_ttend_from_convect_deep
