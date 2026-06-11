! Zero stub for aerosol optical properties.
! Sets all aerosol optics arrays to zero (no aerosol), for use in
! test SDFs where full aerosol optics calculations are not needed.
module aerosol_optics_zero_stub
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: aerosol_optics_zero_stub_run

contains

!> \section arg_table_aerosol_optics_zero_stub_run Argument Table
!! \htmlinclude aerosol_optics_zero_stub_run.html
!!
  subroutine aerosol_optics_zero_stub_run( &
    aer_tau, aer_tau_w, aer_tau_w_g, aer_lw_abs, &
    errmsg, errflg)

    ! Outputs
    real(kind_phys),  intent(out) :: aer_tau(:,:,:)     !< SW extinction OD
    real(kind_phys),  intent(out) :: aer_tau_w(:,:,:)   !< SW ssa*tau
    real(kind_phys),  intent(out) :: aer_tau_w_g(:,:,:) !< SW asy*ssa*tau
    real(kind_phys),  intent(out) :: aer_lw_abs(:,:,:)  !< LW absorption OD

    character(len=*), intent(out) :: errmsg             !< CCPP error message
    integer,          intent(out) :: errflg             !< CCPP error flag

    errmsg = ''
    errflg = 0

    aer_tau(:,:,:)     = 0._kind_phys
    aer_tau_w(:,:,:)   = 0._kind_phys
    aer_tau_w_g(:,:,:) = 0._kind_phys
    aer_lw_abs(:,:,:)  = 0._kind_phys

  end subroutine aerosol_optics_zero_stub_run

end module aerosol_optics_zero_stub
