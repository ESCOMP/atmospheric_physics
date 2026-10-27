! Copy ZM tendencies (core, i.e., before evap and momentum transport)
! to capture T and Q tendencies that MCSP redistributes.
module save_core_tend_from_convect_deep
  implicit none
  private

  public :: save_core_tend_from_convect_deep_run

contains

!> \section arg_table_save_core_tend_from_convect_deep_run Argument Table
!! \htmlinclude save_core_tend_from_convect_deep_run.html
  subroutine save_core_tend_from_convect_deep_run(ncol, pver, cpair, tend_s, qtend, &
                                                  ttend_dp_core, qtend_dp_core, errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: cpair                ! specific heat of dry air at constant pressure [J kg-1 K-1]
    real(kind_phys),    intent(in)  :: tend_s(:,:)          ! deep convective dry static energy tendency [J kg-1 s-1]
    real(kind_phys),    intent(in)  :: qtend(:,:)           ! deep convective water vapor tendency [kg kg-1 s-1]
    real(kind_phys),    intent(out) :: ttend_dp_core(:,:)   ! core deep convective temperature tendency [K s-1]
    real(kind_phys),    intent(out) :: qtend_dp_core(:,:)   ! core deep convective water vapor tendency [kg kg-1 s-1]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: i, k

    errmsg = ''
    errflg = 0

    do k = 1, pver
       do i = 1, ncol
          ttend_dp_core(i,k) = tend_s(i,k) / cpair
          qtend_dp_core(i,k) = qtend(i,k)
       end do
    end do

  end subroutine save_core_tend_from_convect_deep_run

end module save_core_tend_from_convect_deep
