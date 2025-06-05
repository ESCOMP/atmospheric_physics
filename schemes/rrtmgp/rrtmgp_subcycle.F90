!> This module contains the control (timestep init and timestep final)
!!  for the subcycle iteration
module rrtmgp_subcycle
  implicit none
  private

  public rrtmgp_subcycle_run
contains

!> \section arg_table_rrtmgp_subcycle_run Argument Table
!! \htmlinclude rrtmgp_subcycle_run.html
!!
   subroutine rrtmgp_subcycle_run(diag_cur, errmsg, errcode)
      integer,          intent(inout) :: diag_cur
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errcode

      diag_cur = diag_cur + 1

  end subroutine rrtmgp_subcycle_run
end module rrtmgp_subcycle
