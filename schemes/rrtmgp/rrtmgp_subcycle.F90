!> This module contains the control (timestep init and timestep final)
!!  for the subcycle iteration
module rrtmgp_subcycle
  implicit none
  private

  public rrtmgp_subcycle_run
contains

!> \section arg_table_rrtmgp_subcycle_init Argument Table
!! \htmlinclude rrtmgp_subcycle_init
!!
   subroutine rrtmgp_subcycle_init(diag_cur, errmsg, errcode)
      integer,            intent(out) :: diag_cur
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errcode

      errmsg = ''
      errflg = 0
      diag_cur = 1
   end subroutine rrtmgp_subcycle_init

!> \section arg_table_rrtmgp_subcycle_run Argument Table
!! \htmlinclude rrtmgp_subcycle_run.html
!!
   subroutine rrtmgp_subcycle_run(diag_cur, num_diag_cycles, errmsg, errcode)
      integer,             intent(in) :: num_diag_cycles
      integer,          intent(inout) :: diag_cur
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errcode

      errmsg = ''
      errflg = 0
      diag_cur = diag_cur + 1
      if (diag_cur > num_diag_cycles) then
         diag_cur = 1
      end if

  end subroutine rrtmgp_subcycle_run
end module rrtmgp_subcycle
