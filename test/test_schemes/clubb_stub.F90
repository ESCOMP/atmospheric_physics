module clubb_stub
! This stub version of CLUBB is for testing purposes to initialize values needed by other routines

implicit none
private

public :: clubb_stub_register

contains

!> \section arg_table_clubb_stub_register  Argument Table
!! \htmlinclude clubb_stub_register.html
subroutine clubb_stub_register(do_clubb_sgs, errmsg, errcode)
    ! Dummy variables
    logical, intent(out)            :: do_clubb_sgs
    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errcode

    do_clubb_sgs = .true.
    errcode = 0
    errmsg = ''
  end subroutine clubb_stub_register
end module clubb_stub
