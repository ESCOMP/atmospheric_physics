module to_be_ccppized_temporary
! This module is a TEMPORARY place to put calls to initialization routines which have not yet
! been CCPP'ized, and the run methods are being called directly in CCPP'ized routines.

! Once a module has been CCPP'ized, then the call in this routine needs to be removed

implicit none

contains

!> \section arg_table_to_be_ccppized_temporary_init Argument Table
!! \htmlinclude to_be_ccppized_temporary_init.html
!!
subroutine to_be_ccppized_temporary_init(errmsg, errflg)

   use wv_saturation, only: wv_sat_init

   character(len=512), intent(out)      :: errmsg
   integer, intent(out)                 :: errflg

   errmsg = ' '
   errflg = 0

   call wv_sat_init()

end subroutine to_be_ccppized_temporary_init

end module to_be_ccppized_temporary

