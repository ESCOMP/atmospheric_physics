module read_scheme_namelists
! This routine will be auto-generated

implicit none
private
save

public read_scheme_nl

subroutine read_scheme_nl(nlfile, mpi_communicator, rootprocid, isrootproc)
   use qneg_readnl_mod, only: qneg_readnl

   character(len=*), intent(in) :: nlfile
   integer,          intent(in) :: mpi_communicator
   integer,          intent(in) :: rootprocid
   logical,          intent(in) :: isrootproc

   call qneg_readnl(nlfile, mpi_communicator, rootprocid, isrootproc)

end subroutine read_nl

end module read_scheme_namelists
