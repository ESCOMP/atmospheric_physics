module qneg_readnl_mod

! Module for reading qneg namelist variables
! Will be auto-generated in the future

  use cam_logfile,         only: iulog
  use cam_abortutils,      only: endrun

  implicit none
  private
  save

  ! Public interface.

  public :: qneg_readnl

  ! Protected module variables
  character(len=8), public, protected :: print_qneg_warn

contains
  subroutine qneg_readnl(nlfile, mpi_communicator, rootprocid, isrootproc)
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    ! File containing namelist input.
    character(len=*), intent(in) :: nlfile
    integer,          intent(in) :: mpi_communicator
    integer,          intent(in) :: rootprocid
    logical,          intent(in) :: isrootproc

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: sub = 'qneg_readnl'

    !> \section qneg_readnl Argument Table
    !! \htmlinclude qneg_readnl.html
    namelist /qneg_nl/ print_qneg_warn

    print_qneg_warn = ''

    if (isrootproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'qneg_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, qneg_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(sub // ':: ERROR reading namelist qneg_nl')
          end if
       end if
       close(unitn)
    end if

    call mpi_bcast(print_qneg_warn, len(print_qneg_warn), mpi_character, mstrid, mpi_communicator, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: print_qneg_warn")

    if (isrootproc) then
       if (trim(print_qneg_warn) == "timestep") then
          write(iulog, *) sub, ": QNEG statistics will be collected and printed for each timestep"
       else if (trim(print_qneg_warn) == "summary") then
          write(iulog, *) sub, ": QNEG statistics will be collected and printed at the end of the run"
       else if (trim(print_qneg_warn) == "off") then
          write(iulog, *) sub, ": QNEG statistics will not be collected"
       else
          write(iulog,*) sub, ": invalid value for print_qneg_warn"
       end if
    end if

  end subroutine qneg_readnl

end module qneg_readnl_mod
