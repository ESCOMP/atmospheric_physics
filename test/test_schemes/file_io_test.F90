module file_io_test

!This test scheme exists to exercise the
!abstract File I/O interface designed for
!use by host-agnostic CCPP schemes.

!The test code below use input files
!needed by the RRTMGP radiation scheme.

implicit none
private

public :: file_io_test_init
public :: file_io_test_run

contains

!> \section arg_table_file_io_test_init  Argument Table
!! \htmlinclude file_io_test_init.html
  subroutine file_io_test_init(file_path, errmsg, errcode)

     !Portable use statements:
     use ccpp_kinds, only: kind_phys

     !Non-portable (CAM-SIMA specific) use statements:
     use ioFileMod,     only: getfil
     use cam_pio_utils, only: cam_pio_openfile
     use pio,           only: file_desc_t
     use pio,           only: PIO_NOWRITE
     use pio,           only: PIO_BCAST_ERROR
     use pio,           only: PIO_NOERR
     use pio,           only: pio_inq_dimid
     use pio,           only: pio_inq_dimlen
     use pio,           only: pio_inq_varid
     use pio,           only: pio_get_var

     !Input variables:
     character(len=*),   intent(in)  :: file_path

     !Output variables:
     character(len=512), intent(out) :: errmsg
     integer,            intent(out) :: errcode

     !Variables read from file:
     real(kind_phys), dimension(:), allocatable :: press_ref

     !NetCDF dimensions:
     integer :: pressure

     !Local variables:
     type(file_desc_t)  :: fh              ! pio file handle
     character(len=512) :: local_file_path ! path to file on local storage
     integer :: dim_id                     ! NetCDF dimension ID
     integer :: var_id                     ! NetCDF variable ID
     integer :: ierr

     !Initialize output variables:
     errcode = 0
     errmsg  = ''

     ! Open file
     call getfil(file_path, local_file_path, 0)
     call cam_pio_openfile(fh, local_file_path, PIO_NOWRITE)

     call pio_seterrorhandling(fh, PIO_BCAST_ERROR)

     !Get relevant dimensions:
     ierr = pio_inq_dimid(fh, 'pressure', dim_id)
     if(ierr /= PIO_NOERR) then
        errcode = 1
        errmsg = "Failed to find 'pressure' in '"//local_file_path//"'"
     end if
     ierr = pio_inq_dimlen(fh, dim_id, pressure)
     if(ierr /= PIO_NOERR) then
        errcode = 1
        errmsg = "Failed to read 'pressure' dimension in '"//local_file_path//"'"
     end if

     !Allocate RRTMGP reference pressure:
     allocate(press_ref(pressure), stat=errcode)
     if(errcode /= 0) then
        errmsg = "Failed to allocate 'pres_ref(pressure)'"
        return
     end if

     !Look for reference pressure on file:
     ierr = pio_inq_varid(fh, 'press_ref', var_id)
     if (ierr /= PIO_NOERR) then
        errcode = 1
        errmsg  = "Failed to find 'press_ref' in '"//local_file_path//"'"
        return
     end if

     !Read reference pressure from file:
     ierr = pio_get_var(fh, var_id, press_ref)
     if (ierr /= PIO_NOERR) then
        errcode = 1
        errmsg  = "Failed to read 'press_ref' in '"//local_file_path//"'"
     end if

     !Write max pressure value to stdout:
     write(*,*) 'Max RRTMGP reference pressure value = ', max(press_ref)

  end subroutine file_io_test_init

!> \section arg_table_file_io_test_run  Argument Table
!! \htmlinclude file_io_test_run.html
  subroutine file_io_test_run(errmsg, errcode)
     character(len=512), intent(out) :: errmsg
     integer,            intent(out) :: errcode

     errcode = 0
     errmsg = ''

  end subroutine file_io_test_run
end module file_io_test
