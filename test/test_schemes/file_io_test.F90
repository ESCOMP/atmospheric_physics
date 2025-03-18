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
     use sima_ccpp_FileIO, only: sima_get_netcdf_var
     use sima_ccpp_FileIO, only: sima_get_netcdf_dim
     use ioFileMod,        only: cam_get_file
     use cam_pio_utils,    only: cam_pio_openfile
     use pio,              only: file_desc_t
     use pio,              only: PIO_NOWRITE
     use pio,              only: PIO_BCAST_ERROR
     use pio,              only: PIO_NOERR
     use pio,              only: pio_inq_dimid
     use pio,              only: pio_inq_dimlen
     use pio,              only: pio_closefile
     use pio,              only: pio_seterrorhandling

     !Input variables:
     character(len=*),   intent(in)  :: file_path

     !Output variables:
     character(len=512), intent(out) :: errmsg
     integer,            intent(out) :: errcode

     !Variables read from file:
     character(32), dimension(:),  allocatable    :: gas_names
     integer,  dimension(:,:), allocatable        :: band2gpt
     real(kind_phys), dimension(:), allocatable   :: press_ref
     real(kind_phys), dimension(:,:), allocatable :: band_lims_wavenum

     !NetCDF dimensions:
     integer :: bnd
     integer :: absorber
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
     call cam_get_file(file_path, local_file_path)
     call cam_pio_openfile(fh, local_file_path, PIO_NOWRITE)

     call pio_seterrorhandling(fh, PIO_BCAST_ERROR)

     !Get relevant dimensions:
     !-----------------------

     !Bands:
     call sima_get_netcdf_dim(fh, local_file_path, 'bnd', errcode, errmsg, dimlen=bnd)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Absorbers:
     call sima_get_netcdf_dim(fh, local_file_path, 'absorber', errcode, errmsg, dimlen=absorber)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Pressure:
     call sima_get_netcdf_dim(fh, local_file_path, 'pressure', errcode, errmsg, dimlen=pressure)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Allocate relevant variables:
     !-----------------------

     !Allocate RRTMGP absorbing gas names:
     allocate(gas_names(absorber), stat=errcode, errmsg=errmsg)
     if(errcode /= 0) then
        return
     end if

     !Allocate RRTMGP reference pressure:
     allocate(press_ref(pressure), stat=errcode, errmsg=errmsg)
     if(errcode /= 0) then
        return
     end if

     !Allocate RRTMGP starting and ending wavenumber band indices:
     allocate(band2gpt(2,bnd), stat=errcode, errmsg=errmsg)
     if(errcode /= 0) then
        return
     end if

     !Allocate RRTMGP starting and ending wavenumber for each band:
     allocate(band_lims_wavenum(2,bnd), stat=errcode, errmsg=errmsg)
     if(errcode /= 0) then
        return
     end if

     !Read variables from NetCDF file:
     !-----------------------

     !Attempt to get absorbing gas names from file:
     call sima_get_netcdf_var(fh, local_file_path, 'gas_names', len(gas_names), gas_names, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Attempt to get reference pressure from file:
     call sima_get_netcdf_var(fh, local_file_path, 'press_ref', press_ref, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Attempt to get wavenumber band grid start/end indices from file:
     call sima_get_netcdf_var(fh, local_file_path, 'bnd_limits_gpt', band2gpt, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Attempt to get wavenumber band start/end values from file:
     call sima_get_netcdf_var(fh, local_file_path, 'bnd_limits_wavenumber', band_lims_wavenum, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     ! Close file
     call pio_closefile(fh)

     !Write dimension lengths to stdout:
     write(*,*) 'Number of gas absorbers = ', absorber
     write(*,*) 'Pressure dimension length = ', pressure
     write(*,*) 'Band (bnd) dimension length = ', bnd

     !Write max values to stdout:
     write(*,*) 'First absorbing gas name = ', gas_names(1)
     write(*,*) 'Max RRTMGP reference pressure value = ', maxval(press_ref)
     write(*,*) 'Max RRTMGP band starting wave idx   = ', maxval(band2gpt(1,:))
     write(*,*) 'Max RRTMGP band starting wavenumber = ', maxval(band_lims_wavenum(1,:))

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
