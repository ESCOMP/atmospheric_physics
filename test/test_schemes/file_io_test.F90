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
     use sima_ccpp_FileIO, only: sima_open_netcdf_file
     use sima_ccpp_FileIO, only: sima_close_netcdf_file
     use sima_ccpp_FileIO, only: sima_get_netcdf_var
     use sima_ccpp_FileIO, only: sima_get_netcdf_dim

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
     integer :: file_id                    ! NetCDF file ID provided by host

     !Initialize output variables:
     errcode = 0
     errmsg  = ''

     ! Open file
     call sima_open_netcdf_file(file_path, file_id, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Get relevant dimensions:
     !-----------------------

     !Bands:
     call sima_get_netcdf_dim(file_id, 'bnd', errcode, errmsg, dimlen=bnd)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Absorbers:
     call sima_get_netcdf_dim(file_id, 'absorber', errcode, errmsg, dimlen=absorber)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Pressure:
     call sima_get_netcdf_dim(file_id, 'pressure', errcode, errmsg, dimlen=pressure)
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
     call sima_get_netcdf_var(file_id, 'gas_names', len(gas_names), gas_names, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Attempt to get reference pressure from file:
     call sima_get_netcdf_var(file_id, 'press_ref', press_ref, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Attempt to get wavenumber band grid start/end indices from file:
     call sima_get_netcdf_var(file_id, 'bnd_limits_gpt', band2gpt, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Attempt to get wavenumber band start/end values from file:
     call sima_get_netcdf_var(file_id, 'bnd_limits_wavenumber', band_lims_wavenum, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     ! Close file
     call sima_close_netcdf_file(file_id, errcode, errmsg)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

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
