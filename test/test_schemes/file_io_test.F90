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
     use ccpp_io_reader, only: abstract_netcdf_reader_t, create_netcdf_reader_t

     !Input variables:
     character(len=*),   intent(in)  :: file_path

     !Output variables:
     character(len=512), intent(out) :: errmsg
     integer,            intent(out) :: errcode

     !Variables read from file:
     character(len=:), allocatable :: gas_names(:)
     integer,          allocatable :: band2gpt(:,:)
     real(kind_phys),  allocatable :: press_ref(:)
     real(kind_phys),  allocatable :: band_lims_wavenum(:,:)

     class(abstract_netcdf_reader_t), pointer :: reader

     !NetCDF dimensions:
     integer :: bnd
     integer :: absorber

     reader => create_netcdf_reader_t()

     !Initialize output variables:
     errcode = 0
     errmsg  = ''

     ! Open file
     call reader%open_file(file_path, errmsg, errcode)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if


     !Read variables from NetCDF file:
     !-----------------------

     !Attempt to get absorbing gas names from file:
     call reader%get_var('gas_names', gas_names, errmsg, errcode)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Attempt to get reference pressure from file:
     call reader%get_var('press_ref', press_ref, errmsg, errcode)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Attempt to get wavenumber band grid start/end indices from file:
     call reader%get_var('bnd_limits_gpt', band2gpt, errmsg, errcode)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Attempt to get wavenumber band start/end values from file:
     call reader%get_var('bnd_limits_wavenumber', band_lims_wavenum, errmsg, errcode)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     ! Close file
     call reader%close_file(errmsg, errcode)
     if (errcode /= 0) then
        return !Error has occurred, so exit scheme
     end if

     !Write array shape information:
     write(*,*) 'gas_names length', len(gas_names)
     write(*,*) 'gas_names shape ', shape(gas_names)
     write(*,*) 'press_ref shape ', shape(press_ref)
     write(*,*) 'band2gpt shape ', shape(band2gpt)
     write(*,*) 'band_lims_wavenum ', shape(band_lims_wavenum)

     !Write max values to stdout:
     write(*,*) 'First absorbing gas name = ', gas_names(1)
     write(*,*) 'Max RRTMGP reference pressure value = ', maxval(press_ref)
     write(*,*) 'Max RRTMGP band starting wave idx   = ', maxval(band2gpt(1,:))
     write(*,*) 'Max RRTMGP band starting wavenumber = ', maxval(band_lims_wavenum(1,:))

     !Reset file reader:
     deallocate(reader)
     nullify(reader)

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
