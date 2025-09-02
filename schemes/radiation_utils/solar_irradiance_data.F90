!-------------------------------------------------------------------------------
! This module uses the solar irradiance data 
! to provide a spectral scaling factor
! to approximate the spectral distribution of irradiance
! when the radiation scheme might use a different solar source function
!-------------------------------------------------------------------------------
! peverwhee - dependencies = time_coordinate
module solar_irradiance_data
  use cam_time_coord, only: time_coordinate
  use ccpp_kinds,     only: kind_phys

  implicit none
  save

  private
  public :: solar_irradiance_data_register
  public :: solar_irradiance_data_init
  public :: solar_irradiance_data_run

  type(time_coordinate) :: time_coord
  real(kind_phys), allocatable :: ref_tsi
  real(kind_phys), public, protected, allocatable :: sol_etf(:)
  real(kind_phys), public, protected, allocatable :: ssi_ref(:)  ! a reference spectrum constructed from 3 solar cycles of data
  real(kind_phys), allocatable :: irradi(:,:)
  real(kind_phys), allocatable :: irrad_fac(:)
  real(kind_phys), allocatable :: etf_fac(:)
  logical, protected :: has_ref_spectrum = .false.
  logical, protected :: has_tsi = .false.
  logical, protected :: initialized = .false.
  logical, protected :: fixed_scon = .false.

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> \section arg_table_solar_irradiance_data_register Argument Table
!! \htmlinclude solar_irradiance_data_register.html
!!
  subroutine solar_irradiance_data_register(irrad_file_path, nbins, nbinsp, errmsg, errflg)
    use ccpp_io_reader,   only: abstract_netcdf_reader_t, create_netcdf_reader_t
    ! Arguments
    character(len=*), intent(in)  :: irrad_file_path
    integer,          intent(out) :: nbins
    integer,          intent(out) :: nbinsp
    character(len=512), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    real(kind_phys), allocatable :: lambda(:)
    class(abstract_netcdf_reader_t), pointer :: file_reader
    integer, parameter :: missing_variable_error_code = 3

    ! Set error variables
    errmsg = ''
    errflg = 0

    nbins = 0
    nbinsp = 0

    file_reader => create_netcdf_reader_t()

    ! Open the solar irradiance data file
    call file_reader%open_file(irrad_file_path, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Read the wavelengths variable
    call file_reader%get_var('wavelength', lambda, errmsg, errflg)
    if (errflg /= 0 .and. errflg /= missing_variable_error_code) then
       return
    else if (errflg == missing_variable_error_code) then
       ! Check old name (for backward compatibility
       call file_reader%get_var('wvl', lambda, errmsg, errflg)
    end if

    ! Close the solar irradiance file
    call file_reader%close_file(errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    deallocate(file_reader)
    nullify(file_reader)

    if (errflg /= 0) then
       ! Override the errflg, it's ok if there is no wavelength info on file in some scenarios
       errflg = 0
       return
    end if

    ! Set output variables (dimensions)
    nbins = size(lambda)
    nbinsp = nbins + 1

  end subroutine solar_irradiance_data_register
!> \section arg_table_solar_irradiance_data_init Argument Table
!! \htmlinclude solar_irradiance_data_init.html
!!
  subroutine solar_irradiance_data_init(irrad_file_path, solar_data_type, solar_data_ymd, solar_data_tod, solar_const, &
                  solar_heating_spectral_scl, speed_of_light, planck_const, nbins, nbinsp, do_spectral_scaling, has_spectrum, sol_tsi, &
                  wave_end, sol_irrad, errmsg, errflg)
    use ccpp_io_reader,   only: abstract_netcdf_reader_t, create_netcdf_reader_t
    ! Arguments
    character(len=*), intent(in) :: irrad_file_path
    character(len=*), intent(in) :: solar_data_type
    integer, intent(in) :: solar_data_ymd
    integer, intent(in) :: solar_data_tod
    real(kind_phys), intent(in) :: solar_const
    logical, intent(in) :: solar_heating_spectral_scl
    real(kind_phys), intent(in) :: speed_of_light
    real(kind_phys), intent(in) :: planck_const
    logical, intent(out) :: do_spectral_scaling   ! flag to do spectral scaling
    logical, intent(out) :: has_spectrum        ! flag for whether solar input file has irradiance spectrum
    real(kind_phys), intent(out) :: sol_tsi
    real(kind_phys), allocatable, intent(out) :: wave_end(:)
    real(kind_phys), allocatable, intent(out) :: sol_irrad(:)
    integer, intent(in) :: nbins
    integer, intent(in) :: nbinsp
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    logical :: fixed
    integer :: idx
    real(kind_phys), allocatable :: ssi(:,:)
    real(kind_phys), allocatable :: ssi_ref(:)
    real(kind_phys), allocatable :: tsi(:)
    real(kind_phys), allocatable :: lambda(:)
    real(kind_phys), allocatable :: dellam(:)
    integer, allocatable         :: wvl_vid
    class(abstract_netcdf_reader_t), pointer :: file_reader
    integer, parameter :: missing_variable_error_code = 3
    character(len=256) :: alloc_errmsg
    real(kind_phys) :: fac

    ! Set error variables
    errmsg = ''
    errflg = 0

    sol_tsi = -1.0_kind_phys
    fac = 1._kind_phys/(planck_const*speed_of_light)

    has_spectrum = .false.

    if (irrad_file_path /= 'NONE') then
       fixed_scon = .false.
    else
       fixed_scon = .true.
    end if

    if (solar_const>0._kind_phys) then
       sol_tsi = solar_const
    end if

    if ( fixed_scon ) return

    fixed = trim(solar_data_type) == 'FIXED'

    call time_coord%initialize(irrad_file_path, fixed=fixed, fixed_ymd=solar_data_ymd, fixed_tod=solar_data_tod, &
                                force_time_interp=.true., try_dates=.true.)

    file_reader => create_netcdf_reader_t()

    ! Open the solar irradiance data file
    call file_reader%open_file(irrad_file_path, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Check what the file contains
    call file_reader%get_var('ssi', ssi, errmsg, errflg)
    if (errflg /= 0 .and. errflg /= missing_variable_error_code) then
       return
    else if (errflg /= missing_variable_error_code) then
       has_spectrum = .true.
    end if

    call file_reader%get_var('tsi', tsi, errmsg, errflg)
    if (errflg /= 0 .and. errflg /= missing_variable_error_code) then
       return
    else if (errflg /= missing_variable_error_code .and. solar_const < 0._kind_phys) then
       has_tsi = .true.
    end if

    call file_reader%get_var('ssi_ref', ssi_ref, errmsg, errflg)
    if (errflg /= 0 .and. errflg /= missing_variable_error_code) then
       return
    else if (errflg /= missing_variable_error_code) then
       has_ref_spectrum = .true.
    end if

    if (has_ref_spectrum) then
       call file_reader%get_var('tsi_ref', ref_tsi, errmsg, errflg)
       if (errflg /= 0) then
          return
       end if
    end if

    do_spectral_scaling = has_spectrum .and. solar_heating_spectral_scl

    ! Read in data
    if (has_spectrum) then
       call file_reader%get_var('wavelength', lambda, errmsg, errflg)
       if (errflg /= 0 .and. errflg /= missing_variable_error_code) then
          return
       else if (errflg == missing_variable_error_code) then
          ! Check old name (for backward compatibility
          call file_reader%get_var('wvl', lambda, errmsg, errflg)
          if (errflg /= 0) then
             return
          end if
       end if
       call file_reader%get_var('band_width', dellam, errmsg, errflg)
       if (errflg /= 0) then
          return
       end if
    end if

    ! Close the solar irradiance file
    call file_reader%close_file(errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    deallocate(file_reader)
    nullify(file_reader)

    allocate(irrad_fac(nbins), stat=errflg, errmsg=alloc_errmsg)
    if( errflg /= 0 ) then
       write(errmsg,*) 'solar_data_init: failed to allocate irrad_fac; error = ', alloc_errmsg
       return
    end if
    allocate(etf_fac(nbins), stat=errflg, errmsg=alloc_errmsg)
    if( errflg /= 0 ) then
       write(errmsg,*) 'solar_data_init: failed to allocate etf_fac; error = ', alloc_errmsg
       return
    end if

    ! Calculate wavelength ends and convert units
    if ( has_spectrum ) then
       allocate(wave_end(nbins+1), stat=errflg, errmsg=alloc_errmsg)
       if( errflg /= 0 ) then
          write(errmsg,*) 'solar_data_init: failed to allocate wave_end; error = ', alloc_errmsg
          return
       end if
       allocate(sol_etf(nbins), stat=errflg, errmsg=alloc_errmsg)
       if( errflg /= 0 ) then
          write(errmsg,*) 'solar_data_init: failed to allocate sol_etf; error = ', alloc_errmsg
          return
       end if

       wave_end(:nbins)  = lambda(:nbins) - 0.5_kind_phys*dellam(:nbins)
       wave_end(nbins+1) = lambda(nbins)  + 0.5_kind_phys*dellam(nbins)
       do idx = 1,nbins
          irrad_fac(idx) = 1.e-3_kind_phys                  ! mW/m2/nm --> W/m2/nm
          etf_fac(idx)   = 1.e-16_kind_phys*lambda(idx)*fac ! mW/m2/nm --> photons/cm2/sec/nm
       enddo
       if(has_ref_spectrum) then
          ssi_ref = ssi_ref * 1.e-3_kind_phys           ! mW/m2/nm --> W/m2/nm
       endif
    endif

    deallocate(lambda)
    deallocate(dellam)

    ! need to force data loading when the model starts at a time =/ 00:00:00.000
    ! -- may occur in restarts also
    call solar_irradiance_data_run(irrad_file_path, nbins, nbinsp, has_spectrum, do_spectral_scaling, &
            sol_irrad, wave_end, sol_tsi, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    initialized = .true.

  end subroutine solar_irradiance_data_init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!> \section arg_table_solar_irradiance_data_run Argument Table
!! \htmlinclude solar_irradiance_data_run.html
!!
  subroutine solar_irradiance_data_run(irrad_file_path, nbins, nbinsp, has_spectrum, do_spectral_scaling, &
                sol_irrad, sol_tsi, errmsg, errflg)
     use ccpp_io_reader,   only: abstract_netcdf_reader_t, create_netcdf_reader_t
     ! Arguments 
     character(len=*),   intent(in)    :: irrad_file_path
     integer,            intent(in)    :: nbins                 ! number of bins
     integer,            intent(in)    :: nbinsp                ! number of bins plus one
     logical,            intent(in)    :: has_spectrum
     logical,            intent(in)    :: do_spectral_scaling     ! flag to do spectral scaling
     real(kind_phys),    intent(out)   :: sol_tsi               ! total solar irradiance
     real(kind_phys),    intent(out)   :: sol_irrad(:)          ! solar irradiance
     character(len=512), intent(out)   :: errmsg
     integer,            intent(out)   :: errflg

     ! Local variables 
     integer  :: idx, index, nt
     integer  :: offset(2), count(2)
     integer, allocatable :: itsi(:)
     logical  :: read_data
     real(kind_phys) :: data(nbins)
     integer  :: ierr
     real(kind_phys) :: delt
    class(abstract_netcdf_reader_t), pointer :: file_reader

     ! Initialize error variables
     errflg = 0
     errmsg = ''

     if ( fixed_scon ) return
     if ( time_coord%fixed .and. initialized ) return

     index = -1

     read_data = time_coord%read_more() .or. .not.initialized
     call time_coord%advance()

     if ( read_data ) then
        file_reader => create_netcdf_reader_t()

        ! Open the solar irradiance data file
        call file_reader%open_file(irrad_file_path, errmsg, errflg)
        if (errflg /= 0) then
           return
        end if
        nt = 2
        index = time_coord%indxs(1)

        ! get the surrounding time slices
        offset = (/ 1, index /)
        count =  (/ nbins, nt /)

        if (has_spectrum) then
           call file_reader%get_var('ssi', irradi, errmsg, errflg, offset, count)
           if (errflg /= 0) then
              return
           end if
        end if
        if (has_tsi .and. (.not. do_spectral_scaling)) then
           call file_reader%get_var('tsi', itsi, errmsg, errflg, (/index/), (/nt/))
           if (errflg /= 0) then
              return
           end if
           if ( any(itsi(:nt) < 0._kind_phys) ) then
              write(errmsg,*) 'solar_data_advance: invalid or missing tsi data'
              errflg = 1
              return
           end if
        end if
        ! Close the solar irradiance file
        call file_reader%close_file(errmsg, errflg)
        if (errflg /= 0) then
           return
        end if
        deallocate(file_reader)
        nullify(file_reader)
     end if

     delt = time_coord%wghts(2)

     if (has_spectrum) then
        data(:) = irradi(:,1) + delt*( irradi(:,2) - irradi(:,1) )

        do idx = 1,nbins
           sol_irrad(idx) = data(idx)*irrad_fac(idx) ! W/m2/nm
           sol_etf(idx)   = data(idx)*etf_fac(idx)   ! photons/cm2/sec/nm
        end do
     end if
     if (has_tsi .and. (.not.do_spectral_scaling)) then
        sol_tsi = itsi(1) + delt*( itsi(2) - itsi(1) )
     end if

     if (has_spectrum) then
        deallocate(irradi)
     end if
    
  end subroutine solar_irradiance_data_run

end module solar_irradiance_data
