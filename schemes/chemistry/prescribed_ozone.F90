! read and provide prescribed ozone for radiation
! Original Author: Francis Vitt
! CCPP version: Haipeng Lin, October 2025
module prescribed_ozone
  use ccpp_kinds,  only: kind_phys

  ! CAM-SIMA host model dependency to read chemistry data.
  use tracer_data, only: trfile     ! data information and file read state.
  use tracer_data, only: trfld      ! tracer data container.

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: prescribed_ozone_init
  public :: prescribed_ozone_run

  ! fields to store tracer_data state and information.
  type(trfld), pointer :: tracer_data_fields(:)
  type(trfile)         :: tracer_data_file

  ! parameters for ozone
  ! TODO - this should use a centralized chemistry species database
  real(kind_phys),  parameter :: ozone_mw = 47.9981995_kind_phys

  ! namelist options
  logical                     :: has_prescribed_ozone = .false.
  character(len=8), parameter :: ozone_name = 'ozone' ! name of the output field

  character(len=16)    :: fld_name = 'ozone'
  character(len=256)   :: filename = ' '
  character(len=256)   :: filelist = ' '
  character(len=256)   :: datapath = ' '
  character(len=32)    :: data_type = 'SERIAL'
  integer              :: cycle_yr  = 0
  integer              :: fixed_ymd = 0
  integer              :: fixed_tod = 0

contains

!> \section arg_table_prescribed_ozone_init  Argument Table
!! \htmlinclude prescribed_ozone_init.html
  subroutine prescribed_ozone_init( &
    amIRoot, iulog, &
    fld_name_nl, filename_nl, filelist_nl, datapath_nl, &
    data_type_nl, &
    cycle_yr_nl, fixed_ymd_nl, fixed_tod_nl, &
    errmsg, errflg)

    use cam_history, only: history_add_field
    use tracer_data, only: trcdata_init

    logical,            intent(in)  :: amIRoot                   ! MPI root flag
    integer,            intent(in)  :: iulog                     ! log output unit
    character(len=*),   intent(in)  :: fld_name_nl               ! field name from namelist
    character(len=*),   intent(in)  :: filename_nl               ! input filename from namelist
    character(len=*),   intent(in)  :: filelist_nl               ! input filelist from namelist
    character(len=*),   intent(in)  :: datapath_nl               ! input datapath from namelist
    character(len=*),   intent(in)  :: data_type_nl               ! data type from namelist
    integer,            intent(in)  :: cycle_yr_nl               ! cycle year from namelist
    integer,            intent(in)  :: fixed_ymd_nl              ! fixed year-month-day from namelist (YYYYMMDD) [1]
    integer,            intent(in)  :: fixed_tod_nl              ! fixed time of day from namelist [s]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    character(len=32) :: tracer_data_specifier(1)

    errmsg = ''
    errflg = 0

    fld_name = fld_name_nl
    filename = filename_nl
    filelist = filelist_nl
    datapath = datapath_nl
    data_type = data_type_nl
    cycle_yr  = cycle_yr_nl
    fixed_ymd = fixed_ymd_nl
    fixed_tod = fixed_tod_nl

    ! check if user has specified an input dataset
    if(filename /= 'UNSET' .and. len_trim(filename) > 0) then
      has_prescribed_ozone = .true.

      if(amIRoot) then
        write(iulog,*) 'prescribed_ozone_init: ozone is prescribed in: ' // trim(filename)
      endif
    else
      return
    endif

    ! initialize dataset in tracer_data module.
    ! construct field specifier - one field
    ! format is (internal field name):(netCDF name)
    ! the latter is namelist-configurable together with the file source.
    tracer_data_specifier(1) = trim(ozone_name)//':'//trim(fld_name)

    call trcdata_init( &
      specifier      = tracer_data_specifier(:), &
      filename       = filename, &
      filelist       = filelist, &
      datapath       = datapath, &
      flds           = tracer_data_fields, & ! ptr
      file           = tracer_data_file, &
      data_cycle_yr  = cycle_yr, &
      data_fixed_ymd = fixed_ymd, &
      data_fixed_tod = fixed_tod, &
      data_type      = data_type)

    ! add history field for diagnostic purposes
    call history_add_field('ozone', 'prescribed_ozone', 'lev', 'inst', 'mol mol-1')

  end subroutine prescribed_ozone_init

!> \section arg_table_prescribed_ozone_run  Argument Table
!! \htmlinclude prescribed_ozone_run.html
  subroutine prescribed_ozone_run( &
    ncol, pver, &
    mwdry, boltz, &
    t, pmiddry, &
    prescribed_ozone, &
    errmsg, errflg)

    use tracer_data, only: advance_trcdata
    use cam_history, only: history_out_field

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: mwdry                     ! molecular_weight_of_dry_air [g mol-1]
    real(kind_phys),    intent(in)  :: boltz                     ! boltzmann_constant [J K-1]
    real(kind_phys),    intent(in)  :: t(:,:)                    ! temperature [K]
    real(kind_phys),    intent(in)  :: pmiddry(:,:)              ! dry air pressure [Pa]
    real(kind_phys),    intent(out) :: prescribed_ozone(:,:)     ! prescribed ozone mass mixing ratio [kg kg-1 dry]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! conversion factor to mass mixing ratio (kg kg-1 dry)
    real(kind_phys) :: to_mmr(ncol, pver)

    ! units from file
    character(len=32) :: units_str

    errmsg = ''
    errflg = 0

    if(.not. has_prescribed_ozone) then
      return
    endif

    ! advance data in tracer_data to current time.
    call advance_trcdata(tracer_data_fields, tracer_data_file)
    units_str = trim(tracer_data_fields(1)%units)

    ! copy field from tracer_data container.
    prescribed_ozone(:ncol,:pver) = tracer_data_fields(1)%data(:ncol, :pver)

    ! unit conversions needed for this field to be in kg kg-1 dry?
    select case(units_str)
      case('kg/kg', 'mmr')
        to_mmr = 1._kind_phys
      case('mol/mol', 'mole/mole', 'vmr', 'fraction')
        to_mmr = ozone_mw/mwdry
      case('molec/cm3', '/cm3', 'molecules/cm3', 'cm^-3', 'cm**-3')
        to_mmr(:ncol,:pver) = (ozone_mw*1e6_kind_phys*boltz*t(:ncol,:pver))/&
                              (mwdry*pmiddry(:ncol,:pver))
      case default
        errflg = 1
        errmsg = 'prescribed_ozone_run: unit' // units_str //' are not recognized'
    end select

    ! convert to mol mol-1 (dry) only for diagnostic output
    call history_out_field('ozone', prescribed_ozone(:ncol,:pver)*(mwdry/ozone_mw))

  end subroutine prescribed_ozone_run

end module prescribed_ozone
