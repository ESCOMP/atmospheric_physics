! read and provide prescribed ozone for radiation
! this is a simple example of a CCPP scheme using the CAM-SIMA tracer_data utility.
!
! Based on original CAM version from: Francis Vitt
module prescribed_ozone
  use ccpp_kinds,  only: kind_phys

  ! CAM-SIMA host model dependency to read chemistry data.
  use tracer_data, only: trfile     ! data information and file read state.
  use tracer_data, only: trfld      ! tracer data container.

  implicit none
  private

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
  character(len=8), parameter :: ozone_name = 'O3' ! standard name of the output field

contains

!> \section arg_table_prescribed_ozone_init  Argument Table
!! \htmlinclude prescribed_ozone_init.html
  subroutine prescribed_ozone_init( &
    amIRoot, iulog, &
    fld_name, filename, filelist, datapath, &
    data_type, &
    cycle_yr, fixed_ymd, fixed_tod, &
    errmsg, errflg)

    use cam_history, only: history_add_field
    use tracer_data, only: trcdata_init

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog

    ! input fields from namelist to initialize tracer_data
    character(len=*),   intent(in)  :: fld_name   ! field name
    character(len=*),   intent(in)  :: filename   ! input filename
    character(len=*),   intent(in)  :: filelist   ! input filelist
    character(len=*),   intent(in)  :: datapath   ! input datapath
    character(len=*),   intent(in)  :: data_type  ! data type
    integer,            intent(in)  :: cycle_yr   ! cycle year
    integer,            intent(in)  :: fixed_ymd  ! fixed year-month-day (YYYYMMDD) [1]
    integer,            intent(in)  :: fixed_tod  ! fixed time of day [s]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    character(len=32) :: tracer_data_specifier(1)

    errmsg = ''
    errflg = 0

    ! check if user has specified an input dataset
    if(filename /= 'UNSET' .and. len_trim(filename) > 0) then
      has_prescribed_ozone = .true.

      if(amIRoot) then
        write(iulog,*) 'prescribed_ozone_init: ozone is prescribed in: ' // trim(filename)
      end if
    else
      return
    end if

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
    const_props, &
    mwdry, boltz, &
    t, pmiddry, &
    pmid, pint, phis, zi, & ! necessary fields for trcdata read.
    constituents, &
    errmsg, errflg)

    ! host model dependency for tracer_data
    use tracer_data, only: advance_trcdata

    ! host model dependency for history output
    use cam_history, only: history_out_field

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! dependency to get constituent index
    use ccpp_const_utils,          only: ccpp_const_get_idx

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    type(ccpp_constituent_prop_ptr_t), &
                        intent(in)    :: const_props(:)      ! CCPP constituent properties pointer
    real(kind_phys),    intent(in)    :: mwdry               ! molecular_weight_of_dry_air [g mol-1]
    real(kind_phys),    intent(in)    :: boltz               ! boltzmann_constant [J K-1]
    real(kind_phys),    intent(in)    :: t(:,:)              ! air temperature [K]
    real(kind_phys),    intent(in)    :: pmiddry(:,:)        ! dry air pressure [Pa]
    real(kind_phys),    intent(in)    :: pmid(:,:)           ! air pressure [Pa]
    real(kind_phys),    intent(in)    :: pint(:,:)           ! air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)    :: phis(:)             ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zi(:,:)             ! geopotential height above surface, interfaces [m]

    real(kind_phys),    intent(inout) :: constituents(:,:,:) ! constituent array (ncol, pver, pcnst)

    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! conversion factor to mass mixing ratio (kg kg-1 dry)
    real(kind_phys) :: to_mmr(ncol, pver)

    ! prescribed ozone mass mixing ratio [kg kg-1 dry]
    real(kind_phys) :: prescribed_ozone(:,:)

    ! units from file
    character(len=32) :: units_str

    integer :: id_o3

    errmsg = ''
    errflg = 0

    if(.not. has_prescribed_ozone) then
      return
    end if

    ! check for 'O3' constituent where prescribed ozone will be written to
    ! which will be read by radiation.
    call ccpp_const_get_idx(const_props, &
         trim(ozone_name), &
         id_o3, errmsg, errflg)
    if (errflg /= 0) return

    ! could not find the constituent.
    if (id_o3 < 0) then
      return
    end if

    ! advance data in tracer_data to current time.
    call advance_trcdata(tracer_data_fields, tracer_data_file, &
                         pmid, pint, phis, zi)
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
        return
    end select

    ! convert to kg kg-1 (dry)
    prescribed_ozone = to_mmr * prescribed_ozone

    ! write to constituent array
    constituents(:ncol, :pver, id_o3) = prescribed_ozone

    ! convert to mol mol-1 (dry) only for diagnostic output
    call history_out_field('ozone', prescribed_ozone(:ncol,:pver)*(mwdry/ozone_mw))

  end subroutine prescribed_ozone_run

end module prescribed_ozone
