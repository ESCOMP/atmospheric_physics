! Manages reading and interpolation of prescribed volcanic aerosol concentrations.
!
! This module uses CCPP constituents (non-advected) to store prescribed volcanic aero
! fields:
! 1) volcanic aerosol mass mixing ratio (from prescribed dataset)
! 2) geometric-mean wet aerosol radius (derived from mass)
!
! Based on original CAM version from: Francis Vitt
module prescribed_volcanic_aerosol
  use ccpp_kinds, only: kind_phys

  ! CAM-SIMA host model dependency to read aerosol data
  use tracer_data, only: trfile     ! data information and file read state
  use tracer_data, only: trfld      ! tracer data container

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: prescribed_volcanic_aerosol_register
  public :: prescribed_volcanic_aerosol_init
  public :: prescribed_volcanic_aerosol_run

  ! fields to store tracer_data state and information
  type(trfld), pointer :: tracer_data_fields(:)
  type(trfile)         :: tracer_data_file

  ! module state variables
  logical :: has_prescribed_volcaero = .false.

  ! Constituent names
  character(len=*), parameter :: volcaero_const_name = 'VOLC_MMR'
  character(len=*), parameter :: volcrad_const_name  = 'VOLC_RAD_GEOM'

  ! Molecular weight of volcanic aerosol species (sulfate) [g mol-1]
  real(kind_phys), parameter :: molmass_volcaero = 47.9981995_kind_phys

  ! WACCM-derived empirical coefficient relating mass concentration
  ! to wet aerosol geometric-mean radius [m (kg m-3)^(-1/3)]
  real(kind_phys), parameter :: radius_conversion = 1.9e-4_kind_phys

  ! TODO: infrastructure for writing (and reading) tracer_data restart information.
  ! see CAM/prescribed_volcrad_aero::{init,read,write}_prescribed_volcaero_restart
  ! !!! Restarts will not be bit-for-bit without this !!!
  ! TODO when SIMA implements restarts.

contains

  ! Register prescribed volcanic aerosol constituents.
!> \section arg_table_prescribed_volcanic_aerosol_register Argument Table
!! \htmlinclude prescribed_volcanic_aerosol_register.html
  subroutine prescribed_volcanic_aerosol_register( &
    amIRoot, iulog, &
    prescribed_volcaero_file, &
    volcaero_constituents, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

    ! Input arguments:
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_volcaero_file          ! input filename from namelist

    ! Output arguments:
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: volcaero_constituents(:)
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    character(len=*), parameter :: subname = 'prescribed_volcanic_aerosol_register'

    errmsg = ''
    errflg = 0

    ! Check if prescribed volcanic aerosols are enabled
    if (prescribed_volcaero_file == 'UNSET' .or. &
        len_trim(prescribed_volcaero_file) == 0) then
      has_prescribed_volcaero = .false.
      if (amIRoot) then
        write(iulog,*) subname//': No prescribed volcanic aerosols specified'
      end if
      return
    end if

    has_prescribed_volcaero = .true.

    ! Register two constituents: aerosol MMR and geometric-mean radius
    allocate(volcaero_constituents(2), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ": " // trim(errmsg)
      return
    end if

    ! (1) Volcanic aerosol dry mass mixing ratio
    call volcaero_constituents(1)%instantiate( &
         std_name          = volcaero_const_name, &
         diag_name         = volcaero_const_name, &
         long_name         = 'prescribed volcanic aerosol dry mass mixing ratio', &
         units             = 'kg kg-1', &
         vertical_dim      = 'vertical_layer_dimension', &
         min_value         = 0.0_kind_phys, &
         advected          = .false., &
         water_species     = .false., &
         mixing_ratio_type = 'dry', &
         errcode           = errflg, &
         errmsg            = errmsg)
    if (errflg /= 0) return

    ! (2) Volcanic aerosol geometric-mean radius (derived quantity)
    call volcaero_constituents(2)%instantiate( &
         std_name          = volcrad_const_name, &
         diag_name         = volcrad_const_name, &
         long_name         = 'prescribed volcanic aerosol geometric-mean radius derived from mass', &
         units             = 'm', &
         vertical_dim      = 'vertical_layer_dimension', &
         min_value         = 0.0_kind_phys, &
         advected          = .false., &
         water_species     = .false., &
         mixing_ratio_type = 'dry', &
         errcode           = errflg, &
         errmsg            = errmsg)
    if (errflg /= 0) return

    if (amIRoot) then
      write(iulog,*) trim(subname)//': Registered 2 prescribed volcanic aerosol constituents'
    end if

  end subroutine prescribed_volcanic_aerosol_register

  ! Initialize prescribed volcanic aerosol reading via tracer_data.
!> \section arg_table_prescribed_volcanic_aerosol_init Argument Table
!! \htmlinclude prescribed_volcanic_aerosol_init.html
  subroutine prescribed_volcanic_aerosol_init( &
    amIRoot, iulog, &
    prescribed_volcaero_name, &
    prescribed_volcaero_file, &
    prescribed_volcaero_filelist, &
    prescribed_volcaero_datapath, &
    prescribed_volcaero_type, &
    prescribed_volcaero_cycle_yr, &
    prescribed_volcaero_fixed_ymd, &
    prescribed_volcaero_fixed_tod, &
    errmsg, errflg)

    ! host model dependency for tracer_data read utility
    use tracer_data, only: trcdata_init

    ! host model dependency for history output
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    ! Input arguments:
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_volcaero_name          ! netCDF field name for volcanic aerosol
    character(len=*),   intent(in)  :: prescribed_volcaero_file          ! input filename from namelist
    character(len=*),   intent(in)  :: prescribed_volcaero_filelist      ! input filelist from namelist
    character(len=*),   intent(in)  :: prescribed_volcaero_datapath      ! input datapath from namelist
    character(len=*),   intent(in)  :: prescribed_volcaero_type          ! data type from namelist
    integer,            intent(in)  :: prescribed_volcaero_cycle_yr      ! cycle year from namelist [1]
    integer,            intent(in)  :: prescribed_volcaero_fixed_ymd     ! fixed year-month-day from namelist (YYYYMMDD) [1]
    integer,            intent(in)  :: prescribed_volcaero_fixed_tod     ! fixed time of day from namelist [s]

    ! Output arguments:
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables:
    character(len=64) :: specifier(1)

    character(len=*), parameter :: subname = 'prescribed_volcanic_aerosol_init'

    errmsg = ''
    errflg = 0

    if (.not. has_prescribed_volcaero) return

    if (amIRoot) then
      write(iulog,*) trim(subname)//': volcanic aerosol is prescribed in: '// &
           trim(prescribed_volcaero_file)
    end if

    ! Build specifier: constituent_name:ncdf_field_name
    specifier(1) = trim(volcaero_const_name) // ':' // trim(prescribed_volcaero_name)

    ! Initialize tracer_data module with file and field information
    call trcdata_init( &
      specifier      = specifier, &
      filename       = prescribed_volcaero_file, &
      filelist       = prescribed_volcaero_filelist, &
      datapath       = prescribed_volcaero_datapath, &
      flds           = tracer_data_fields, &
      file           = tracer_data_file, &
      data_cycle_yr  = prescribed_volcaero_cycle_yr, &
      data_fixed_ymd = prescribed_volcaero_fixed_ymd, &
      data_fixed_tod = prescribed_volcaero_fixed_tod, &
      data_type      = prescribed_volcaero_type)

    ! Verify tracer_data is correctly initialized
    if (.not. associated(tracer_data_fields)) then
      errflg = 1
      errmsg = subname//': tracer_data_fields not associated after trcdata_init'
      return
    end if

    ! Register history fields
    call history_add_field(volcaero_const_name, &
         'prescribed volcanic aerosol dry mass mixing ratio', &
         'lev', 'inst', 'kg kg-1')
    call history_add_field(volcrad_const_name, &
         'volcanic aerosol geometric-mean radius', &
         'lev', 'inst', 'm')
    call history_add_field('VOLC_MASS', &
         'volcanic aerosol vertical mass path in layer', &
         'lev', 'inst', 'kg m-2')
    call history_add_field('VOLC_MASS_C', &
         'volcanic aerosol column mass', &
         horiz_only, 'inst', 'kg m-2')

    if (amIRoot) then
      write(iulog,*) trim(subname)//': Initialized volcanic aerosol field from tracer_data'
    end if

  end subroutine prescribed_volcanic_aerosol_init

  ! Advance prescribed volcanic aerosol data, convert units, apply tropopause
  ! masking, and compute geometric-mean radius.
!> \section arg_table_prescribed_volcanic_aerosol_run Argument Table
!! \htmlinclude prescribed_volcanic_aerosol_run.html
  subroutine prescribed_volcanic_aerosol_run( &
    ncol, pver, pcnst, &
    const_props, &
    mwdry, boltz, gravit, &
    T, pmiddry, pdel, zm, &
    pmid, pint, phis, zi, &
    tropLev, &
    constituents, &
    errmsg, errflg)

    ! host model dependency for tracer_data
    use tracer_data,               only: advance_trcdata

    ! host model dependency for history output
    use cam_history,               only: history_out_field

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! dependency to get constituent index
    use ccpp_const_utils,          only: ccpp_const_get_idx

    ! dependency for unit string handling
    use string_utils,              only: to_lower, get_last_significant_char

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: pcnst
    type(ccpp_constituent_prop_ptr_t), &
                        intent(in)    :: const_props(:)
    real(kind_phys),    intent(in)    :: mwdry             ! molecular weight of dry air [g mol-1]
    real(kind_phys),    intent(in)    :: boltz             ! Boltzmann constant [J K-1 molecule-1]
    real(kind_phys),    intent(in)    :: gravit            ! gravitational acceleration [m s-2]
    real(kind_phys),    intent(in)    :: T(:,:)            ! air temperature [K] (layer centers)
    real(kind_phys),    intent(in)    :: pmiddry(:,:)      ! dry air pressure [Pa] (layer centers)
    real(kind_phys),    intent(in)    :: pdel(:,:)         ! air pressure thickness [Pa] (layer centers)
    real(kind_phys),    intent(in)    :: zm(:,:)           ! geopotential height wrt surface [m] (layer centers)
    real(kind_phys),    intent(in)    :: pmid(:,:)         ! air pressure [Pa] (layer centers)
    real(kind_phys),    intent(in)    :: pint(:,:)         ! air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)    :: phis(:)           ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zi(:,:)           ! geopotential height wrt surface at interfaces [m]
    integer,            intent(in)    :: tropLev(:)        ! tropopause vertical layer index [index]

    real(kind_phys),    intent(inout) :: constituents(:,:,:)

    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local variables
    integer         :: i, k
    integer         :: mmr_idx, rad_idx
    real(kind_phys) :: to_mmr(ncol, pver)    ! unit conversion factor to MMR [1]
    real(kind_phys) :: mmrvolc               ! volcanic aerosol MMR [kg kg-1]
    real(kind_phys) :: concvolc              ! mass concentration of volcanic aerosol [kg m-3]
    real(kind_phys) :: volcmass(ncol, pver)  ! volcanic aerosol mass path in layer [kg m-2]
    real(kind_phys) :: columnmass(ncol)      ! volcanic aerosol column mass [kg m-2]

    character(len=*), parameter :: subname = 'prescribed_volcanic_aerosol_run'

    errmsg = ''
    errflg = 0

    if (.not. has_prescribed_volcaero) return

    ! Advance tracer_data to current time
    call advance_trcdata(tracer_data_fields, tracer_data_file, &
                         pmid, pint, phis, zi)

    ! Get constituent indices for MMR and radius
    call ccpp_const_get_idx(const_props, volcaero_const_name, &
         mmr_idx, errmsg, errflg)
    if (errflg /= 0) return

    call ccpp_const_get_idx(const_props, volcrad_const_name, &
         rad_idx, errmsg, errflg)
    if (errflg /= 0) return

    ! Determine unit conversion factor based on units in the input file
    select case ( to_lower(trim(tracer_data_fields(1)%units(:get_last_significant_char(tracer_data_fields(1)%units)))) )
    case ("molec/cm3", "/cm3", "molecules/cm3", "cm^-3", "cm**-3")
      ! Number density [molecules cm-3] -> MMR [kg kg-1]
      ! mmr = (M * 1e6 * k_B * T) / (M_air * p_dry)
      to_mmr(:ncol,:) = (molmass_volcaero * 1.0e6_kind_phys * boltz * T(:ncol,:)) &
                       / (mwdry * pmiddry(:ncol,:))
    case ('kg/kg', 'mmr', 'kg kg-1')
      to_mmr(:ncol,:) = 1.0_kind_phys
    case ('mol/mol', 'mole/mole', 'vmr', 'fraction')
      to_mmr(:ncol,:) = molmass_volcaero / mwdry
    case default
      errflg = 1
      errmsg = subname//': unrecognized units: '//trim(tracer_data_fields(1)%units)
      return
    end select

    ! Convert tracer_data field to MMR and store in constituent array.
    ! Apply tropopause masking: zero below tropopause.
    ! Compute geometric-mean radius where MMR > 0 above tropopause.
    constituents(:ncol, :pver, rad_idx) = 0.0_kind_phys

    do k = 1, pver
      do i = 1, ncol
        ! Apply unit conversion
        mmrvolc = to_mmr(i, k) * tracer_data_fields(1)%data(i, k)

        ! Zero below tropopause
        if (k >= tropLev(i)) then
          mmrvolc = 0.0_kind_phys
        end if

        constituents(i, k, mmr_idx) = mmrvolc

        ! Compute geometric-mean wet aerosol radius from mass concentration
        if (mmrvolc > 0.0_kind_phys) then
          ! concvolc [kg m-3] = mmr [kg kg-1] * pdel [Pa] / (g [m s-2] * zm [m])
          concvolc = (mmrvolc * pdel(i, k)) / (gravit * zm(i, k))
          constituents(i, k, rad_idx) = radius_conversion * (concvolc ** (1.0_kind_phys / 3.0_kind_phys))
        end if
      end do
    end do

    ! Compute volcanic aerosol mass path in each layer [kg m-2]
    volcmass(:ncol, :pver)  = constituents(:ncol, :pver, mmr_idx) * pdel(:ncol, :pver) / gravit
    columnmass(:ncol) = sum(volcmass(:ncol, :pver), 2)

    ! History output
    call history_out_field(volcaero_const_name, constituents(:, :, mmr_idx))
    call history_out_field(volcrad_const_name,  constituents(:, :, rad_idx))
    call history_out_field('VOLC_MASS',   volcmass(:, :))
    call history_out_field('VOLC_MASS_C', columnmass(:))

  end subroutine prescribed_volcanic_aerosol_run

end module prescribed_volcanic_aerosol
