! Prescribe time-varying chemical lower boundary (or total column) conditions.
!
! Reads a CHEM_LBC_FILE dataset of surface mole fractions
! on a zonal-mean (lat x time) or longitude-latitude (lon x lat x time) grid
! and advances its read window every timestep.
!
! See fill mechanics below in the flbc struct.
!
! This scheme is inactive unless flbc_file is set.
!
! Based on original CAM version from Francis Vitt et al.
module prescribe_lower_boundary_conditions
  use ccpp_kinds,     only: kind_phys
  use ccpp_io_reader, only: abstract_netcdf_reader_t

  ! Host model dependency for interpolation:
  use interpolate_data, only: interp_type

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: prescribe_lower_boundary_conditions_init
  public :: prescribe_lower_boundary_conditions_timestep_init
  public :: prescribe_lower_boundary_conditions_final

  ! interfaces for other CCPP schemes to retrieve LBC information read by this module
  ! at timestep_init:
  public :: lbc_initialized
  public :: lbc_species_count
  public :: lbc_species_name
  public :: lbc_claim_species
  public :: lbc_column_vmr

  ! per-species lower boundary condition state
  type :: flbc
    character(len=16)            :: species = ' '         ! species name as given in flbc_list
    character(len=32)            :: fldname = ' '         ! netCDF variable name in flbc_file
    integer                      :: const_idx = -1        ! constituent array index the species fills

    ! Each species in flbc_list is either "claimed" (and applied by another scheme)
    ! or filled by this scheme.
    !
    !  - Advected species for prognostic chemistry (claimed = .true.).
    !    In chemistry, the bottom-layer vmr-space is overwritten with the LBC data.
    !    This is handled inside the chemistry step and thus this scheme only reads
    !    the data and provides per-column values via lbc_column_vmr.
    !    We do not write to constituents from this scheme in this case.
    !
    !  - Non-advected GHGs (claimed = .false.)
    !    For GHGs only (since only GHG molar masses are known in radiation_utils)
    !    Time-interpolated *global mean* vmrs are written as a whole-column uniform mmr
    !    (after conversion).
    logical                      :: is_advected = .false. ! advected constituent: must be claimed, never filled here
    logical                      :: claimed = .false.     ! applied by a chemistry scheme, not filled here

    ! greenhouse gas with a molar mass known to radiation_utils?
    logical                      :: is_ghg = .false.
    ! conversion factor: vmr -> mass mixing ratio w.r.t. dry air
    real(kind_phys)              :: mmr_factor = -1._kind_phys

    real(kind_phys), allocatable :: vmr(:, :)             ! surface vmr on model columns (columns, time window)
    logical                      :: has_mean = .false.    ! file provides <species>_LBC_mean global means
    real(kind_phys), allocatable :: vmr_mean(:)           ! global mean vmr from file (time window)
  end type flbc

  ! number of time samples to read ahead of the current time (SERIAL)
  integer, parameter :: time_span = 1

  ! dataset time coordinate and current read window
  integer                       :: ntimes = 0
  integer                       :: tim_ndx(2) = 0
  integer,         allocatable  :: dates(:)
  real(kind_phys), allocatable  :: times(:)

  ! namelist-derived module state
  character(len=256) :: filename = ' '     ! resolved local path of flbc_file
  character(len=32)  :: lbc_type = 'SERIAL'
  integer            :: lbc_cycle_yr = 0
  integer            :: lbc_fixed_ymd = 0
  integer            :: lbc_fixed_tod = 0

  integer                  :: flbc_cnt = 0 ! number of active flbc species
  type(flbc), allocatable  :: flbcs(:)

  ! module state resolved at init
  logical, protected :: lbc_initialized = .false. ! the init phase has run (idle or not)
  logical :: is_root = .false.
  integer :: log_unit = -1
  real(kind_phys) :: pi_const = 0._kind_phys

  ! time interpolation of the current step, set by timestep_init and read by
  ! lbc_column_vmr in the run phase (CAM computes it from the date its
  ! timestep_init stored, so both see the same weights)
  real(kind_phys) :: dels_cur = 0._kind_phys
  integer         :: last_cur = -1
  integer         :: next_cur = -1
  logical         :: dels_valid = .false.

  ! dataset grid sizes and horizontal interpolation weights to the physics columns.
  ! The columns and the dataset grid do not change, so the weights are built once.
  integer           :: nlat = 0
  integer           :: nlon = 0
  type(interp_type) :: lat_wgts
  type(interp_type) :: lon_wgts
  logical           :: lat_wgts_built = .false.
  logical           :: lon_wgts_built = .false. ! only (lon, lat, time) fields need lon weights

  ! netCDF reader for the flbc_file, created at init
  class(abstract_netcdf_reader_t), pointer :: file_reader => null()

  ! ccpp_io_reader error codes tolerated by probing reads:
  integer, parameter :: missing_variable_error_code = 3
  integer, parameter :: wrong_rank_error_code = 5

contains

!> \section arg_table_prescribe_lower_boundary_conditions_init  Argument Table
!! \htmlinclude prescribe_lower_boundary_conditions_init.html
  subroutine prescribe_lower_boundary_conditions_init( &
    amIRoot, iulog, pi, ncol, lat, lon, &
    curr_date, curr_sec, &
    flbc_file, flbc_list, flbc_type, &
    flbc_cycle_yr, flbc_fixed_ymd, flbc_fixed_tod, &
    const_props, &
    errmsg, errflg)

    ! portable netCDF reader for the dataset
    use ccpp_io_reader, only: create_netcdf_reader_t

    ! Host model dependency for interpolation:
    use interpolate_data, only: lininterp_init

    ! framework dependencies to look up and describe constituents
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t, int_unassigned
    use ccpp_scheme_utils,         only: ccpp_constituent_index

    ! portable dependencies
    use atmos_phys_string_utils, only: to_upper
    use radiation_utils,         only: get_molar_mass_ratio

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    real(kind_phys),    intent(in)  :: pi     ! pi constant [1]
    integer,            intent(in)  :: ncol   ! number of columns [count]
    real(kind_phys),    intent(in)  :: lat(:) ! latitude of columns [rad]
    real(kind_phys),    intent(in)  :: lon(:) ! longitude of columns [rad]
    integer,            intent(in)  :: curr_date ! current model date [YYYYMMDD]
    integer,            intent(in)  :: curr_sec  ! elapsed seconds of current date [s]

    ! input fields from namelist
    character(len=*),   intent(in)  :: flbc_file      ! dataset filename ('UNSET' disables the scheme)
    character(len=*),   intent(in)  :: flbc_list(:)   ! species to prescribe
    character(len=*),   intent(in)  :: flbc_type      ! time interpolation: 'SERIAL', 'CYCLICAL', or 'FIXED'
    integer,            intent(in)  :: flbc_cycle_yr  ! cycle year (CYCLICAL only)
    integer,            intent(in)  :: flbc_fixed_ymd ! fixed date YYYYMMDD (FIXED only)
    integer,            intent(in)  :: flbc_fixed_tod ! fixed time of day [s] (FIXED only)

    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:) ! constituent properties

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! local variables
    integer            :: m, n
    integer            :: yr, mon, day
    integer            :: wrk_date, wrk_sec
    real(kind_phys)    :: wrk_time
    character(len=32)  :: target_constituent_name
    logical            :: is_mass_mmr, is_dry_mmr
    real(kind_phys), allocatable :: file_lat(:)
    real(kind_phys)    :: d2r
    character(len=512) :: read_errmsg
    character(len=*), parameter :: subname = 'prescribe_lower_boundary_conditions_init'

    ! greenhouse gas species accepted in flbc_list.
    ! CFC11eq is the CMIP CFC11-equivalent (all other halogens expressed as
    ! CFC11) and fills the CFC11 constituent.
    character(len=8), parameter :: ghg_names(*) = &
         [character(len=8) :: 'CO2', 'CH4', 'N2O', 'CFC11', 'CFC12', 'CFC11eq']

    errmsg = ''
    errflg = 0

    is_root = amIRoot
    log_unit = iulog
    pi_const = pi
    lbc_initialized = .true.
    flbc_cnt = 0
    dels_valid = .false.

    ! check if user has specified an input dataset
    if (flbc_file == 'UNSET' .or. len_trim(flbc_file) == 0) then
      return
    end if

    mon = mod(curr_date, 10000)/100
    day = mod(curr_date, 100)

    ! Check timing
    lbc_type = to_upper(flbc_type)
    lbc_cycle_yr = flbc_cycle_yr
    lbc_fixed_ymd = flbc_fixed_ymd
    lbc_fixed_tod = flbc_fixed_tod

    if (lbc_type /= 'SERIAL' .and. lbc_type /= 'CYCLICAL' .and. lbc_type /= 'FIXED') then
      errmsg = subname // ': flbc_type "' // trim(flbc_type) // &
               '" is not SERIAL, CYCLICAL, or FIXED'
      errflg = 1
      return
    end if
    if (lbc_cycle_yr > 0 .and. lbc_type /= 'CYCLICAL') then
      errmsg = subname // ': cannot specify flbc_cycle_yr if flbc_type is not CYCLICAL'
      errflg = 1
      return
    end if
    if ((lbc_fixed_ymd > 0 .or. lbc_fixed_tod > 0) .and. lbc_type /= 'FIXED') then
      errmsg = subname // ': cannot specify flbc_fixed_ymd or ' // &
               'flbc_fixed_tod if flbc_type is not FIXED'
      errflg = 1
      return
    end if

    wrk_sec = curr_sec
    if (lbc_type == 'SERIAL') then
      wrk_date = curr_date
    else if (lbc_type == 'CYCLICAL') then
      wrk_date = lbc_cycle_yr*10000 + mod(curr_date, 10000)
      ! If this is a leap-day, we have to avoid asking for a non-leap-year
      ! on a cyclical dataset. When this happens, just use Feb 28 instead
      if ((mon == 2) .and. (day == 29)) then
        wrk_date = wrk_date - 1
        if (is_root) then
          write(iulog, *) subname // ' WARNING: using ' // &
                          'Feb 28 instead of Feb 29 for cyclical dataset'
        end if
      end if
    else
      wrk_date = lbc_fixed_ymd
      wrk_sec = lbc_fixed_tod
    end if
    call flt_date(wrk_date, wrk_sec, wrk_time)

    ! Species with fixed lbc: map each to constituent
    flbc_cnt = 0
    count_loop: do m = 1, size(flbc_list)
      if (len_trim(flbc_list(m)) == 0 .or. trim(flbc_list(m)) == 'UNSET') then
        exit count_loop
      end if
      flbc_cnt = flbc_cnt + 1
    end do count_loop

    if (flbc_cnt == 0) then
      return
    end if

    allocate (flbcs(flbc_cnt), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ': failed to allocate flbcs: ' // trim(errmsg)
      return
    end if

    species_loop: do m = 1, flbc_cnt

      flbcs(m)%species = trim(flbc_list(m))
      flbcs(m)%is_ghg = any(ghg_names == flbc_list(m))

      ! netCDF variable naming: CFC11 and CFC12 are stored under their chemical formulas.
      !REMOVECAM: in the distant future where we have the opportunity to restructure input files
      ! we should remove this semi-permanent quirk.
      if (trim(flbcs(m)%species) == 'CFC11') then
        flbcs(m)%fldname = 'CFCL3_LBC'
      else if (trim(flbcs(m)%species) == 'CFC12') then
        flbcs(m)%fldname = 'CF2CL2_LBC'
      else
        flbcs(m)%fldname = trim(flbcs(m)%species)//'_LBC'
      end if
      !REMOVECAM_END

      ! CFC11eq fills the CFC11 constituent.
      ! Every other species fills the constituent matching its own name.
      ! The molar mass of the filled constituent is used for the vmr -> mmr conversion.
      target_constituent_name = trim(flbcs(m)%species)
      if (trim(flbcs(m)%species) == 'CFC11eq') then
        target_constituent_name = 'CFC11'
      end if

      call ccpp_constituent_index(trim(target_constituent_name), flbcs(m)%const_idx, errflg, errmsg)
      if (errflg /= 0) then
        return
      end if
      if (flbcs(m)%const_idx == int_unassigned) then
        errmsg = subname // ': no constituent found for flbc_list member "' // &
                 trim(flbcs(m)%species) // '"'
        errflg = 1
        return
      end if

      call const_props(flbcs(m)%const_idx)%is_advected(flbcs(m)%is_advected, errflg, errmsg)
      if (errflg /= 0) then
        return
      end if
      call const_props(flbcs(m)%const_idx)%is_mass_mixing_ratio(is_mass_mmr, errflg, errmsg)
      if (errflg /= 0) then
        return
      end if
      call const_props(flbcs(m)%const_idx)%is_dry(is_dry_mmr, errflg, errmsg)
      if (errflg /= 0) then
        return
      end if

      ! Only the greenhouse gases have a molar mass here, so only they can be
      ! filled by this scheme; anything else must be claimed by a chemistry
      ! scheme before the first timestep, which is checked in timestep_init.
      if (flbcs(m)%is_ghg) then
        ! conversion factor from vmr to mass mixing ratio w.r.t. dry air:
        ! get_molar_mass_ratio returns (molar mass of dry air)/(molar mass of species)
        call get_molar_mass_ratio(trim(target_constituent_name), flbcs(m)%mmr_factor, errmsg, errflg)
        if (errflg /= 0) then
          return
        end if
        flbcs(m)%mmr_factor = 1._kind_phys/flbcs(m)%mmr_factor
      end if

      ! The global-mean fill writes the constituent as a dry mass mixing ratio (MMR).
      ! LBCs are read as VMR, so we need to be able to convert them (i.e., mmr factor is available).
      ! Consistency check:
      if(.not. flbcs(m)%is_advected .and. &           ! global-mean fill (needs VMR -> MMR conversion)
               flbcs(m)%is_ghg      .and. &           ! is a GHG (will be converted to MMR)
          (.not. is_mass_mmr .or. .not. is_dry_mmr))  then  ! but constituent not initd as dry and MMR
        errmsg = subname // ': non-advected constituent for flbc_list member "' // &
                 trim(flbcs(m)%species) // '" is not a dry mass mixing ratio;' // &
                 ' the global-mean fill is only implemented for that convention'
        errflg = 1
        return
      end if

      ! Consistency check: ensure we are not filling the same constituent from two LBC entries
      if (any(flbcs(1:m - 1)%const_idx == flbcs(m)%const_idx)) then
        errmsg = subname // ': flbc_list members "' // trim(flbcs(m)%species) // &
                 '" and another entry (CFC11 and CFC11eq?) fill the same constituent'
        errflg = 1
        return
      end if

    end do species_loop

    if (is_root) then
      write(iulog, *) subname // ': species with prescribed lower boundary values:'
      do m = 1, flbc_cnt
        write(iulog, *) '  ', trim(flbcs(m)%species), ' (file variable ', trim(flbcs(m)%fldname), ')'
      end do
      write(iulog, *) subname // ': lower boundary timing specifications:'
      write(iulog, *) '  type = ', trim(lbc_type)
      if (lbc_type == 'CYCLICAL') then
        write(iulog, *) '  cycle year = ', lbc_cycle_yr
      else if (lbc_type == 'FIXED') then
        write(iulog, *) '  fixed date = ', lbc_fixed_ymd
        write(iulog, *) '  fixed time = ', lbc_fixed_tod
      end if
    end if

    ! Get timing information, allocate arrays, and read in dates
    filename = trim(flbc_file)

    file_reader => create_netcdf_reader_t()
    call file_reader%open_file(trim(filename), errmsg, errflg)
    if (errflg /= 0) then
      return
    end if

    call file_reader%get_var('date', dates, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if
    ntimes = size(dates)

    allocate (times(ntimes), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ': failed to allocate times array: ' // trim(errmsg)
      return
    end if

    do n = 1, ntimes
      call flt_date(dates(n), 0, times(n))
    end do

    if (lbc_type /= 'CYCLICAL') then
      if (wrk_time < times(1) .or. wrk_time > times(ntimes)) then
        errmsg = subname // ': time out of bounds for dataset "' // trim(filename) // '"'
        errflg = 1
        return
      end if
      do n = 2, ntimes
        if (wrk_time <= times(n)) then
          exit
        end if
      end do
      tim_ndx(1) = n - 1
    else
      yr = lbc_cycle_yr
      do n = 1, ntimes
        if (yr == dates(n)/10000) then
          exit
        end if
      end do
      if (n >= ntimes) then
        errmsg = subname // ': cycle year out of bounds for dataset "' // trim(filename) // '"'
        errflg = 1
        return
      end if
      tim_ndx(1) = n
    end if

    select case (lbc_type)
    case ('FIXED')
      tim_ndx(2) = tim_ndx(1) + 1
    case ('CYCLICAL')
      do n = tim_ndx(1), ntimes
        if (yr /= dates(n)/10000) then
          exit
        end if
      end do
      tim_ndx(2) = n - 1
      if ((tim_ndx(2) - tim_ndx(1)) < 2) then
        errmsg = subname // ': cyclical lower-boundary conditions require at least two time points'
        errflg = 1
        return
      end if
    case ('SERIAL')
      tim_ndx(2) = min(ntimes, tim_ndx(1) + time_span)
    end select

    ! Get the latitude coordinate from the file and build the latitude weights
    call file_reader%get_var('lat', file_lat, read_errmsg, errflg)
    if (errflg /= 0) then
      errmsg = subname // ': cannot read lat coordinate of "' // &
               trim(filename) // '": ' // trim(read_errmsg)
      return
    end if
    nlat = size(file_lat)
    d2r = pi_const/180._kind_phys
    file_lat(:nlat) = file_lat(:nlat)*d2r

    call lininterp_init(file_lat, nlat, lat, ncol, 1, lat_wgts)
    lat_wgts_built = .true.

    ! Read in the flbc vmr for the current time window
    do m = 1, flbc_cnt
      call flbc_get(flbcs(m), ncol, lat, lon, errmsg, errflg)
      if (errflg /= 0) then
        return
      end if
    end do

    call file_reader%close_file(errmsg, errflg)
    if (errflg /= 0) then
      return
    end if

  end subroutine prescribe_lower_boundary_conditions_init

!> \section arg_table_prescribe_lower_boundary_conditions_timestep_init  Argument Table
!! \htmlinclude prescribe_lower_boundary_conditions_timestep_init.html
  subroutine prescribe_lower_boundary_conditions_timestep_init( &
    ncol, lat, lon, curr_date, curr_sec, constituents, &
    errmsg, errflg)

    integer,            intent(in)    :: ncol   ! number of columns [count]
    real(kind_phys),    intent(in)    :: lat(:) ! latitude of columns [rad]
    real(kind_phys),    intent(in)    :: lon(:) ! longitude of columns [rad]
    integer,            intent(in)    :: curr_date ! current model date [YYYYMMDD]
    integer,            intent(in)    :: curr_sec  ! elapsed seconds of current date [s]
    real(kind_phys),    intent(inout) :: constituents(:,:,:) ! constituent array (ncol, pver, pcnst)

    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! local variables
    integer           :: m
    integer           :: last, next
    real(kind_phys)   :: dels
    real(kind_phys)   :: wrk_time
    real(kind_phys)   :: vmr_gmean

    errmsg = ''
    errflg = 0

    if (flbc_cnt < 1) then
      return
    end if

    ! Advance the SERIAL read window when current time passes it
    if (lbc_type == 'SERIAL') then
      call flt_date(curr_date, curr_sec, wrk_time)
      if (wrk_time > times(tim_ndx(2))) then
        tim_ndx(1) = tim_ndx(2)
        tim_ndx(2) = min(ntimes, tim_ndx(1) + time_span)

        call file_reader%open_file(trim(filename), errmsg, errflg)
        if (errflg /= 0) then
          return
        end if
        do m = 1, flbc_cnt
          call flbc_get(flbcs(m), ncol, lat, lon, errmsg, errflg)
          if (errflg /= 0) then
            return
          end if
        end do
        call file_reader%close_file(errmsg, errflg)
        if (errflg /= 0) then
          return
        end if
      end if
    end if

    call get_dels(curr_date, curr_sec, dels, last, next, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if
    dels_cur = dels
    last_cur = last
    next_cur = next
    dels_valid = .true.

    do m = 1, flbc_cnt
      ! If this species is handled externally then we stop here.
      if (flbcs(m)%claimed) then
        cycle
      end if

      ! Fill each unclaimed species' non-advected constituent with the global
      ! mean vmr converted to a whole-column uniform dry mass mixing ratio.
      ! Species claimed by a chemistry scheme are applied there from the same weights.
      if (.not. flbcs(m)%is_ghg) then
        errmsg = 'prescribe_lower_boundary_conditions_timestep_init: flbc_list member "' // &
                 trim(flbcs(m)%species) // '" is not a greenhouse gas and no chemistry' // &
                 ' scheme claimed it (lbc_claim_species) during init'
        errflg = 1
        return
      end if
      ! CAM pins an advected species inside the chemistry step (flbc_set), so
      ! the scheme applying it must claim it; the core never writes one.
      if (flbcs(m)%is_advected) then
        errmsg = 'prescribe_lower_boundary_conditions_timestep_init: flbc_list member "' // &
                 trim(flbcs(m)%species) // '" fills an advected constituent, which must be' // &
                 ' applied by a chemistry scheme (lbc_claim_species) instead of this fill'
        errflg = 1
        return
      end if
      call global_mean_vmr(flbcs(m), ncol, dels, last, next, vmr_gmean)
      constituents(:, :, flbcs(m)%const_idx) = vmr_gmean*flbcs(m)%mmr_factor
    end do

  end subroutine prescribe_lower_boundary_conditions_timestep_init

!> \section arg_table_prescribe_lower_boundary_conditions_final  Argument Table
!! \htmlinclude prescribe_lower_boundary_conditions_final.html
  subroutine prescribe_lower_boundary_conditions_final(errmsg, errflg)

    ! Host model dependency for interpolation:
    use interpolate_data, only: lininterp_finish

    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    if (lat_wgts_built) then
      call lininterp_finish(lat_wgts)
      lat_wgts_built = .false.
    end if
    if (lon_wgts_built) then
      call lininterp_finish(lon_wgts)
      lon_wgts_built = .false.
    end if

    if (allocated(flbcs)) then
      deallocate (flbcs)
    end if
    flbc_cnt = 0
    dels_valid = .false.
    lbc_initialized = .false.
    if (allocated(dates)) then
      deallocate (dates)
    end if
    if (allocated(times)) then
      deallocate (times)
    end if

    ! Deallocate the module-level file reader object
    if (associated(file_reader)) then
      deallocate (file_reader)
      nullify (file_reader)
    end if

  end subroutine prescribe_lower_boundary_conditions_final

  ! Private helpers:

  ! Read one species' lower bndy values for the current time window and
  ! interpolate horizontally to the physics columns.
  subroutine flbc_get(lbcs, ncol, to_lats, to_lons, errmsg, errflg)

    ! Host model dependency for interpolation:
    use interpolate_data, only: lininterp_init, lininterp

    type(flbc),        intent(inout) :: lbcs
    integer,           intent(in)    :: ncol       ! number of columns [count]
    real(kind_phys),   intent(in)    :: to_lats(:) ! latitude of columns [rad]
    real(kind_phys),   intent(in)    :: to_lons(:) ! longitude of columns [rad]
    character(len=*),  intent(out)   :: errmsg
    integer,           intent(out)   :: errflg

    ! local variables
    integer :: m
    integer :: t1, t2, tcnt
    real(kind_phys), allocatable :: lon(:)
    real(kind_phys), allocatable :: wrk(:, :, :), wrk_zonal(:, :)
    logical            :: zonal_field
    character(len=512) :: read_errmsg
    real(kind_phys)    :: d2r, twopi
    real(kind_phys), parameter :: zero = 0._kind_phys

    errmsg = ''
    errflg = 0

    d2r = pi_const/180._kind_phys
    twopi = 2._kind_phys*pi_const

    t1 = tim_ndx(1)
    t2 = tim_ndx(2)
    tcnt = t2 - t1 + 1

    ! Reallocate only when the read window size changes (e.g. at the end of a SERIAL dataset)
    if (allocated(lbcs%vmr)) then
      if (size(lbcs%vmr, 2) /= tcnt) then
        deallocate (lbcs%vmr)
      end if
    end if
    if (.not. allocated(lbcs%vmr)) then
      allocate (lbcs%vmr(ncol, tcnt), stat=errflg, errmsg=errmsg)
      if (errflg /= 0) then
        errmsg = 'prescribe_lower_boundary_conditions (flbc_get): failed to allocate vmr: ' // trim(errmsg)
        return
      end if
    end if
    lbcs%vmr(:, :) = 0._kind_phys

    !-----------------------------------------------------------------------
    ! ... read the current time window of the field.  Try the zonal mean
    !     (lat, time) shape first, and fall back to (lon, lat, time) when
    !     the variable turns out to have a different rank.
    !-----------------------------------------------------------------------
    call file_reader%get_var(trim(lbcs%fldname), wrk_zonal, read_errmsg, errflg, &
         start=[1, t1], count=[nlat, tcnt])
    zonal_field = (errflg == 0)

    if (errflg == wrong_rank_error_code) then
      errflg = 0
      ! Build the longitude weights on the first (lon, lat, time) field;
      ! this happens during init, where every species is first read.
      if (.not. lon_wgts_built) then
        call file_reader%get_var('lon', lon, read_errmsg, errflg)
        if (errflg /= 0) then
          errmsg = 'prescribe_lower_boundary_conditions (flbc_get): cannot read lon coordinate of "' // &
                   trim(filename) // '": ' // trim(read_errmsg)
          return
        end if
        nlon = size(lon)
        lon(:nlon) = lon(:nlon)*d2r

        call lininterp_init(lon, nlon, to_lons, ncol, 2, lon_wgts, zero, twopi)
        lon_wgts_built = .true.
      end if

      call file_reader%get_var(trim(lbcs%fldname), wrk, read_errmsg, errflg, &
           start=[1, 1, t1], count=[nlon, nlat, tcnt])
    end if

    if (errflg /= 0) then
      errmsg = 'prescribe_lower_boundary_conditions (flbc_get): failed to read "' // &
               trim(lbcs%fldname) // '" from file "' // trim(filename) // '": ' // trim(read_errmsg)
      return
    end if

    !-----------------------------------------------------------------------
    ! ... interpolate to the physics columns
    !-----------------------------------------------------------------------
    if (zonal_field) then
      do m = 1, tcnt
        call lininterp(wrk_zonal(:, m), nlat, lbcs%vmr(:, m), ncol, lat_wgts)
      end do
    else
      do m = 1, tcnt
        call lininterp(wrk(:, :, m), nlon, nlat, lbcs%vmr(:, m), ncol, lon_wgts, lat_wgts)
      end do
    end if

    !-----------------------------------------------------------------------
    ! ... read the global mean directly if the file provides it
    !-----------------------------------------------------------------------
    call file_reader%get_var(trim(lbcs%fldname)//'_mean', lbcs%vmr_mean, read_errmsg, errflg, &
         start=[t1], count=[tcnt])
    if (errflg == 0) then
      lbcs%has_mean = .true.
    else if (errflg == missing_variable_error_code) then
      lbcs%has_mean = .false.
      errflg = 0
    else
      errmsg = 'prescribe_lower_boundary_conditions (flbc_get): failed to read "' // &
               trim(lbcs%fldname) // '_mean" from file "' // trim(filename) // '": ' // trim(read_errmsg)
      return
    end if

    if (is_root) then
      write(log_unit, *) 'prescribe_lower_boundary_conditions: read ', trim(lbcs%fldname), &
           ' dates ', dates(t1), ' - ', dates(t2)
    end if

  end subroutine flbc_get

  ! Compute the time interpolation factor and bracketing window-relative
  ! time indices for the current date:
  subroutine get_dels(ncdate, ncsec, dels, last, next, errmsg, errflg)

    ! Host model dependency
    use tracer_data, only: findplb

    integer,          intent(in)  :: ncdate ! current date YYYYMMDD
    integer,          intent(in)  :: ncsec  ! current time of day [s]
    real(kind_phys),  intent(out) :: dels
    integer,          intent(out) :: last   ! index into the read window (1-based)
    integer,          intent(out) :: next
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! local variables
    integer         :: wrk_date, wrk_sec
    integer         :: tcnt, n
    real(kind_phys) :: wrk_time

    errmsg = ''
    errflg = 0
    last = -1
    next = -1

    !--------------------------------------------------------
    ! ... setup the time interpolation
    !--------------------------------------------------------
    wrk_sec = ncsec
    select case (lbc_type)
    case ('SERIAL')
      wrk_date = ncdate
    case ('CYCLICAL')
      wrk_date = lbc_cycle_yr*10000 + mod(ncdate, 10000)
    case ('FIXED')
      wrk_date = lbc_fixed_ymd
      wrk_sec = lbc_fixed_tod
    end select

    call flt_date(wrk_date, wrk_sec, wrk_time)

    !--------------------------------------------------------
    ! ... set time interpolation factor
    !--------------------------------------------------------
    if (lbc_type /= 'CYCLICAL') then
      do n = tim_ndx(1) + 1, tim_ndx(2)
        if (wrk_time <= times(n)) then
          last = n - 1
          next = n
          exit
        end if
      end do
      if (last < 0) then
        errmsg = 'prescribe_lower_boundary_conditions (get_dels): interp time is out of bounds for dataset ' // &
                 trim(filename)
        errflg = 1
        return
      end if
      dels = (wrk_time - times(last))/(times(next) - times(last))
    else
      tcnt = tim_ndx(2) - tim_ndx(1) + 1
      call findplb(times(tim_ndx(1):tim_ndx(2)), tcnt, wrk_time, n)
      if (n < tcnt) then
        last = tim_ndx(1) + n - 1
        next = last + 1
        dels = (wrk_time - times(last))/(times(next) - times(last))
      else
        next = tim_ndx(1)
        last = tim_ndx(2)
        dels = wrk_time - times(last)
        if (dels < 0._kind_phys) then
          dels = 365._kind_phys + dels
        end if
        dels = dels/(365._kind_phys + times(next) - times(last))
      end if
    end if

    dels = max(min(1._kind_phys, dels), 0._kind_phys)

    ! convert absolute time indices to indices into the read window
    last = last - tim_ndx(1) + 1
    next = next - tim_ndx(1) + 1

  end subroutine get_dels

  ! Time-interpolated global mean vmr of one species:
  ! from the file's global mean variable when present, otherwise a global
  ! mean over the horizontally interpolated columns.
  subroutine global_mean_vmr(lbcs, ncol, dels, last, next, vmr_out)

    ! This scheme is non-portable due to dependency: Global mean module gmean from src/utils
    use gmean_mod,    only: gmean

    type(flbc),      intent(in)  :: lbcs
    integer,         intent(in)  :: ncol   ! number of columns [count]
    real(kind_phys), intent(in)  :: dels
    integer,         intent(in)  :: last
    integer,         intent(in)  :: next
    real(kind_phys), intent(out) :: vmr_out

    ! local variables
    real(kind_phys) :: vmr_arr(ncol)

    if (lbcs%has_mean) then
      vmr_out = lbcs%vmr_mean(last) &
                + dels*(lbcs%vmr_mean(next) - lbcs%vmr_mean(last))
    else
      vmr_arr(:ncol) = lbcs%vmr(:ncol, last) &
                       + dels*(lbcs%vmr(:ncol, next) - lbcs%vmr(:ncol, last))
      call gmean(vmr_arr, vmr_out)
    end if

  end subroutine global_mean_vmr

  ! Interface for a chemistry scheme that applies the lower boundary itself:

  ! Number of listed species (zero while the scheme is idle).
  pure integer function lbc_species_count()
    lbc_species_count = flbc_cnt
  end function lbc_species_count

  ! Name of listed species m as given in flbc_list.
  pure function lbc_species_name(m) result(name)
    integer, intent(in) :: m
    character(len=16)   :: name
    name = flbcs(m)%species
  end function lbc_species_name

  ! Hand a listed species to the caller: timestep_init no longer touches its
  ! constituent. Called from an init phase that runs after this scheme's.
  ! handle is the species' index for lbc_column_vmr, or 0 when the name is
  ! not in flbc_list (the caller decides whether that is an error).
  ! Claiming an already-claimed species is an error: every species has
  ! exactly one applier.
  subroutine lbc_claim_species(name, handle, errmsg, errflg)
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: handle
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: m

    errmsg = ''
    errflg = 0
    handle = 0
    do m = 1, flbc_cnt
      if (trim(flbcs(m)%species) == trim(name)) then
        if (flbcs(m)%claimed) then
          errmsg = 'lbc_claim_species: species "' // trim(name) // &
                   '" is already claimed by another scheme'
          errflg = 1
          return
        end if
        flbcs(m)%claimed = .true.
        handle = m
        return
      end if
    end do
  end subroutine lbc_claim_species

  ! Per-column vmr of a claimed species at the current step, interpolated in
  ! time with the weights of this step's timestep_init: CAM's flbc_set
  ! expression, so a caller writing it into its bottom layer matches CAM.
  ! handle comes from lbc_claim_species.
  subroutine lbc_column_vmr(handle, ncol, vmr_out, errmsg, errflg)
    integer,          intent(in)  :: handle
    integer,          intent(in)  :: ncol       ! number of columns [count]
    real(kind_phys),  intent(out) :: vmr_out(:) ! (ncol) [mol mol-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: m

    errmsg = ''
    errflg = 0
    m = handle
    if (m < 1 .or. m > flbc_cnt) then
      write(errmsg, '(a,i0,a,i0,a)') 'lbc_column_vmr: species handle ', m, &
           ' is outside the flbc_list of ', flbc_cnt, ' species'
      errflg = 1
      return
    end if
    if (.not. dels_valid) then
      errmsg = 'lbc_column_vmr: prescribe_lower_boundary_conditions_timestep_init has not run this step'
      errflg = 1
      return
    end if
    vmr_out(:ncol) = flbcs(m)%vmr(:ncol, last_cur) &
                     + dels_cur*(flbcs(m)%vmr(:ncol, next_cur) - flbcs(m)%vmr(:ncol, last_cur))
  end subroutine lbc_column_vmr

  ! Convert an integer date and seconds of day to a float time coordinate
  ! consistent with the model calendar.
  subroutine flt_date(ncdate, ncsec, time)

    ! Host model dependency for calendar math, applied to dataset dates as well
    ! as the current date, so it cannot be threaded as data by standard name.
    ! A portable replacement needs the calendar type and the reference date
    ! threaded as arguments to stay bit-for-bit; until such a utility exists
    ! this stays a host dependency like tracer_data and interpolate_data.
    use time_manager, only: timemgr_datediff

    integer,         intent(in)  :: ncdate ! date YYYYMMDD
    integer,         intent(in)  :: ncsec  ! seconds of day
    real(kind_phys), intent(out) :: time

    ! Days since 0001-01-01, CAM's time_utils flt_date, which mo_flbc and
    ! chem_surfvals use. The fixed reference keeps the fractional-day rounding,
    ! and so the interpolation weights, bit-for-bit with CAM and independent
    ! of the reference date of the run's clock.
    integer, parameter :: refymd = 00010101
    integer, parameter :: reftod = 0

    call timemgr_datediff(refymd, reftod, ncdate, ncsec, time)

  end subroutine flt_date

end module prescribe_lower_boundary_conditions
