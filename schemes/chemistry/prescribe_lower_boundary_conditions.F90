! Read time-varying chemical lower boundary conditions and serve them to the
! schemes that apply them.
!
! Reads a CHEM_LBC_FILE dataset of surface volume mixing ratios (VMRs)
! on a zonal-mean (lat x time) or longitude-latitude (lon x lat x time) grid
! and advances its read window every timestep.
!
! This scheme is a provider of LBC data and does not write to
! the constituents array. Instead, it provides data to other CCPP schemes
! following this lifecycle:
!
! init: use lbc_activate_species to "activate" this species.
! timestep_init/run:
!   use lbc_column_vmr to get per-column VMRs (generally used for prognostic chem) -or-
!   use lbc_global_mean_vmr for global mean VMRs (generally used for GHGs)
!
! for consistency, all species must be "active" at first timestep (LBC read)
! this ensures that ALL flbc entries in the namelist are applied to the model
! by another scheme and not silently unused.
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

  ! interfaces for the consuming CCPP schemes that apply the data this module reads:
  public :: lbc_initialized
  public :: lbc_species_count
  public :: lbc_species_name
  public :: lbc_activate_species
  public :: lbc_column_vmr
  public :: lbc_global_mean_vmr

  ! per-species lower boundary condition state
  type :: flbc
    character(len=16)            :: species = ' '         ! species name as given in flbc_list
    character(len=32)            :: fldname = ' '         ! netCDF variable name in flbc_file
    logical                      :: active = .false.      ! some scheme looked this species up at init
    real(kind_phys), allocatable :: vmr(:, :)             ! surface vmr on model columns (columns, time window)
    logical                      :: has_mean = .false.    ! file provides <species>_LBC_mean global means
    real(kind_phys), allocatable :: vmr_mean(:)           ! global mean vmr from file (time window)
  end type flbc

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

  ! every listed species has been looked up by a consumer: checked once at
  ! the first timestep_init, after every consumer's init has had its chance
  logical :: lbc_active_check = .false.

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

  ! number of time samples to read ahead of the current time (SERIAL)
  integer, parameter :: time_span = 1

contains

!> \section arg_table_prescribe_lower_boundary_conditions_init  Argument Table
!! \htmlinclude prescribe_lower_boundary_conditions_init.html
  subroutine prescribe_lower_boundary_conditions_init( &
    amIRoot, iulog, pi, ncol, lat, lon, &
    curr_date, curr_sec, &
    flbc_file, flbc_list, flbc_type, &
    flbc_cycle_yr, flbc_fixed_ymd, flbc_fixed_tod, &
    errmsg, errflg)

    ! portable netCDF reader for the dataset
    use ccpp_io_reader, only: create_netcdf_reader_t

    ! Host model dependency for interpolation:
    use interpolate_data, only: lininterp_init

    ! portable dependencies
    use atmos_phys_string_utils, only: to_upper

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

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! local variables
    integer            :: m, n
    integer            :: yr, mon, day
    integer            :: wrk_date, wrk_sec
    real(kind_phys)    :: wrk_time
    real(kind_phys), allocatable :: file_lat(:)
    real(kind_phys)    :: d2r
    character(len=512) :: read_errmsg
    character(len=*), parameter :: subname = 'prescribe_lower_boundary_conditions_init'

    errmsg = ''
    errflg = 0

    ! CCPP init phases are idempotent:
    if(lbc_initialized) return

    is_root = amIRoot
    log_unit = iulog
    pi_const = pi
    lbc_initialized = .true.
    flbc_cnt = 0
    dels_valid = .false.
    lbc_active_check = .false.

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

      ! Disallow duplicate entries:
      if (any(flbcs(1:m - 1)%species == flbcs(m)%species)) then
        errmsg = subname // ': flbc_list member "' // trim(flbcs(m)%species) // &
                 '" is listed twice'
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
    ncol, lat, lon, curr_date, curr_sec, &
    errmsg, errflg)

    integer,            intent(in)    :: ncol   ! number of columns [count]
    real(kind_phys),    intent(in)    :: lat(:) ! latitude of columns [rad]
    real(kind_phys),    intent(in)    :: lon(:) ! longitude of columns [rad]
    integer,            intent(in)    :: curr_date ! current model date [YYYYMMDD]
    integer,            intent(in)    :: curr_sec  ! elapsed seconds of current date [s]

    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! local variables
    integer           :: m
    integer           :: last, next
    real(kind_phys)   :: dels
    real(kind_phys)   :: wrk_time

    errmsg = ''
    errflg = 0

    if (flbc_cnt < 1) then
      return
    end if

    ! Verify that all LBC entries have been activated by at least one CCPP scheme.
    if (.not. lbc_active_check) then
      do m = 1, flbc_cnt
        if (.not. flbcs(m)%active) then
          ! abort here: all flbc_list entries MUST be claimed
          ! this is to prevent silent failures in the flbc_list where an entry is added
          ! but it actually has no effect on the model:
          errmsg = 'prescribe_lower_boundary_conditions_timestep_init: no scheme in the suite' // &
                   ' applies flbc_list member "' // trim(flbcs(m)%species) // &
                   '" (no lbc_activate_species lookup during init)'
          errflg = 1
          return
        end if
      end do
    end if

    lbc_active_check = .true.

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

    ! Cache this step's time interpolation for the consumers' reads
    call get_dels(curr_date, curr_sec, dels, last, next, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if
    dels_cur = dels
    last_cur = last
    next_cur = next
    dels_valid = .true.

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
    lbc_active_check = .false.
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

  ! Read one species' lower boundary values for the current time window and
  ! interpolate horizontally to the physics columns.
  subroutine flbc_get(lbcs, ncol, to_lats, to_lons, errmsg, errflg)

    ! Host model dependency for interpolation:
    use interpolate_data, only: lininterp_init, lininterp

    type(flbc),        intent(inout) :: lbcs
    integer,           intent(in)    :: ncol
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

  ! Interface for the consuming schemes that apply the lower boundary data:

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

  ! Mark a species in the LBC list as "activated",
  ! i.e., another CCPP scheme promises (by contract) that it will retrieve the LBC value
  ! of that species and apply it to the constituents array.
  ! This is a one (LBC) to many (>= 1) relationship.
  ! If a species is not matched, the handle returns 0.
  subroutine lbc_activate_species(name, handle)
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: handle

    integer :: m

    handle = 0
    do m = 1, flbc_cnt
      if (trim(flbcs(m)%species) == trim(name)) then
        flbcs(m)%active = .true.
        handle = m
        return
      end if
    end do
  end subroutine lbc_activate_species

  ! Get per-column time-interpolated VMR of a LBC species at the current timestep.
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

  ! Get global mean time-interpolated VMR of a LBC species at the current timestep.
  subroutine lbc_global_mean_vmr(handle, ncol, vmr_out, errmsg, errflg)

    ! Non-portable due to host dependency: Global mean module gmean from src/utils
    use gmean_mod,    only: gmean

    integer,          intent(in)  :: handle
    integer,          intent(in)  :: ncol     ! number of columns [count]
    real(kind_phys),  intent(out) :: vmr_out  ! [mol mol-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! local variables
    real(kind_phys) :: vmr_arr(ncol)

    errmsg = ''
    errflg = 0
    if (handle < 1 .or. handle > flbc_cnt) then
      write(errmsg, '(a,i0,a,i0,a)') 'lbc_global_mean_vmr: species handle ', handle, &
           ' is outside the flbc_list of ', flbc_cnt, ' species'
      errflg = 1
      return
    end if
    if (.not. dels_valid) then
      errmsg = 'lbc_global_mean_vmr: prescribe_lower_boundary_conditions_timestep_init has not run this step'
      errflg = 1
      return
    end if

    if (flbcs(handle)%has_mean) then
      ! mean is from file
      vmr_out = flbcs(handle)%vmr_mean(last_cur) &
                + dels_cur*(flbcs(handle)%vmr_mean(next_cur) - flbcs(handle)%vmr_mean(last_cur))
    else
      ! mean computed from host
      vmr_arr(:ncol) = flbcs(handle)%vmr(:ncol, last_cur) &
                       + dels_cur*(flbcs(handle)%vmr(:ncol, next_cur) - flbcs(handle)%vmr(:ncol, last_cur))
      call gmean(vmr_arr, vmr_out)
    end if

  end subroutine lbc_global_mean_vmr

  ! Convert an integer date and seconds of day to a float time coordinate
  ! consistent with the model calendar.
  subroutine flt_date(ncdate, ncsec, time)

    ! Host model dependency for calendar math on dataset and model times.
    use time_manager, only: timemgr_datediff

    integer,         intent(in)  :: ncdate ! date YYYYMMDD
    integer,         intent(in)  :: ncsec  ! seconds of day
    real(kind_phys), intent(out) :: time

    ! Reference datetime matching CAM time_utils::flt_date for bfb:
    integer, parameter :: refymd = 00010101
    integer, parameter :: reftod = 0

    call timemgr_datediff(refymd, reftod, ncdate, ncsec, time)

  end subroutine flt_date

end module prescribe_lower_boundary_conditions
