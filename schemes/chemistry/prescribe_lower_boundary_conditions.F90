! Prescribe time-varying chemical lower boundary (or total column) conditions.
!
! Reads a CHEM_LBC_FILE dataset of surface mole fractions
! on a zonal-mean (lat x time) or longitude-latitude (lon x lat x time) grid
! and, every timestep, updates each species in flbc_list based on its
! constituent's advected flag:
!  - non-advected: the time-interpolated global-mean volume mixing ratio (vmr)
!    is written as a whole-column uniform mass mixing ratio.
!
!  - advected (prognostic chemistry): the time-interpolated per-column value
!    is pinned into the bottom two vertical levels as a dry mass mixing ratio.
!
! The scheme is inactive unless flbc_file is set.
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

  ! per-species lower boundary condition state
  type :: flbc
     character(len=16)            :: species = ' '   ! species name as given in flbc_list
     character(len=32)            :: fldname = ' '   ! netCDF variable name in flbc_file
     integer                      :: const_idx = -1  ! constituent array index the species fills
     logical                      :: is_advected = .false. ! advected constituent: pin bottom levels instead of global-mean fill
     real(kind_phys)              :: mmr_factor = 1._kind_phys ! conversion factor: vmr -> mass mixing ratio w.r.t. dry air
     real(kind_phys), allocatable :: vmr(:, :)       ! surface vmr on model columns (columns, time window)
     logical                      :: has_mean = .false. ! file provides <species>_LBC_mean global means
     real(kind_phys), allocatable :: vmr_mean(:)     ! global mean vmr from file (time window)
  end type flbc

  ! number of time samples to read ahead of the current time (SERIAL)
  integer, parameter :: time_span = 1

  ! greenhouse gas species accepted in flbc_list.
  ! CFC11eq is the CMIP CFC11-equivalent (all other halogens expressed as
  ! CFC11) and fills the CFC11 constituent.
  integer, parameter :: nghg = 6
  character(len=8), parameter :: ghg_names(nghg) = &
       [character(len=8) :: 'CO2', 'CH4', 'N2O', 'CFC11', 'CFC12', 'CFC11eq']

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
  logical :: is_root = .false.
  integer :: log_unit = -1
  real(kind_phys) :: pi_const = 0._kind_phys

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
    flbc_file, flbc_list, flbc_type, &
    flbc_cycle_yr, flbc_fixed_ymd, flbc_fixed_tod, &
    const_props, &
    errmsg, errflg)

    ! portable netCDF reader for the dataset
    use ccpp_io_reader, only: create_netcdf_reader_t

    ! Host model dependency for interpolation:
    use interpolate_data, only: lininterp_init

    ! CAM-SIMA host model dependency for the model date
    use time_manager,   only: get_curr_date, set_time_float_from_date

    ! framework dependency to describe constituents
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! portable dependencies
    use atmos_phys_string_utils, only: to_upper
    use radiation_utils,         only: get_molar_mass_ratio

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    real(kind_phys),    intent(in)  :: pi     ! pi constant [1]
    integer,            intent(in)  :: ncol   ! number of columns [count]
    real(kind_phys),    intent(in)  :: lat(:) ! latitude of columns [rad]
    real(kind_phys),    intent(in)  :: lon(:) ! longitude of columns [rad]

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
    integer            :: yr, mon, day, ncsec, ncdate
    integer            :: wrk_date, wrk_sec
    real(kind_phys)    :: wrk_time
    character(len=32)  :: mw_species
    logical            :: is_mass_mmr, is_dry_mmr
    character(len=256) :: diag_name
    real(kind_phys), allocatable :: file_lat(:)
    real(kind_phys)    :: d2r
    character(len=512) :: read_errmsg
    character(len=*), parameter :: subname = 'prescribe_lower_boundary_conditions_init'

    errmsg = ''
    errflg = 0

    is_root = amIRoot
    log_unit = iulog
    pi_const = pi

    ! check if user has specified an input dataset
    if (flbc_file == 'UNSET' .or. len_trim(flbc_file) == 0) then
      return
    end if

    call get_curr_date(yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day

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

    wrk_sec = ncsec
    if (lbc_type == 'SERIAL') then
      wrk_date = ncdate
    else if (lbc_type == 'CYCLICAL') then
      ! If this is a leap-day, we have to avoid asking for a non-leap-year
      ! on a cyclical dataset. When this happens, just use Feb 28 instead
      if ((mon == 2) .and. (day == 29)) then
        ncdate = yr*10000 + mon*100 + (day - 1)
        if (is_root) then
          write(iulog, *) subname // ' WARNING: using ' // &
                          'Feb 28 instead of Feb 29 for cyclical dataset'
        end if
      end if
      wrk_date = lbc_cycle_yr*10000 + mod(ncdate, 10000)
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

      if (.not. any(ghg_names == flbc_list(m))) then
        errmsg = subname // ': flbc_list member "' // trim(flbc_list(m)) // &
                 '" is not allowed (only greenhouse gas species are supported: molar masses' // &
                 ' are resolved through radiation_utils)'
        errflg = 1
        return
      end if

      flbcs(m)%species = trim(flbc_list(m))

      ! netCDF variable naming: CFC11 and CFC12 are
      ! stored under their chemical formulas.
      if (trim(flbcs(m)%species) == 'CFC11') then
        flbcs(m)%fldname = 'CFCL3_LBC'
      else if (trim(flbcs(m)%species) == 'CFC12') then
        flbcs(m)%fldname = 'CF2CL2_LBC'
      else
        flbcs(m)%fldname = trim(flbcs(m)%species)//'_LBC'
      end if

      ! CFC11eq fills the CFC11 constituent; every other species fills the
      ! constituent matching its own name.  The molar mass of the filled
      ! constituent is used for the vmr -> mmr conversion.
      if (trim(flbcs(m)%species) == 'CFC11eq') then
        mw_species = 'CFC11'
      else
        mw_species = trim(flbcs(m)%species)
      end if

      ! find the constituent whose diagnostic name matches the species
      flbcs(m)%const_idx = -1
      const_loop: do n = 1, size(const_props)
        call const_props(n)%diagnostic_name(diag_name, errcode=errflg, errmsg=errmsg)
        if (errflg /= 0) then
          return
        end if
        if (trim(diag_name) == trim(mw_species)) then
          call const_props(n)%const_index(flbcs(m)%const_idx, errflg, errmsg)
          if (errflg /= 0) then
            return
          end if
          call const_props(n)%is_advected(flbcs(m)%is_advected, errflg, errmsg)
          if (errflg /= 0) then
            return
          end if
          call const_props(n)%is_mass_mixing_ratio(is_mass_mmr, errflg, errmsg)
          if (errflg /= 0) then
            return
          end if
          call const_props(n)%is_dry(is_dry_mmr, errflg, errmsg)
          if (errflg /= 0) then
            return
          end if
          exit const_loop
        end if
      end do const_loop

      if (flbcs(m)%const_idx < 0) then
        errmsg = subname // ': no constituent found for flbc_list member "' // &
                 trim(flbcs(m)%species) // '"'
        errflg = 1
        return
      end if

      ! Bottom-boundary pinning writes the constituent as a dry mass mixing ratio;
      ! require consistency:
      if (flbcs(m)%is_advected .and. .not. (is_mass_mmr .and. is_dry_mmr)) then
        errmsg = subname // ': advected constituent for flbc_list member "' // &
                 trim(flbcs(m)%species) // '" is not a dry mass mixing ratio;' // &
                 ' lower-boundary pinning is only implemented for that convention'
        errflg = 1
        return
      end if

      if (any(flbcs(1:m - 1)%const_idx == flbcs(m)%const_idx)) then
        errmsg = subname // ': flbc_list members "' // trim(flbcs(m)%species) // &
                 '" and another entry (CFC11 and CFC11eq?) fill the same constituent'
        errflg = 1
        return
      end if

      ! conversion factor from vmr to mass mixing ratio w.r.t. dry air:
      ! get_molar_mass_ratio returns (molar mass of dry air)/(molar mass of species)
      call get_molar_mass_ratio(trim(mw_species), flbcs(m)%mmr_factor, errmsg, errflg)
      if (errflg /= 0) then
        return
      end if
      flbcs(m)%mmr_factor = 1._kind_phys/flbcs(m)%mmr_factor

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
    ncol, pver, lat, lon, constituents, &
    errmsg, errflg)

    ! CAM-SIMA host model dependency for the model date
    use time_manager,  only: get_curr_date

    integer,            intent(in)    :: ncol   ! number of columns [count]
    integer,            intent(in)    :: pver   ! number of vertical layers [count]
    real(kind_phys),    intent(in)    :: lat(:) ! latitude of columns [rad]
    real(kind_phys),    intent(in)    :: lon(:) ! longitude of columns [rad]
    real(kind_phys),    intent(inout) :: constituents(:,:,:) ! constituent array (ncol, pver, pcnst)

    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! local variables
    integer           :: m
    integer           :: yr, mon, day, ncsec, ncdate
    integer           :: last, next
    real(kind_phys)   :: dels
    real(kind_phys)   :: wrk_time
    real(kind_phys)   :: vmr_gmean

    errmsg = ''
    errflg = 0

    if (flbc_cnt < 1) then
      return
    end if

    call get_curr_date(yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day

    ! Advance the SERIAL read window when current time passes it
    if (lbc_type == 'SERIAL') then
      call flt_date(ncdate, ncsec, wrk_time)
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

    ! Update each species' constituent converted to a dry mass mixing ratio:
    ! non-advected constituents get a whole-column uniform fill with the global mean vmr;
    ! advected constituents get the per-column value pinned into the bottom two levels.
    call get_dels(ncdate, ncsec, dels, last, next, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if

    do m = 1, flbc_cnt
      if (flbcs(m)%is_advected) then
        constituents(:, pver, flbcs(m)%const_idx) = &
             (flbcs(m)%vmr(:, last) &
              + dels*(flbcs(m)%vmr(:, next) - flbcs(m)%vmr(:, last)))*flbcs(m)%mmr_factor
        constituents(:, pver - 1, flbcs(m)%const_idx) = constituents(:, pver, flbcs(m)%const_idx)
      else
        call global_mean_vmr(flbcs(m), ncol, dels, last, next, vmr_gmean)
        constituents(:, :, flbcs(m)%const_idx) = vmr_gmean*flbcs(m)%mmr_factor
      end if
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

  ! Convert an integer date and seconds of day to a float time coordinate
  ! consistent with the model calendar.
  subroutine flt_date(ncdate, ncsec, time)

    ! Host model time dependency...
    use time_manager, only: set_time_float_from_date

    integer,         intent(in)  :: ncdate ! date YYYYMMDD
    integer,         intent(in)  :: ncsec  ! seconds of day
    real(kind_phys), intent(out) :: time

    call set_time_float_from_date(time, ncdate/10000, mod(ncdate, 10000)/100, mod(ncdate, 100), ncsec)

  end subroutine flt_date

end module prescribe_lower_boundary_conditions
