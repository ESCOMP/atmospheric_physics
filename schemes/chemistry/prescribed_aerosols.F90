! Manages reading and interpolation of prescribed aerosol concentrations.
!
! This module uses CCPP constituents (non-advected) to store prescribed aerosol fields
! for access by radiation and other schemes.
!
! The prescribed_aero_specifier namelist variable specifies a list of
! variable names of the concentration fields in the netCDF dataset (ncdf_fld_name)
! and the corresponding constituent names, both of which are arbitrary:
!
! prescribed_aero_specifier = 'constituent_name1:ncdf_fld_name1','constituent_name2:ncdf_fld_name2', ...
!
! If there is no ":" then the specified name is used as both the
! constituent_name and ncdf_fld_name
!
! Based on original CAM version from: Francis Vitt
module prescribed_aerosols

  use ccpp_kinds, only: kind_phys

  ! CAM-SIMA host model dependency to read aerosol data
  use tracer_data, only: trfile     ! data information and file read state
  use tracer_data, only: trfld      ! tracer data container

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: prescribed_aerosols_register
  public :: prescribed_aerosols_init
  public :: prescribed_aerosols_run

  ! maximum of specified aerosol fields
  ! if need to be expanded, need to change namelist definition as well.
  integer,   parameter :: N_AERO_MAX = 50

  ! fields to store tracer_data state and information
  type(trfld), pointer :: tracer_data_fields(:)
  type(trfile)         :: tracer_data_file

  ! local derived type to store constituent name to
  ! netCDF field name and tracer_data read index mapping:
  type :: aero_constituent_map
    character(len=64)  :: constituent_name       ! CCPP constituent standard name

    ! if modal aerosols and is an interstitial (_a) species, then the species
    ! is diagnosed from log-mean (_logm) and log-variance (_logv) with a randomized
    ! log-normal distribution;
    ! in this case, field_index is used for log-mean, and a separate field for variance
    ! is stored here.
    logical            :: is_modal_aero_interstitial
    character(len=64)  :: trcdata_field_name_logv! tracer_data field name (fldnam, not netCDF variable name) (log-variance)
    integer            :: field_index_logv       ! index into tracer_data_fields array (log-variance)

    ! all other cases, including modal log-mean:
    character(len=64)  :: trcdata_field_name     ! tracer_data field name (fldnam, not netCDF variable name)
    integer            :: field_index            ! index into tracer_data_fields array
  end type aero_constituent_map

  type(aero_constituent_map), allocatable :: aero_map_list(:)

  ! module state variables
  logical :: has_prescribed_aerosols = .false.
  logical :: clim_modal_aero         = .false.
  integer :: aero_cnt                       ! # of aerosol constituents
  integer :: aero_cnt_c                     ! # of cloud-borne species (for modal aerosols only)

  ! Normal random number which persists from one timestep to the next
  ! (used for modal aerosol sampling)
  real(kind_phys) :: randn_persists = 0.0_kind_phys

  ! TODO: infrastructure for writing (and reading) randn_persists to persist during restart runs.
  ! see CAM/prescribed_aero::{read,write}_prescribed_aero_restart
  ! !!! Restarts will not be bit-for-bit without this !!!
  ! TODO when SIMA implements restarts.

contains

  ! Register prescribed aerosols in constituents object.
!> \section arg_table_prescribed_aerosols_register Argument Table
!! \htmlinclude prescribed_aerosols_register.html
  subroutine prescribed_aerosols_register( &
    amIRoot, iulog, &
    prescribed_aero_specifier, &
    prescribed_aero_file, &
    prescribed_aero_filelist, &
    prescribed_aero_datapath, &
    prescribed_aero_type, &
    prescribed_aero_cycle_yr, &
    prescribed_aero_fixed_ymd, &
    prescribed_aero_fixed_tod, &
    prescribed_aero_model, &
    aerosol_constituents, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

    ! Input arguments:
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_aero_specifier(:)     ! aerosol specifiers from namelist
    character(len=*),   intent(in)  :: prescribed_aero_file             ! input filename from namelist
    character(len=*),   intent(in)  :: prescribed_aero_filelist         ! input filelist from namelist
    character(len=*),   intent(in)  :: prescribed_aero_datapath         ! input datapath from namelist
    character(len=*),   intent(in)  :: prescribed_aero_type             ! data type from namelist
    integer,            intent(in)  :: prescribed_aero_cycle_yr         ! cycle year from namelist [1]
    integer,            intent(in)  :: prescribed_aero_fixed_ymd        ! fixed year-month-day from namelist (YYYYMMDD) [1]
    integer,            intent(in)  :: prescribed_aero_fixed_tod        ! fixed time of day from namelist [s]
    character(len=*),   intent(in)  :: prescribed_aero_model            ! type of aerosol representation [none]

    ! Output arguments:
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: aerosol_constituents(:) ! prescribed aero runtime CCPP constituents
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables:
    character(len=256) :: tmpstr
    character(len=64)  :: constituent_name
    character(len=64)  :: ncdf_fld_name
    character(len=64)  :: unit_name
    integer            :: strlen
    integer            :: idx, idx2
    integer            :: aero_idx
    integer            :: i, j
    logical            :: is_modal_aero_interstitial
    logical            :: skip_spec

    character(len=*), parameter :: subname = 'prescribed_aerosols_register'

    errmsg = ''
    errflg = 0

    ! Store module-level settings
    clim_modal_aero = (trim(prescribed_aero_model) == "modal")

    ! Check if prescribed aerosols are enabled
    if (prescribed_aero_file == 'UNSET' .or. &
        len_trim(prescribed_aero_file) == 0) then
      has_prescribed_aerosols = .false.
      if (amIRoot) then
        write(iulog,*) subname//': No prescribed aerosols specified'
      end if
      return
    end if

    ! Parse the aerosol format specifier from namelist into mapping ddt.
    ! We need two scans. First to determine the count, the second to populate
    ! the information into the ddt.
    aero_cnt   = 0
    aero_cnt_c = 0 ! cloud borne species count
    cnt_loop: do i = 1, N_AERO_MAX
      ! FIXME: should I be responsible for handling this? I feel like I should not handle this
      if(prescribed_aero_specifier(i) == 'UNSET' .or. &
         len_trim(prescribed_aero_specifier(i)) == 0) exit cnt_loop

      skip_spec = .false.
      if(clim_modal_aero) then
        ! For modal aerosols, interstitial species (*_a) are diagnosed from
        ! their *_logm and *_logv counterparts (e.g. soa_a1 is diagnosed from
        ! soa_a1_logm and soa_a1_logv). Therefore, only *_logm and *_logv and cloud
        ! borne (*_c) species are specified in the build-namelist.

        ! In the following cnt_loop, we will count the cloud borne species and *_logm species
        ! (in lieu of *_a species). We will skip *_logv species.
        ! This will ensure that aero_cnt variable is the sum of cloud borne and
        ! interstitial species (which we will manually add in the names to the ddt later).
        ! We are also counting cloud borne (*_c) species which will help
        ! adding the same number of interstitial species to the ddt.
        !
        ! For modal aerosols, skip counting species ending with *_logv
        if(index(prescribed_aero_specifier(i),'_c') >= 1)    aero_cnt_c = aero_cnt_c + 1
        if(index(prescribed_aero_specifier(i),'_logv') >= 1) skip_spec = .true.
      endif

      if(.not. skip_spec) aero_cnt = aero_cnt + 1
    end do cnt_loop

    if(aero_cnt == 0) then
      has_prescribed_aerosols = .false.
      return
    endif

    has_prescribed_aerosols = .true.

    ! Allocate mapping list of ddt
    allocate(aero_map_list(aero_cnt), stat=errflg, errmsg=errmsg)
    if(errflg /= 0) then
      return
    end if

    ! Now populate the information into the ddt.
    !
    ! In CAM modal aerosols, the interstitial species (_a) are added at the end of the cloud borne (_c) species.
    ! TODO hplin 10/22/25 We might not need this? Fields were identified by their pbuf name anyway

    aero_idx = 1 ! pointer for current aero_map_list index.
    ddt_loop: do i = 1, N_AERO_MAX
      ! Parse specifier
      tmpstr = trim(adjustl(prescribed_aero_specifier(i)))

      ! FIXME: should I be responsible for handling this? I feel like I should not handle this
      if(tmpstr == 'UNSET' .or. len_trim(tmpstr) == 0) exit ddt_loop

      idx = index(tmpstr, ':')
      if(idx > 0) then
        ! format as constituent_name:ncdf_name
        constituent_name = tmpstr(:idx-1)
        ncdf_fld_name    = tmpstr(idx+1:)
      else
        ! format as single_name_for_both
        constituent_name = tmpstr
        ncdf_fld_name    = tmpstr
      end if

      is_modal_aero_interstitial = .false.

      ! For modal aerosols only,
      ! skip the _logv species; remove the _logm suffix to directly get the interstitial (_a) name
      ! e.g. in loop: num_c1, num_a1_logm, num_a1_logv, ... --> num_c1, num_a1, (skip), ...
      if(clim_modal_aero) then
        if(index(constituent_name,'_logv') >= 1) cycle ddt_loop
        idx = index(constituent_name,'_logm')
        if(idx > 0) then
          ! use this as the proxy for interstitial species and initialize it here
          is_modal_aero_interstitial = .true.
          aero_map_list(aero_idx)%is_modal_aero_interstitial = .true.

          aero_map_list(aero_idx)%constituent_name           = constituent_name(:idx-1) ! removed _logm
          aero_map_list(aero_idx)%trcdata_field_name         = constituent_name  ! log-mean specifier name

          ! construct log-variance specifier name by changing _logm with _logv
          aero_map_list(aero_idx)%trcdata_field_name_logv    = constituent_name(:idx-1) // '_logv'

          ! The tracer data specifier index will be rescanned once tracer_data_fields is initialized
          ! at the init phase.
          aero_map_list(aero_idx)%field_index_logv     = -1
        end if ! logm
      end if ! clim_modal_aero

      if(.not. is_modal_aero_interstitial) then
        ! Store mapping for all other cases except modal aerosol interstitial species.
        aero_map_list(aero_idx)%is_modal_aero_interstitial = .false.

        aero_map_list(aero_idx)%constituent_name = constituent_name
        aero_map_list(aero_idx)%trcdata_field_name = constituent_name

        ! The tracer data specifier index will be rescanned once tracer_data_fields is initialized
        ! at the init phase.
        aero_map_list(aero_idx)%field_index      = -1
      end if

      ! We added a new aero_map_list entry, so advance the pointer.
      aero_idx = aero_idx + 1
    end do ddt_loop

    ! Sanity check
    if(aero_idx /= aero_cnt+1) then
      errflg = 1
      write(errmsg,*) subname//': consistency check 1 failure; at the end of ddt allocation, aero_idx is not aero_cnt+1', aero_idx, aero_cnt
      return
    end if

    ! Allocate CCPP dynamic constituents object for prescribed aerosols.
    allocate(aerosol_constituents(aero_cnt), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) return

    ! Now register constituents in the CCPP constituent properties object.
    reg_loop: do i = 1, aero_cnt
      ! check units. at this point, we do not know the units from file
      ! because tracer_data has not read any data yet.
      ! number concentrations are units of 1 kg-1; all others are kg kg-1
      if(index(aero_map_list(i)%constituent_name, 'num_') == 1) then
        unit_name = '1 kg-1'
      else
        unit_name = 'kg kg-1'
      end if

      call aerosol_constituents(i)%instantiate( &
           std_name          = trim(aero_map_list(i)%constituent_name), &
           long_name         = 'prescribed aerosol '//trim(aero_map_list(i)%constituent_name), &
           units             = unit_name, &
           vertical_dim      = 'vertical_layer_dimension', &
           min_value         = 0.0_kind_phys, &
           advected          = .false., &
           water_species     = .false., &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if(errflg /= 0) return
    end do reg_loop

    if (amIRoot) then
      write(iulog,*) trim(subname)//': Registered ', aero_cnt, ' prescribed aerosol constituents'
    end if

  end subroutine prescribed_aerosols_register

  ! Initialize prescribed aerosol reading via tracer_data.
!> \section arg_table_prescribed_aerosols_init Argument Table
!! \htmlinclude prescribed_aerosols_init.html
  subroutine prescribed_aerosols_init( &
    amIRoot, iulog, &
    prescribed_aero_specifier, &
    prescribed_aero_file, &
    prescribed_aero_filelist, &
    prescribed_aero_datapath, &
    prescribed_aero_type, &
    prescribed_aero_cycle_yr, &
    prescribed_aero_fixed_ymd, &
    prescribed_aero_fixed_tod, &
    errmsg, errflg)

    ! host model dependency for tracer_data read utility
    use tracer_data, only: trcdata_init

    ! host model dependency for history output
    use cam_history, only: history_add_field

    ! Input arguments:
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_aero_specifier(:)     ! aerosol specifier from namelist
    character(len=*),   intent(in)  :: prescribed_aero_file             ! input filename from namelist
    character(len=*),   intent(in)  :: prescribed_aero_filelist         ! input filelist from namelist
    character(len=*),   intent(in)  :: prescribed_aero_datapath         ! input datapath from namelist
    character(len=*),   intent(in)  :: prescribed_aero_type             ! data type from namelist
    integer,            intent(in)  :: prescribed_aero_cycle_yr         ! cycle year from namelist [1]
    integer,            intent(in)  :: prescribed_aero_fixed_ymd        ! fixed year-month-day from namelist (YYYYMMDD) [1]
    integer,            intent(in)  :: prescribed_aero_fixed_tod        ! fixed time of day from namelist [s]

    ! Output arguments:
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables:
    integer            :: i
    integer            :: idx
    character(len=64)  :: unit_name

    ! Local parameters:
    character(len=*), parameter :: subname = 'prescribed_aerosols_init'

    ! Initialize output arguments:
    errmsg = ''
    errflg = 0

    if (.not. has_prescribed_aerosols) return

    if (amIRoot) then
      write(iulog,*) trim(subname)//': aerosols are prescribed in: '//trim(prescribed_aero_file)
    end if

    ! Initialize tracer_data module with file and field information
    call trcdata_init( &
      specifier      = prescribed_aero_specifier(:aero_cnt), &
      filename       = prescribed_aero_file, &
      filelist       = prescribed_aero_filelist, &
      datapath       = prescribed_aero_datapath, &
      flds           = tracer_data_fields, &
      file           = tracer_data_file, &
      data_cycle_yr  = prescribed_aero_cycle_yr, &
      data_fixed_ymd = prescribed_aero_fixed_ymd, &
      data_fixed_tod = prescribed_aero_fixed_tod, &
      data_type      = prescribed_aero_type)

    ! Verify tracer_data is correctly initialized
    if (.not. associated(tracer_data_fields)) then
      errflg = 1
      errmsg = subname//': tracer_data_fields not associated after trcdata_init'
      return
    end if

    ! Note: because in modal aerosols, interstitial fields are derived from the
    ! log-mean and log-variance (logm, logv) fields, the number of tracer_data
    ! fields is not equal to aero_cnt(= a + c), but (logm + logv) + c.

    ! Based on aero_map_list, scan the tracer data fields to populate correct field_index
    do i = 1, aero_cnt
      ! Find the matching field in tracer_data_fields for the primary field
      do idx = 1, size(tracer_data_fields)
        if (trim(tracer_data_fields(idx)%fldnam) == trim(aero_map_list(i)%trcdata_field_name)) then
          aero_map_list(i)%field_index = idx
          exit
        end if
      end do

      ! For modal aerosol interstitial species, also find the logv field
      if (aero_map_list(i)%is_modal_aero_interstitial) then
        do idx = 1, size(tracer_data_fields)
          if (trim(tracer_data_fields(idx)%fldnam) == trim(aero_map_list(i)%trcdata_field_name_logv)) then
            aero_map_list(i)%field_index_logv = idx
            exit
          end if
        end do
      end if
    end do

    ! Check aero_map_list for any unpopulated field indices (consistency check),
    ! Register history field, and
    ! Print out each aero_map_list field information
    do i = 1, aero_cnt
      if (aero_map_list(i)%field_index <= 0) then
        errflg = 1
        write(errmsg, '(3a)') trim(subname), ': Field not found in tracer_data for constituent: ', &
             trim(aero_map_list(i)%constituent_name)
        return
      end if

      ! Add history field
      ! Check units. at this point, we do not know the units from file
      ! because tracer_data has not read any data yet.
      ! number concentrations are units of 1 kg-1; all others are kg kg-1
      if(index(aero_map_list(i)%constituent_name, 'num_') == 1) then
        unit_name = '1 kg-1'
      else
        unit_name = 'kg kg-1'
      end if
      call history_add_field(trim(aero_map_list(i)%constituent_name) // '_D', &
            'prescribed aero ' // trim(aero_map_list(i)%constituent_name), &
            'lev', 'avg', unit_name)

      ! Informational printout
      if (amIRoot) then
        if (aero_map_list(i)%is_modal_aero_interstitial) then
          if (aero_map_list(i)%field_index_logv <= 0) then
            errflg = 1
            write(errmsg, '(3a)') trim(subname), ': logv field not found for interstitial constituent: ', &
                 trim(aero_map_list(i)%constituent_name)
            return
          end if

          write(iulog, '(a,i3,2a)') '  ', i, ': ', trim(aero_map_list(i)%constituent_name)
          write(iulog, '(3a,i3)') ' log-mean field: ', trim(aero_map_list(i)%trcdata_field_name), &
               ' (trcdata index ', aero_map_list(i)%field_index, ')'
          write(iulog, '(3a,i3)') ' log-variance field: ', trim(aero_map_list(i)%trcdata_field_name_logv), &
               ' (trcdata index ', aero_map_list(i)%field_index_logv, ')'
        else
          write(iulog, '(a,i3,5a,i3,a)') '  ', i, ': ', trim(aero_map_list(i)%constituent_name), &
               ' from ', trim(aero_map_list(i)%trcdata_field_name), ' (trcdata index ', aero_map_list(i)%field_index, ')'
        end if
      end if
    end do

    if (amIRoot) then
      write(iulog,*) trim(subname)//': Initialized ', aero_cnt, ' aerosol fields'
    end if

  end subroutine prescribed_aerosols_init

!> \section arg_table_prescribed_aerosols_run Argument Table
!! \htmlinclude prescribed_aerosols_run.html
  subroutine prescribed_aerosols_run( &
    ncol, pver, pcnst, &
    const_props, &
    pi, &
    pmid, pint, phis, zi, & ! necessary fields for trcdata read.
    constituents, &
    errmsg, errflg)

    !
    use tracer_data,               only: advance_trcdata

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! dependency to get constituent index
    use ccpp_const_utils,          only: ccpp_const_get_idx

    ! host model dependency for history output
    use cam_history,               only: history_out_field

    ! Input arguments:
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: pcnst            ! # of CCPP constituents [count]
    type(ccpp_constituent_prop_ptr_t), &
                        intent(in)    :: const_props(:)   ! CCPP constituent properties pointer
    real(kind_phys),    intent(in)    :: pi
    real(kind_phys),    intent(in)    :: pmid(:,:)        ! air pressure [Pa]
    real(kind_phys),    intent(in)    :: pint(:,:)        ! air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)    :: phis(:)          ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zi(:,:)          ! height above surface, interfaces [m]

    ! Input/Output arguments:
    real(kind_phys),    intent(inout) :: constituents(:,:,:) ! constituent array (ncol, pver, pcnst) [kg kg-1 dry]

    ! Output arguments:
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local variables:
    integer            :: i
    integer            :: const_idx

    ! Local parameters:
    character(len=*), parameter :: subname = 'prescribed_aerosols_run'

    ! Initialize output arguments:
    errmsg = ''
    errflg = 0

    if (.not. has_prescribed_aerosols) return

    ! Advance tracer_data to current time
    call advance_trcdata(tracer_data_fields, tracer_data_file, &
                         pmid, pint, phis, zi)

    ! Copy the prescribed aerosol data to constituent array
    ! Loop over aero_map_list to retrieve data from tracer_data_fields.
    ! For most species (non-modal; cloud borne) just retrieve data and save to constituent.
    ! For interstitial species (is_modal_aero_interstitial) construct mixing ratio based on
    ! the logv and logm values read from tracer_data (port of rand_sample_prescribed_aero)
    do i = 1, aero_cnt
      ! Get constituent index
      call ccpp_const_get_idx(const_props, &
           trim(aero_map_list(i)%constituent_name), &
           const_idx, errmsg, errflg)
      if (errflg /= 0) return

      if (.not. aero_map_list(i)%is_modal_aero_interstitial) then
        ! For non-interstitial species (cloud-borne or bulk aerosols),
        ! directly copy field data from tracer_data container
        constituents(:ncol, :pver, const_idx) = &
             tracer_data_fields(aero_map_list(i)%field_index)%data(:ncol,:pver)
      else
        ! For modal aerosol interstitial species (_a), compute from log-normal distribution
        ! using log-mean (logm) and log-variance (logv)
        call compute_modal_aero_interstitial( &
             ncol, pver, &
             pi, &
             tracer_data_fields(aero_map_list(i)%field_index)%data(:ncol,:pver), &
             tracer_data_fields(aero_map_list(i)%field_index_logv)%data(:ncol,:pver), &
             constituents(:ncol, :pver, const_idx))
      end if

      ! History output
      call history_out_field(trim(aero_map_list(i)%constituent_name) // '_D', &
                             constituents(:ncol, :pver, const_idx))
    end do

  end subroutine prescribed_aerosols_run

  ! Compute modal aerosol interstitial species from log-normal distribution
  ! Port of rand_sample_prescribed_aero from CAM prescribed_aero.F90
  !
  ! Original authors: Balwinder Singh (12/14/2012) adapted from Jin-Ho Yoon
  subroutine compute_modal_aero_interstitial(ncol, pver, pi, logm_data, logv_data, mixing_ratio)

    ! Input arguments
    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    real(kind_phys), intent(in)  :: pi
    real(kind_phys), intent(in)  :: logm_data(:, :)           ! Log-mean of aerosol distribution [ln(kg kg-1)]
    real(kind_phys), intent(in)  :: logv_data(:, :)           ! Log-variance of aerosol distribution [ln(kg kg-1)^2]

    ! Output arguments
    real(kind_phys), intent(out) :: mixing_ratio(:, :)        ! Computed aerosol mixing ratio [kg kg-1]

    ! Local variables
    integer                      :: i, k
    real(kind_phys)              :: logm2                     ! Log-mean squared [ln(kg kg-1)^2]
    real(kind_phys)              :: variance                  ! Variance [(ln(kg kg-1))^2]
    real(kind_phys)              :: std                       ! Standard deviation [ln(kg kg-1)]
    real(kind_phys)              :: mean_max                  ! Maximum mean value [kg kg-1]
    real(kind_phys)              :: std_max                   ! Maximum std value [kg kg-1]
    real(kind_phys)              :: randn                     ! Normal random number [1]

    ! Local parameters
    real(kind_phys), parameter   :: mean_max_val = 5.0_kind_phys
    real(kind_phys), parameter   :: std_max_val  = 3.0_kind_phys

    ! Generate/use persisting random number
    randn = randn_prescribed_aero(pi)

    ! This loop logic is a good candidate for making into a pure elemental function
    ! despite the impurity of the RNG dependency above.
    do k = 1, pver
      do i = 1, ncol
        logm2 = logm_data(i, k) * logm_data(i, k)

        ! Compute (non-negative) variance
        variance = max(0.0_kind_phys, (logv_data(i, k) - logm2))

        ! Standard deviation
        std = sqrt(variance)

        ! Bounds to keep mixing ratios from going unphysical
        mean_max = exp(logm_data(i, k)) * mean_max_val
        std_max  = exp(logm_data(i, k) + std_max_val * std)

        ! Compute mixing ratio with random sampling
        mixing_ratio(i, k) = min(exp(logm_data(i, k) + randn * std), mean_max, std_max)
      end do
    end do
  end subroutine compute_modal_aero_interstitial

  ! Generate normally distributed random number that persists for entire day
  ! Port of randn_prescribed_aero and boxMuller from CAM prescribed_aero.F90
  !
  ! NOTE: I think this can be moved into the host model to be used as a
  ! common RNG for physics schemes to avoid having dependency on time_manager.
  function randn_prescribed_aero(pi) result(randn)
    use time_manager, only: is_end_curr_day, is_first_step, get_nstep

    real(kind_phys), intent(in)  :: pi

    ! Local variables
    integer                    :: i, seed_size, nstep
    integer,  allocatable      :: seed(:)
    real(kind_phys)            :: randu1, randu2          ! Uniform random numbers [1]
    real(kind_phys)            :: randn                   ! Normal random number [1]
    real(kind_phys)            :: ur, theta               ! Box-Muller variables [1]

    ! Local parameters for random seed generation
    integer, parameter         :: rconst1_1 = 5000000
    integer, parameter         :: rconst1_2 = 50
    integer, parameter         :: rconst2_1 = 10000000
    integer, parameter         :: rconst2_2 = 10

    ! Use same random number for the entire day; generate new one at start of new day
    if (is_first_step() .or. is_end_curr_day()) then
      ! Generate two uniformly distributed random numbers (between 0 and 1)
      call random_seed(size=seed_size)
      allocate(seed(seed_size))

      ! Using nstep as a seed to generate same sequence
      nstep = get_nstep()
      do i = 1, seed_size
        seed(i) = rconst1_1 * nstep + rconst1_2 * (i - 1)
      end do
      call random_seed(put=seed)
      call random_number(randu1)

      do i = 1, seed_size
        seed(i) = rconst2_1 * nstep + rconst2_2 * (i - 1)
      end do
      call random_seed(put=seed)
      call random_number(randu2)
      deallocate(seed)

      ! Box-Muller method: convert uniform to normal distribution (mean=0, std=1)
      ur    = sqrt(-2.0_kind_phys * log(randu1))
      theta = 2.0_kind_phys * pi * randu2
      randn = ur * cos(theta)

      ! Store for use throughout the day
      randn_persists = randn
    else
      ! Use the previously generated random number
      randn = randn_persists
    end if

  end function randn_prescribed_aero

end module prescribed_aerosols
