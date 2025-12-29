! Manages reading and interpolation of prescribed aerosol deposition
! fluxes. These are the deposition fluxes sent to the surface model.
!
! Note: this module currently only implements bulk aerosol deposition fluxes.
!
! This scheme requires the concurrent use of prescribed_aerosols scheme
! to provide prescribed aerosol concentrations.
!
! Based on original CAM version from: Francis Vitt
module prescribed_aerosol_deposition_flux
  use ccpp_kinds, only: kind_phys

  ! CAM-SIMA host model dependency to read aerosol data.
  use tracer_data, only: trfile      ! data information and file read state.
  use tracer_data, only: trfld       ! tracer data container.

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: prescribed_aerosol_deposition_flux_init
  public :: prescribed_aerosol_deposition_flux_run

  ! fields to store tracer_data state and information.
  type(trfld), pointer :: tracer_data_fields(:)
  type(trfile)         :: tracer_data_file

  ! is this module active?
  logical :: has_aerodep_flx = .false.

  ! for bulk aerosol fluxes
  integer, parameter :: N_BULK = 14
  logical :: is_bulk = .false.

  ! list of bulk aerosol flux field names
  ! these and the index fields below are explicitly enumerated
  ! because they are explicitly enumerated in the cam_out coupler data structure.
  ! no runtime configurability is achievable here if cam_out hardcodes these.
  character(len=12), parameter :: bulk_names(N_BULK) = [&
                                  'BCDEPWET    ', 'BCPHODRY    ', 'BCPHIDRY    ', &
                                  'OCDEPWET    ', 'OCPHODRY    ', 'OCPHIDRY    ', &
                                  'DSTX01DD    ', 'DSTX02DD    ', 'DSTX03DD    ', 'DSTX04DD    ', &
                                  'DSTX01WD    ', 'DSTX02WD    ', 'DSTX03WD    ', 'DSTX04WD    ']
  integer :: index_bulk_map(N_BULK)
  integer :: ibcphiwet, ibcphidry, ibcphodry
  integer :: iocphiwet, iocphidry, iocphodry
  integer :: idst1dry, idst2dry, idst3dry, idst4dry
  integer :: idst1wet, idst2wet, idst3wet, idst4wet

  ! for modal aerosol fluxes (unavailable)
  integer, parameter :: N_MODAL = 22
  logical :: is_modal = .false.

contains
!> \section arg_table_prescribed_aerosol_deposition_flux_init  Argument Table
!! \htmlinclude prescribed_aerosol_deposition_flux_init.html
  subroutine prescribed_aerosol_deposition_flux_init( &
    amIRoot, iulog, &
    specifier, &
    filename, filelist, datapath, &
    data_type, &
    cycle_yr, fixed_ymd, fixed_tod, &
    prescribed_aero_model, &
    errmsg, errflg)

    !
    use tracer_data, only: trcdata_init

    ! host model dependency for history diagnostics
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: specifier(:)              ! field specifiers for tracer data
    character(len=*),   intent(in)  :: filename                  ! input filename
    character(len=*),   intent(in)  :: filelist                  ! input filelist
    character(len=*),   intent(in)  :: datapath                  ! input datapath
    character(len=*),   intent(in)  :: data_type                 ! data type
    integer,            intent(in)  :: cycle_yr                  ! cycle year
    integer,            intent(in)  :: fixed_ymd                 ! fixed year-month-day (YYYYMMDD) [1]
    integer,            intent(in)  :: fixed_tod                 ! fixed time of day [s]
    character(len=*),   intent(in)  :: prescribed_aero_model     ! type of aerosol representation
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: i, ndx, number_flds
    character(len=*), parameter :: subname = 'prescribed_aerosol_deposition_flux_init'

    errmsg = ''
    errflg = 0

    ! check if user has specified an input dataset
    if (filename /= 'UNSET' .and. len_trim(filename) > 0) then
      has_aerodep_flx = .true.

      if (amIRoot) then
        write(iulog,*) subname//': aerosol deposition fluxes prescribed in: '//trim(filename)
        write(iulog,*) subname//': aerosol representation is: '//trim(prescribed_aero_model)
      end if
    else
      return
    end if

    select case(prescribed_aero_model)
      case('bulk')
        is_bulk = .true.
      case('modal')
        is_modal = .true.
        errflg = 1
        errmsg = subname//': modal aerosol mode not yet implemented'
        return
      case default
        errflg = 1
        errmsg = subname//': ERROR: unknown aerosol representation: "'//trim(prescribed_aero_model)//'"'
    end select

    ! count number of specifiers
    number_flds = 0
    specifier_count_loop: do i = 1, size(specifier)
      ! remove empty or unset specifiers
      if(len_trim(specifier(i)) > 0 .and. trim(specifier(i)) /= 'UNSET') then
        number_flds = number_flds + 1
      else
        ! Assume all remaining specifier entries are also not set.
        exit specifier_count_loop
      end if
    end do specifier_count_loop

    ! initialize dataset in tracer_data module.
    call trcdata_init( &
         specifier      = specifier(:number_flds), &
         filename       = filename, &
         filelist       = filelist, &
         datapath       = datapath, &
         flds           = tracer_data_fields, & ! ptr
         file           = tracer_data_file, &
         data_cycle_yr  = cycle_yr, &
         data_fixed_ymd = fixed_ymd, &
         data_fixed_tod = fixed_tod, &
         data_type      = data_type)

    number_flds = 0
    if (associated(tracer_data_fields)) number_flds = size(tracer_data_fields)

    if (number_flds < 1) then
      has_aerodep_flx = .false.
      if (amIRoot) then
        write (iulog, *) subname//': no aerosol deposition fluxes have been specified'
      end if
      return
    end if

    ! map field names to indices
    index_bulk_map(:) = -1

    do i = 1, number_flds
      ! match specifier field name to known bulk aerosol names
      ndx = get_ndx(tracer_data_fields(i)%fldnam, bulk_names)
      if (ndx > 0) then
        index_bulk_map(ndx) = i

        ! add history field for diagnostics
        call history_add_field(trim(tracer_data_fields(i)%fldnam)//'_D', &
                               'prescribed aero deposition flux ' // trim(tracer_data_fields(i)%fldnam), &
                               horiz_only, 'avg', &
                               tracer_data_fields(i)%units)
      else
        ! when modal aerosols are available can match against modal names.

        ! if all else fails
        errflg = 1
        errmsg = 'prescribed_aerosol_deposition_flux_init: aerosol flux name not recognized: ' &
                 //trim(tracer_data_fields(i)%fldnam)
        return
      end if
    end do

    ! set up index pointers for bulk fluxes
    ! these correspond exactly to the order in bulk_names
    ibcphiwet = index_bulk_map(1)
    ibcphodry = index_bulk_map(2)
    ibcphidry = index_bulk_map(3)
    iocphiwet = index_bulk_map(4)
    iocphodry = index_bulk_map(5)
    iocphidry = index_bulk_map(6)
    idst1dry  = index_bulk_map(7)
    idst2dry  = index_bulk_map(8)
    idst3dry  = index_bulk_map(9)
    idst4dry  = index_bulk_map(10)
    idst1wet  = index_bulk_map(11)
    idst2wet  = index_bulk_map(12)
    idst3wet  = index_bulk_map(13)
    idst4wet  = index_bulk_map(14)

  end subroutine prescribed_aerosol_deposition_flux_init

!> \section arg_table_prescribed_aerosol_deposition_flux_run  Argument Table
!! \htmlinclude prescribed_aerosol_deposition_flux_run.html
  subroutine prescribed_aerosol_deposition_flux_run( &
    ncol, &
    pmid, pint, phis, zi, & ! necessary fields for trcdata read.
    bcphiwet, bcphidry, bcphodry, &
    ocphiwet, ocphidry, ocphodry, &
    dst1dry, dst2dry, dst3dry, dst4dry, &
    dst1wet, dst2wet, dst3wet, dst4wet, &
    errmsg, errflg)

    use tracer_data, only: advance_trcdata
    use cam_history, only: history_out_field

    ! Input arguments:
    integer,            intent(in)  :: ncol
    real(kind_phys),    intent(in)  :: pmid(:, :)  ! air pressure at centers [Pa]
    real(kind_phys),    intent(in)  :: pint(:, :)  ! air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)  :: phis(:)     ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)  :: zi(:, :)    ! height above surface, interfaces [m]

    ! Output arguments to coupler:
    real(kind_phys),    intent(out) :: bcphiwet(:) ! wet deposition of hydrophilic black carbon [kg m-2 s-1]
    real(kind_phys),    intent(out) :: bcphidry(:) ! dry deposition of hydrophilic black carbon [kg m-2 s-1]
    real(kind_phys),    intent(out) :: bcphodry(:) ! dry deposition of hydrophobic black carbon [kg m-2 s-1]
    real(kind_phys),    intent(out) :: ocphiwet(:) ! wet deposition of hydrophilic organic carbon [kg m-2 s-1]
    real(kind_phys),    intent(out) :: ocphidry(:) ! dry deposition of hydrophilic organic carbon [kg m-2 s-1]
    real(kind_phys),    intent(out) :: ocphodry(:) ! dry deposition of hydrophobic organic carbon [kg m-2 s-1]
    real(kind_phys),    intent(out) :: dst1dry(:)  ! dry deposition of dust bin 1 [kg m-2 s-1]
    real(kind_phys),    intent(out) :: dst2dry(:)  ! dry deposition of dust bin 2 [kg m-2 s-1]
    real(kind_phys),    intent(out) :: dst3dry(:)  ! dry deposition of dust bin 3 [kg m-2 s-1]
    real(kind_phys),    intent(out) :: dst4dry(:)  ! dry deposition of dust bin 4 [kg m-2 s-1]
    real(kind_phys),    intent(out) :: dst1wet(:)  ! wet deposition of dust bin 1 [kg m-2 s-1]
    real(kind_phys),    intent(out) :: dst2wet(:)  ! wet deposition of dust bin 2 [kg m-2 s-1]
    real(kind_phys),    intent(out) :: dst3wet(:)  ! wet deposition of dust bin 3 [kg m-2 s-1]
    real(kind_phys),    intent(out) :: dst4wet(:)  ! wet deposition of dust bin 4 [kg m-2 s-1]

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! initialize outputs to zero
    bcphiwet(:) = 0._kind_phys
    bcphidry(:) = 0._kind_phys
    bcphodry(:) = 0._kind_phys
    ocphiwet(:) = 0._kind_phys
    ocphidry(:) = 0._kind_phys
    ocphodry(:) = 0._kind_phys
    dst1dry (:) = 0._kind_phys
    dst2dry (:) = 0._kind_phys
    dst3dry (:) = 0._kind_phys
    dst4dry (:) = 0._kind_phys
    dst1wet (:) = 0._kind_phys
    dst2wet (:) = 0._kind_phys
    dst3wet (:) = 0._kind_phys
    dst4wet (:) = 0._kind_phys

    if (.not. has_aerodep_flx) return

    if (is_modal) then
      ! modal aerosol mode not yet implemented
      errflg = 1
      errmsg = 'prescribed_aerosol_deposition_flux_run: modal aerosol mode not yet implemented'
      return
    end if

    ! advance data in tracer_data to current time.
    call advance_trcdata(tracer_data_fields, tracer_data_file, &
                         pmid, pint, phis, zi)

    if(is_bulk) then
      ! set bulk aerosol fluxes
      call set_fluxes(bcphiwet, ibcphiwet)
      call set_fluxes(bcphidry, ibcphidry)
      call set_fluxes(bcphodry, ibcphodry)
      call set_fluxes(ocphiwet, iocphiwet)
      call set_fluxes(ocphidry, iocphidry)
      call set_fluxes(ocphodry, iocphodry)
      call set_fluxes(dst1dry,  idst1dry)
      call set_fluxes(dst2dry,  idst2dry)
      call set_fluxes(dst3dry,  idst3dry)
      call set_fluxes(dst4dry,  idst4dry)
      call set_fluxes(dst1wet,  idst1wet)
      call set_fluxes(dst2wet,  idst2wet)
      call set_fluxes(dst3wet,  idst3wet)
      call set_fluxes(dst4wet,  idst4wet)
    end if

  end subroutine prescribed_aerosol_deposition_flux_run

  ! helper subroutine to set fluxes based on target (hardcoded) array and tracer data index
  subroutine set_fluxes(fluxes, fld_idx)
    ! host model dependency for history output
    use cam_history, only: history_out_field

    real(kind_phys), intent(inout) :: fluxes(:)
    integer,         intent(in)    :: fld_idx

    integer :: i

    if (fld_idx < 1) return

    ! copy from tracer_data container (surface field only)
    fluxes(:) = tracer_data_fields(fld_idx)%data(:, 1)

    ! output diagnostic
    call history_out_field(trim(tracer_data_fields(fld_idx)%fldnam)//'_D', &
                           fluxes(:))

  end subroutine set_fluxes

  pure integer function get_ndx(name, list)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: list(:)

    integer :: i
    integer :: maxnum

    maxnum = size(list)

    get_ndx = -1
    do i = 1, maxnum
      if (trim(name) == trim(list(i))) then
        get_ndx = i
        return
      end if
    end do

  end function get_ndx
end module prescribed_aerosol_deposition_flux
