! Initialize prescribed aerosols from file.
! This is the "pbuf" equivalent of initialize_constituents:
! when aerosols are prescribed (in a "-chem none" configuration in CAM)
! they are not constituents, but pbuf_ fields instead (rad_climate N: specifier)
!
! In CAM-SIMA, the "provider" of chemistry/aerosol needs to register these species,
! thus a scheme needs to register these constituents based on pbuf_<aerosol>
! entries in the ncdata input file.
!
! This scheme should ONLY be used for snapshot runs.
! For non-snapshot runs, use the prescribed_aerosols and associated schemes.
module initialize_pbuf_aerosols

  implicit none
  private

  public :: initialize_pbuf_aerosols_register

  integer, parameter :: num_bam_species = 13
  character(len=5), parameter :: bam_species(num_bam_species) = &
    (/'sulf ', 'bcar1', 'bcar2', 'ocar1', 'ocar2', &
      'sslt1', 'sslt2', 'sslt3', 'sslt4', &
      'dust1', 'dust2', 'dust3', 'dust4'/)

  integer, parameter :: num_mam_species = 30
  character(len=6), parameter :: mam_species(num_mam_species) = &
    (/'num_a1', 'bc_a1 ', 'dst_a1', 'ncl_a1', 'pom_a1', 'so4_a1', 'soa_a1', &
      'num_a2', 'ncl_a2', 'so4_a2', 'soa_a2', &
      'num_a3', 'dst_a3', 'ncl_a3', 'so4_a3', &
      'num_c1', 'bc_c1 ', 'dst_c1', 'ncl_c1', 'pom_c1', 'so4_c1', 'soa_c1', &
      'num_c2', 'ncl_c2', 'so4_c2', 'soa_c2', &
      'num_c3', 'dst_c3', 'ncl_c3', 'so4_c3'/)

contains

!> \section arg_table_initialize_pbuf_aerosols_register  Argument Table
!! \htmlinclude initialize_pbuf_aerosols_register.html
  subroutine initialize_pbuf_aerosols_register(constituents, errmsg, errflg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use cam_initfiles,             only: initial_file_get_id
    use pio,                       only: file_desc_t, var_desc_t, &
                                         pio_inq_varid, pio_seterrorhandling, &
                                         PIO_BCAST_ERROR, PIO_NOERR
    use ccpp_kinds,                only: kind_phys
    use ccpp_chem_utils,           only: chem_constituent_qmin

    ! Dummy variables
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituents(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    type(file_desc_t), pointer :: ncdata
    type(var_desc_t)           :: varid
    integer                    :: ierr
    integer                    :: old_err_method
    integer                    :: i, nspecies
    logical                    :: is_bam, is_mam
    character(len=6), allocatable :: species_list(:)
    character(len=7)           :: unit_name

    character(len=*), parameter :: subname = 'initialize_pbuf_aerosols_register'

    errflg = 0
    errmsg = ''

    ncdata => initial_file_get_id()

    ! Suppress PIO fatal errors for variable probing
    call pio_seterrorhandling(ncdata, PIO_BCAST_ERROR, old_err_method)

    ! Detect aerosol mode by probing for indicator variables
    is_mam = .false.
    is_bam = .false.

    ierr = pio_inq_varid(ncdata, 'pbuf_num_a1', varid)
    if (ierr == PIO_NOERR) is_mam = .true.

    ierr = pio_inq_varid(ncdata, 'pbuf_sulf', varid)
    if (ierr == PIO_NOERR) is_bam = .true.

    ! Restore PIO error handling
    call pio_seterrorhandling(ncdata, old_err_method)

    if (.not. is_bam .and. .not. is_mam) then
      errflg = 1
      errmsg = subname//': No BAM (pbuf_sulf) or MAM (pbuf_num_a1) aerosol fields found on initial file. This scheme should NOT be used for non-snapshot runs or physics tests that do not use aerosol'
      return
    end if

    ! Build species list based on detected mode (MAM takes priority)
    if (is_mam) then
      nspecies = num_mam_species
      allocate(species_list(nspecies))
      species_list(:) = mam_species(:)
    else
      nspecies = num_bam_species
      allocate(species_list(nspecies))
      species_list(:) = bam_species(:)
    end if

    ! Register constituents
    allocate(constituents(nspecies), stat=ierr, errmsg=errmsg)
    if (ierr /= 0) then
      errflg = 1
      errmsg = subname // ': Failed to allocate constituents: ' // trim(errmsg)
      return
    end if

    do i = 1, nspecies
      if (index(species_list(i), 'num_') == 1) then
        unit_name = 'kg-1'
      else
        unit_name = 'kg kg-1'
      end if

      call constituents(i)%instantiate( &
        std_name          = trim(species_list(i)), &
        diag_name         = trim(species_list(i)), &
        long_name         = 'prescribed aerosol ' // trim(species_list(i)), &
        units             = trim(unit_name), &
        vertical_dim      = 'vertical_layer_dimension', &
        min_value         = chem_constituent_qmin(trim(species_list(i))), &
        advected          = .false., &
        water_species     = .false., &
        mixing_ratio_type = 'dry', &
        errcode           = errflg, &
        errmsg            = errmsg)
      if (errflg /= 0) return
    end do

    deallocate(species_list)

  end subroutine initialize_pbuf_aerosols_register

end module initialize_pbuf_aerosols
