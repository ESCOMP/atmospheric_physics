!> \file rrtmgp_sw_gas_optics.F90
!!

!> This module contains an init routine to initialize the shortwave gas optics object
!>  with data read in from file on the host side
module rrtmgp_sw_gas_optics

  implicit none
  private
  public :: rrtmgp_sw_gas_optics_init
  public :: rrtmgp_sw_gas_optics_run

contains
!> \section arg_table_rrtmgp_sw_gas_optics_init Argument Table
!! \htmlinclude rrtmgp_sw_gas_optics_init.html
!!
  subroutine rrtmgp_sw_gas_optics_init(kdist, sw_filename, available_gases, &
                  errmsg, errcode)
    use machine,                 only: kind_phys
    use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
    use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
    use radiation_tools,         only: check_error_msg
    use ccpp_io_reader,          only: abstract_netcdf_reader_t, create_netcdf_reader_t
    use mo_rte_kind,             only: wl

    ! Inputs
    character(len=*),                    intent(in) :: sw_filename                   ! Full path to RRTMGP shortwave coefficients file
    class(ty_gas_concs_ccpp),            intent(in) :: available_gases               ! Gas concentrations object
 
    ! Outputs
    class(ty_gas_optics_rrtmgp_ccpp),    intent(out) :: kdist                        ! RRTMGP gas optics object
    character(len=512),                  intent(out) :: errmsg                       ! CCPP error message
    integer,                             intent(out) :: errcode                      ! CCPP error code

    ! Local variables
    class(abstract_netcdf_reader_t),   allocatable :: file_reader
    character(len=:),  dimension(:),   pointer :: gas_names                          ! Names of absorbing gases
    character(len=:),  dimension(:),   pointer :: gas_minor                          ! Name of absorbing minor gas
    character(len=:),  dimension(:),   pointer :: identifier_minor                   ! Unique string identifying minor gas
    character(len=:),  dimension(:),   pointer :: minor_gases_lower                  ! Names of minor absorbing gases in lower atmosphere
    character(len=:),  dimension(:),   pointer :: minor_gases_upper                  ! Names of minor absorbing gases in upper atmosphere
    character(len=:),  dimension(:),   pointer :: scaling_gas_lower                  ! Absorption also depends on the concentration of this gas in the lower atmosphere
    character(len=:),  dimension(:),   pointer :: scaling_gas_upper                  ! Absorption also depends on the concentration of this gas in the upper atmosphere
    integer,        dimension(:,:,:),  pointer :: key_species                        ! Key species pair for each band
    integer,        dimension(:,:),    pointer :: band2gpt                           ! Array for converting shortwave band limits to g-points
    integer,        dimension(:,:),    pointer :: minor_limits_gpt_lower             ! Beginning and ending gpoint for each minor interval in lower atmosphere
    integer,        dimension(:,:),    pointer :: minor_limits_gpt_upper             ! Beginning and ending gpoint for each minor interval in upper atmosphere
    integer,        dimension(:),      pointer :: kminor_start_lower                 ! Starting index in the [1,nContributors] vector for a contributor given by "minor_gases_lower"
    integer,        dimension(:),      pointer :: kminor_start_upper                 ! Starting index in the [1,nContributors] vector for a contributor given by "minor_gases_upper"
    logical(wl),     dimension(:),     pointer :: minor_scales_with_density_lower    ! Density scaling is applied to minor absorption coefficients in the lower atmosphere
    logical(wl),     dimension(:),     pointer :: scale_by_complement_lower          ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T) in the lower atmosphere
    logical(wl),     dimension(:),     pointer :: minor_scales_with_density_upper    ! Density scaling is applied to minor absorption coefficients in the upper atmosphere
    logical(wl),     dimension(:),     pointer :: scale_by_complement_upper          ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T) in the upper atmosphere
    real(kind_phys), dimension(:,:,:,:), pointer :: kmajor                           ! Stored absorption coefficients due to major absorbing gases
    real(kind_phys), dimension(:,:,:), pointer :: kminor_lower                       ! Transformed from [nTemp x nEta x nGpt x nAbsorbers] array to [nTemp x nEta x nContributors] array
    real(kind_phys), dimension(:,:,:), pointer :: kminor_upper                       ! Transformed from [nTemp x nEta x nGpt x nAbsorbers] array to [nTemp x nEta x nContributors] array
    real(kind_phys), dimension(:,:,:), pointer :: vmr_ref                            ! Volume mixing ratios for reference atmosphere [mol mol-1]
    real(kind_phys), dimension(:,:),   pointer :: band_lims_wavenum                  ! Beginning and ending wavenumber for each band [cm-1]
    real(kind_phys), dimension(:),     pointer :: press_ref                          ! Pressures for reference atmosphere [Pa]
    real(kind_phys), dimension(:),     pointer :: temp_ref                           ! Temperatures for reference atmosphere [K]
    real(kind_phys), dimension(:),     pointer :: solar_src_quiet                    ! Quiet sun term of incoming solar irradiance [W m-2]
    real(kind_phys), dimension(:),     pointer :: solar_src_facular                  ! Facular brightening term of incoming solar irradiance [W m-2]
    real(kind_phys), dimension(:),     pointer :: solar_src_sunspot                  ! Sunspot dimming term of incoming solar irradiance [W m-2]
    real(kind_phys),                   pointer :: mg_default                         ! Mean value of Mg2 solar activity index [1]
    real(kind_phys),                   pointer :: sb_default                         ! Mean value of sunspot index [1]
    real(kind_phys),                   pointer :: tsi_default                        ! Default total solar irradiance [W m-2]
    real(kind_phys),                   pointer :: press_ref_trop                     ! Reference pressure separating the lower and upper atmosphere [Pa]
    real(kind_phys),                   pointer :: temp_ref_p                         ! Standard spectroscopic reference pressure [Pa]
    real(kind_phys),                   pointer :: temp_ref_t                         ! Standard spectroscopic reference temperature [K]
    real(kind_phys), dimension(:,:,:), pointer :: rayl_lower                         ! Stored coefficients due to rayleigh scattering contribution in lower part of atmosphere
    real(kind_phys), dimension(:,:,:), pointer :: rayl_upper                         ! Stored coefficients due to rayleigh scattering contribution in upper part of atmosphere
    real(kind_phys), dimension(:,:,:), allocatable :: rayl_lower_allocatable         ! Stored coefficients due to rayleigh scattering contribution in lower part of atmosphere
    real(kind_phys), dimension(:,:,:), allocatable :: rayl_upper_allocatable         ! Stored coefficients due to rayleigh scattering contribution in upper part of atmosphere
    integer,             dimension(:), pointer :: int2log                            ! use this to convert integer-to-logical.
    integer,                         parameter :: missing_variable_error_code = 3
    character(len=256)                         :: alloc_errmsg
    integer                                    :: idx

    ! Initialize error variables
    errmsg = ''
    errcode = 0

    file_reader = create_netcdf_reader_t()

    ! Open the shortwave coefficients file
    call file_reader%open_file(sw_filename, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if

    ! Read the coefficients from the file
    call file_reader%get_var('gas_names', gas_names, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('key_species', key_species, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('bnd_limits_gpt', band2gpt, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('bnd_limits_wavenumber', band_lims_wavenum, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('press_ref', press_ref, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('press_ref_trop', press_ref_trop, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('temp_ref', temp_ref, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('absorption_coefficient_ref_T', temp_ref_t, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('absorption_coefficient_ref_P', temp_ref_p, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('vmr_ref', vmr_ref, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('kmajor', kmajor, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('kminor_lower', kminor_lower, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('kminor_upper', kminor_upper, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('solar_source_quiet', solar_src_quiet, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('solar_source_facular', solar_src_facular, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('solar_source_sunspot', solar_src_sunspot, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('tsi_default', tsi_default, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('mg_default', mg_default, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('sb_default', sb_default, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('rayl_lower', rayl_lower, errmsg, errcode)
    ! OK if variable is not on file
    if (errcode /= 0 .and. errcode /= missing_variable_error_code) then
       return
    end if
    if (errcode /= missing_variable_error_code) then
       allocate(rayl_lower_allocatable(size(rayl_lower,1), size(rayl_lower,2), size(rayl_lower,3)))
       rayl_lower_allocatable = rayl_lower
    end if
    call file_reader%get_var('rayl_upper', rayl_upper, errmsg, errcode)
    ! OK if variable is not on file
    if (errcode /= 0 .and. errcode /= missing_variable_error_code) then
       return
    end if
    if (errcode /= missing_variable_error_code) then
       allocate(rayl_upper_allocatable(size(rayl_upper,1), size(rayl_upper,2), size(rayl_upper,3)))
       rayl_upper_allocatable = rayl_upper
    end if
    call file_reader%get_var('gas_minor', gas_minor, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('identifier_minor', identifier_minor, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('minor_gases_lower', minor_gases_lower, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('minor_gases_upper', minor_gases_upper, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('minor_limits_gpt_lower', minor_limits_gpt_lower, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('minor_limits_gpt_upper', minor_limits_gpt_upper, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('minor_scales_with_density_lower', int2log, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    allocate(minor_scales_with_density_lower(size(int2log)), stat=errcode, errmsg=alloc_errmsg)
    if (errcode /= 0) then
       write(errmsg, '(a,a)') 'Error allocating "minor_scales_with_density_lower" - message: ', alloc_errmsg
       return
    end if
    do idx = 1, size(int2log)
       if (int2log(idx) == 0) then
           minor_scales_with_density_lower(idx) = .false.
       else
           minor_scales_with_density_lower(idx) = .true.
       end if
    end do
    deallocate(int2log)
    nullify(int2log)
    call file_reader%get_var('scale_by_complement_lower', int2log, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    allocate(scale_by_complement_lower(size(int2log)), stat=errcode, errmsg=alloc_errmsg)
    if (errcode /= 0) then
       write(errmsg, '(a,a)') 'Error allocating "scale_by_complement_lower" - message: ', alloc_errmsg
       return
    end if
    do idx = 1, size(int2log)
       if (int2log(idx) == 0) then
           scale_by_complement_lower(idx) = .false.
       else
           scale_by_complement_lower(idx) = .true.
       end if
    end do
    deallocate(int2log)
    nullify(int2log)
    call file_reader%get_var('minor_scales_with_density_upper', int2log, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    allocate(minor_scales_with_density_upper(size(int2log)), stat=errcode, errmsg=alloc_errmsg)
    if (errcode /= 0) then
       write(errmsg, '(a,a)') 'Error allocating "minor_scales_with_density_upper" - message: ', alloc_errmsg
       return
    end if
    do idx = 1, size(int2log)
       if (int2log(idx) == 0) then
           minor_scales_with_density_upper(idx) = .false.
       else
           minor_scales_with_density_upper(idx) = .true.
       end if
    end do
    deallocate(int2log)
    nullify(int2log)
    call file_reader%get_var('scale_by_complement_upper', int2log, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    allocate(scale_by_complement_upper(size(int2log)), stat=errcode, errmsg=alloc_errmsg)
    if (errcode /= 0) then
       write(errmsg, '(a,a)') 'Error allocating "scale_by_complement_upper" - message: ', alloc_errmsg
       return
    end if
    do idx = 1, size(int2log)
       if (int2log(idx) == 0) then
           scale_by_complement_upper(idx) = .false.
       else
           scale_by_complement_upper(idx) = .true.
       end if
    end do
    deallocate(int2log)
    nullify(int2log)
    call file_reader%get_var('scaling_gas_lower', scaling_gas_lower, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('scaling_gas_upper', scaling_gas_upper, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('kminor_start_lower', kminor_start_lower, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if
    call file_reader%get_var('kminor_start_upper', kminor_start_upper, errmsg, errcode)
    if (errcode /= 0) then
       return
    end if

    ! Close the shortwave coefficients file
    call file_reader%close_file(errmsg, errcode)
    if (errcode /= 0) then
       return
    end if

    ! Initialize the gas optics object with data.
    errmsg = kdist%gas_props%load( &
         available_gases%gas_concs, gas_names, key_species,     &
         band2gpt, band_lims_wavenum,                           &
         press_ref, press_ref_trop, temp_ref,                   &
         temp_ref_p, temp_ref_t, vmr_ref,                       &
         kmajor, kminor_lower, kminor_upper,                    &
         gas_minor, identifier_minor,                           &
         minor_gases_lower, minor_gases_upper,                  &
         minor_limits_gpt_lower, minor_limits_gpt_upper,        &
         minor_scales_with_density_lower,                       &
         minor_scales_with_density_upper,                       &
         scaling_gas_lower, scaling_gas_upper,                  &
         scale_by_complement_lower, scale_by_complement_upper,  &
         kminor_start_lower, kminor_start_upper,                &
         solar_src_quiet, solar_src_facular, solar_src_sunspot, &
         tsi_default, mg_default, sb_default,                   &
         rayl_lower_allocatable, rayl_upper_allocatable)

    if (len_trim(errmsg) > 0) then
       errcode = 1
    end if
    call check_error_msg('rrtmgp_sw_gas_optics_init_load', errmsg)

    ! Deallocate pointer variables
    deallocate(gas_names, gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, &
            scaling_gas_lower, scaling_gas_upper, key_species, band2gpt, minor_limits_gpt_lower, &
            minor_limits_gpt_upper, kminor_start_lower, kminor_start_upper, minor_scales_with_density_lower, &
            minor_scales_with_density_upper, scale_by_complement_lower, scale_by_complement_upper, &
            kmajor, kminor_lower, kminor_upper, vmr_ref, band_lims_wavenum, solar_src_quiet, &
            solar_src_facular, solar_src_sunspot, mg_default, sb_default, tsi_default, press_ref, temp_ref, &
            press_ref_trop, temp_ref_p, temp_ref_t)
    nullify(gas_names)
    nullify(gas_minor)
    nullify(identifier_minor)
    nullify(minor_gases_lower)
    nullify(minor_gases_upper)
    nullify(scaling_gas_lower)
    nullify(scaling_gas_upper)
    nullify(key_species)
    nullify(band2gpt)
    nullify(minor_limits_gpt_lower)
    nullify(minor_limits_gpt_upper)
    nullify(kminor_start_lower)
    nullify(kminor_start_upper)
    nullify(minor_scales_with_density_lower)
    nullify(minor_scales_with_density_upper)
    nullify(scale_by_complement_lower)
    nullify(scale_by_complement_upper)
    nullify(kmajor)
    nullify(solar_src_quiet)
    nullify(solar_src_facular)
    nullify(solar_src_sunspot)
    nullify(mg_default)
    nullify(sb_default)
    nullify(tsi_default)
    nullify(kminor_lower)
    nullify(kminor_upper)
    nullify(vmr_ref)
    nullify(band_lims_wavenum)
    nullify(press_ref)
    nullify(temp_ref)
    nullify(press_ref_trop)
    nullify(temp_ref_p)
    nullify(temp_ref_t)
    if (associated(rayl_lower)) then
       deallocate(rayl_lower)
       nullify(rayl_lower)
    end if
    if (associated(rayl_upper)) then
       deallocate(rayl_upper)
       nullify(rayl_upper)
    end if

  end subroutine rrtmgp_sw_gas_optics_init

!> \section arg_table_rrtmgp_sw_gas_optics_run Argument Table
!! \htmlinclude rrtmgp_sw_gas_optics_run.html
!!
  subroutine rrtmgp_sw_gas_optics_run(dosw, iter_num, ncol, rrtmgp_phys_blksz, p_lay, p_lev, t_lay,  &
             gas_concs, sw_optical_props, sw_gas_props, toa_src_sw, errmsg, errflg)
   use machine,                 only: kind_phys
   use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
   use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
   use ccpp_optical_props,      only: ty_optical_props_2str_ccpp
   use radiation_tools,         only: check_error_msg
   ! Inputs
   logical,                           intent(in) :: dosw                        !< Flag for whether to perform longwave calculation
   integer,                           intent(in) :: iter_num                    !< Subcycle iteration number
   integer,                           intent(in) :: ncol                        !< Total number of columns
   integer,                           intent(in) :: rrtmgp_phys_blksz           !< Number of horizontal points to process at once
   real(kind_phys), dimension(:,:),   intent(in) :: p_lay                       !< Air pressure at midpoints [Pa]
   real(kind_phys), dimension(:,:),   intent(in) :: p_lev                       !< Air pressure at interfaces [Pa]
   real(kind_phys), dimension(:,:),   intent(in) :: t_lay                       !< Air temperature at midpoints [K]
   type(ty_gas_concs_ccpp),           intent(in) :: gas_concs                   !< RRTMGP gas concentrations object

   ! Outputs
   type(ty_optical_props_2str_ccpp),  intent(inout) :: sw_optical_props         !< Clearsky optical properties
   type(ty_gas_optics_rrtmgp_ccpp),   intent(inout) :: sw_gas_props             !< RRTMGP gas optics object
   real(kind_phys),                   intent(out)   :: toa_src_sw(:,:)
   character(len=512),                intent(out)   :: errmsg
   integer,                           intent(out)   :: errflg

   ! Local variables
   integer :: iCol, iCol2

   ! Set error variables
   errmsg = ''
   errflg = 0

   if (.not. dosw) then
      return
   end if

   iCol = ((iter_num - 1) * rrtmgp_phys_blksz) + 1
   iCol2= min(iCol + rrtmgp_phys_blksz - 1, ncol)

   errmsg = sw_gas_props%gas_props%gas_optics(&
         p_lay(iCol:iCol2,:),             & ! IN  - Pressure @ layer-centers (Pa)
         p_lev(iCol:iCol2,:),             & ! IN  - Pressure @ layer-interfaces (Pa)
         t_lay(iCol:iCol2,:),             & ! IN  - Temperature @ layer-centers (K)
         gas_concs%gas_concs,             & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
         sw_optical_props%optical_props,  & ! OUT - RRTMGP DDT: Shortwave optical properties, by
                                            !                   spectral point (tau,ssa,g)
         toa_src_sw)                        ! OUT - TOA incident shortwave radiation (spectral)

   if (len_trim(errmsg) /= 0) then
      errflg = 1
   end if

  end subroutine rrtmgp_sw_gas_optics_run

end module rrtmgp_sw_gas_optics
