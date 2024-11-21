module musica_ccpp_tuvx

  ! Note: "tuvx_t" is included in an external pre-built TUV-x library that the host
  ! model is responsible for linking to during compilation
  use ccpp_kinds,           only: kind_phys
  use musica_ccpp_namelist, only: filename_of_tuvx_configuration
  use musica_ccpp_util,     only: has_error_occurred
  use musica_tuvx,          only: tuvx_t, grid_t, profile_t
  use musica_util,          only: mappings_t, index_mappings_t

  implicit none
  private

  public :: tuvx_init, tuvx_run, tuvx_final

  type(tuvx_t),           pointer :: tuvx => null()
  type(grid_t),           pointer :: height_grid => null()
  type(grid_t),           pointer :: wavelength_grid => null()
  type(profile_t),        pointer :: temperature_profile => null()
  type(profile_t),        pointer :: surface_albedo_profile => null()
  type(profile_t),        pointer :: extraterrestrial_flux_profile => null()
  type(index_mappings_t), pointer :: photolysis_rate_constants_mapping => null( )
  integer                         :: number_of_photolysis_rate_constants = 0

contains

  subroutine reset_tuvx_map_state( grids, profiles, radiators )
    use musica_tuvx, only: grid_map_t, profile_map_t, radiator_map_t

    type(grid_map_t),     pointer :: grids
    type(profile_map_t),  pointer :: profiles
    type(radiator_map_t), pointer :: radiators

    if (associated( grids )) deallocate( grids )
    if (associated( profiles )) deallocate( profiles )
    if (associated( radiators )) deallocate( radiators )

  end subroutine reset_tuvx_map_state

  !> This is a helper subroutine created to deallocate objects associated with TUV-x
  subroutine cleanup_tuvx_resources()

    if (associated( height_grid )) then
      deallocate( height_grid )
      height_grid => null()
    end if

    if (associated( wavelength_grid )) then
      deallocate( wavelength_grid )
      wavelength_grid => null()
    end if

    if (associated( temperature_profile )) then
      deallocate( temperature_profile )
      temperature_profile => null()
    end if

    if (associated( surface_albedo_profile )) then
      deallocate( surface_albedo_profile )
      surface_albedo_profile => null()
    end if

    if (associated( extraterrestrial_flux_profile )) then
      deallocate( extraterrestrial_flux_profile )
      extraterrestrial_flux_profile => null()
    end if

    if (associated( photolysis_rate_constants_mapping )) then
      deallocate( photolysis_rate_constants_mapping )
      photolysis_rate_constants_mapping => null()
    end if

  end subroutine cleanup_tuvx_resources

  !> Initializes TUV-x
  subroutine tuvx_init(vertical_layer_dimension, vertical_interface_dimension, &
                       wavelength_grid_interfaces, micm_rate_parameter_ordering, &
                       errmsg, errcode)
    use musica_tuvx, only: grid_map_t, profile_map_t, radiator_map_t
    use musica_util, only: error_t, configuration_t
    use musica_ccpp_namelist, only: filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_height_grid, &
      only: create_height_grid, height_grid_label, height_grid_unit
    use musica_ccpp_tuvx_wavelength_grid, &
      only: create_wavelength_grid, wavelength_grid_label, wavelength_grid_unit
    use musica_ccpp_tuvx_temperature, &
      only: create_temperature_profile, temperature_label, temperature_unit
    use musica_ccpp_tuvx_surface_albedo, &
      only: create_surface_albedo_profile, surface_albedo_label, surface_albedo_unit
    use musica_ccpp_tuvx_extraterrestrial_flux, &
      only: create_extraterrestrial_flux_profile, extraterrestrial_flux_label, &
            extraterrestrial_flux_unit

    integer,            intent(in)  :: vertical_layer_dimension      ! (count)
    integer,            intent(in)  :: vertical_interface_dimension  ! (count)
    real(kind_phys),    intent(in)  :: wavelength_grid_interfaces(:) ! m
    type(mappings_t),   intent(in)  :: micm_rate_parameter_ordering  ! index mappings for MICM rate parameters
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! local variables
    type(grid_map_t),      pointer :: grids
    type(profile_map_t),   pointer :: profiles
    type(radiator_map_t),  pointer :: radiators
    type(configuration_t)          :: config
    type(mappings_t),      pointer :: photolysis_rate_constants_ordering
    type(error_t)                  :: error

    grids => grid_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) return

    height_grid => create_height_grid( vertical_layer_dimension, &
        vertical_interface_dimension, errmsg, errcode )
    if (errcode /= 0) then
      call reset_tuvx_map_state( grids, null(), null() )
      return
    endif

    call grids%add( height_grid, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call reset_tuvx_map_state( grids, null(), null() )
      call cleanup_tuvx_resources()
      return
    end if

    wavelength_grid => create_wavelength_grid( wavelength_grid_interfaces, &
                                               errmsg, errcode )
    if (errcode /= 0) then
      call reset_tuvx_map_state( grids, null(), null() )
      call cleanup_tuvx_resources()
      return
    endif

    call grids%add( wavelength_grid, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call reset_tuvx_map_state( grids, null(), null() )
      call cleanup_tuvx_resources()
      return
    end if

    profiles => profile_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call reset_tuvx_map_state( grids, null(), null() )
      call cleanup_tuvx_resources()
      return
    end if

    temperature_profile => create_temperature_profile( height_grid, errmsg, errcode )
    if (errcode /= 0) then
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    endif

    call profiles%add( temperature_profile, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    end if

    surface_albedo_profile => create_surface_albedo_profile( wavelength_grid, &
                                                             errmsg, errcode )
    if (errcode /= 0) then
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    endif

    call profiles%add( surface_albedo_profile, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    end if

    extraterrestrial_flux_profile => create_extraterrestrial_flux_profile( &
                wavelength_grid, wavelength_grid_interfaces, errmsg, errcode )
    if (errcode /= 0) then
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    endif

    call profiles%add( extraterrestrial_flux_profile, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    end if

    radiators => radiator_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    end if

    tuvx => tuvx_t( trim(filename_of_tuvx_configuration), grids, profiles, &
                    radiators, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call reset_tuvx_map_state( grids, profiles, radiators )
      call cleanup_tuvx_resources()
      return
    end if

    call reset_tuvx_map_state( grids, profiles, radiators )
    call cleanup_tuvx_resources()

    ! Gets resources associated with TUV-x from 'tuvx' pointer
    grids => tuvx%get_grids( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      return
    end if

    height_grid => grids%get( height_grid_label, height_grid_unit, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      call reset_tuvx_map_state( grids, null(), null() )
      return
    end if

    wavelength_grid => grids%get( wavelength_grid_label, wavelength_grid_unit, &
                                  error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      call reset_tuvx_map_state( grids, null(), null() )
      call cleanup_tuvx_resources()
      return
    end if

    profiles => tuvx%get_profiles( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      call reset_tuvx_map_state( grids, null(), null() )
      call cleanup_tuvx_resources()
      return
    end if

    temperature_profile => profiles%get( temperature_label, temperature_unit, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    end if

    surface_albedo_profile => profiles%get( surface_albedo_label, surface_albedo_unit, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    end if

    extraterrestrial_flux_profile => &
      profiles%get( extraterrestrial_flux_label, extraterrestrial_flux_unit, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      call reset_tuvx_map_state( grids, profiles, null() )
      call cleanup_tuvx_resources()
      return
    end if

    call reset_tuvx_map_state( grids, profiles, null() )

    ! 'photolysis_rate_constants_ordering' is a local variable
    photolysis_rate_constants_ordering => &
        tuvx%get_photolysis_rate_constants_ordering( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      return
    end if
    number_of_photolysis_rate_constants = photolysis_rate_constants_ordering%size()

    call config%load_from_file( trim(filename_of_tuvx_micm_mapping_configuration), error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      deallocate( photolysis_rate_constants_ordering )
      return
    end if

    photolysis_rate_constants_mapping => &
        index_mappings_t( config, photolysis_rate_constants_ordering, &
                          micm_rate_parameter_ordering, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      deallocate( photolysis_rate_constants_ordering )
      return
    end if

    deallocate( photolysis_rate_constants_ordering )

  end subroutine tuvx_init

  !> Calculates photolysis rate constants for the current model conditions
  subroutine tuvx_run(temperature, dry_air_density,                  &
                      geopotential_height_wrt_surface_at_midpoint,   &
                      geopotential_height_wrt_surface_at_interface,  &
                      surface_geopotential, surface_temperature,     &
                      surface_albedo,                                &
                      number_of_photolysis_wavelength_grid_sections, &
                      photolysis_wavelength_grid_interfaces,         &
                      extraterrestrial_flux,                         &
                      standard_gravitational_acceleration,           &
                      rate_parameters, errmsg, errcode)
    use musica_util,                            only: error_t
    use musica_ccpp_tuvx_height_grid,           only: set_height_grid_values, calculate_heights
    use musica_ccpp_tuvx_temperature,           only: set_temperature_values
    use musica_ccpp_util,                       only: has_error_occurred
    use musica_ccpp_tuvx_surface_albedo,        only: set_surface_albedo_values
    use musica_ccpp_tuvx_extraterrestrial_flux, only: set_extraterrestrial_flux_values

    real(kind_phys),    intent(in)    :: temperature(:,:)                                  ! K (column, layer)
    real(kind_phys),    intent(in)    :: dry_air_density(:,:)                              ! kg m-3 (column, layer)
    real(kind_phys),    intent(in)    :: geopotential_height_wrt_surface_at_midpoint(:,:)  ! m (column, layer)
    real(kind_phys),    intent(in)    :: geopotential_height_wrt_surface_at_interface(:,:) ! m (column, interface)
    real(kind_phys),    intent(in)    :: surface_geopotential(:)                           ! m2 s-2
    real(kind_phys),    intent(in)    :: surface_temperature(:)                            ! K
    real(kind_phys),    intent(in)    :: surface_albedo                                    ! unitless
    integer,            intent(in)    :: number_of_photolysis_wavelength_grid_sections     ! (count)
    real(kind_phys),    intent(in)    :: photolysis_wavelength_grid_interfaces(:)          ! nm
    real(kind_phys),    intent(in)    :: extraterrestrial_flux(:)                          ! photons cm-2 s-1 nm-1
    real(kind_phys),    intent(in)    :: standard_gravitational_acceleration               ! m s-2
    real(kind_phys),    intent(inout) :: rate_parameters(:,:,:)                            ! various units (column, layer, reaction)
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errcode

    ! local variables
    real(kind_phys), dimension(size(geopotential_height_wrt_surface_at_midpoint, dim = 2))  :: height_midpoints
    real(kind_phys), dimension(size(geopotential_height_wrt_surface_at_interface, dim = 2)) :: height_interfaces
    real(kind_phys), dimension(size(rate_parameters, dim=2)+2, &
                               number_of_photolysis_rate_constants) :: photolysis_rate_constants, & ! s-1
                                                                       heating_rates                ! K s-1 (TODO: check units)
    real(kind_phys) :: reciprocal_of_gravitational_acceleration ! s2 m-1
    real(kind_phys) :: solar_zenith_angle ! degrees
    real(kind_phys) :: earth_sun_distance ! AU
    type(error_t)   :: error
    integer         :: i_col, i_level

    reciprocal_of_gravitational_acceleration = 1.0_kind_phys / standard_gravitational_acceleration

    ! surface albedo with respect to direct UV/visible radiation
    call set_surface_albedo_values( surface_albedo_profile, surface_albedo, errmsg, errcode )
    if (errcode /= 0) return

    call set_extraterrestrial_flux_values( extraterrestrial_flux_profile,                 &
                                           number_of_photolysis_wavelength_grid_sections, &
                                           photolysis_wavelength_grid_interfaces,         &
                                           extraterrestrial_flux, errmsg, errcode )
    if (errcode /= 0) return

    do i_col = 1, size(temperature, dim=1)
      call calculate_heights( geopotential_height_wrt_surface_at_midpoint(i_col,:),  &
                              geopotential_height_wrt_surface_at_interface(i_col,:), &
                              surface_geopotential(i_col),                           &
                              reciprocal_of_gravitational_acceleration,              &
                              height_midpoints, height_interfaces )
      call set_height_grid_values( height_grid, height_midpoints, height_interfaces, &
                                   errmsg, errcode )
      if (errcode /= 0) return

      call set_temperature_values( temperature_profile, temperature(i_col,:), &
                                   surface_temperature(i_col), errmsg, errcode )
      if (errcode /= 0) return

      ! temporary values until these are available from the host model
      solar_zenith_angle = 0.0_kind_phys
      earth_sun_distance = 1.0_kind_phys

      ! calculate photolysis rate constants and heating rates
      call tuvx%run( solar_zenith_angle, earth_sun_distance, &
                     photolysis_rate_constants(:,:), heating_rates(:,:), &
                     error )
      if (has_error_occurred( error, errmsg, errcode )) return

      ! filter out negative photolysis rate constants
      photolysis_rate_constants(:,:) = &
          max( photolysis_rate_constants(:,:), 0.0_kind_phys )

      ! map photolysis rate constants to the host model's rate parameters and vertical grid
      do i_level = 1, size(rate_parameters, dim=2)
        call photolysis_rate_constants_mapping%copy_data( &
            photolysis_rate_constants(size(rate_parameters, dim=2)-i_level+2,:), &
            rate_parameters(i_col,i_level,:) )
      end do
    end do

  end subroutine tuvx_run

  !> Finalizes TUV-x
  subroutine tuvx_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

    if (associated( tuvx )) then
      deallocate( tuvx )
      tuvx => null()
    end if

    call cleanup_tuvx_resources()

  end subroutine tuvx_final

end module musica_ccpp_tuvx