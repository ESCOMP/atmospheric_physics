module musica_ccpp_tuvx

  ! Note: "tuvx_t" is included in an external pre-built tuvx library that the host
  ! model is responsible for linking to during compilation
  use musica_tuvx,          only: tuvx_t, grid_t, profile_t
  use musica_ccpp_util,     only: has_error_occurred
  use ccpp_kinds,           only: kind_phys
  use musica_ccpp_namelist, only: filename_of_tuvx_configuration

  implicit none
  private

  public :: tuvx_init, tuvx_run, tuvx_final

  type(tuvx_t),    pointer :: tuvx => null()
  type(grid_t),    pointer :: height_grid => null()
  type(grid_t),    pointer :: wavelength_grid => null()
  type(profile_t), pointer :: temperature_profile => null()
  type(profile_t), pointer :: surface_albedo_profile => null()

contains

  !> Intitialize TUV-x
  subroutine tuvx_init(vertical_layer_dimension, vertical_interface_dimension, &
                       wavelength_grid_interfaces, errmsg, errcode)
    use musica_tuvx, only: grid_map_t, profile_map_t, radiator_map_t
    use musica_util, only: error_t
    use musica_ccpp_tuvx_util, only: tuvx_deallocate
    use musica_ccpp_tuvx_height_grid, &
      only: create_height_grid, height_grid_label, height_grid_unit
    use musica_ccpp_tuvx_wavelength_grid, &
      only: create_wavelength_grid, wavelength_grid_label, wavelength_grid_unit
    use musica_ccpp_tuvx_temperature, &
      only: create_temperature_profile, temperature_label, temperature_unit
    use musica_ccpp_tuvx_surface_albedo, &
      only: create_surface_albedo_profile, surface_albedo_label, surface_albedo_unit

    integer,            intent(in)  :: vertical_layer_dimension      ! (count)
    integer,            intent(in)  :: vertical_interface_dimension  ! (count)
    real(kind_phys),    intent(in)  :: wavelength_grid_interfaces(:) ! m
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! local variables
    type(grid_map_t),     pointer :: grids
    type(profile_map_t),  pointer :: profiles
    type(radiator_map_t), pointer :: radiators
    type(error_t)                 :: error

    grids => grid_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) return

    height_grid => create_height_grid( vertical_layer_dimension, &
        vertical_interface_dimension, errmsg, errcode )
    if (errcode /= 0) then
      deallocate( grids )
      return
    endif

    call grids%add( height_grid, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, null(), null(), null(), height_grid, null(), &
                            null(), null() )
      return
    end if

    wavelength_grid => create_wavelength_grid( wavelength_grid_interfaces, &
                                               errmsg, errcode )
    if (errcode /= 0) then
      call tuvx_deallocate( grids, null(), null(), null(), height_grid, null(), &
                            null(), null() )
      return
    endif

    call grids%add( wavelength_grid, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, null(), null(), null(), height_grid, &
                            wavelength_grid, null(), null() )
      return
    end if

    profiles => profile_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, null(), null(), null(), height_grid, &
                            wavelength_grid, null(), null() )
      return
    end if

    temperature_profile => create_temperature_profile( height_grid, errmsg, errcode )
    if (errcode /= 0) then
      call tuvx_deallocate( grids, profiles, null(), null(), height_grid, &
                            wavelength_grid, null(), null() )
      return
    endif

    call profiles%add( temperature_profile, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, profiles, null(), null(), height_grid, &
                            wavelength_grid, temperature_profile, null() )
      return
    end if

    surface_albedo_profile => create_surface_albedo_profile( wavelength_grid, &
                                                             errmsg, errcode )
    if (errcode /= 0) then
      call tuvx_deallocate( grids, profiles, null(), null(), height_grid, &
                            wavelength_grid, temperature_profile, null() )
      return
    endif

    call profiles%add( surface_albedo_profile, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, profiles, null(), null(), height_grid, &
            wavelength_grid, temperature_profile, surface_albedo_profile )
      return
    end if

    radiators => radiator_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, profiles, null(), null(), height_grid, &
            wavelength_grid, temperature_profile, surface_albedo_profile )
      return
    end if

    tuvx => tuvx_t( filename_of_tuvx_configuration, grids, profiles, &
                    radiators, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, profiles, radiators, null(), height_grid, &
            wavelength_grid, temperature_profile, surface_albedo_profile )
      return
    end if

    call tuvx_deallocate( grids, profiles, radiators, null(), height_grid, &
          wavelength_grid, temperature_profile, surface_albedo_profile )

    grids => tuvx%get_grids( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      tuvx => null()
      return
    end if

    height_grid => grids%get( height_grid_label, height_grid_unit, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, null(), null(), tuvx, null(), null(), &
                            null(), null() )
      return
    end if

    wavelength_grid => grids%get( wavelength_grid_label, wavelength_grid_unit, &
                                  error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, null(), null(), tuvx, height_grid, null(), &
                            null(), null() )
      return
    end if

    profiles => tuvx%get_profiles( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, null(), null(), tuvx, height_grid, &
            wavelength_grid, null(), null() )
      return
    end if

    temperature_profile => profiles%get( temperature_label, temperature_unit, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, profiles, null(), tuvx, height_grid, &
            wavelength_grid, null(), null() )
      return
    end if

    surface_albedo_profile => profiles%get( surface_albedo_label, surface_albedo_unit, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      call tuvx_deallocate( grids, profiles, null(), tuvx, height_grid, &
            wavelength_grid, temperature_profile, null() )
      return
    end if

    call tuvx_deallocate( grids, profiles, null(), null(), null(), null(), &
                          null(), null() )

  end subroutine tuvx_init

  !> Calculates photolysis rate constants for the current model conditions
  !!
  !! The variable name of 'surface_temperature' is set as 'blackbody_temperature_at_surface'
  !! according to CCPP standards, assuming that TUV-x needs the temperature everywhere.
  !! This is a topic that CGD needs to discuss with NOAA.
  subroutine tuvx_run(temperature, dry_air_density,                 &
                      geopotential_height_wrt_surface_at_midpoint,  &
                      geopotential_height_wrt_surface_at_interface, &
                      surface_temperature, surface_geopotential,    &
                      surface_albedo,                               &
                      standard_gravitational_acceleration,          &
                      photolysis_rate_constants, errmsg, errcode)
    use musica_util,                     only: error_t
    use musica_ccpp_tuvx_height_grid,    only: set_height_grid_values, calculate_heights
    use musica_ccpp_tuvx_temperature,    only: set_temperature_values
    use musica_ccpp_tuvx_surface_albedo, only: set_surface_albedo_values

    real(kind_phys),    intent(in)  :: temperature(:,:)                                  ! K (column, layer)
    real(kind_phys),    intent(in)  :: dry_air_density(:,:)                              ! kg m-3 (column, layer)
    real(kind_phys),    intent(in)  :: geopotential_height_wrt_surface_at_midpoint(:,:)  ! m (column, layer)
    real(kind_phys),    intent(in)  :: geopotential_height_wrt_surface_at_interface(:,:) ! m (column, interface)
    real(kind_phys),    intent(in)  :: surface_temperature(:)                            ! K
    real(kind_phys),    intent(in)  :: surface_geopotential(:)                           ! m2 s-2
    real(kind_phys),    intent(in)  :: surface_albedo                                    ! unitless
    real(kind_phys),    intent(in)  :: standard_gravitational_acceleration               ! m s-2
    ! temporarily set to Chapman mechanism and 1 dimension
    ! until mapping between MICM and TUV-x is implemented
    real(kind_phys),    intent(out) :: photolysis_rate_constants(:) ! s-1 (column, reaction)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! local variables
    real(kind_phys), dimension(size(geopotential_height_wrt_surface_at_midpoint, dim = 2))  :: height_midpoints
    real(kind_phys), dimension(size(geopotential_height_wrt_surface_at_interface, dim = 2)) :: height_interfaces
    real(kind_phys) :: reciprocal_of_gravitational_acceleration ! s2 m-1
    integer         :: i_col

    reciprocal_of_gravitational_acceleration = 1.0_kind_phys / standard_gravitational_acceleration

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
    end do

    call set_surface_albedo_values( surface_albedo_profile, surface_albedo, errmsg, errcode )
    if (errcode /= 0) return

    ! stand-in until actual photolysis rate constants are calculated
    photolysis_rate_constants(:) = 1.0e-6_kind_phys

  end subroutine tuvx_run

  !> Finalize tuvx
  subroutine tuvx_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

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

    if (associated( tuvx )) then
      deallocate( tuvx )
      tuvx => null()
    end if

  end subroutine tuvx_final

end module musica_ccpp_tuvx