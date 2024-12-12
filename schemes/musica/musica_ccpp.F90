!> Top-level wrapper for MUSICA chemistry components
module musica_ccpp
  use musica_ccpp_micm,     only: micm_register, micm_init, micm_run, micm_final
  use musica_ccpp_namelist, only: filename_of_tuvx_micm_mapping_configuration
  use musica_ccpp_tuvx,     only: tuvx_register, tuvx_init, tuvx_run, tuvx_final
  use musica_util,          only: index_mappings_t

  implicit none
  private

  public :: musica_ccpp_register, musica_ccpp_init, musica_ccpp_run, musica_ccpp_final

contains

  !> \section arg_table_musica_ccpp_register Argument Table
  !! \htmlinclude musica_ccpp_register.html
  subroutine musica_ccpp_register(constituent_props, errmsg, &
                                  errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_namelist,      only: micm_solver_type

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
    character(len=512),                               intent(out) :: errmsg
    integer,                                          intent(out) :: errcode

    type(ccpp_constituent_properties_t), allocatable :: constituent_props_subset(:)
    integer :: number_of_grid_cells

    ! Temporary fix until the number of grid cells is only needed to create a MICM state
    ! instead of when the solver is created.
    ! The number of grid cells is not know at this point, so we set it to 1 and recreate
    ! the solver when the number of grid cells is known at the init stage.
    number_of_grid_cells = 1
    call micm_register(micm_solver_type, number_of_grid_cells, constituent_props_subset, &
                       errmsg, errcode)
    if (errcode /= 0) return
    constituent_props = constituent_props_subset
    deallocate(constituent_props_subset)

    call tuvx_register(constituent_props_subset, errmsg, errcode)
    if (errcode /= 0) return
    constituent_props = [ constituent_props, constituent_props_subset ]

  end subroutine musica_ccpp_register

  !> \section arg_table_musica_ccpp_init Argument Table
  !! \htmlinclude musica_ccpp_init.html
  subroutine musica_ccpp_init(horizontal_dimension, vertical_layer_dimension, &
                              vertical_interface_dimension, &
                              photolysis_wavelength_grid_interfaces, &
                              constituent_props, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t, ccpp_constituent_prop_ptr_t
    use ccpp_kinds, only : kind_phys
    use musica_ccpp_micm, only: micm
    use musica_ccpp_namelist, only: micm_solver_type
    use musica_ccpp_util, only: has_error_occurred
    integer,                           intent(in)  :: horizontal_dimension                     ! (count)
    integer,                           intent(in)  :: vertical_layer_dimension                 ! (count)
    integer,                           intent(in)  :: vertical_interface_dimension             ! (count)
    real(kind_phys),                   intent(in)  :: photolysis_wavelength_grid_interfaces(:) ! m
    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props(:)
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errcode

    integer :: number_of_grid_cells
    type(ccpp_constituent_properties_t), allocatable :: micm_species_props(:)

    ! Temporary fix until the number of grid cells is only needed to create a MICM state
    ! instead of when the solver is created.
    ! Re-create the MICM solver with the correct number of grid cells
    number_of_grid_cells = horizontal_dimension * vertical_layer_dimension
    call micm_register(micm_solver_type, number_of_grid_cells, micm_species_props, errmsg, errcode)
    call micm_init(errmsg, errcode)
    if (errcode /= 0) return
    call tuvx_init(vertical_layer_dimension, vertical_interface_dimension, &
                   photolysis_wavelength_grid_interfaces, &
                   micm%user_defined_reaction_rates, &
                   constituent_props, errmsg, errcode)
    if (errcode /= 0) return

  end subroutine musica_ccpp_init

  !> \section arg_table_musica_ccpp_run Argument Table
  !! \htmlinclude musica_ccpp_run.html
  !!
  !! The standard name for the variable 'surface_temperature' is
  !! 'blackbody_temperature_at_surface' because this is what we have as
  !! the standard name for 'cam_in%ts', whcih represents the same quantity.
  subroutine musica_ccpp_run(time_step, temperature, pressure, dry_air_density, constituent_props,   &
                             constituents, geopotential_height_wrt_surface_at_midpoint,              &
                             geopotential_height_wrt_surface_at_interface, surface_geopotential,     &
                             surface_temperature, surface_albedo,                                    &
                             photolysis_wavelength_grid_interfaces, extraterrestrial_flux,           &
                             standard_gravitational_acceleration, cloud_area_fraction,               &
                             air_pressure_thickness, solar_zenith_angle,                             &
                             earth_sun_distance, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_kinds,                only: kind_phys
    use musica_ccpp_micm,          only: number_of_rate_parameters
    use musica_ccpp_micm_util,     only: convert_to_mol_per_cubic_meter, convert_to_mass_mixing_ratio

    real(kind_phys),         intent(in)    :: time_step                                         ! s
    real(kind_phys), target, intent(in)    :: temperature(:,:)                                  ! K
    real(kind_phys), target, intent(in)    :: pressure(:,:)                                     ! Pa
    real(kind_phys), target, intent(in)    :: dry_air_density(:,:)                              ! kg m-3
    type(ccpp_constituent_prop_ptr_t), &
                             intent(in)    :: constituent_props(:)
    real(kind_phys), target, intent(inout) :: constituents(:,:,:)                               ! kg kg-1
    real(kind_phys),         intent(in)    :: geopotential_height_wrt_surface_at_midpoint(:,:)  ! m (column, layer)
    real(kind_phys),         intent(in)    :: geopotential_height_wrt_surface_at_interface(:,:) ! m (column, interface)
    real(kind_phys),         intent(in)    :: surface_geopotential(:)                           ! m2 s-2
    real(kind_phys),         intent(in)    :: surface_temperature(:)                            ! K
    real(kind_phys),         intent(in)    :: surface_albedo                                    ! unitless
    real(kind_phys),         intent(in)    :: photolysis_wavelength_grid_interfaces(:)          ! nm
    real(kind_phys),         intent(in)    :: extraterrestrial_flux(:)                          ! photons cm-2 s-1 nm-1
    real(kind_phys),         intent(in)    :: standard_gravitational_acceleration               ! m s-2
    real(kind_phys),         intent(in)    :: cloud_area_fraction(:,:)                          ! unitless (column, level)
    real(kind_phys),         intent(in)    :: air_pressure_thickness(:,:)                       ! Pa (column, level)
    real(kind_phys),         intent(in)    :: solar_zenith_angle(:)                             ! radians (column)
    real(kind_phys),         intent(in)    :: earth_sun_distance                                ! AU
    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errcode

    ! local variables
    real(kind_phys), dimension(size(constituents, dim=3))    :: molar_mass_arr    ! kg mol-1
    real(kind_phys), dimension(size(constituents, dim=1), &
                               size(constituents, dim=2), &
                               number_of_rate_parameters)    :: rate_parameters ! various units
    integer :: i_elem

    ! Calculate photolysis rate constants using TUV-x
    call tuvx_run(temperature, dry_air_density,                  &
                  geopotential_height_wrt_surface_at_midpoint,   &
                  geopotential_height_wrt_surface_at_interface,  &
                  surface_geopotential, surface_temperature,     &
                  surface_albedo,                                &
                  photolysis_wavelength_grid_interfaces,         &
                  extraterrestrial_flux,                         &
                  standard_gravitational_acceleration,           &
                  cloud_area_fraction,                           &
                  solar_zenith_angle,                            &
                  earth_sun_distance,                            &
                  constituents,                                  &
                  air_pressure_thickness, rate_parameters,       &
                  errmsg, errcode)

    ! Get the molar mass that is set in the call to instantiate()
    do i_elem = 1, size(molar_mass_arr)
      call constituent_props(i_elem)%molar_mass(molar_mass_arr(i_elem), errcode, errmsg)
      if (errcode /= 0) then
        errmsg = "[MUSICA Error] Unable to get molar mass."
        return
      end if
    end do

    ! TODO(jiwon) Check molar mass is non zero as it becomes a denominator for unit converison
    ! this code will be deleted when the framework does the check
    do i_elem = 1, size(molar_mass_arr)
      if (molar_mass_arr(i_elem) <= 0) then
        errcode = 1
        errmsg = "[MUSICA Error] Molar mass must be greater than zero."
        return
      end if
    end do

    ! Convert CAM-SIMA unit to MICM unit (kg kg-1  ->  mol m-3)
    call convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)

    ! Solve chemistry at the current time step
    call micm_run(time_step, temperature, pressure, dry_air_density, rate_parameters, &
                  constituents, errmsg, errcode)

    ! Convert MICM unit back to CAM-SIMA unit (mol m-3  ->  kg kg-1)
    call convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, constituents)

  end subroutine musica_ccpp_run

  !> \section arg_table_musica_ccpp_final Argument Table
  !! \htmlinclude musica_ccpp_final.html
  subroutine musica_ccpp_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    call tuvx_final(errmsg, errcode)
    call micm_final(errmsg, errcode)

  end subroutine musica_ccpp_final

end module musica_ccpp
