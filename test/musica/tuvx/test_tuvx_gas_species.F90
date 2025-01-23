program test_tuvx_gas_species_profile

  use musica_ccpp_tuvx_gas_species

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_create_gas_species_profile()
  call test_initialize_tuvx_species()

contains

  !> Returns the wavelength edges used in the CAM-Chem photolysis rate constant lookup table
  !! These are the values that will be used in CAM-SIMA and correspond to the wavelength
  !! bins used in the CAM-Chem photolysis rate constant lookup table.
  !!
  !! We're using the actual values here because several of the TS1/TSMLT photolysis
  !! rate constant configurations are sensitive to the wavelength grid.
  subroutine get_wavelength_edges(edges)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), dimension(:), intent(out) :: edges

    edges = (/ &
      120.0e-9_kind_phys, &
      121.4e-9_kind_phys, &
      121.9e-9_kind_phys, &
      123.5e-9_kind_phys, &
      124.3e-9_kind_phys, &
      125.5e-9_kind_phys, &
      126.3e-9_kind_phys, &
      127.1e-9_kind_phys, &
      130.1e-9_kind_phys, &
      131.1e-9_kind_phys, &
      135.0e-9_kind_phys, &
      140.0e-9_kind_phys, &
      145.0e-9_kind_phys, &
      150.0e-9_kind_phys, &
      155.0e-9_kind_phys, &
      160.0e-9_kind_phys, &
      165.0e-9_kind_phys, &
      168.0e-9_kind_phys, &
      171.0e-9_kind_phys, &
      173.0e-9_kind_phys, &
      174.4e-9_kind_phys, &
      175.4e-9_kind_phys, &
      177.0e-9_kind_phys, &
      178.6e-9_kind_phys, &
      180.2e-9_kind_phys, &
      181.8e-9_kind_phys, &
      183.5e-9_kind_phys, &
      185.2e-9_kind_phys, &
      186.9e-9_kind_phys, &
      188.7e-9_kind_phys, &
      190.5e-9_kind_phys, &
      192.3e-9_kind_phys, &
      194.2e-9_kind_phys, &
      196.1e-9_kind_phys, &
      198.0e-9_kind_phys, &
      200.0e-9_kind_phys, &
      202.0e-9_kind_phys, &
      204.1e-9_kind_phys, &
      206.2e-9_kind_phys, &
      208.0e-9_kind_phys, &
      211.0e-9_kind_phys, &
      214.0e-9_kind_phys, &
      217.0e-9_kind_phys, &
      220.0e-9_kind_phys, &
      223.0e-9_kind_phys, &
      226.0e-9_kind_phys, &
      229.0e-9_kind_phys, &
      232.0e-9_kind_phys, &
      235.0e-9_kind_phys, &
      238.0e-9_kind_phys, &
      241.0e-9_kind_phys, &
      244.0e-9_kind_phys, &
      247.0e-9_kind_phys, &
      250.0e-9_kind_phys, &
      253.0e-9_kind_phys, &
      256.0e-9_kind_phys, &
      259.0e-9_kind_phys, &
      263.0e-9_kind_phys, &
      267.0e-9_kind_phys, &
      271.0e-9_kind_phys, &
      275.0e-9_kind_phys, &
      279.0e-9_kind_phys, &
      283.0e-9_kind_phys, &
      287.0e-9_kind_phys, &
      291.0e-9_kind_phys, &
      295.0e-9_kind_phys, &
      298.5e-9_kind_phys, &
      302.5e-9_kind_phys, &
      305.5e-9_kind_phys, &
      308.5e-9_kind_phys, &
      311.5e-9_kind_phys, &
      314.5e-9_kind_phys, &
      317.5e-9_kind_phys, &
      322.5e-9_kind_phys, &
      327.5e-9_kind_phys, &
      332.5e-9_kind_phys, &
      337.5e-9_kind_phys, &
      342.5e-9_kind_phys, &
      347.5e-9_kind_phys, &
      350.0e-9_kind_phys, &
      355.0e-9_kind_phys, &
      360.0e-9_kind_phys, &
      365.0e-9_kind_phys, &
      370.0e-9_kind_phys, &
      375.0e-9_kind_phys, &
      380.0e-9_kind_phys, &
      385.0e-9_kind_phys, &
      390.0e-9_kind_phys, &
      395.0e-9_kind_phys, &
      400.0e-9_kind_phys, &
      405.0e-9_kind_phys, &
      410.0e-9_kind_phys, &
      415.0e-9_kind_phys, &
      420.0e-9_kind_phys, &
      430.0e-9_kind_phys, &
      440.0e-9_kind_phys, &
      450.0e-9_kind_phys, &
      500.0e-9_kind_phys, &
      550.0e-9_kind_phys, &
      600.0e-9_kind_phys, &
      650.0e-9_kind_phys, &
      700.0e-9_kind_phys, &
      750.0e-9_kind_phys &
    /)

  end subroutine get_wavelength_edges

  ! Calculates the expected values for comparison with the answer
  subroutine calculate_gas_species_interfaces_and_densities( &
      molar_mass, dry_air_density, constituents, height_deltas, &
      is_O3, interfaces, densities)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in)    :: molar_mass
    real(kind_phys), intent(in)    :: dry_air_density(:)
    real(kind_phys), intent(in)    :: constituents(:)
    real(kind_phys), intent(in)    :: height_deltas(:)
    logical,         intent(in)    :: is_O3
    real(kind_phys), intent(inout) :: interfaces(size(constituents) + 2)
    real(kind_phys), intent(inout) :: densities(size(constituents) + 1)

    ! local variables
    real(kind_phys) :: constituent_mol_per_cm_3(size(constituents)) ! mol cm-3
    integer         :: num_vertical_levels

    constituent_mol_per_cm_3(:) = constituents(:) * dry_air_density(:) / molar_mass / m_3_to_cm_3

    num_vertical_levels = size(constituents)
    interfaces(1) = constituent_mol_per_cm_3(num_vertical_levels)
    interfaces(2:num_vertical_levels+1) = constituent_mol_per_cm_3(num_vertical_levels:1:-1)
    interfaces(num_vertical_levels+2) = constituent_mol_per_cm_3(1)

    if ( is_O3 ) then
      densities(:) = height_deltas(:) * km_to_cm           &
                  * ( interfaces(1:num_vertical_levels+1)  &
                  + interfaces(2:num_vertical_levels+2) ) * 0.5_kind_phys
    else
      densities(:) = height_deltas(:) * km_to_cm             &
                 * sqrt(interfaces(1:num_vertical_levels+1)) &
                 * sqrt(interfaces(2:num_vertical_levels+2))
    end if

  end subroutine calculate_gas_species_interfaces_and_densities

  subroutine test_create_gas_species_profile()
    use ccpp_kinds,                    only: kind_phys
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_properties_t
    use musica_tuvx,                   only: grid_t, profile_t
    use musica_util,                   only: error_t
    use musica_ccpp,                   only: musica_ccpp_register, musica_ccpp_final
    use musica_ccpp_tuvx_height_grid,  only: create_height_grid
    use musica_ccpp_namelist,          only: filename_of_micm_configuration, &
                                             filename_of_tuvx_configuration, &
                                             filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_load_species, only: index_dry_air, index_O2, index_O3, &
                                             SCALE_HEIGHT_DRY_AIR, SCALE_HEIGHT_O2, SCALE_HEIGHT_O3, &
                                             MOLAR_MASS_DRY_AIR, MOLAR_MASS_O2, MOLAR_MASS_O3
    use musica_ccpp_species,           only: tuvx_species_set, MUSICA_INT_UNASSIGNED

    integer,         parameter          :: NUM_COLUMNS = 2
    integer,         parameter          :: NUM_LAYERS = 2
    integer,         parameter          :: NUM_TUVX_SPECIES = 4
    real,            parameter          :: ABS_ERROR = 10.0
    real,            parameter          :: m_to_cm = 100.0_kind_phys
    type(ccpp_constituent_properties_t), &
                    allocatable, target :: constituent_props(:)
    type(grid_t),    pointer            :: height_grid => null()
    type(profile_t), pointer            :: dry_air_profile => null()
    type(profile_t), pointer            :: O2_profile => null()
    type(profile_t), pointer            :: O3_profile => null()
    real(kind_phys)                     :: dry_air_density(NUM_COLUMNS,NUM_LAYERS) ! kg m-3
    real(kind_phys)                     :: constituents(NUM_COLUMNS,NUM_LAYERS,NUM_TUVX_SPECIES)
    real(kind_phys)                     :: height_deltas(NUM_LAYERS+1) ! km
    integer                             :: errcode
    character(len=512)                  :: errmsg
    type(error_t)                       :: error
    character(len=50)                   :: name, unit
    real(kind_phys)                     :: molar_mass   ! kg mol-1
    real(kind_phys)                     :: scale_height ! km
    integer                             :: index_musica
    logical                             :: tmp_bool
    integer                             :: i_elem, i_col, i, j, k
    real(kind_phys)                     :: dry_air_interfaces(NUM_LAYERS+2), expected_dry_air_interfaces(NUM_COLUMNS,NUM_LAYERS+2)
    real(kind_phys)                     :: O2_interfaces(NUM_LAYERS+2), expected_O2_interfaces(NUM_COLUMNS,NUM_LAYERS+2)
    real(kind_phys)                     :: O3_interfaces(NUM_LAYERS+2), expected_O3_interfaces(NUM_COLUMNS,NUM_LAYERS+2)
    real(kind_phys)                     :: dry_air_densities(NUM_LAYERS+1), expected_dry_air_densities(NUM_COLUMNS,NUM_LAYERS+1)
    real(kind_phys)                     :: O2_densities(NUM_LAYERS+1), expected_O2_densities(NUM_COLUMNS,NUM_LAYERS+1)
    real(kind_phys)                     :: O3_densities(NUM_LAYERS+1), expected_O3_densities(NUM_COLUMNS,NUM_LAYERS+1)
    real(kind_phys)                     :: dry_air_exo_layer_density
    real(kind_phys)                     :: expected_dry_air_exo_layer_density = SCALE_HEIGHT_DRY_AIR
    real(kind_phys)                     :: O2_exo_layer_density
    real(kind_phys)                     :: expected_O2_exo_layer_density = SCALE_HEIGHT_O2
    real(kind_phys)                     :: O3_exo_layer_density
    real(kind_phys)                     :: expected_O3_exo_layer_density = SCALE_HEIGHT_O3

    filename_of_micm_configuration = 'musica_configurations/chapman/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/chapman/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/chapman/tuvx_micm_mapping.json'

    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
    height_deltas(:) = (/ 0.5_kind_phys, 1.5_kind_phys, 1.0_kind_phys /)

    ! Initialize the constituents array with real values between 0 and 1
    do k = 1, NUM_TUVX_SPECIES
      do j = 1, NUM_LAYERS
        do i = 1, NUM_COLUMNS
          constituents(i, j, k) = (i + j + k) * 0.1_kind_phys
        end do
      end do
    end do

    call musica_ccpp_register(constituent_props, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    do i_col = 1, NUM_COLUMNS
      call calculate_gas_species_interfaces_and_densities( &
        MOLAR_MASS_DRY_AIR, dry_air_density(i_col,:), constituents(i_col,:,index_dry_air), &
        height_deltas, .false., expected_dry_air_interfaces(i_col,:), expected_dry_air_densities(i_col,:) )

      call calculate_gas_species_interfaces_and_densities( &
        MOLAR_MASS_O2, dry_air_density(i_col,:), constituents(i_col,:,index_O2), &
        height_deltas, .false., expected_O2_interfaces(i_col,:), expected_O2_densities(i_col,:) )

      call calculate_gas_species_interfaces_and_densities( &
        MOLAR_MASS_O3, dry_air_density(i_col,:), constituents(i_col,:,index_O3), &
        height_deltas, .true., expected_O3_interfaces(i_col,:), expected_O3_densities(i_col,:) )
    end do

    height_grid => create_height_grid( NUM_LAYERS, NUM_LAYERS + 1 , errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(height_grid))

    dry_air_profile => create_dry_air_profile( height_grid, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(dry_air_profile))

    O2_profile => create_O2_profile( height_grid, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(O2_profile))

    O3_profile => create_O3_profile( height_grid, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(O3_profile))

    do i_col = 1, NUM_COLUMNS
      call set_gas_species_values( dry_air_profile, dry_air_density(i_col,:), &
            constituents(i_col,:,index_dry_air), height_deltas, index_dry_air, &
            errmsg, errcode )
      ASSERT(errcode == 0)

      call set_gas_species_values( O2_profile, dry_air_density(i_col,:), &
            constituents(i_col,:,index_O2), height_deltas, index_O2, &
            errmsg, errcode )
      ASSERT(errcode == 0)

      call set_gas_species_values( O3_profile, dry_air_density(i_col,:), &
            constituents(i_col,:,index_O3), height_deltas, index_O3, &
            errmsg, errcode )
      ASSERT(errcode == 0)

      ! Retrieve the values set by the profile
      ! Dry air
      call dry_air_profile%get_edge_values( dry_air_interfaces, error )
      ASSERT(error%is_success())
      do i = 1, size(dry_air_interfaces)
        ASSERT(dry_air_interfaces(i) == expected_dry_air_interfaces(i_col, i))
      end do
      call dry_air_profile%get_layer_densities( dry_air_densities, error )
      ASSERT(error%is_success())
      do i = 1, size(dry_air_densities) - 1
        ASSERT(dry_air_densities(i) == expected_dry_air_densities(i_col, i))
      end do
      ! the calculate_exo_layer_density call uses the scale height and the density of the uppermost
      ! layer to estimate the density above the column top, which affects the uppermost top value
      ASSERT_NEAR(dry_air_densities(size(dry_air_densities)), expected_dry_air_densities(i_col, size(dry_air_densities)) * SCALE_HEIGHT_DRY_AIR * m_to_cm, ABS_ERROR)
      dry_air_exo_layer_density = dry_air_profile%get_exo_layer_density( error )
      ASSERT(error%is_success())
      expected_dry_air_exo_layer_density = expected_dry_air_densities(i_col, size(dry_air_densities)) * SCALE_HEIGHT_DRY_AIR * m_to_cm
      ASSERT_NEAR(dry_air_exo_layer_density, expected_dry_air_exo_layer_density, ABS_ERROR)

      ! O2
      call O2_profile%get_edge_values( O2_interfaces, error )
      ASSERT(error%is_success())
      do i = 1, size(O2_interfaces)
        ASSERT(O2_interfaces(i) == expected_O2_interfaces(i_col, i))
      end do
      call O2_profile%get_layer_densities( O2_densities, error )
      ASSERT(error%is_success())
      do i = 1, size(O2_densities) - 1
        ASSERT(O2_densities(i) == expected_O2_densities(i_col, i))
      end do
      ASSERT_NEAR(O2_densities(size(O2_densities)), expected_O2_densities(i_col, size(O2_densities)) * SCALE_HEIGHT_O2 * m_to_cm, ABS_ERROR)
      O2_exo_layer_density = O2_profile%get_exo_layer_density( error )
      ASSERT(error%is_success())
      expected_O2_exo_layer_density = expected_O2_densities(i_col, size(O2_densities)) * SCALE_HEIGHT_O2 * m_to_cm
      ASSERT_NEAR(O2_exo_layer_density, expected_O2_exo_layer_density, ABS_ERROR)

      ! O3
      call O3_profile%get_edge_values( O3_interfaces, error )
      ASSERT(error%is_success())
      do i = 1, size(O3_interfaces)
        ASSERT(O3_interfaces(i) == expected_O3_interfaces(i_col, i))
      end do
      call O3_profile%get_layer_densities( O3_densities, error )
      ASSERT(error%is_success())
      do i = 1, size(O3_densities) - 1
        ASSERT(O3_densities(i) == expected_O3_densities(i_col, i))
      end do
      ASSERT_NEAR(O3_densities(size(O3_densities)), expected_O3_densities(i_col, size(O3_densities)) * SCALE_HEIGHT_O3 * m_to_cm, ABS_ERROR)
      O3_exo_layer_density = O3_profile%get_exo_layer_density( error )
      ASSERT(error%is_success())
      expected_O3_exo_layer_density = expected_O3_densities(i_col, size(O3_densities)) * SCALE_HEIGHT_O3 * m_to_cm
      ASSERT_NEAR(O3_exo_layer_density, expected_O3_exo_layer_density, ABS_ERROR)
    end do

    ! The gas species index starts at 2 because index 1 is reserved for cloud liquid
    do i_elem = 2, size(tuvx_species_set)
      name = tuvx_species_set(i_elem)%name
      unit = tuvx_species_set(i_elem)%unit
      molar_mass = tuvx_species_set(i_elem)%molar_mass
      scale_height = tuvx_species_set(i_elem)%scale_height
      index_musica = tuvx_species_set(i_elem)%index_musica_species
      tmp_bool = (trim(name) == 'air' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0289644_kind_phys .and. &
                  scale_height == 8.01_kind_phys .and. index_musica == 2) .or.  &
                 (trim(name) == 'O2' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0319988_kind_phys .and. &
                  scale_height == 7.0_kind_phys .and. index_musica == 3) .or.  &
                 (trim(name) == 'O3' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0479982_kind_phys .and. &
                  scale_height == 7.0_kind_phys .and. index_musica == 4)
      ASSERT(tmp_bool)
    end do

    deallocate( dry_air_profile )
    deallocate( O2_profile )
    deallocate( O3_profile )
    deallocate( height_grid )

    call musica_ccpp_final(errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

  end subroutine test_create_gas_species_profile

  subroutine test_initialize_tuvx_species()
    use ccpp_kinds,                    only: kind_phys
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_properties_t
    use musica_tuvx,                   only: grid_t, profile_t
    use musica_ccpp,                   only: musica_ccpp_register, musica_ccpp_init, musica_ccpp_final
    use musica_ccpp_tuvx_height_grid,  only: create_height_grid
    use musica_ccpp_namelist,          only: filename_of_micm_configuration, &
                                             filename_of_tuvx_configuration, &
                                             filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_load_species, only: index_dry_air, index_O2, index_O3
    use musica_ccpp_species,           only: tuvx_species_set, MUSICA_INT_UNASSIGNED

    integer,         parameter                       :: NUM_COLUMNS = 2
    integer,         parameter                       :: NUM_LAYERS = 2
    integer,         parameter                       :: NUM_WAVELENGTH_BINS = 102
    integer,         parameter                       :: NUM_TUVX_SPECIES = 4
    real(kind_phys)                                  :: photolysis_wavelength_grid_interfaces(NUM_WAVELENGTH_BINS+1) ! m
    type(ccpp_constituent_properties_t), allocatable, &
                                              target :: constituent_props(:)
    type(ccpp_constituent_prop_ptr_t),   allocatable :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), pointer     :: const_prop
    real(kind_phys)                                  :: dry_air_density(NUM_COLUMNS,NUM_LAYERS) ! kg m-3
    real(kind_phys)                                  :: constituents(NUM_COLUMNS,NUM_LAYERS, NUM_TUVX_SPECIES)
    integer                                          :: errcode
    character(len=512)                               :: errmsg
    character(len=50)                                :: name, unit
    real(kind_phys)                                  :: molar_mass   ! kg mol-1
    real(kind_phys)                                  :: scale_height ! km
    integer                                          :: index_musica, index_constituent_props
    logical                                          :: tmp_bool, has_profile
    integer                                          :: i_elem, i_col, i, j, k

    filename_of_micm_configuration = 'musica_configurations/chapman/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/chapman/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/chapman/tuvx_micm_mapping.json'

    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
    call get_wavelength_edges( photolysis_wavelength_grid_interfaces )
    do k = 1, 4
      do j = 1, 2
        do i = 1, 2
          constituents(i, j, k) = (i + j + k) * 0.1_kind_phys
        end do
      end do
    end do

    call musica_ccpp_register( constituent_props, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set( const_prop, errcode, errmsg )
    end do

    call musica_ccpp_init( NUM_COLUMNS, NUM_LAYERS, NUM_LAYERS+1, photolysis_wavelength_grid_interfaces, &
                          constituent_props_ptr, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    ! The gas species index starts at 2 because index 1 is reserved for cloud liquid
    do i_elem = 2, size(tuvx_species_set)
      name = tuvx_species_set(i_elem)%name
      unit = tuvx_species_set(i_elem)%unit
      molar_mass = tuvx_species_set(i_elem)%molar_mass
      scale_height = tuvx_species_set(i_elem)%scale_height
      index_musica = tuvx_species_set(i_elem)%index_musica_species
      index_constituent_props = tuvx_species_set(i_elem)%index_constituent_props
      has_profile = tuvx_species_set(i_elem)%profiled
      tmp_bool = (trim(name) == 'air' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0289644_kind_phys .and. &
                 scale_height == 8.01_kind_phys .and. index_musica == 2 .and. index_constituent_props /= MUSICA_INT_UNASSIGNED &
                 .and. has_profile) .or.  &
                (trim(name) == 'O2' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0319988_kind_phys  .and. &
                 scale_height == 7.0_kind_phys .and. index_musica == 3 .and. index_constituent_props /= MUSICA_INT_UNASSIGNED &
                 .and. has_profile) .or.  &
                (trim(name) == 'O3' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0479982_kind_phys .and. &
                 scale_height == 7.0_kind_phys .and. index_musica == 4 .and. index_constituent_props /= MUSICA_INT_UNASSIGNED &
                 .and. has_profile)
      ASSERT(tmp_bool)
    end do

    call musica_ccpp_final(errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

  end subroutine  test_initialize_tuvx_species

end program test_tuvx_gas_species_profile