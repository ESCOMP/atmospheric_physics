program test_tuvx_gas_species_profile

  use musica_ccpp_tuvx_gas_species_profile

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_create_gas_species_profile()

contains

  subroutine test_create_gas_species_profile()
    use musica_micm,                  only: Rosenbrock, RosenbrockStandardOrder
    use musica_tuvx_grid,             only: grid_t
    use ccpp_kinds,                   only: kind_phys
    use ccpp_constituent_prop_mod,    only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod,    only: ccpp_constituent_properties_t
    use musica_ccpp,                  only: musica_ccpp_register
    use musica_ccpp_tuvx_height_grid, only: create_height_grid
    use musica_ccpp_namelist,         only: filename_of_micm_configuration, &
                                            filename_of_tuvx_configuration, &
                                            filename_of_tuvx_micm_mapping_configuration

    integer,         parameter                       :: NUM_COLUMNS = 2
    integer,         parameter                       :: NUM_LAYERS = 2
    integer                                          :: NUM_GRID_CELLS = NUM_COLUMNS * NUM_LAYERS
    integer                                          :: solver_type = Rosenbrock
    ! Molar mass for O2, O3 is read from 'cam-sima-chemistry-data/mechanisms/chapman/micm/species.json'
    real(kind_phys), parameter                       :: MOLAR_MASS_O2 = 0.0319988_kind_phys
    real(kind_phys), parameter                       :: MOLAR_MASS_O3 = 0.0479982_kind_phys
    integer,         parameter                       :: NUMBER_GAS_SPECIES = 3
    integer,         parameter                       :: INDEX_NOT_KNOWN = -9999
    type(ccpp_constituent_prop_ptr_t),   allocatable :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, &
                                              target :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer     :: const_prop
    type(gas_species_t),                 allocatable :: gas_species_group(:)
    type(profile_group_t),               allocatable :: profile_gas_species_group(:)
    type(grid_t),                        pointer     :: height_grid => null()
    real(kind_phys)                                  :: dry_air_density(NUM_COLUMNS,NUM_LAYERS) ! kg m-3
    real(kind_phys)                                  :: gas_species_constituents(NUM_COLUMNS,NUM_LAYERS)
    real(kind_phys)                                  :: height_deltas(NUM_LAYERS) ! km
    integer                                          :: errcode
    character(len=512)                               :: errmsg
    character(len=50)                                :: label
    character(len=50)                                :: unit
    real(kind_phys)                                  :: molar_mass   ! kg mol-1
    real(kind_phys)                                  :: scale_height ! km
    integer                                          :: index
    logical                                          :: tmp_bool
    integer                                          :: i, i_gas_species, i_col

    filename_of_micm_configuration = 'musica_configurations/chapman/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/chapman/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/chapman/tuvx_micm_mapping.json'

    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
    gas_species_constituents(:,1) = (/ 0.35_kind_phys, 0.45_kind_phys /)
    gas_species_constituents(:,2) = (/ 0.55_kind_phys, 0.65_kind_phys /)
    height_deltas(:) = (/ 0.5_kind_phys, 1.5_kind_phys /)

    call musica_ccpp_register(solver_type, NUM_GRID_CELLS, constituent_props, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set( const_prop, errcode, errmsg )
    end do

    height_grid => create_height_grid( NUM_LAYERS, NUM_LAYERS + 1 , errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(associated(height_grid))

    call configure_gas_species( constituent_props_ptr, gas_species_group, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(size(gas_species_group) == NUMBER_GAS_SPECIES)

    allocate( profile_gas_species_group( size(gas_species_group) ) )

    call create_gas_species_profile_group( height_grid, gas_species_group, &
                                profile_gas_species_group, errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(size(profile_gas_species_group) == NUMBER_GAS_SPECIES)

    do i_gas_species = 1, size(gas_species_group)
      label = gas_species_group(i_gas_species)%label
      unit = gas_species_group(i_gas_species)%unit
      molar_mass = gas_species_group(i_gas_species)%molar_mass
      scale_height = gas_species_group(i_gas_species)%scale_height
      index = gas_species_group(i_gas_species)%index_constituent_props

      tmp_bool = (trim(label) == "air" .and. trim(unit) == "molecule cm-3" .and. molar_mass == MOLAR_MASS_DRY_AIR .and. &
                      scale_height == 8.01_kind_phys .and. index == INDEX_NOT_KNOWN) .or.  &
                 (trim(label) == "O2" .and. trim(unit) == "molecule cm-3" .and. molar_mass == MOLAR_MASS_O2 .and. &
                      scale_height == 7.0_kind_phys .and. index /= INDEX_NOT_KNOWN .and. index /= 0) .or.  &
                 (trim(label) == "O3" .and. trim(unit) == "molecule cm-3" .and. molar_mass == MOLAR_MASS_O3 .and. &
                      scale_height == 7.0_kind_phys .and. index /= INDEX_NOT_KNOWN .and. index /= 0)
      ASSERT(tmp_bool)

      ASSERT(associated(profile_gas_species_group(i_gas_species)%profile))
    end do

    do i_col = 1, NUM_COLUMNS
      do i_gas_species = 1, size(gas_species_group)
        call set_gas_species_values(profile_gas_species_group(i_gas_species)%profile, &
                                    gas_species_group(i_gas_species),                 &
                                    gas_species_constituents(i_col, :),               & ! mol m-3
                                    height_deltas, errmsg, errcode)
        if (errcode /= 0) return
      end do
    end do

    call deallocate_gas_species(gas_species_group)
    call deallocate_gas_species_profile_group(profile_gas_species_group)

  end subroutine test_create_gas_species_profile

end program test_tuvx_gas_species_profile