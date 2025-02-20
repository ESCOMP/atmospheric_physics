program test_tuvx_no_photolysis

  use musica_ccpp_tuvx_no_photolysis

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_calculate_NO_photolysis_rate()

contains

  subroutine test_calculate_NO_photolysis_rate()
    use ccpp_kinds,                          only: kind_phys
    use ccpp_constituent_prop_mod,           only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod,           only: ccpp_constituent_properties_t
    use musica_ccpp,                         only: musica_ccpp_register, musica_ccpp_init, musica_ccpp_final
    use musica_ccpp_namelist,                only: filename_of_micm_configuration, &
                                                   filename_of_tuvx_configuration, &
                                                   filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_load_species,       only: index_cloud_liquid_water_content, index_dry_air, &
                                                   tuvx_index_O2 => index_O2, &
                                                   tuvx_index_O3 => index_O3, &
                                                   tuvx_index_NO => index_NO, &
                                                   tuvx_index_N2 => index_N2
    use musica_ccpp_species,                 only: tuvx_indices_constituent_props, extract_subset_constituents, &
                                                   initialize_musica_species_indices, initialize_molar_mass_array, &
                                                   check_initialization
    use test_musica_data,                    only: get_wavelength_edges
    use musica_ccpp_util,                    only: PI

    implicit none

    integer, parameter                                                  :: NUM_SPECIES = 10
    integer, parameter                                                  :: NUM_TUVX_CONSTITUENTS = 1     ! cloud liquid water content
    integer, parameter                                                  :: NUM_TUVX_ONLY_GAS_SPECIES = 1 ! dry air
    integer, parameter                                                  :: NUM_TUVX_CONSTITUENTS_INCLUDE_SHARED_SPECIES = 4
    integer, parameter                                                  :: NUM_NO_CONSTITUENTS = 2
    ! This test requires that the number of grid cells = 4, which is the default
    ! vector dimension for MICM. This restriction will be removed once
    ! https://github.com/NCAR/musica/issues/217 is finished.
    integer, parameter                                                  :: NUM_COLUMNS = 2
    integer, parameter                                                  :: NUM_LAYERS = 2
    integer, parameter                                                  :: NUM_WAVELENGTH_BINS = 102
    integer, parameter                                                  :: NUM_MICM_RATE_PARAMETERS = 4
    integer, parameter                                                  :: INDEX_MICM_JNO = 3
    integer                                                             :: errcode
    character(len=512)                                                  :: errmsg
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS+1)                   :: photolysis_wavelength_grid_interfaces        ! m
    integer, parameter                                                  :: num_photolysis_wavelength_grid_sections = 8  ! (count)
    real(kind_phys), dimension(num_photolysis_wavelength_grid_sections) :: extraterrestrial_flux                        ! photons cm-2 s-1 nm-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                  :: dry_air_density                              ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
           NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES) :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS)                             :: solar_zenith_angle                           ! radians
    real(kind_phys), dimension(NUM_COLUMNS, NUM_LAYERS+1)               :: height_at_interfaces                         ! km
    real(kind_phys), dimension(NUM_COLUMNS, NUM_LAYERS, &
           NUM_MICM_RATE_PARAMETERS)                                    :: rate_parameters
    type(ccpp_constituent_prop_ptr_t),   allocatable                    :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target            :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer                        :: const_prop
    real(kind_phys)                                                     :: base_conc
    character(len=256)                                                  :: constituent_name
    integer                                                             :: i, j, k, i_col
    integer                                                             :: Ar_index, CO2_index, H2O_index, N_index
    integer                                                             :: O_index, O1D_index, O2_index, O3_index
    integer                                                             :: N2_index, NO_index, cloud_index, air_index

    filename_of_micm_configuration = 'musica_configurations/NO_photolysis/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/NO_photolysis/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/NO_photolysis/tuvx_micm_mapping.json'

    call get_wavelength_edges(photolysis_wavelength_grid_interfaces)
    dry_air_density(:,1) = (/ 0.95_kind_phys, 0.93_kind_phys /)
    dry_air_density(:,2) = (/ 1.2255_kind_phys, 1.15_kind_phys /)
    extraterrestrial_flux(:) = &
      (/ 1.5e13_kind_phys, 1.5e13_kind_phys, 1.4e13_kind_phys, 1.4e13_kind_phys, &
        1.3e13_kind_phys, 1.2e13_kind_phys, 1.1e13_kind_phys, 1.0e13_kind_phys /)
    solar_zenith_angle(:) = (/ 0.0_kind_phys, 2.1_kind_phys /)
    height_at_interfaces(1,:) = (/ 3.01_kind_phys, 1.01_kind_phys, 0.01_kind_phys /)
    height_at_interfaces(2,:) = (/ 6.01_kind_phys, 1.01_kind_phys, 0.51_kind_phys /)

    call musica_ccpp_register( constituent_props, errmsg, errcode )
    ASSERT(errcode == 0)

    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set( const_prop, errcode, errmsg )
      ASSERT(errcode == 0)
    end do

    do i = 1, size(constituent_props)
      call constituent_props(i)%standard_name(constituent_name)
      if (constituent_name == "O2") then
        O2_index = i
        base_conc = 0.21_kind_phys
      else if (constituent_name == "Ar") then
        Ar_index = i
        base_conc = 1.0e-9_kind_phys
      else if (constituent_name == "CO2") then
        CO2_index = i
        base_conc = 1.0e-9_kind_phys
      else if (constituent_name == "H2O") then
        H2O_index = i
        base_conc = 1.0e-9_kind_phys
      else if (constituent_name == "N") then
        N_index = i
        base_conc = 1.0e-9_kind_phys
      else if (constituent_name == "NO") then
        NO_index = i
        base_conc = 1.0e-9_kind_phys
      else if (constituent_name == "O") then
        O_index = i
        base_conc = 1.0e-9_kind_phys
      else if (constituent_name == "O1D") then
        O1D_index = i
        base_conc = 1.0e-9_kind_phys
      else if (constituent_name == "O3") then
        O3_index = i
        base_conc = 1.0e-4_kind_phys
      else if (constituent_name == "N2") then
        N2_index = i
        base_conc = 0.79_kind_phys
      else if (constituent_name == "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water") then
        cloud_index = i
        base_conc = 1.0e-3_kind_phys
      else if (constituent_name == "air") then
        air_index = i
        base_conc = 1.2_kind_phys
      else
        write(*,*) "Unknown species: ", constituent_name
        stop 3
      endif
    end do

    call musica_ccpp_init( NUM_COLUMNS, NUM_LAYERS, NUM_LAYERS+1, photolysis_wavelength_grid_interfaces, &
                          constituent_props_ptr, errmsg, errcode )
    ASSERT(errcode == 0)

    ! check species indices
    ASSERT( O2_index    == tuvx_indices_constituent_props( tuvx_index_O2 ) )
    ASSERT( O3_index    == tuvx_indices_constituent_props( tuvx_index_O3 ) )
    ASSERT( NO_index    == tuvx_indices_constituent_props( tuvx_index_NO ) )
    ASSERT( N2_index    == tuvx_indices_constituent_props( tuvx_index_N2 ) )
    ASSERT( cloud_index == tuvx_indices_constituent_props( index_cloud_liquid_water_content ) )
    ASSERT( air_index   == tuvx_indices_constituent_props( index_dry_air ) )

    ! Set initial species concentrations
    do i = 1, size(constituent_props)
      do j = 1, NUM_COLUMNS
        do k = 1, NUM_LAYERS
          constituents(j,k,i) = base_conc * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
        end do
      end do
    end do

    rate_parameters = 0.0_kind_phys
    call calculate_NO_photolysis_rate_constants( &
                solar_zenith_angle, extraterrestrial_flux, &
                height_at_interfaces, dry_air_density,     &
                constituents, rate_parameters )
    do i_col = 1, size(dry_air_density, dim=1)
      ! TODO: figure out why the uppermost layer has a jNO of 0.0
      do i = 1, NUM_LAYERS-1
        do j = 1, NUM_MICM_RATE_PARAMETERS
          if (j == INDEX_MICM_JNO) then
            if (solar_zenith_angle(i_col) <= 110.0_kind_phys / 180.0_kind_phys * PI) then
              ASSERT( rate_parameters(i_col,i,j) > 0.0_kind_phys )
            else
              ASSERT( rate_parameters(i_col,i,j) == 0.0_kind_phys )
            end if
          else
            ASSERT( rate_parameters(i_col,i,j) == 0.0_kind_phys )
          end if
        end do
      end do
    end do

    call musica_ccpp_final(errmsg, errcode)
    ASSERT(errcode == 0)

  end subroutine test_calculate_NO_photolysis_rate

end program test_tuvx_no_photolysis