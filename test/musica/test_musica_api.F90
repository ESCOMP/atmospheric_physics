program run_test_musica_ccpp

  use ccpp_kinds,       only: kind_phys
  use test_musica_data, only: get_wavelength_edges
  use musica_ccpp

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  real(kind_phys), parameter :: DEGREE_TO_RADIAN = 3.14159265358979323846_kind_phys / 180.0_kind_phys

  write(*,*) "[MUSICA Test] Running the Chapman test"
  call test_chapman()
  write(*,*) "[MUSICA Test] Ends the Chapman test"

  write(*,*) "[MUSICA Test] Running the Terminator test"
  call test_terminator()
  write(*,*) "[MUSICA Test] Ends the Terminator test"

  write(*,*) "[MUSICA Test] Running the TUVX with NO photolysis test"
  call test_chapman_with_no_photolysis()
  write(*,*) "[MUSICA Test] Ends the TUVX with NO photolysis test"

contains

  !> Tests the Chapman chemistry scheme
  subroutine test_chapman()
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_properties_t
    use musica_ccpp_micm,              only: micm
    use musica_ccpp_namelist,          only: filename_of_micm_configuration, &
                                             filename_of_tuvx_configuration, &
                                             filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_load_species, only: index_dry_air, index_O2, index_O3

    implicit none

    integer, parameter                                                    :: NUM_SPECIES = 5
    integer, parameter                                                    :: NUM_TUVX_CONSTITUENTS = 1
    integer, parameter                                                    :: NUM_TUVX_ONLY_GAS_SPECIES = 1
    ! This test requires that the number of grid cells = 4, which is the default
    ! vector dimension for MICM. This restriction will be removed once
    ! https://github.com/NCAR/musica/issues/217 is finished.
    integer, parameter                                                    :: NUM_COLUMNS = 2
    integer, parameter                                                    :: NUM_LAYERS = 2
    integer, parameter                                                    :: NUM_WAVELENGTH_BINS = 102
    integer                                                               :: NUM_GRID_CELLS = NUM_COLUMNS * NUM_LAYERS
    integer                                                               :: errcode
    character(len=512)                                                    :: errmsg
    real(kind_phys)                                                       :: time_step = 60._kind_phys                    ! s
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS+1)                     :: photolysis_wavelength_grid_interfaces        ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: geopotential_height_wrt_surface_at_midpoint  ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS+1)                  :: geopotential_height_wrt_surface_at_interface ! m
    real(kind_phys), dimension(NUM_COLUMNS)                               :: surface_geopotential                         ! m2 s-2
    real(kind_phys), dimension(NUM_COLUMNS)                               :: surface_temperature                          ! K
    real(kind_phys), dimension(NUM_COLUMNS)                               :: surface_albedo                               ! fraction
    integer, parameter                                                    :: num_photolysis_wavelength_grid_sections = 8  ! (count)
    real(kind_phys), dimension(num_photolysis_wavelength_grid_sections+1) :: flux_data_photolysis_wavelength_interfaces   ! nm
    real(kind_phys), dimension(num_photolysis_wavelength_grid_sections)   :: extraterrestrial_flux                        ! photons cm-2 s-1 nm-1
    real(kind_phys)                                                       :: standard_gravitational_acceleration          ! s2 m-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: temperature                                  ! K
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: pressure                                     ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: dry_air_density                              ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: cloud_area_fraction                          ! fraction
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: air_pressure_thickness                       ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
      NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)        :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
        NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)      :: initial_constituents                         ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS)                               :: solar_zenith_angle                           ! radians
    real(kind_phys)                                                       :: earth_sun_distance                           ! AU
    type(ccpp_constituent_prop_ptr_t),   allocatable                      :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target              :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer                          :: const_prop
    real(kind_phys)                                                       :: molar_mass, base_conc
    character(len=512)                                                    :: species_name, units
    character(len=:), allocatable                                         :: micm_species_name
    logical                                                               :: tmp_bool, is_advected
    integer                                                               :: i, j, k
    integer                                                               :: N2_index, O2_index, O_index, O1D_index, O3_index
    real(kind_phys)                                                       :: total_O, total_O_init

    call get_wavelength_edges(photolysis_wavelength_grid_interfaces)
    time_step = 60._kind_phys
    geopotential_height_wrt_surface_at_midpoint(1,:) = (/ 2000.0_kind_phys, 500.0_kind_phys /)
    geopotential_height_wrt_surface_at_midpoint(2,:) = (/ 2000.0_kind_phys, -500.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(1,:) = (/ 3000.0_kind_phys, 1000.0_kind_phys, 0.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(2,:) = (/ 3000.0_kind_phys, 500.0_kind_phys, -1500.0_kind_phys /)
    surface_temperature = (/ 300.0_kind_phys, 300.0_kind_phys /)
    surface_geopotential = (/ 100.0_kind_phys, 200.0_kind_phys /)
    surface_albedo(:) = 0.10_kind_phys
    standard_gravitational_acceleration = 10.0_kind_phys
    temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    pressure(:,1) = (/ 6000.04_kind_phys, 7000.04_kind_phys /)
    pressure(:,2) = (/ 8000.04_kind_phys, 9000.04_kind_phys /)
    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
    flux_data_photolysis_wavelength_interfaces(:) = &
      (/ 200.0_kind_phys, 210.0_kind_phys, 220.0_kind_phys, 230.0_kind_phys, &
        240.0_kind_phys, 250.0_kind_phys, 260.0_kind_phys, 270.0_kind_phys, 280.0_kind_phys /)
    extraterrestrial_flux(:) = &
      (/ 1.5e13_kind_phys, 1.5e13_kind_phys, 1.4e13_kind_phys, 1.4e13_kind_phys, &
        1.3e13_kind_phys, 1.2e13_kind_phys, 1.1e13_kind_phys, 1.0e13_kind_phys /)
    cloud_area_fraction(:,1) = (/ 0.1_kind_phys, 0.2_kind_phys /)
    cloud_area_fraction(:,2) = (/ 0.3_kind_phys, 0.4_kind_phys /)
    air_pressure_thickness(:,1) = (/ 900.0_kind_phys, 905.0_kind_phys /)
    air_pressure_thickness(:,2) = (/ 910.0_kind_phys, 915.0_kind_phys /)
    solar_zenith_angle = (/ 0.0_kind_phys, 2.1_kind_phys /)
    earth_sun_distance = 1.04_kind_phys

    filename_of_micm_configuration = 'musica_configurations/chapman/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/chapman/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/chapman/tuvx_micm_mapping.json'

    call musica_ccpp_register(constituent_props, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    ASSERT(allocated(constituent_props))
    ASSERT(size(constituent_props) == NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)
    do i = 1, size(constituent_props)
      ASSERT(constituent_props(i)%is_instantiated(errcode, errmsg))
      ASSERT(errcode == 0)
      call constituent_props(i)%standard_name(species_name, errcode, errmsg)
      ASSERT(errcode == 0)
      call constituent_props(i)%molar_mass(molar_mass, errcode, errmsg)
      ASSERT(errcode == 0)
      call constituent_props(i)%is_advected(is_advected, errcode, errmsg)
      ASSERT(errcode == 0)
      tmp_bool = (trim(species_name) == "O2" .and. molar_mass == 0.0319988_kind_phys .and. is_advected) .or.  &
                (trim(species_name) == "O" .and. molar_mass == 0.0159994_kind_phys .and. .not. is_advected) .or.   &
                (trim(species_name) == "O1D" .and. molar_mass == 0.0159994_kind_phys .and. .not. is_advected) .or. &
                (trim(species_name) == "O3" .and. molar_mass == 0.0479982_kind_phys .and. is_advected) .or. &
                (trim(species_name) == "N2" .and. molar_mass == 0.0280134_kind_phys .and. is_advected) .or. &
                (trim(species_name) == "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water" .and. &
                 molar_mass == 0.018_kind_phys .and. is_advected) .or. &
                (trim(species_name) == "air" .and. molar_mass == 0.0289644_kind_phys .and. .not. is_advected)
      ASSERT(tmp_bool)
      call constituent_props(i)%units(units, errcode, errmsg)
      if (errcode /= 0) then
        write(*,*) errcode, trim(errmsg)
        stop 3
      endif
      ASSERT(trim(units) == 'kg kg-1')
    end do
    if (errcode /= 0) then
      write(*,*) errcode, trim(errmsg)
      stop 3
    end if

    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set(const_prop, errcode, errmsg)
    end do

    call musica_ccpp_init(NUM_COLUMNS, NUM_LAYERS, NUM_LAYERS+1, photolysis_wavelength_grid_interfaces, &
                          constituent_props_ptr, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    do i = 1, micm%species_ordering%size()
      micm_species_name = micm%species_ordering%name(i)
      if (micm_species_name == "O2") then
        O2_index = i
        base_conc = 0.21_kind_phys
      else if (micm_species_name == "O") then
        O_index = i
        base_conc = 1.0e-9_kind_phys
      else if (micm_species_name == "O1D") then
        O1D_index = i
        base_conc = 1.0e-9_kind_phys
      else if (micm_species_name == "O3") then
        O3_index = i
        base_conc = 1.0e-4_kind_phys
      else if (micm_species_name == "N2") then
        N2_index = i
        base_conc = 0.79_kind_phys
      else
        write(*,*) "Unknown species: ", micm_species_name
        stop 3
      endif
      do j = 1, NUM_COLUMNS
        do k = 1, NUM_LAYERS
          constituents(j,k,i) = base_conc * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
        end do
      end do
    end do
    ! set initial cloud liquid water mixing ratio to ~1e-3 kg kg-1
    do j = 1, NUM_COLUMNS
      do k = 1, NUM_LAYERS
        constituents(j,k,NUM_SPECIES+1) = 1.0e-3_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
      end do
    end do
    do j = 1, NUM_COLUMNS
      do k = 1, NUM_LAYERS
        constituents(j,k,NUM_SPECIES+index_dry_air) = 1.0e-3_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
      end do
    end do
    initial_constituents(:,:,:) = constituents(:,:,:)

    write(*,*) "[MUSICA INFO] Initial Time Step"
    write(*,fmt="(1x,f10.2)") time_step
    write(*,*) "[MUSICA INFO] Initial Temperature"
    write(*,fmt="(4(1x,f10.4))") temperature
    write(*,*) "[MUSICA INFO] Initial Pressure"
    write(*,fmt="(4(1x,f10.4))") pressure
    write(*,*) "[MUSICA INFO] Initial Concentrations"
    write(*,fmt="(4(3x,e13.6))") constituents

    call musica_ccpp_run( time_step, temperature, pressure, dry_air_density, constituent_props_ptr,     &
                          constituents, geopotential_height_wrt_surface_at_midpoint,                    &
                          geopotential_height_wrt_surface_at_interface, surface_geopotential,           &
                          surface_temperature, surface_albedo,                                          &
                          flux_data_photolysis_wavelength_interfaces, extraterrestrial_flux,            &
                          standard_gravitational_acceleration, cloud_area_fraction,                     &
                          air_pressure_thickness, solar_zenith_angle, earth_sun_distance, errmsg,       &
                          errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    write(*,*) "[MUSICA INFO] Solved Concentrations"
    write(*,fmt="(4(3x,e13.6))") constituents

    call musica_ccpp_final(errmsg, errcode)

    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    ! Check the results
    do i = 1, NUM_COLUMNS
      do j = 1, NUM_LAYERS
        ! N2 should be unchanged
        ASSERT_NEAR(constituents(i,j,N2_index), initial_constituents(i,j,N2_index), 1.0e-13)
        ! O3 and O2 should be relatively unchanged
        ASSERT_NEAR(constituents(i,j,O3_index), initial_constituents(i,j,O3_index), 1.0e-4)
        ASSERT_NEAR(constituents(i,j,O2_index), initial_constituents(i,j,O2_index), 1.0e-4)
        ! O and O1D should be < 1e-10 kg kg-1
        ASSERT(constituents(i,j,O_index) < 1.0e-10)
        ASSERT(constituents(i,j,O1D_index) < 1.0e-10)
        ! total O mass should be conserved
        total_O = constituents(i,j,O_index) + constituents(i,j,O1D_index) + &
                  constituents(i,j,O2_index) + constituents(i,j,O3_index)
        total_O_init = initial_constituents(i,j,O_index) + initial_constituents(i,j,O1D_index) + &
                       initial_constituents(i,j,O2_index) + initial_constituents(i,j,O3_index)
        ASSERT_NEAR(total_O, total_O_init, 1.0e-13)
        ! cloud liquid water mixing ratio should be unchanged
        ASSERT_NEAR(constituents(i,j,NUM_SPECIES+1), initial_constituents(i,j,NUM_SPECIES+1), 1.0e-13)
      end do
    end do
    do j = 1, NUM_LAYERS
      ! O and O1D should be lower in the nighttime column
      ASSERT(constituents(2,j,O_index) < constituents(1,j,O_index))
      ASSERT(constituents(2,j,O1D_index) < constituents(1,j,O1D_index))
    end do

    deallocate(constituent_props_ptr)

  end subroutine test_chapman

  !> Tests the simple Terminator chemistry scheme
  subroutine test_terminator()
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_properties_t
    use musica_ccpp_micm,              only: micm
    use musica_ccpp_namelist,          only: filename_of_micm_configuration, &
                                             filename_of_tuvx_configuration, &
                                             filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_load_species, only: index_dry_air, index_O2, index_O3
    implicit none

    integer, parameter                                                    :: NUM_SPECIES = 2
    integer, parameter                                                    :: NUM_TUVX_CONSTITUENTS = 1
    integer, parameter                                                    :: NUM_TUVX_ONLY_GAS_SPECIES = 3
    ! This test requires that the number of grid cells = 4, which is the default
    ! vector dimension for MICM. This restriction will be removed once
    ! https://github.com/NCAR/musica/issues/217 is finished.
    integer, parameter                                                    :: NUM_COLUMNS = 2
    integer, parameter                                                    :: NUM_LAYERS = 2
    integer, parameter                                                    :: NUM_WAVELENGTH_BINS = 102
    integer                                                               :: NUM_GRID_CELLS = NUM_COLUMNS * NUM_LAYERS
    integer                                                               :: errcode
    character(len=512)                                                    :: errmsg
    real(kind_phys)                                                       :: time_step = 60._kind_phys                    ! s
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS+1)                     :: photolysis_wavelength_grid_interfaces        ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: geopotential_height_wrt_surface_at_midpoint  ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS+1)                  :: geopotential_height_wrt_surface_at_interface ! m
    real(kind_phys), dimension(NUM_COLUMNS)                               :: surface_geopotential                         ! m2 s-2
    real(kind_phys), dimension(NUM_COLUMNS)                               :: surface_temperature                          ! K
    real(kind_phys), dimension(NUM_COLUMNS)                               :: surface_albedo                               ! fraction
    integer, parameter                                                    :: num_photolysis_wavelength_grid_sections = 8  ! (count)
    real(kind_phys), dimension(num_photolysis_wavelength_grid_sections+1) :: flux_data_photolysis_wavelength_interfaces   ! nm
    real(kind_phys), dimension(num_photolysis_wavelength_grid_sections)   :: extraterrestrial_flux                        ! photons cm-2 s-1 nm-1
    real(kind_phys)                                                       :: standard_gravitational_acceleration          ! s2 m-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: temperature                                  ! K
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: pressure                                     ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: dry_air_density                              ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: cloud_area_fraction                          ! fraction
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: air_pressure_thickness                       ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
      NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)        :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
      NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)        :: initial_constituents                         ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS)                               :: solar_zenith_angle                           ! radians
    real(kind_phys)                                                       :: earth_sun_distance                           ! AU
    type(ccpp_constituent_prop_ptr_t),   allocatable                      :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target              :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer                          :: const_prop
    real(kind_phys)                                                       :: molar_mass, base_conc
    character(len=512)                                                    :: species_name, units
    character(len=:), allocatable                                         :: micm_species_name
    logical                                                               :: tmp_bool, is_advected
    integer                                                               :: i, j, k
    integer                                                               :: Cl_index, Cl2_index
    real(kind_phys)                                                       :: total_Cl, total_Cl_init

    call get_wavelength_edges(photolysis_wavelength_grid_interfaces)
    time_step = 60._kind_phys
    geopotential_height_wrt_surface_at_midpoint(1,:) = (/ 2000.0_kind_phys, 500.0_kind_phys /)
    geopotential_height_wrt_surface_at_midpoint(2,:) = (/ 2000.0_kind_phys, -500.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(1,:) = (/ 3000.0_kind_phys, 1000.0_kind_phys, 0.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(2,:) = (/ 3000.0_kind_phys, 500.0_kind_phys, -1500.0_kind_phys /)
    surface_temperature = (/ 300.0_kind_phys, 300.0_kind_phys /)
    surface_geopotential = (/ 100.0_kind_phys, 200.0_kind_phys /)
    surface_albedo(:) = 0.10_kind_phys
    standard_gravitational_acceleration = 10.0_kind_phys
    temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    pressure(:,1) = (/ 6000.04_kind_phys, 7000.04_kind_phys /)
    pressure(:,2) = (/ 8000.04_kind_phys, 9000.04_kind_phys /)
    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
    flux_data_photolysis_wavelength_interfaces(:) = &
      (/ 200.0_kind_phys, 210.0_kind_phys, 220.0_kind_phys, 230.0_kind_phys, &
        240.0_kind_phys, 250.0_kind_phys, 260.0_kind_phys, 270.0_kind_phys, 280.0_kind_phys /)
    extraterrestrial_flux(:) = &
      (/ 1.5e13_kind_phys, 1.5e13_kind_phys, 1.4e13_kind_phys, 1.4e13_kind_phys, &
        1.3e13_kind_phys, 1.2e13_kind_phys, 1.1e13_kind_phys, 1.0e13_kind_phys /)
    cloud_area_fraction(:,1) = (/ 0.1_kind_phys, 0.2_kind_phys /)
    cloud_area_fraction(:,2) = (/ 0.3_kind_phys, 0.4_kind_phys /)
    air_pressure_thickness(:,1) = (/ 900.0_kind_phys, 905.0_kind_phys /)
    air_pressure_thickness(:,2) = (/ 910.0_kind_phys, 915.0_kind_phys /)
    solar_zenith_angle = (/ 0.0_kind_phys, 2.1_kind_phys /)
    earth_sun_distance = 1.04_kind_phys

    filename_of_micm_configuration = 'musica_configurations/terminator/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/terminator/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/terminator/tuvx_micm_mapping.json'

    call musica_ccpp_register(constituent_props, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    ASSERT(allocated(constituent_props))
    ASSERT(size(constituent_props) == NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)
    do i = 1, size(constituent_props)
      ASSERT(constituent_props(i)%is_instantiated(errcode, errmsg))
      ASSERT(errcode == 0)
      call constituent_props(i)%standard_name(species_name, errcode, errmsg)
      ASSERT(errcode == 0)
      call constituent_props(i)%molar_mass(molar_mass, errcode, errmsg)
      ASSERT(errcode == 0)
      call constituent_props(i)%is_advected(is_advected, errcode, errmsg)
      ASSERT(errcode == 0)
      tmp_bool = (trim(species_name) == "Cl" .and. molar_mass == 0.035453_kind_phys .and. is_advected) .or.  &
                 (trim(species_name) == "Cl2" .and. molar_mass == 0.070906_kind_phys .and. is_advected) .or. &
                 (trim(species_name) == "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water" &
                                        .and. molar_mass == 0.018_kind_phys .and. is_advected) .or. &
                 (trim(species_name) == "air" .and. molar_mass == 0.0289644_kind_phys .and. .not. is_advected) .or. &
                 (trim(species_name) == "O2" .and. molar_mass == 0.0319988_kind_phys .and. .not. is_advected) .or. &
                 (trim(species_name) == "O3" .and. molar_mass == 0.0479982_kind_phys .and. .not. is_advected)
      ASSERT(tmp_bool)
      call constituent_props(i)%units(units, errcode, errmsg)
      if (errcode /= 0) then
        write(*,*) errcode, trim(errmsg)
        stop 3
      endif
      ASSERT(trim(units) == 'kg kg-1')
    end do
    if (errcode /= 0) then
      write(*,*) errcode, trim(errmsg)
      stop 3
    end if

    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set(const_prop, errcode, errmsg)
    end do

    call musica_ccpp_init(NUM_COLUMNS, NUM_LAYERS, NUM_LAYERS+1, photolysis_wavelength_grid_interfaces, &
                          constituent_props_ptr, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    do i = 1, micm%species_ordering%size()
      micm_species_name = micm%species_ordering%name(i)
      if (micm_species_name == "Cl") then
        Cl_index = i
        base_conc = 1.0e-10_kind_phys
      else if (micm_species_name == "Cl2") then
        Cl2_index = i
        base_conc = 1.0e-6_kind_phys
      else
        write(*,*) "Unknown species: ", micm_species_name
        stop 3
      endif
      do j = 1, NUM_COLUMNS
        do k = 1, NUM_LAYERS
          constituents(j,k,i) = base_conc * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
        end do
      end do
    end do
    ! set initial cloud liquid water mixing ratio to ~1e-3 kg kg-1
    do j = 1, NUM_COLUMNS
      do k = 1, NUM_LAYERS
        constituents(j,k,NUM_SPECIES+1) = 1.0e-3_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
      end do
    end do
    ! set concentrations for TUV-x gas species
    do j = 1, NUM_COLUMNS
      do k = 1, NUM_LAYERS
        constituents(j,k,NUM_SPECIES+index_dry_air) = 1.0e-3_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
        constituents(j,k,NUM_SPECIES+index_O2) = 1.0e-3_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
        constituents(j,k,NUM_SPECIES+index_O3) = 1.0e-3_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
      end do
    end do
    initial_constituents(:,:,:) = constituents(:,:,:)

    write(*,*) "[MUSICA INFO] Initial Time Step"
    write(*,fmt="(1x,f10.2)") time_step
    write(*,*) "[MUSICA INFO] Initial Temperature"
    write(*,fmt="(4(1x,f10.4))") temperature
    write(*,*) "[MUSICA INFO] Initial Pressure"
    write(*,fmt="(4(1x,f10.4))") pressure
    write(*,*) "[MUSICA INFO] Initial Concentrations"
    write(*,fmt="(4(3x,e13.6))") constituents

    call musica_ccpp_run( time_step, temperature, pressure, dry_air_density, constituent_props_ptr,     &
                          constituents, geopotential_height_wrt_surface_at_midpoint,                    &
                          geopotential_height_wrt_surface_at_interface, surface_geopotential,           &
                          surface_temperature, surface_albedo,                                          &
                          flux_data_photolysis_wavelength_interfaces, extraterrestrial_flux,            &
                          standard_gravitational_acceleration, cloud_area_fraction,                     &
                          air_pressure_thickness, solar_zenith_angle, earth_sun_distance, errmsg,       &
                          errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    write(*,*) "[MUSICA INFO] Solved Concentrations"
    write(*,fmt="(4(3x,e13.6))") constituents

    call musica_ccpp_final(errmsg, errcode)

    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    ! Check the results
    do i = 1, NUM_COLUMNS
      do j = 1, NUM_LAYERS
        ! Cl2 should not be unchanged
        ASSERT(constituents(i,j,Cl2_index) .ne. initial_constituents(i,j,Cl2_index))
        ! Cl should not be unchanged
        ASSERT(constituents(i,j,Cl_index) .ne. initial_constituents(i,j,Cl_index))
        ! total Cl mass should be conserved
        total_Cl = constituents(i,j,Cl_index) + constituents(i,j,Cl2_index)
        total_Cl_init = initial_constituents(i,j,Cl_index) + initial_constituents(i,j,Cl2_index)
        ASSERT_NEAR(total_Cl, total_Cl_init, 1.0e-13)
        ! cloud liquid water should be unchanged
        ASSERT_NEAR(constituents(i,j,NUM_SPECIES+1), initial_constituents(i,j,NUM_SPECIES+1), 1.0e-13)
      end do
    end do
    do j = 1, NUM_LAYERS
      ! Cl should be lower in the nighttime column
      ASSERT(constituents(2,j,Cl_index) < constituents(1,j,Cl_index))
      ! Cl2 should be higher in the nighttime column
      ASSERT(constituents(2,j,Cl2_index) > constituents(1,j,Cl2_index))
    end do

    deallocate(constituent_props_ptr)

  end subroutine test_terminator

  !> Tests the Chapman chemistry scheme but with an NO photolysis reaction
  subroutine test_chapman_with_no_photolysis()
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_micm,          only: micm
    use musica_ccpp_namelist,      only: filename_of_micm_configuration, &
                                         filename_of_tuvx_configuration, &
                                         filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_util,          only: PI

    implicit none

    integer, parameter                                                    :: NUM_SPECIES = 10
    integer, parameter                                                    :: NUM_TUVX_CONSTITUENTS = 1
    integer, parameter                                                    :: NUM_TUVX_ONLY_GAS_SPECIES = 1
    ! This test requires that the number of grid cells = 4, which is the default
    ! vector dimension for MICM. This restriction will be removed once
    ! https://github.com/NCAR/musica/issues/217 is finished.
    integer, parameter                                                    :: NUM_COLUMNS = 2
    integer, parameter                                                    :: NUM_LAYERS = 2
    integer, parameter                                                    :: NUM_WAVELENGTH_BINS = 102
    integer                                                               :: NUM_GRID_CELLS = NUM_COLUMNS * NUM_LAYERS
    integer                                                               :: errcode
    character(len=512)                                                    :: errmsg
    real(kind_phys)                                                       :: time_step = 60._kind_phys                    ! s
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS+1)                     :: photolysis_wavelength_grid_interfaces        ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: geopotential_height_wrt_surface_at_midpoint  ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS+1)                  :: geopotential_height_wrt_surface_at_interface ! m
    real(kind_phys), dimension(NUM_COLUMNS)                               :: surface_geopotential                         ! m2 s-2
    real(kind_phys), dimension(NUM_COLUMNS)                               :: surface_temperature                          ! K
    real(kind_phys), dimension(NUM_COLUMNS)                               :: surface_albedo                               ! unitless
    integer, parameter                                                    :: num_photolysis_wavelength_grid_sections = 8  ! (count)
    real(kind_phys), dimension(num_photolysis_wavelength_grid_sections+1) :: flux_data_photolysis_wavelength_interfaces   ! nm
    real(kind_phys), dimension(num_photolysis_wavelength_grid_sections)   :: extraterrestrial_flux                        ! photons cm-2 s-1 nm-1
    real(kind_phys)                                                       :: standard_gravitational_acceleration          ! s2 m-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: temperature                                  ! K
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: pressure                                     ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: dry_air_density                              ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: cloud_area_fraction                          ! unitless
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                    :: air_pressure_thickness                       ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
        NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)      :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
        NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)      :: initial_constituents                         ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS)                               :: solar_zenith_angle                           ! radians
    real(kind_phys)                                                       :: earth_sun_distance                           ! AU
    type(ccpp_constituent_prop_ptr_t),   allocatable                      :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target              :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer                          :: const_prop
    real(kind_phys)                                                       :: molar_mass, base_conc
    character(len=512)                                                    :: species_name, units
    character(len=:), allocatable                                         :: micm_species_name
    logical                                                               :: tmp_bool, is_advected
    integer                                                               :: i, j, k
    integer                                                               :: NO_index

    call get_wavelength_edges(photolysis_wavelength_grid_interfaces)
    time_step = 60._kind_phys
    geopotential_height_wrt_surface_at_midpoint(1,:) = (/ 2000.0_kind_phys, 500.0_kind_phys /)
    geopotential_height_wrt_surface_at_midpoint(2,:) = (/ 2000.0_kind_phys, -500.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(1,:) = (/ 3000.0_kind_phys, 1000.0_kind_phys, 0.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(2,:) = (/ 3000.0_kind_phys, 500.0_kind_phys, -1500.0_kind_phys /)
    surface_temperature = (/ 300.0_kind_phys, 300.0_kind_phys /)
    surface_geopotential = (/ 100.0_kind_phys, 200.0_kind_phys /)
    surface_albedo(:) = 0.10_kind_phys
    standard_gravitational_acceleration = 10.0_kind_phys
    temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    pressure(:,1) = (/ 6000.04_kind_phys, 7000.04_kind_phys /)
    pressure(:,2) = (/ 8000.04_kind_phys, 9000.04_kind_phys /)
    dry_air_density(:,1) = (/ 0.95_kind_phys, 0.93_kind_phys /)
    dry_air_density(:,2) = (/ 1.225_kind_phys, 1.15_kind_phys /)
    flux_data_photolysis_wavelength_interfaces(:) = &
      (/ 200.0_kind_phys, 210.0_kind_phys, 220.0_kind_phys, 230.0_kind_phys, &
        240.0_kind_phys, 250.0_kind_phys, 260.0_kind_phys, 270.0_kind_phys, 280.0_kind_phys /)
    extraterrestrial_flux(:) = &
      (/ 1.5e13_kind_phys, 1.5e13_kind_phys, 1.4e13_kind_phys, 1.4e13_kind_phys, &
        1.3e13_kind_phys, 1.2e13_kind_phys, 1.1e13_kind_phys, 1.0e13_kind_phys /)
    cloud_area_fraction(:,1) = (/ 0.1_kind_phys, 0.2_kind_phys /)
    cloud_area_fraction(:,2) = (/ 0.3_kind_phys, 0.4_kind_phys /)
    air_pressure_thickness(:,1) = (/ 900.0_kind_phys, 905.0_kind_phys /)
    air_pressure_thickness(:,2) = (/ 910.0_kind_phys, 915.0_kind_phys /)
    solar_zenith_angle = (/ 0.0_kind_phys, 2.1_kind_phys /)
    earth_sun_distance = 1.04_kind_phys

    filename_of_micm_configuration = 'musica_configurations/NO_photolysis/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/NO_photolysis/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/NO_photolysis/tuvx_micm_mapping.json'

    call musica_ccpp_register(constituent_props, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    ASSERT(allocated(constituent_props))
    ASSERT(size(constituent_props) == NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)

    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set(const_prop, errcode, errmsg)
    end do

    call musica_ccpp_init(NUM_COLUMNS, NUM_LAYERS, NUM_LAYERS+1, photolysis_wavelength_grid_interfaces, &
                          constituent_props_ptr, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    ASSERT(NUM_SPECIES == micm%species_ordering%size())

    constituents = 0.0_kind_phys
    do i = 1, micm%species_ordering%size()
      micm_species_name = micm%species_ordering%name(i)
      if (micm_species_name == "NO") then
        NO_index = i
      endif
      do j = 1, NUM_COLUMNS
        do k = 1, NUM_LAYERS
          constituents(j,k,i) = 0.21_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
        end do
      end do
    end do
    ! set initial cloud liquid water mixing ratio to ~1e-3 kg kg-1
    do j = 1, NUM_COLUMNS
      do k = 1, NUM_LAYERS
        constituents(j,k,NUM_SPECIES+1) = 1.0e-3_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
      end do
    end do
    initial_constituents(:,:,:) = constituents(:,:,:)

    write(*,*) "[MUSICA INFO] Initial Time Step"
    write(*,fmt="(1x,f10.2)") time_step
    write(*,*) "[MUSICA INFO] Initial Temperature"
    write(*,fmt="(4(1x,f10.4))") temperature
    write(*,*) "[MUSICA INFO] Initial Pressure"
    write(*,fmt="(4(1x,f10.4))") pressure
    write(*,*) "[MUSICA INFO] Initial Concentrations"
    write(*,fmt="(4(3x,e13.6))") constituents

    call musica_ccpp_run( time_step, temperature, pressure, dry_air_density, constituent_props_ptr,     &
                          constituents, geopotential_height_wrt_surface_at_midpoint,                    &
                          geopotential_height_wrt_surface_at_interface, surface_geopotential,           &
                          surface_temperature, surface_albedo,                                          &
                          flux_data_photolysis_wavelength_interfaces, extraterrestrial_flux,            &
                          standard_gravitational_acceleration, cloud_area_fraction,                     &
                          air_pressure_thickness, solar_zenith_angle, earth_sun_distance, errmsg,       &
                          errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    write(*,*) "[MUSICA INFO] Solved Concentrations"
    write(*,fmt="(4(3x,e13.6))") constituents

    call musica_ccpp_final(errmsg, errcode)

    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    ! Check the results
    do i = 1, NUM_COLUMNS
      do j = 1, NUM_LAYERS
        ! NO should be different
        ! TODO: find a more accurate way to check this
        ASSERT(constituents(i,j,NO_index) .ne. initial_constituents(i,j,NO_index))
      end do
    end do

    deallocate(constituent_props_ptr)

  end subroutine test_chapman_with_no_photolysis

end program run_test_musica_ccpp