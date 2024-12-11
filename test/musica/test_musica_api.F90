program run_test_musica_ccpp

  use ccpp_kinds, only: kind_phys
  use musica_ccpp

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  real(kind_phys), parameter :: DEGREE_TO_RADIAN = 3.14159265358979323846_kind_phys / 180.0_kind_phys

  call test_chapman()
  call test_terminator()

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

  !> Tests the Chapman chemistry scheme
  subroutine test_chapman()
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_micm,          only: micm
    use musica_ccpp_namelist,      only: filename_of_micm_configuration, &
                                         filename_of_tuvx_configuration, &
                                         filename_of_tuvx_micm_mapping_configuration

    implicit none

    integer, parameter                                             :: NUM_SPECIES = 5
    integer, parameter                                             :: NUM_TUVX_CONSTITUENTS = 1
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
    real(kind_phys)                                                       :: surface_albedo                               ! unitless
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
                               NUM_SPECIES+NUM_TUVX_CONSTITUENTS)         :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
                               NUM_SPECIES+NUM_TUVX_CONSTITUENTS)         :: initial_constituents                         ! kg kg-1
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
    surface_albedo = 0.10_kind_phys
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
    ASSERT(size(constituent_props) == NUM_SPECIES+NUM_TUVX_CONSTITUENTS)
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
                 molar_mass == 0.018_kind_phys .and. is_advected)
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
                          surface_temperature, surface_albedo, num_photolysis_wavelength_grid_sections, &
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
        ! cloud liquid water mixing ratio should be unchanged
        ASSERT_NEAR(constituents(i,j,NUM_SPECIES+1), initial_constituents(i,j,NUM_SPECIES+1), 1.0e-13)
        ASSERT_NEAR(total_O, total_O_init, 1.0e-13)
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
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_micm,          only: micm
    use musica_ccpp_namelist,      only: filename_of_micm_configuration, &
                                         filename_of_tuvx_configuration, &
                                         filename_of_tuvx_micm_mapping_configuration

    implicit none

    integer, parameter                                             :: NUM_SPECIES = 2
    integer, parameter                                             :: NUM_TUVX_CONSTITUENTS = 1
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
    real(kind_phys)                                                       :: surface_albedo                               ! unitless
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
                               NUM_SPECIES+NUM_TUVX_CONSTITUENTS)         :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
                               NUM_SPECIES+NUM_TUVX_CONSTITUENTS)         :: initial_constituents                         ! kg kg-1
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
    surface_albedo = 0.10_kind_phys
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
    ASSERT(size(constituent_props) == NUM_SPECIES+NUM_TUVX_CONSTITUENTS)
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
                    .and. molar_mass == 0.018_kind_phys .and. is_advected)
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
                          surface_temperature, surface_albedo, num_photolysis_wavelength_grid_sections, &
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

end program run_test_musica_ccpp