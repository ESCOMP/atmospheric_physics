! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program run_test_musica_ccpp

  use musica_ccpp
  use musica_test_data, only: get_wavelength_edges, get_extrarterrestrial_fluxes
  use ccpp_kinds,       only: kind_phys

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) \
  if ( (abs(a - b) > abs_error) .and. (abs(a - b) /= 0.0) .and. (a /= 0.0) .and. (b /= 0.0) ) then; \
    write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a = ", a, ", b = ", b; stop 1; \
  endif

  real(kind_phys), parameter :: DEGREE_TO_RADIAN = 3.14159265358979323846_kind_phys / 180.0_kind_phys
  real(kind_phys), parameter :: AVOGADRO = 6.02214179e23_kind_phys        ! mol-1
  real(kind_phys), parameter :: GAS_CONSTANT = 8.31446261815324_kind_phys ! J K-1 mol-1
  real(kind_phys), parameter :: MOLAR_MASS_DRY_AIR = 0.0289644_kind_phys  ! kg mol-1
  real(kind_phys), parameter :: MOLAR_MASS_DRY_AIR__G_MOL = MOLAR_MASS_DRY_AIR * 1.0e3_kind_phys ! g mol-1

  integer, parameter :: STDOUT = 6

  type :: ArrheniusReaction
    real(kind_phys) :: A_ = 1.0
    real(kind_phys) :: B_ = 0.0
    real(kind_phys) :: C_ = 0.0
    real(kind_phys) :: D_ = 300.0
    real(kind_phys) :: E_ = 0.0
  end type ArrheniusReaction

  write(*,*) "[MUSICA Test] Running the Chapman test"
  call test_chapman()
  write(*,*) "[MUSICA Test] Ends the Chapman test"

  write(*,*) "[MUSICA Test] Running the Terminator test"
  call test_terminator()
  write(*,*) "[MUSICA Test] Ends the Terminator test"

  write(*,*) "[MUSICA Test] Running the Analytical test with Rosenbrock solver"
  call test_rosenbrock()
  write(*,*) "[MUSICA Test] Ends the Analytical test with Rosenbrock solver"

  write(*,*) "[MUSICA Test] Running the Analytical test with Backward Euler solver"
  call test_backward_euler()
  write(*,*) "[MUSICA Test] Ends the Analytical test with Backward Euler solver"

contains

  !> Calculate the rate constant for an Arrhenius reaction
  function calculate_arrhenius( reaction, temperature, pressure ) result( rate )
    type(ArrheniusReaction), intent(in) :: reaction
    real(kind_phys), intent(in) :: temperature
    real(kind_phys), intent(in) :: pressure
    real(kind_phys) :: rate
    rate = reaction%A_ * exp( reaction%C_ / temperature ) &
           * (temperature / reaction%D_) ** reaction%B_ &
           * (1.0 + reaction%E_ * pressure)
  end function calculate_arrhenius

  !> Tests the Chapman chemistry scheme
  subroutine test_chapman()
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_properties_t
    use ccpp_const_utils,              only: ccpp_const_get_idx
    use musica_ccpp_namelist,          only: filename_of_micm_configuration, &
                                             filename_of_tuvx_configuration, &
                                             filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_load_species, only: index_dry_air, index_O2, index_O3

    implicit none

    ! This test requires that the number of grid cells = 4, which is the default
    ! vector dimension for MICM. This restriction will be removed once
    ! https://github.com/NCAR/musica/issues/217 is finished.
    integer, parameter                                               :: NUM_COLUMNS = 2
    integer, parameter                                               :: NUM_LAYERS = 2
    integer, parameter                                               :: NUM_WAVELENGTH_BINS = 102
    integer                                                          :: NUM_GRID_CELLS = NUM_COLUMNS * NUM_LAYERS
    integer, parameter                                               :: NUM_SPECIES = 5
    integer, parameter                                               :: NUM_TUVX_CONSTITUENTS = 1
    integer, parameter                                               :: NUM_TUVX_ONLY_GAS_SPECIES = 1
    integer                                                          :: errcode
    character(len=512)                                               :: errmsg
    real(kind_phys)                                                  :: time_step = 60._kind_phys                    ! s
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS+1)                :: photolysis_wavelength_grid_interfaces        ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)               :: geopotential_height_wrt_surface_at_midpoint  ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS+1)             :: geopotential_height_wrt_surface_at_interface ! m
    real(kind_phys), dimension(NUM_COLUMNS)                          :: surface_geopotential                         ! m2 s-2
    real(kind_phys), dimension(NUM_COLUMNS)                          :: surface_temperature                          ! K
    real(kind_phys), dimension(NUM_COLUMNS)                          :: surface_albedo                               ! fraction
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS)                  :: extraterrestrial_radiation_flux              ! photons cm-2 s-1 nm-1
    real(kind_phys)                                                  :: standard_gravitational_acceleration          ! s2 m-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)               :: temperature                                  ! K
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)               :: pressure                                     ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)               :: dry_air_density                              ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)               :: cloud_area_fraction                          ! fraction
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)               :: air_pressure_thickness                       ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
      NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES)   :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
        NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES) :: initial_constituents                         ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS)                          :: solar_zenith_angle                           ! radians
    real(kind_phys)                                                  :: earth_sun_distance                           ! AU
    type(ccpp_constituent_prop_ptr_t),   allocatable                 :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target         :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer                     :: const_prop
    real(kind_phys)                                                  :: molar_mass, base_conc
    character(len=512)                                               :: species_name, units
    character(len=:), allocatable                                    :: micm_species_name
    logical                                                          :: tmp_bool, is_advected
    integer                                                          :: i, j
    integer                                                          :: N2_index, O_index, O1D_index, O2_index, O3_index
    integer                                                          :: cloud_index, air_index
    real(kind_phys)                                                  :: total_O, total_O_init

    call get_wavelength_edges(photolysis_wavelength_grid_interfaces)
    geopotential_height_wrt_surface_at_midpoint(1,:) = (/ 2000.0_kind_phys, 500.0_kind_phys /)
    geopotential_height_wrt_surface_at_midpoint(2,:) = (/ 2000.0_kind_phys, -500.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(1,:) = (/ 3000.0_kind_phys, 1000.0_kind_phys, 0.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(2,:) = (/ 3000.0_kind_phys, 500.0_kind_phys, -1500.0_kind_phys /)
    surface_temperature = (/ 300.0_kind_phys, 300.0_kind_phys /)
    surface_geopotential = (/ 100.0_kind_phys, 200.0_kind_phys /)
    surface_albedo(:) = 0.10_kind_phys
    call get_extrarterrestrial_fluxes(extraterrestrial_radiation_flux)
    standard_gravitational_acceleration = 10.0_kind_phys
    temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    pressure(:,1) = (/ 6000.04_kind_phys, 7000.04_kind_phys /)
    pressure(:,2) = (/ 8000.04_kind_phys, 9000.04_kind_phys /)
    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
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
                          constituent_props_ptr, MOLAR_MASS_DRY_AIR__G_MOL, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    call ccpp_const_get_idx( constituent_props_ptr, "N2", N2_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "O", O_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "O1D", O1D_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "O2", O2_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "O3", O3_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water", &
                             cloud_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "air", air_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    do i = 1, NUM_COLUMNS
      do j = 1, NUM_LAYERS
        constituents(i,j,N2_index) = 0.79_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,O_index) = 1.0e-9_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,O1D_index) = 1.0e-9_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,O2_index) = 0.21_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,O3_index) = 1.0e-4_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,cloud_index) = 1.0e-3_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,air_index) = 1.0e-3_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
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

    call musica_ccpp_run( time_step, temperature, pressure, dry_air_density, constituent_props_ptr, &
                          constituents, geopotential_height_wrt_surface_at_midpoint,                &
                          geopotential_height_wrt_surface_at_interface, surface_geopotential,       &
                          surface_temperature, surface_albedo, extraterrestrial_radiation_flux,     &
                          standard_gravitational_acceleration, cloud_area_fraction,                 &
                          air_pressure_thickness, solar_zenith_angle, earth_sun_distance, STDOUT,   &
                          errmsg, errcode )
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
        ASSERT_NEAR(constituents(i,j,N2_index), initial_constituents(i,j,N2_index), 1.0e-13_kind_phys)
        ! O, O1D, O2, and O3 should be changed
        ASSERT(constituents(i,j,O_index) .ne. initial_constituents(i,j,O_index))
        ASSERT(constituents(i,j,O1D_index) .ne. initial_constituents(i,j,O1D_index))
        ASSERT(constituents(i,j,O2_index) .ne. initial_constituents(i,j,O2_index))
        ASSERT(constituents(i,j,O3_index) .ne. initial_constituents(i,j,O3_index))
        ! total O mass should be conserved
        total_O = constituents(i,j,O_index) + constituents(i,j,O1D_index) + &
                  constituents(i,j,O2_index) + constituents(i,j,O3_index)
        total_O_init = initial_constituents(i,j,O_index) + initial_constituents(i,j,O1D_index) + &
                       initial_constituents(i,j,O2_index) + initial_constituents(i,j,O3_index)
        ASSERT_NEAR(total_O, total_O_init, 1.0e-13_kind_phys)
        ! cloud liquid water mixing ratio should be unchanged
        ASSERT_NEAR(constituents(i,j,cloud_index), initial_constituents(i,j,cloud_index), 1.0e-13_kind_phys)
        ! air should be unchanged
        ASSERT_NEAR(constituents(i,j,air_index), initial_constituents(i,j,air_index), 1.0e-13_kind_phys)

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
    use ccpp_const_utils,              only: ccpp_const_get_idx
    use musica_ccpp_namelist,          only: filename_of_micm_configuration, &
                                             filename_of_tuvx_configuration, &
                                             filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_load_species, only: index_dry_air, index_O2, index_O3

    implicit none

    ! This test requires that the number of grid cells = 4, which is the default
    ! vector dimension for MICM. This restriction will be removed once
    ! https://github.com/NCAR/musica/issues/217 is finished.
    integer, parameter                                             :: NUM_COLUMNS = 2
    integer, parameter                                             :: NUM_LAYERS = 2
    integer, parameter                                             :: NUM_WAVELENGTH_BINS = 102
    integer                                                        :: NUM_GRID_CELLS = NUM_COLUMNS * NUM_LAYERS
    integer, parameter                                             :: NUM_SPECIES = 2
    integer, parameter                                             :: NUM_TUVX_CONSTITUENTS = 1
    integer, parameter                                             :: NUM_TUVX_ONLY_GAS_SPECIES = 3
    integer                                                        :: errcode
    character(len=512)                                             :: errmsg
    real(kind_phys)                                                :: time_step = 60._kind_phys                    ! s
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS+1)              :: photolysis_wavelength_grid_interfaces        ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: geopotential_height_wrt_surface_at_midpoint  ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS+1)           :: geopotential_height_wrt_surface_at_interface ! m
    real(kind_phys), dimension(NUM_COLUMNS)                        :: surface_geopotential                         ! m2 s-2
    real(kind_phys), dimension(NUM_COLUMNS)                        :: surface_temperature                          ! K
    real(kind_phys), dimension(NUM_COLUMNS)                        :: surface_albedo                               ! fraction
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS)                :: extraterrestrial_radiation_flux              ! photons cm-2 s-1 nm-1
    real(kind_phys)                                                :: standard_gravitational_acceleration          ! s2 m-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: temperature                                  ! K
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: pressure                                     ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: dry_air_density                              ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: cloud_area_fraction                          ! fraction
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: air_pressure_thickness                       ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
      NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES) :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
      NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES) :: initial_constituents                         ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS)                        :: solar_zenith_angle                           ! radians
    real(kind_phys)                                                :: earth_sun_distance                           ! AU
    type(ccpp_constituent_prop_ptr_t),   allocatable               :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target       :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer                   :: const_prop
    real(kind_phys)                                                :: molar_mass, base_conc
    character(len=512)                                             :: species_name, units
    character(len=:), allocatable                                  :: micm_species_name
    logical                                                        :: tmp_bool, is_advected
    integer                                                        :: i, j
    integer                                                        :: Cl_index, Cl2_index
    integer                                                        :: cloud_index, air_index, O2_index, O3_index
    real(kind_phys)                                                :: total_Cl, total_Cl_init

    call get_wavelength_edges(photolysis_wavelength_grid_interfaces)
    geopotential_height_wrt_surface_at_midpoint(1,:) = (/ 2000.0_kind_phys, 500.0_kind_phys /)
    geopotential_height_wrt_surface_at_midpoint(2,:) = (/ 2000.0_kind_phys, -500.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(1,:) = (/ 3000.0_kind_phys, 1000.0_kind_phys, 0.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(2,:) = (/ 3000.0_kind_phys, 500.0_kind_phys, -1500.0_kind_phys /)
    surface_temperature = (/ 300.0_kind_phys, 300.0_kind_phys /)
    surface_geopotential = (/ 100.0_kind_phys, 200.0_kind_phys /)
    surface_albedo(:) = 0.10_kind_phys
    call get_extrarterrestrial_fluxes(extraterrestrial_radiation_flux)
    standard_gravitational_acceleration = 10.0_kind_phys
    temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    pressure(:,1) = (/ 6000.04_kind_phys, 7000.04_kind_phys /)
    pressure(:,2) = (/ 8000.04_kind_phys, 9000.04_kind_phys /)
    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
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
                          constituent_props_ptr, MOLAR_MASS_DRY_AIR__G_MOL, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    call ccpp_const_get_idx( constituent_props_ptr, "Cl", Cl_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "Cl2", Cl2_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water", &
                             cloud_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "air", air_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "O2", O2_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call ccpp_const_get_idx( constituent_props_ptr, "O3", O3_index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    do i = 1, NUM_COLUMNS
      do j = 1, NUM_LAYERS
        constituents(i,j,Cl_index) =  1.0e-10_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,Cl2_index) = 1.0e-6_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,cloud_index) = 1.0e-3_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,air_index) = 1.0e-3_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,O2_index) = 0.21_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
        constituents(i,j,O3_index) = 1.0e-4_kind_phys * (1.0 + 0.1 * (i-1) + 0.01 * (j-1))
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

    call musica_ccpp_run( time_step, temperature, pressure, dry_air_density, constituent_props_ptr, &
                          constituents, geopotential_height_wrt_surface_at_midpoint,                &
                          geopotential_height_wrt_surface_at_interface, surface_geopotential,       &
                          surface_temperature, surface_albedo, extraterrestrial_radiation_flux,     &
                          standard_gravitational_acceleration, cloud_area_fraction,                 &
                          air_pressure_thickness, solar_zenith_angle, earth_sun_distance, STDOUT,   &
                          errmsg, errcode )
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
        ASSERT_NEAR(total_Cl, total_Cl_init, 1.0e-13_kind_phys)
        ! cloud liquid water should be unchanged
        ASSERT_NEAR(constituents(i,j,cloud_index), initial_constituents(i,j,cloud_index), 1.0e-13_kind_phys)
        ! dry air should be unchanged
        ASSERT_NEAR(constituents(i,j,air_index), initial_constituents(i,j,air_index), 1.0e-13_kind_phys)
        ! O2 should be unchanged
        ASSERT_NEAR(constituents(i,j,O2_index), initial_constituents(i,j,O2_index), 1.0e-13_kind_phys)
        ! O3 should be unchanged
        ASSERT_NEAR(constituents(i,j,O3_index), initial_constituents(i,j,O3_index), 1.0e-13_kind_phys)
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

  subroutine get_index_and_molar_mass(constituent_props, species_name, index, molar_mass)
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_properties_t
    use ccpp_const_utils,              only: ccpp_const_get_idx

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props(:)
    character(len=*),                  intent(in)  :: species_name
    integer,                           intent(out) :: index
    real(kind_phys),                   intent(out) :: molar_mass

    character(len=512) :: errmsg
    integer            :: errcode

    call ccpp_const_get_idx( constituent_props, species_name, index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call constituent_props(index)%molar_mass(molar_mass, errcode, errmsg)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

  end subroutine get_index_and_molar_mass

  subroutine test_analytical(number_of_columns, number_of_layers, test_accuracy)
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_prop_ptr_t, &
                                             ccpp_constituent_properties_t
    use ccpp_const_utils,              only: ccpp_const_get_idx
    use musica_ccpp_namelist,          only: filename_of_micm_configuration, &
                                             filename_of_tuvx_configuration, &
                                             filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_load_species, only: index_dry_air, index_O2, index_O3

    integer,         intent(in) :: number_of_columns
    integer,         intent(in) :: number_of_layers
    real(kind_phys), intent(in) :: test_accuracy     ! Relative tolrance for checking results

    integer, parameter                                             :: NUM_SPECIES = 6
    integer                                                        :: NUM_GRID_CELLS

    integer                                                        :: errcode
    character(len=512)                                             :: errmsg
    real(kind_phys)                                                :: time_step = 60._kind_phys ! s
    real(kind_phys), dimension(number_of_columns,number_of_layers) :: temperature               ! K
    real(kind_phys), dimension(number_of_columns,number_of_layers) :: pressure                  ! Pa
    real(kind_phys), dimension(number_of_columns,number_of_layers) :: air_density               ! mol m-3
    real(kind_phys), dimension(number_of_columns,number_of_layers) :: dry_air_mass_density      ! kg m-3
    real(kind_phys), dimension(number_of_columns,number_of_layers,NUM_SPECIES) :: constituents         ! kg kg-1
    real(kind_phys), dimension(number_of_columns,number_of_layers,NUM_SPECIES) :: initial_constituents ! kg kg-1
    type(ccpp_constituent_prop_ptr_t),   allocatable               :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target       :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer                   :: const_prop
    character(len=512)                                             :: species_name, units
    character(len=:), allocatable                                  :: micm_species_name
    integer                                                        :: i, j
    integer                                                        :: A_index, B_index, C_index, D_index, E_index, F_index
    real(kind_phys)                                                :: A_MW, B_MW, C_MW, D_MW, E_MW, F_MW
    type(ArrheniusReaction)                                        :: arr1, arr2, arr3, arr4
    real(kind_phys)                                                :: initial_A, initial_C, initial_D, initial_F
    real(kind_phys)                                                :: A, B, C, D, E, F
    real(kind_phys)                                                :: k1, k2, k3, k4
    real(kind_phys)                                                :: dummy_array_1D(1), dummy_array_2D(1,1)

    NUM_GRID_CELLS = number_of_columns * number_of_layers
    dummy_array_1D = -HUGE(0.0_kind_phys)
    dummy_array_2D = -HUGE(0.0_kind_phys)

    filename_of_micm_configuration = 'test/musica/configuration/analytical/config.json'
    filename_of_tuvx_configuration = 'none'
    filename_of_tuvx_micm_mapping_configuration = 'none'

    ! MUSICA registration
    call musica_ccpp_register(constituent_props, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    ASSERT(allocated(constituent_props))
    ASSERT(size(constituent_props) == NUM_SPECIES)
    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set(const_prop, errcode, errmsg)
    end do

    ! Get indices and molar masses for chemical species
    call get_index_and_molar_mass(constituent_props_ptr, "A", A_index, A_MW)
    call get_index_and_molar_mass(constituent_props_ptr, "B", B_index, B_MW)
    call get_index_and_molar_mass(constituent_props_ptr, "C", C_index, C_MW)
    call get_index_and_molar_mass(constituent_props_ptr, "D", D_index, D_MW)
    call get_index_and_molar_mass(constituent_props_ptr, "E", E_index, E_MW)
    call get_index_and_molar_mass(constituent_props_ptr, "F", F_index, F_MW)

    ! Check that the molar masses are correct
    ASSERT_NEAR(A_MW, 0.01802_kind_phys,   1.0e-5_kind_phys)
    ASSERT_NEAR(B_MW, 0.02897_kind_phys,   1.0e-5_kind_phys)
    ASSERT_NEAR(C_MW, 0.0319988_kind_phys, 1.0e-5_kind_phys)
    ASSERT_NEAR(D_MW, 0.0479982_kind_phys, 1.0e-5_kind_phys)
    ASSERT_NEAR(E_MW, 0.07254_kind_phys,   1.0e-5_kind_phys)
    ASSERT_NEAR(F_MW, 0.082356_kind_phys,  1.0e-5_kind_phys)

    ! MUSICA initialization
    call musica_ccpp_init(number_of_columns, number_of_layers, number_of_layers+1, &
                          dummy_array_1D, constituent_props_ptr, MOLAR_MASS_DRY_AIR__G_MOL, &
                          errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    do i = 1, number_of_layers
      do j = 1, number_of_columns
        temperature(j,i) = 250.0_kind_phys + 100.0_kind_phys * (i-1) / number_of_columns + 10.0_kind_phys * (j-1) / number_of_layers
        pressure(j,i) = 10000.0_kind_phys + 100000.0_kind_phys * (i-1) / number_of_columns + 10000.0_kind_phys * (j-1) / number_of_layers
        air_density(j,i) = pressure(j,i) / (GAS_CONSTANT * temperature(j,i))
        dry_air_mass_density(j,i) = air_density(j,i) * MOLAR_MASS_DRY_AIR
        constituents(j,i,A_index) = 0.75 + 0.1_kind_phys * (i-1) / number_of_columns + 0.01_kind_phys * (j-1) / number_of_layers - 0.05
        constituents(j,i,B_index) = 0.0_kind_phys
        constituents(j,i,C_index) = 0.4_kind_phys + 0.1_kind_phys * (i-1) / number_of_columns + 0.01_kind_phys * (j-1) / number_of_layers - 0.05
        constituents(j,i,D_index) = 0.8_kind_phys + 0.1_kind_phys * (i-1) / number_of_columns + 0.01_kind_phys * (j-1) / number_of_layers - 0.05
        constituents(j,i,E_index) = 0.0_kind_phys
        constituents(j,i,F_index) = 0.1_kind_phys + 0.1_kind_phys * (i-1) / number_of_columns + 0.01_kind_phys * (j-1) / number_of_layers - 0.05
      end do
    end do
    initial_constituents(:,:,:) = constituents(:,:,:)

    ! Convert to kg kg-1 for CAM-SIMA
    do i = 1, number_of_columns
      do j = 1, number_of_layers
        constituents(i,j,A_index) = constituents(i,j,A_index) / dry_air_mass_density(i,j) * A_MW
        constituents(i,j,B_index) = constituents(i,j,B_index) / dry_air_mass_density(i,j) * B_MW
        constituents(i,j,C_index) = constituents(i,j,C_index) / dry_air_mass_density(i,j) * C_MW
        constituents(i,j,D_index) = constituents(i,j,D_index) / dry_air_mass_density(i,j) * D_MW
        constituents(i,j,E_index) = constituents(i,j,E_index) / dry_air_mass_density(i,j) * E_MW
        constituents(i,j,F_index) = constituents(i,j,F_index) / dry_air_mass_density(i,j) * F_MW
      end do
    end do

    ! MUSICA run for one time step
    call musica_ccpp_run( time_step, temperature, pressure, dry_air_mass_density, constituent_props_ptr, &
                          constituents, dummy_array_2D, dummy_array_2D, dummy_array_1D, dummy_array_1D, &
                          dummy_array_1D, dummy_array_1D, -HUGE(0.0_kind_phys), dummy_array_2D, &
                          dummy_array_2D, dummy_array_1D, -HUGE(0.0_kind_phys), STDOUT, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    ! Convert back to mol m-3 for analytical check
    do i = 1, number_of_columns
      do j = 1, number_of_layers
        constituents(i,j,A_index) = constituents(i,j,A_index) * dry_air_mass_density(i,j) / A_MW
        constituents(i,j,B_index) = constituents(i,j,B_index) * dry_air_mass_density(i,j) / B_MW
        constituents(i,j,C_index) = constituents(i,j,C_index) * dry_air_mass_density(i,j) / C_MW
        constituents(i,j,D_index) = constituents(i,j,D_index) * dry_air_mass_density(i,j) / D_MW
        constituents(i,j,E_index) = constituents(i,j,E_index) * dry_air_mass_density(i,j) / E_MW
        constituents(i,j,F_index) = constituents(i,j,F_index) * dry_air_mass_density(i,j) / F_MW
      end do
    end do

    ! Check the results
    arr1%A_ = 0.004_kind_phys
    arr1%C_ = 50.0_kind_phys
    arr2%A_ = 0.012_kind_phys
    arr2%B_ = -2.0_kind_phys
    arr2%C_ = 75.0_kind_phys
    arr2%D_ = 50.0_kind_phys
    arr2%E_ = 1.0e-6_kind_phys
    arr3%A_ = 0.001_kind_phys
    arr4%A_ = 0.002_kind_phys

    do i = 1, number_of_columns
      do j = 1, number_of_layers
        initial_A = initial_constituents(i,j,A_index)
        initial_C = initial_constituents(i,j,C_index)
        initial_D = initial_constituents(i,j,D_index)
        initial_F = initial_constituents(i,j,F_index)
        k1 = calculate_arrhenius( arr1, temperature(i,j), pressure(i,j) )
        k2 = calculate_arrhenius( arr2, temperature(i,j), pressure(i,j) )
        k3 = calculate_arrhenius( arr3, temperature(i,j), pressure(i,j) )
        k4 = calculate_arrhenius( arr4, temperature(i,j), pressure(i,j) )
        A = initial_A * exp( -k1 * time_step )
        B = initial_A * (k1 / (k2 - k1)) * &
            (exp( -k1 * time_step ) - exp( -k2 * time_step ))
        C = initial_C + initial_A * &
            (1.0 + (k1 * exp(-k2 * time_step) - k2 * exp(-k1 * time_step)) / (k2 - k1))
        D = initial_D * exp( -k3 * time_step )
        E = initial_D * (k3 / (k4 - k3)) * &
            (exp( -k3 * time_step ) - exp(-k4 * time_step ))
        F = initial_F + initial_D * &
            (1.0 + (k3 * exp(-k4 * time_step) - k4 * exp(-k3 * time_step)) / (k4 - k3))
        
        ! Check that the results are correct
        ASSERT_NEAR(constituents(i,j,A_index), A, test_accuracy)
        ASSERT_NEAR(constituents(i,j,B_index), B, test_accuracy)
        ASSERT_NEAR(constituents(i,j,C_index), C, test_accuracy)
        ASSERT_NEAR(constituents(i,j,D_index), D, test_accuracy)
        ASSERT_NEAR(constituents(i,j,E_index), E, test_accuracy)
        ASSERT_NEAR(constituents(i,j,F_index), F, test_accuracy)
      end do
    end do

    ! Clean up
    call musica_ccpp_final(errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

  end subroutine test_analytical

  !> Check with Rosenbrock Solver
  subroutine test_rosenbrock()
    use musica_ccpp_namelist, only : micm_solver_type

    real(kind_phys) :: test_accuracy = 5.0e-3_kind_phys

    micm_solver_type = "Rosenbrock"
    call test_analytical(1, 1, test_accuracy)
    call test_analytical(2, 5, test_accuracy)
    call test_analytical(5, 2, test_accuracy)
    call test_analytical(128, 128, test_accuracy)
    call test_analytical(160, 200, test_accuracy)
  end subroutine test_rosenbrock

  !> Check with Backward Euler Solver
  subroutine test_backward_euler()
    use musica_ccpp_namelist, only : micm_solver_type

    real(kind_phys) :: test_accuracy = 0.1_kind_phys

    micm_solver_type = "Backward Euler"
    call test_analytical(1, 1, test_accuracy)
    call test_analytical(2, 5, test_accuracy)
    call test_analytical(5, 2, test_accuracy)
    call test_analytical(128, 128, test_accuracy)
    call test_analytical(160, 200, test_accuracy)
  end subroutine test_backward_euler

end program run_test_musica_ccpp
