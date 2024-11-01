program run_test_musica_ccpp

  use musica_ccpp

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif

  call test_musica_ccpp_api()

contains

  subroutine test_musica_ccpp_api()
    use musica_micm,               only: Rosenbrock, RosenbrockStandardOrder
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

    integer, parameter                                             :: NUM_SPECIES = 4
    integer, parameter                                             :: NUM_COLUMNS = 2
    integer, parameter                                             :: NUM_LAYERS = 2
    integer, parameter                                             :: NUM_WAVELENGTH_BINS = 102
    integer                                                        :: num_grid_cells = NUM_COLUMNS * NUM_LAYERS
    integer                                                        :: solver_type = Rosenbrock
    integer                                                        :: errcode
    character(len=512)                                             :: errmsg
    real(kind_phys)                                                :: time_step = 60._kind_phys                    ! s
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS+1)              :: photolysis_wavelength_grid_interfaces        ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: geopotential_height_wrt_surface_at_midpoint  ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS+1)           :: geopotential_height_wrt_surface_at_interface ! m
    real(kind_phys), dimension(NUM_COLUMNS)                        :: surface_geopotential                         ! m2 s-2
    real(kind_phys), dimension(NUM_COLUMNS)                        :: surface_temperature                          ! K
    real(kind_phys)                                                :: surface_albedo                               ! unitless
    real(kind_phys)                                                :: standard_gravitational_acceleration          ! s2 m-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: temperature                                  ! K
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: pressure                                     ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: dry_air_density                              ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) :: constituents                                 ! kg kg-1
    type(ccpp_constituent_prop_ptr_t), allocatable                 :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target       :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer                   :: const_prop
    real(kind_phys)                                                :: molar_mass
    character(len=512)                                             :: species_name, units
    logical                                                        :: tmp_bool, is_advected
    integer                                                        :: i

    ! These are the values that will be used in CAM-SIMA and correspond to the wavelength
    ! bins used in the CAM-Chem photolysis rate constant lookup table.
    !
    ! We're using the actual values here because several of the TS1/TSMLT photolysis
    ! rate constant configurations are sensitive to the wavelength grid.
    photolysis_wavelength_grid_interfaces = (/ &
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
    geopotential_height_wrt_surface_at_midpoint(1,:) = (/ 2000.0_kind_phys, 500.0_kind_phys /)
    geopotential_height_wrt_surface_at_midpoint(2,:) = (/ 2000.0_kind_phys, -500.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(1,:) = (/ 3000.0_kind_phys, 1000.0_kind_phys, 0.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(2,:) = (/ 3000.0_kind_phys, 500.0_kind_phys, -1500.0_kind_phys /)
    surface_geopotential = (/ 100.0_kind_phys, 200.0_kind_phys /)
    surface_temperature = (/ 300.0_kind_phys, 300.0_kind_phys /)
    surface_albedo = 0.10_kind_phys
    standard_gravitational_acceleration = 10.0_kind_phys
    temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    pressure(:,1) = (/ 6000.04_kind_phys, 7000.04_kind_phys /)
    pressure(:,2) = (/ 8000.04_kind_phys, 9000.04_kind_phys /)
    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
    constituents(1,1,:) = (/ 0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys, 0.4_kind_phys /)
    constituents(1,2,:) = (/ 0.41_kind_phys, 0.42_kind_phys, 0.43_kind_phys, 0.44_kind_phys /)
    constituents(2,1,:) = (/ 0.21_kind_phys, 0.22_kind_phys, 0.23_kind_phys, 0.24_kind_phys /)
    constituents(2,2,:) = (/ 0.31_kind_phys, 0.32_kind_phys, 0.33_kind_phys, 0.34_kind_phys /)

    call musica_ccpp_register(solver_type, num_grid_cells, constituent_props, errmsg, errcode)
    ASSERT(allocated(constituent_props))
    ASSERT(size(constituent_props) == NUM_SPECIES)
    do i = 1, size(constituent_props)
      ASSERT(constituent_props(i)%is_instantiated(errcode, errmsg))
      ASSERT(errcode == 0)
      call constituent_props(i)%standard_name(species_name, errcode, errmsg)
      ASSERT(errcode == 0)
      call constituent_props(i)%molar_mass(molar_mass, errcode, errmsg)
      ASSERT(errcode == 0)
      call constituent_props(i)%is_advected(is_advected, errcode, errmsg)
      ASSERT(errcode == 0)
      tmp_bool = (trim(species_name) == "O2" .and. molar_mass == 0.032_kind_phys .and. .not. is_advected) .or.  &
                (trim(species_name) == "O" .and. molar_mass == 0.016_kind_phys .and. .not. is_advected) .or.   &
                (trim(species_name) == "O1D" .and. molar_mass == 0.016_kind_phys .and. .not. is_advected) .or. &
                (trim(species_name) == "O3" .and. molar_mass == 0.048_kind_phys .and. is_advected)
      ASSERT(tmp_bool)
      call constituent_props(i)%units(units, errcode, errmsg)
      ASSERT(errcode == 0)
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

    call musica_ccpp_init(NUM_LAYERS, NUM_LAYERS+1, photolysis_wavelength_grid_interfaces, &
                          errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    write(*,*) "[MUSICA INFO] Initial Time Step"
    write(*,fmt="(1x,f10.2)") time_step
    write(*,*) "[MUSICA INFO] Initial Temperature"
    write(*,fmt="(4(1x,f10.4))") temperature
    write(*,*) "[MUSICA INFO] Initial Pressure"
    write(*,fmt="(4(1x,f10.4))") pressure
    write(*,*) "[MUSICA INFO] Initial Concentrations"
    write(*,fmt="(4(3x,e13.6))") constituents

    call musica_ccpp_run(time_step, temperature, pressure, dry_air_density, constituent_props_ptr,  &
                         constituents, geopotential_height_wrt_surface_at_midpoint,                 &
                         geopotential_height_wrt_surface_at_interface, surface_geopotential,        &
                         surface_temperature, surface_albedo, standard_gravitational_acceleration,  &
                         errmsg, errcode)
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

    deallocate(constituent_props_ptr)

  end subroutine test_musica_ccpp_api

end program run_test_musica_ccpp