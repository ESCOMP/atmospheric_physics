program run_test_musica_ccpp

  use musica_ccpp

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_musica_ccpp_api()

contains

  subroutine test_musica_ccpp_api()
    use musica_micm,               only: Rosenbrock, RosenbrockStandardOrder
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_micm,          only: micm

    implicit none

    integer, parameter                                             :: NUM_SPECIES = 5
    ! This test requires that the number of grid cells = 4, which is the default
    ! vector dimension for MICM. This restriction will be removed once
    ! https://github.com/NCAR/musica/issues/217 is finished.
    integer, parameter                                             :: NUM_COLUMNS = 2
    integer, parameter                                             :: NUM_LAYERS = 2
    integer                                                        :: solver_type
    integer                                                        :: errcode
    character(len=512)                                             :: errmsg
    real(kind_phys)                                                :: time_step                                    ! s
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: geopotential_height_wrt_surface_at_midpoint  ! m
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS+1)           :: geopotential_height_wrt_surface_at_interface ! m
    real(kind_phys), dimension(NUM_COLUMNS)                        :: surface_temperature                          ! K
    real(kind_phys), dimension(NUM_COLUMNS)                        :: surface_geopotential                         ! m2 s-2
    real(kind_phys)                                                :: standard_gravitational_acceleration          ! s2 m-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: temperature                                  ! K
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: pressure                                     ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: dry_air_density                              ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) :: initial_constituents                         ! kg kg-1
    type(ccpp_constituent_prop_ptr_t), allocatable                 :: constituent_props_ptr(:)

    ! local variables
    type(ccpp_constituent_properties_t), allocatable, target :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer             :: const_prop
    real(kind_phys)                                          :: molar_mass, base_conc
    character(len=512)                                       :: species_name, units
    character(len=:), allocatable                            :: micm_species_name
    logical                                                  :: tmp_bool, is_advected
    integer                                                  :: num_grid_cells
    integer                                                  :: i, j, k
    integer                                                  :: N2_index, O2_index, O_index, O1D_index, O3_index
    real(kind_phys)                                          :: total_O, total_O_init

    solver_type = Rosenbrock
    num_grid_cells = NUM_COLUMNS * NUM_LAYERS
    time_step = 60._kind_phys
    geopotential_height_wrt_surface_at_midpoint(1,:) = (/ 2000.0_kind_phys, 500.0_kind_phys /)
    geopotential_height_wrt_surface_at_midpoint(2,:) = (/ 2000.0_kind_phys, -500.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(1,:) = (/ 3000.0_kind_phys, 1000.0_kind_phys, 0.0_kind_phys /)
    geopotential_height_wrt_surface_at_interface(2,:) = (/ 3000.0_kind_phys, 500.0_kind_phys, -1500.0_kind_phys /)
    surface_temperature = (/ 300.0_kind_phys, 300.0_kind_phys /)
    surface_geopotential = (/ 100.0_kind_phys, 200.0_kind_phys /)
    standard_gravitational_acceleration = 10.0_kind_phys
    temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    pressure(:,1) = (/ 6000.04_kind_phys, 7000.04_kind_phys /)
    pressure(:,2) = (/ 8000.04_kind_phys, 9000.04_kind_phys /)
    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)

    call musica_ccpp_register(solver_type, num_grid_cells, constituent_props, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
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
      tmp_bool = (trim(species_name) == "O2" .and. molar_mass == 0.0319988_kind_phys .and. is_advected) .or.  &
                (trim(species_name) == "O" .and. molar_mass == 0.0159994_kind_phys .and. .not. is_advected) .or.   &
                (trim(species_name) == "O1D" .and. molar_mass == 0.0159994_kind_phys .and. .not. is_advected) .or. &
                (trim(species_name) == "O3" .and. molar_mass == 0.0479982_kind_phys .and. is_advected) .or. &
                (trim(species_name) == "N2" .and. molar_mass == 0.0280134_kind_phys .and. is_advected)
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

    call musica_ccpp_init(NUM_LAYERS, NUM_LAYERS+1, errmsg, errcode)
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
    initial_constituents(:,:,:) = constituents(:,:,:)

    write(*,*) "[MUSICA INFO] Initial Time Step"
    write(*,fmt="(1x,f10.2)") time_step
    write(*,*) "[MUSICA INFO] Initial Temperature"
    write(*,fmt="(4(1x,f10.4))") temperature
    write(*,*) "[MUSICA INFO] Initial Pressure"
    write(*,fmt="(4(1x,f10.4))") pressure
    write(*,*) "[MUSICA INFO] Initial Concentrations"
    write(*,fmt="(4(3x,e13.6))") constituents

    call musica_ccpp_run(time_step, temperature, pressure, dry_air_density, constituent_props_ptr, &
                         constituents, geopotential_height_wrt_surface_at_midpoint,                &
                         geopotential_height_wrt_surface_at_interface, surface_temperature,        &
                         surface_geopotential, standard_gravitational_acceleration, errmsg, errcode)
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
      end do
    end do

  end subroutine test_musica_ccpp_api

end program run_test_musica_ccpp