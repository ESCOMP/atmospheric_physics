program test_micm_util

  use iso_c_binding
  use micm_util
  use ccpp_kinds,   only: kind_phys

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( abs(a - b) > abs_error) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif

  call test_reshape()
  call test_unit_conversion()


contains

  subroutine test_reshape()
    use iso_c_binding, only: c_double
    use ccpp_kinds,    only: kind_phys

    integer, parameter      :: NUM_SPECIES = 4
    integer, parameter      :: NUM_RATES = 3
    integer, parameter      :: NUM_COLUMNS = 2
    integer, parameter      :: NUM_LAYERS = 2
    integer, parameter      :: NUM_GRID_CELLS = 4
    real(kind_phys), target :: temperature(NUM_COLUMNS,NUM_LAYERS)
    real(kind_phys), target :: pressure(NUM_COLUMNS,NUM_LAYERS)
    real(kind_phys), target :: dry_air_density(NUM_COLUMNS,NUM_LAYERS)
    real(kind_phys), target :: constituents(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES)
    real(kind_phys), target :: rate_params(NUM_COLUMNS,NUM_LAYERS,NUM_RATES)
    real(c_double),  target :: m_temperature(NUM_GRID_CELLS)
    real(c_double),  target :: m_pressure(NUM_GRID_CELLS)
    real(c_double),  target :: m_dry_air_density(NUM_GRID_CELLS)
    real(c_double),  target :: m_constituents(NUM_GRID_CELLS*NUM_SPECIES)
    real(c_double),  target :: m_rate_params(NUM_GRID_CELLS*NUM_RATES)

    ! local variables
    real(c_double), dimension(NUM_GRID_CELLS)             :: arr_conditions
    real(c_double), dimension(NUM_GRID_CELLS*NUM_SPECIES) :: arr_constituents
    real(c_double), dimension(NUM_GRID_CELLS*NUM_RATES)   :: arr_rates
    integer                                               :: i_column, i_layer, i_elem, i_arr
    real                                                  :: abs_error = 1e-7

    temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    pressure(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    pressure(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    dry_air_density(:,1) = (/ 100._kind_phys, 200._kind_phys /)
    dry_air_density(:,2) = (/ 300._kind_phys, 400._kind_phys /)
    constituents(1,1,:) = (/ 0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys, 0.4_kind_phys /)
    constituents(1,2,:) = (/ 0.41_kind_phys, 0.42_kind_phys, 0.43_kind_phys, 0.44_kind_phys /)
    constituents(2,1,:) = (/ 0.21_kind_phys, 0.22_kind_phys, 0.23_kind_phys, 0.24_kind_phys /)
    constituents(2,2,:) = (/ 0.31_kind_phys, 0.32_kind_phys, 0.33_kind_phys, 0.34_kind_phys /)
    rate_params(1,1,:) = (/ 100._kind_phys, 200._kind_phys, 300._kind_phys /)
    rate_params(1,2,:) = (/ 400._kind_phys, 500._kind_phys, 600._kind_phys /)
    rate_params(2,1,:) = (/ 700._kind_phys, 800._kind_phys, 900._kind_phys /)
    rate_params(2,2,:) = (/ 1000._kind_phys, 1100._kind_phys, 1200._kind_phys /)

    arr_conditions = (/ 100.0, 200.0, 300.0, 400.0 /)
    arr_constituents = (/ 0.1, 0.2, 0.3, 0.4, 0.21, 0.22, 0.23, 0.24, 0.41, 0.42, 0.43, 0.44, 0.31, 0.32, 0.33, 0.34 /)
    arr_rates = (/ 100., 200., 300., 700., 800., 900., 400., 500., 600., 1000., 1100., 1200. /)
    
    call reshape_into_micm_arr(temperature, pressure, dry_air_density, constituents, &
      rate_params, m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)

    do i_elem = 1, NUM_GRID_CELLS
      ASSERT(m_temperature(i_elem) == arr_conditions(i_elem))
      ASSERT(m_pressure(i_elem) == arr_conditions(i_elem))
      ASSERT(m_dry_air_density(i_elem) == arr_conditions(i_elem))
    end do

    do i_elem = 1, size(m_constituents)
      ASSERT_NEAR(m_constituents(i_elem), arr_constituents(i_elem), abs_error)
    end do

    do i_elem = 1, size(m_rate_params)
      ASSERT_NEAR(m_rate_params(i_elem), arr_rates(i_elem), abs_error)
    end do

    call reshape_into_ccpp_arr(temperature, pressure, dry_air_density, constituents, &
      rate_params, m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)

    i_elem = 1
    do i_layer = 1, NUM_LAYERS
      do i_column = 1, NUM_COLUMNS
        ASSERT(temperature(i_column, i_layer) == arr_conditions(i_elem))
        ASSERT(pressure(i_column, i_layer) == arr_conditions(i_elem))
        ASSERT(dry_air_density(i_column, i_layer) == arr_conditions(i_elem))
        i_elem = i_elem + 1
      end do
    end do

    i_arr = 1
    do i_layer = 1, NUM_LAYERS
      do i_column = 1, NUM_COLUMNS
        do i_elem = 1, NUM_SPECIES
          ASSERT_NEAR(constituents(i_column, i_layer, i_elem), arr_constituents(i_arr), abs_error)
          i_arr = i_arr + 1
        end do
      end do
    end do

    i_arr = 1
    do i_layer = 1, NUM_LAYERS
      do i_column = 1, NUM_COLUMNS
        do i_elem = 1, NUM_RATES
          ASSERT_NEAR(rate_params(i_column, i_layer, i_elem), arr_rates(i_arr), abs_error)
          i_arr = i_arr + 1
        end do
      end do
    end do

  end subroutine test_reshape

  subroutine test_unit_conversion()
    use iso_c_binding, only: c_double

    integer, parameter :: NUM_GRID_CELLS = 4
    integer, parameter :: NUM_SPECIES = 5
    real(c_double)     :: dry_air_density(NUM_GRID_CELLS)     ! kg m-3
    real(c_double)     :: molar_mass_arr(NUM_SPECIES)         ! kg mol-1
    real(c_double)     :: constituents(NUM_SPECIES)           ! kg kg-1
    real(c_double)     :: converted(NUM_SPECIES)              ! mol m-3

    dry_air_density(:) = (/ 3.5, 4.5, 5.5, 6.5 /)
    molar_mass_arr = (/ 200., 200., 200., 200., 200 /)
    constituents(:) = (/ 0.1, 0.2, 0.3, 0.2, 0.1 /)

    converted_constituents(:)
    call convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)

  end subroutine test_unit_conversion






  subroutine convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)
    use iso_c_binding, only: c_double

    integer, parameter      :: NUM_GRID_CELLS = 4
    real(c_double), intent(in)    :: dry_air_density(NUM_GRID_CELLS) ! kg m-3
    real(c_double), intent(in)    :: molar_mass_arr(NUM_GRID_CELLS)  ! kg mol-1
    real(c_double), intent(inout) :: constituents(:)    ! in: kg kg-1 | out: mol m-3


    integer        :: i_elem
    real(c_double) :: val

    dry_air_density(:,1) = (/ 3.5, 4.5, 5.5, 6.5 /)
    
  ! Convert MICM unit to CAM-SIMA unit (mol m-3  ->  kg kg-1)
  ! call convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, constituents)

  ! ! Convert CAM-SIMA unit to MICM unit (kg kg-1  ->  mol m-3)
  ! call convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)


end program test_micm_util