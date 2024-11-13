program test_micm_util

  use musica_ccpp_micm_util

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_reshape()
  call test_unit_conversion()

contains

  subroutine test_reshape()
    use iso_c_binding, only: c_double
    use ccpp_kinds,    only: kind_phys

    integer, parameter :: NUM_SPECIES = 4
    integer, parameter :: NUM_COLUMNS = 2
    integer, parameter :: NUM_LAYERS = 2
    integer, parameter :: NUM_GRID_CELLS = 4
    real(kind_phys)    :: temperature(NUM_COLUMNS,NUM_LAYERS)
    real(kind_phys)    :: pressure(NUM_COLUMNS,NUM_LAYERS)
    real(kind_phys)    :: dry_air_density(NUM_COLUMNS,NUM_LAYERS)
    real(kind_phys)    :: constituents(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES)
    real(c_double)     :: micm_temperature(NUM_GRID_CELLS)
    real(c_double)     :: micm_pressure(NUM_GRID_CELLS)
    real(c_double)     :: micm_dry_air_density(NUM_GRID_CELLS)
    real(c_double)     :: micm_constituents(NUM_GRID_CELLS*NUM_SPECIES)

    ! local variables
    real(c_double), dimension(NUM_GRID_CELLS)             :: arr_conditions
    real(c_double), dimension(NUM_GRID_CELLS*NUM_SPECIES) :: arr_constituents
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
    arr_conditions = (/ 100.0, 200.0, 300.0, 400.0 /)
    arr_constituents = (/ 0.1, 0.2, 0.3, 0.4, 0.21, 0.22, 0.23, 0.24, 0.41, 0.42, 0.43, 0.44, 0.31, 0.32, 0.33, 0.34 /)
    
    call reshape_into_micm_arr(temperature, pressure, dry_air_density, constituents, &
      micm_temperature, micm_pressure, micm_dry_air_density, micm_constituents)

    do i_elem = 1, NUM_GRID_CELLS
      ASSERT(micm_temperature(i_elem) == arr_conditions(i_elem))
      ASSERT(micm_pressure(i_elem) == arr_conditions(i_elem))
      ASSERT(micm_dry_air_density(i_elem) == arr_conditions(i_elem))
    end do

    do i_elem = 1, size(micm_constituents)
      ASSERT_NEAR(micm_constituents(i_elem), arr_constituents(i_elem), abs_error)
    end do

    call reshape_into_ccpp_arr(micm_constituents, constituents)

    i_arr = 1
    do i_layer = 1, NUM_LAYERS
      do i_column = 1, NUM_COLUMNS
        do i_elem = 1, NUM_SPECIES
          ASSERT_NEAR(constituents(i_column, i_layer, i_elem), arr_constituents(i_arr), abs_error)
          i_arr = i_arr + 1
        end do
      end do
    end do

  end subroutine test_reshape

  subroutine test_unit_conversion()
    use ccpp_kinds, only: kind_phys

    integer, parameter                                             :: NUM_COLUMNS = 2
    integer, parameter                                             :: NUM_LAYERS = 2
    integer, parameter                                             :: NUM_SPECIES = 4
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: dry_air_density   ! kg m-3
    real(kind_phys), dimension(NUM_SPECIES)                        :: molar_mass_arr
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) :: constituents
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) :: ccpp_constituents ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) :: micm_constituents ! mol m-3
    integer                                                        :: i_column, i_layer, i_elem
    real                                                           :: abs_error = 1e-3

    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
    molar_mass_arr(:) = (/ 200._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys /)
    constituents(1,1,:) = (/ 0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys, 0.4_kind_phys /)
    constituents(1,2,:) = (/ 0.41_kind_phys, 0.42_kind_phys, 0.43_kind_phys, 0.44_kind_phys /)
    constituents(2,1,:) = (/ 0.21_kind_phys, 0.22_kind_phys, 0.23_kind_phys, 0.24_kind_phys /)
    constituents(2,2,:) = (/ 0.31_kind_phys, 0.32_kind_phys, 0.33_kind_phys, 0.34_kind_phys /)

    ccpp_constituents(1,1,:) = (/ 0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys, 0.4_kind_phys /)
    ccpp_constituents(1,2,:) = (/ 0.41_kind_phys, 0.42_kind_phys, 0.43_kind_phys, 0.44_kind_phys /)
    ccpp_constituents(2,1,:) = (/ 0.21_kind_phys, 0.22_kind_phys, 0.23_kind_phys, 0.24_kind_phys /)
    ccpp_constituents(2,2,:) = (/ 0.31_kind_phys, 0.32_kind_phys, 0.33_kind_phys, 0.34_kind_phys /)

    micm_constituents(1,1,:) = (/ 1.750E-003_kind_phys,   3.500E-003_kind_phys,   5.250E-003_kind_phys,   7.000E-003_kind_phys /)
    micm_constituents(1,2,:) = (/ 1.127E-002_kind_phys,   1.155E-002_kind_phys,   1.182E-002_kind_phys,   1.210E-002_kind_phys /)
    micm_constituents(2,1,:) = (/ 4.725E-003_kind_phys,   4.949E-003_kind_phys,   5.175E-003_kind_phys,   5.400E-003_kind_phys /)
    micm_constituents(2,2,:) = (/ 1.007E-002_kind_phys,   1.040E-002_kind_phys,   1.072E-002_kind_phys,   1.104E-002_kind_phys /)

    call convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)

    do i_column = 1, NUM_COLUMNS
      do i_layer = 1, NUM_LAYERS
        do i_elem = 1, NUM_SPECIES
          ASSERT_NEAR(constituents(i_column, i_layer, i_elem), micm_constituents(i_column, i_layer, i_elem), abs_error)
        end do
      end do
    end do

    call convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, constituents)

    do i_column = 1, NUM_COLUMNS
      do i_layer = 1, NUM_LAYERS
        do i_elem = 1, NUM_SPECIES
          ASSERT_NEAR(constituents(i_column, i_layer, i_elem), ccpp_constituents(i_column, i_layer, i_elem), abs_error)
        end do
      end do
    end do
  end subroutine test_unit_conversion

end program test_micm_util