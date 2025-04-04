! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_micm_util

  use musica_ccpp_micm_util

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_unit_conversion()

contains

  subroutine test_unit_conversion()
    use ccpp_kinds, only: kind_phys

    integer, parameter :: NUM_COLUMNS = 2
    integer, parameter :: NUM_LAYERS = 2
    integer, parameter :: NUM_SPECIES = NUM_COLUMNS * NUM_LAYERS
    real,    parameter :: ABS_ERROR = 1e-3
    real(kind_phys)    :: dry_air_density(NUM_COLUMNS,NUM_LAYERS)
    real(kind_phys)    :: molar_mass_arr(NUM_SPECIES)
    real(kind_phys)    :: constituents(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES)
    real(kind_phys)    :: ccpp_constituents(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) ! kg kg-1
    real(kind_phys)    :: micm_constituents(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) ! mol m-3
    integer            :: i_column, i_layer, i_elem

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
          ASSERT_NEAR(constituents(i_column, i_layer, i_elem), micm_constituents(i_column, i_layer, i_elem), ABS_ERROR)
        end do
      end do
    end do

    call convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, constituents)

    do i_column = 1, NUM_COLUMNS
      do i_layer = 1, NUM_LAYERS
        do i_elem = 1, NUM_SPECIES
          ASSERT_NEAR(constituents(i_column, i_layer, i_elem), ccpp_constituents(i_column, i_layer, i_elem), ABS_ERROR)
        end do
      end do
    end do

  end subroutine test_unit_conversion

end program test_micm_util