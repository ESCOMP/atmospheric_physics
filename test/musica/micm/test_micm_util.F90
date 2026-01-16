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
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_util,    only: set_constants
    use musica_micm,         only: micm_t
    use musica_state,        only: state_t
    use musica_util,         only: error_t
    use musica_ccpp_species, only: micm_molar_mass_array

    character(len=*), parameter :: filename_of_micm_configuration = "test/musica/configuration/micm_util/config.json"
    integer, parameter :: NUM_COLUMNS = 2
    integer, parameter :: NUM_LAYERS = 2
    integer, parameter :: NUM_CELLS = NUM_COLUMNS * NUM_LAYERS
    integer, parameter :: NUM_SPECIES = 4
    integer, parameter :: NUM_RATE_PARAMETERS = 2
    real,    parameter :: ABS_ERROR = 1e-3
    real(kind_phys), parameter :: MOLAR_MASS_DRY_AIR = 0.0289644_kind_phys ! kg mol-1
    real(kind_phys), target :: dry_air_mass_density(NUM_COLUMNS,NUM_LAYERS) ! kg m-3
    real(kind_phys), target :: pressure(NUM_COLUMNS,NUM_LAYERS)        ! Pa
    real(kind_phys), target :: temperature(NUM_COLUMNS,NUM_LAYERS)     ! K
    real(kind_phys), target :: constituents(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES)
    real(kind_phys), target :: ccpp_constituents(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) ! kg kg-1
    real(kind_phys), target :: micm_constituents(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) ! mol m-3
    real(kind_phys), target :: rate_parameters(NUM_COLUMNS,NUM_LAYERS,NUM_RATE_PARAMETERS) ! various units
    integer                 :: i_column, i_layer, i_elem, i_cell
    type(micm_t),  pointer  :: micm
    type(state_t), pointer  :: state
    type(error_t)           :: error

    call set_constants(MOLAR_MASS_DRY_AIR)

    micm => micm_t(trim(filename_of_micm_configuration), NUM_CELLS, error)
    ASSERT(error%is_success( ))
    state => micm%get_state(NUM_CELLS, error)

    temperature(:,1) = (/ 200._kind_phys, 210._kind_phys /)
    temperature(:,2) = (/ 220._kind_phys, 230._kind_phys /)
    pressure(:,1)    = (/ 1000._kind_phys, 1100._kind_phys /)
    pressure(:,2)    = (/ 1200._kind_phys, 1300._kind_phys /)
    dry_air_mass_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_mass_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
    constituents(1,1,:) = (/ 0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys, 0.4_kind_phys /)
    constituents(1,2,:) = (/ 0.41_kind_phys, 0.42_kind_phys, 0.43_kind_phys, 0.44_kind_phys /)
    constituents(2,1,:) = (/ 0.21_kind_phys, 0.22_kind_phys, 0.23_kind_phys, 0.24_kind_phys /)
    constituents(2,2,:) = (/ 0.31_kind_phys, 0.32_kind_phys, 0.33_kind_phys, 0.34_kind_phys /)
    rate_parameters(1,1,:) = (/ 0.1_kind_phys, 0.2_kind_phys /)
    rate_parameters(1,2,:) = (/ 0.3_kind_phys, 0.4_kind_phys /)
    rate_parameters(2,1,:) = (/ 0.5_kind_phys, 0.6_kind_phys /)
    rate_parameters(2,2,:) = (/ 0.7_kind_phys, 0.8_kind_phys /)

    ccpp_constituents(:,:,:) = constituents(:,:,:)
    
    ! Set the expected values for the constituents in mol m-3
    micm_constituents(:,:,1) = ccpp_constituents(:,:,1) * dry_air_mass_density(:,:) / ( MOLAR_MASS_DRY_AIR * micm_molar_mass_array(1) )
    micm_constituents(:,:,2) = ccpp_constituents(:,:,2) * dry_air_mass_density(:,:) / ( MOLAR_MASS_DRY_AIR * micm_molar_mass_array(2) )
    micm_constituents(:,:,3) = ccpp_constituents(:,:,3) * dry_air_mass_density(:,:) / ( MOLAR_MASS_DRY_AIR * micm_molar_mass_array(3) )
    micm_constituents(:,:,4) = ccpp_constituents(:,:,4) * dry_air_mass_density(:,:) / ( MOLAR_MASS_DRY_AIR * micm_molar_mass_array(4) )

    call update_micm_state(state, 0, temperature, pressure, dry_air_mass_density, constituents, rate_parameters)

    do i_column = 1, NUM_COLUMNS
      do i_layer = 1, NUM_LAYERS
        i_cell = i_column + (i_layer - 1) * NUM_COLUMNS
        ASSERT_NEAR(state%conditions(i_cell)%temperature, temperature(i_column, i_layer),     ABS_ERROR)
        ASSERT_NEAR(state%conditions(i_cell)%pressure,    pressure(i_column, i_layer),        ABS_ERROR)
        ASSERT_NEAR(state%conditions(i_cell)%air_density, dry_air_mass_density(i_column, i_layer) / MOLAR_MASS_DRY_AIR, ABS_ERROR)
        do i_elem = 1, NUM_SPECIES
          ASSERT_NEAR(state%concentrations(1 + (i_cell-1)*state%species_strides%grid_cell + (i_elem-1)*state%species_strides%variable), micm_constituents(i_column, i_layer, i_elem), ABS_ERROR)
        end do
        do i_elem = 1, NUM_RATE_PARAMETERS
          ASSERT_NEAR(state%rate_parameters(1 + (i_cell-1)*state%rate_parameters_strides%grid_cell + (i_elem-1)*state%rate_parameters_strides%variable), rate_parameters(i_column, i_layer, i_elem), ABS_ERROR)
        end do
      end do
    end do

    constituents(:,:,:) = -HUGE(1.0_kind_phys)
    call extract_mixing_ratios_from_state(state, 0, constituents)

    do i_column = 1, NUM_COLUMNS
      do i_layer = 1, NUM_LAYERS
        do i_elem = 1, NUM_SPECIES
          ASSERT_NEAR(constituents(i_column, i_layer, i_elem), ccpp_constituents(i_column, i_layer, i_elem), ABS_ERROR)
        end do
      end do
    end do

    deallocate(state)
    deallocate(micm)

  end subroutine test_unit_conversion

end program test_micm_util
