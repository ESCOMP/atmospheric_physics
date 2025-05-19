! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_micm_util

  use ccpp_kinds,           only: kind_phys

  implicit none
  private

  public :: update_micm_state, extract_mixing_ratios_from_state

contains

  !> Populate a MICM state object with conditions from CCPP variables
  !!
  !! The state object is populated with data from the first grid cell that has not
  !! yet been added to a MICM state. Indices for chemical species are mapped from
  !! the CCPP constituent ordering to the MICM species ordering. Mass mixing ratios
  !! are converted to number density (mol m-3) using the dry air density and
  !! molecular weights of the species.
  subroutine update_micm_state(state, state_data_offset, temperature, pressure, &
                               dry_air_density, mixing_ratios, rate_parameters)

    use musica_ccpp_species, only: micm_indices_constituent_props, micm_molar_mass_array
    use musica_state,        only: state_t

    type(state_t),                       intent(inout) :: state
    integer,                             intent(in)    :: state_data_offset      ! number of grid cells already updated
    real(kind_phys), target, contiguous, intent(in)    :: temperature(:,:)       ! K (column, layer)
    real(kind_phys), target, contiguous, intent(in)    :: pressure(:,:)          ! Pa (column, layer)
    real(kind_phys), target, contiguous, intent(in)    :: dry_air_density(:,:)   ! kg m-3 (column, layer)
    real(kind_phys), target, contiguous, intent(in)    :: mixing_ratios(:,:,:)   ! kg kg-1 (column, layer, species)
    real(kind_phys), target, contiguous, intent(in)    :: rate_parameters(:,:,:) ! various units (column, layer, parameter)

    integer :: i_cell, i_var, state_offset, n_cells, n_cells_total
    real(kind_phys), pointer :: temperature_1D(:), pressure_1D(:), air_density_1D(:), &
                                species_1D(:), params_1D(:)

    ! get grid cell dimensions
    n_cells = state%number_of_grid_cells
    n_cells_total = size(temperature, 1) * size(temperature, 2)

    ! Update environmental conditions
    ! collapse 2D arrays to 1D
    ! (column, layer) -> (column*layer)
    temperature_1D(1:n_cells_total) => temperature(:,:)
    pressure_1D(1:n_cells_total) => pressure(:,:)
    air_density_1D(1:n_cells_total) => dry_air_density(:,:)
    do i_cell = 1, n_cells
      state%conditions(i_cell)%temperature = temperature_1D(i_cell + state_data_offset)
      state%conditions(i_cell)%pressure    = pressure_1D(i_cell + state_data_offset)
      state%conditions(i_cell)%air_density = air_density_1D(i_cell + state_data_offset)
    end do

    ! Update species concentrations
    associate(cell_stride => state%species_strides%grid_cell, &
              var_stride  => state%species_strides%variable)
    do i_var = 1, state%number_of_species
      species_1D(1:n_cells_total) => mixing_ratios(:,:,micm_indices_constituent_props(i_var))
      do i_cell = 1, n_cells
        state%concentrations( 1 + ( i_cell - 1 ) * cell_stride + ( i_var - 1 ) * var_stride ) = &
          species_1D(i_cell + state_data_offset) * state%conditions(i_cell)%air_density &
              / micm_molar_mass_array(i_var)
      end do
    end do
    end associate

    ! Update rate parameters
    associate(cell_stride => state%rate_parameters_strides%grid_cell, &
              var_stride  => state%rate_parameters_strides%variable)
    do i_var = 1, state%number_of_rate_parameters
      params_1D(1:n_cells_total) => rate_parameters(:,:,i_var)
      do i_cell = 1, n_cells
        state%rate_parameters( 1 + ( i_cell - 1 ) * cell_stride + ( i_var - 1 ) * var_stride ) = &
          params_1D(i_cell + state_data_offset)
      end do
    end do
    end associate

  end subroutine update_micm_state

  !> Extract mixing ratios from a MICM state object
  !!
  !! Species concentrations are mapped from the MICM species ordering to the
  !! CCPP constituent ordering. The concentrations are converted to mass mixing
  !! ratios (kg kg-1) using the dry air density and molecular weights of the
  !! species.
  subroutine extract_mixing_ratios_from_state(state, state_data_offset, mixing_ratios)

    use musica_ccpp_species, only: micm_indices_constituent_props, micm_molar_mass_array
    use musica_state,        only: state_t

    type(state_t),                       intent(in)    :: state
    integer,                             intent(in)    :: state_data_offset    ! number of grid cells already updated
    real(kind_phys), target, contiguous, intent(inout) :: mixing_ratios(:,:,:) ! kg kg-1 (column, layer, species)

    integer :: i_cell, i_var, state_offset, n_cells, n_cells_total
    real(kind_phys), pointer :: species_1D(:)

    ! get grid cell dimensions
    n_cells = state%number_of_grid_cells
    n_cells_total = size(mixing_ratios, 1) * size(mixing_ratios, 2)

    ! Update species mass mixing ratios
    associate(cell_stride => state%species_strides%grid_cell, &
              var_stride  => state%species_strides%variable)
    do i_var = 1, state%number_of_species
      species_1D(1:n_cells_total) => mixing_ratios(:,:,micm_indices_constituent_props(i_var))
      do i_cell = 1, n_cells
        species_1D(i_cell + state_data_offset) = &
          state%concentrations( 1 + ( i_cell - 1 ) * cell_stride + ( i_var - 1 ) * var_stride ) &
          * micm_molar_mass_array(i_var) / state%conditions(i_cell)%air_density
      end do
    end do
    end associate

  end subroutine extract_mixing_ratios_from_state

end module musica_ccpp_micm_util
