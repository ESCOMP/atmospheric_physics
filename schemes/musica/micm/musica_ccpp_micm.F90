! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_micm

  ! Note: "micm_t" is included in an external pre-built MICM library that the host
  ! model is responsible for linking to during compilation
  use ccpp_kinds,           only: kind_phys
  use musica_ccpp_util,     only: has_error_occurred
  use musica_ccpp_namelist, only: filename_of_micm_configuration
  use musica_micm,          only: micm_t, solver_stats_t, Rosenbrock
  use musica_state,         only: conditions_t, state_t
  use musica_util,          only: mappings_t

  implicit none
  private
  save

  public :: micm_register, micm_init, micm_run, micm_final

  !> MICM solver. The solver will be configured for a specific chemical mechanism.
  !! It then can be used to create and solve MICM states for the mechanism and a
  !! given number of grid cells.
  type(micm_t),  pointer :: micm => null( )
  !> For optimal performance, the grid cells assigned to any particular MPI rank
  !! are solved in sets of a fixed size specified at build time. If the total number
  !! of grid cells is not evenly divisible by the set size, an additional state
  !! is created to handle the residual grid cells.
  !! If the number of grid cells is less than the optimal set size, only the first
  !! state is created and used.
  type(state_t), pointer :: state_1 => null( ) ! state for the optimal set of grid cells
  type(state_t), pointer :: state_2 => null( ) ! state for the residual grid cells
  integer                :: number_of_grid_cells = 0

  type(mappings_t), public, protected :: species_ordering
  type(mappings_t), public, protected :: rate_parameters_ordering

  integer, parameter :: SOLVER_TYPE_ROSENBROCK = 1
  integer, parameter :: SOLVER_TYPE_BACKWARD_EULER = 3
  integer, parameter :: ONE_GRID_CELL = 1

contains

  !> Registers MICM constituent properties with the CCPP
  subroutine micm_register(solver_type, constituent_props, micm_species, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_species,       only: musica_species_t
    use musica_util,               only: error_t
    use iso_c_binding,             only: c_int

    character(len=*),                                 intent(in)  :: solver_type
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
    type(musica_species_t),              allocatable, intent(out) :: micm_species(:)
    character(len=512),                               intent(out) :: errmsg
    integer,                                          intent(out) :: errcode

    ! local variables
    type(error_t)                 :: error
    real(kind=kind_phys)          :: molar_mass
    character(len=:), allocatable :: species_name
    logical                       :: is_advected
    real(kind=kind_phys)          :: default_value
    integer                       :: number_of_species
    integer                       :: i, species_index, solver_type_int
    type(state_t), pointer        :: state
    character(len=:), allocatable :: error_message
    character(len=*), parameter   :: PROPERTY_NOT_FOUND_SUFFIX = 'Property not found'

    if (associated( micm )) then
      deallocate( micm )
      micm => null()
    end if
    if (trim(solver_type) == 'Rosenbrock') then
      solver_type_int = SOLVER_TYPE_ROSENBROCK
    else if (trim(solver_type) == 'Backward Euler') then
      solver_type_int = SOLVER_TYPE_BACKWARD_EULER
    else
      errmsg = "[MUSICA Error] Invalid solver type. Supported types: 'Rosenbrock', 'Backward Euler'." // &
               " Got: '" // trim(solver_type) // "'."
      errcode = 1
      return
    end if
    micm => micm_t(trim(filename_of_micm_configuration), solver_type_int, error)
    if (has_error_occurred(error, errmsg, errcode)) return
    state => micm%get_state(ONE_GRID_CELL, error)
    if (has_error_occurred(error, errmsg, errcode)) return

    number_of_species = state%species_ordering%size()
    allocate(constituent_props(number_of_species), stat=errcode)
    if (errcode /= 0) then
      errmsg = "[MUSICA Error] Failed to allocate memory for constituent properties."
      return
    end if

    allocate(micm_species(number_of_species), stat=errcode)
    if (errcode /= 0) then
      errmsg = "[MUSICA Error] Failed to allocate memory for micm species."
      return
    end if

    do i = 1, number_of_species
    associate( map => state%species_ordering )
      species_name = map%name(i)
      species_index = map%index(i)

      molar_mass = micm%get_species_property_double(species_name, &
                                                    "molecular weight [kg mol-1]", &
                                                    error)
      if (has_error_occurred(error, errmsg, errcode)) return
      is_advected = micm%get_species_property_bool(species_name, &
                                                   "__is advected", &
                                                   error)
      if (has_error_occurred(error, errmsg, errcode)) return
      default_value = micm%get_species_property_double(species_name, &
                                                       "__default mixing ratio [kg kg-1]", &
                                                       error)
      if (.not. error%is_success( )) then
        error_message = error%message( )
        if (len(error_message) >= len(PROPERTY_NOT_FOUND_SUFFIX)) then
          if (error_message(len(error_message)-len(PROPERTY_NOT_FOUND_SUFFIX)+1:len(error_message)) &
              == PROPERTY_NOT_FOUND_SUFFIX) then
            ! If the default mixing ratio is not defined, use zero
            default_value = 0.0_kind_phys
          else
            if (has_error_occurred(error, errmsg, errcode)) return
          end if
        else
          if (has_error_occurred(error, errmsg, errcode)) return
        end if
      end if

      call constituent_props(species_index)%instantiate( &
        std_name = species_name, &
        long_name = species_name, &
        units = 'kg kg-1', &
        vertical_dim = 'vertical_layer_dimension', &
        default_value = default_value, &
        min_value = 0.0_kind_phys, &
        molar_mass = molar_mass, &
        advected = is_advected, &
        errcode = errcode, &
        errmsg = errmsg )
      if (errcode /= 0) return

      ! Species are ordered to match the sequence of the MICM state array
      micm_species(species_index) = musica_species_t( &
        name = species_name, &
        unit = 'kg kg-1', &
        molar_mass = molar_mass, &
        index_musica_species = species_index )
    end associate ! map
    end do
    species_ordering = state%species_ordering
    rate_parameters_ordering = state%rate_parameters_ordering
    deallocate( state )

  end subroutine micm_register

  !> Initializes MICM
  subroutine micm_init(n_grid_cells, errmsg, errcode)
    use musica_util, only: error_t

    integer,            intent(in)  :: n_grid_cells
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    integer :: max_grid_cells, size_1, size_2
    type(error_t) :: error

    errmsg = ''
    errcode = 0

    if (.not. associated( micm )) then
      errmsg = "[MUSICA Error] MICM not registered. Call micm_register first."
      errcode = 1
      return
    end if
    if (n_grid_cells <= 0) then
      errmsg = "[MUSICA Error] Invalid number of grid cells."
      errcode = 1
      return
    end if
    number_of_grid_cells = n_grid_cells
    max_grid_cells = micm%get_maximum_number_of_grid_cells( )
    size_1 = min( number_of_grid_cells, max_grid_cells )
    size_2 = mod( number_of_grid_cells - size_1, max_grid_cells )
    state_1 => micm%get_state( size_1, error )
    if (has_error_occurred(error, errmsg, errcode)) return
    if (size_2 > 0) then
      state_2 => micm%get_state( size_2, error )
      if (has_error_occurred(error, errmsg, errcode)) return
    end if

  end subroutine micm_init

  !> Solves chemistry at the current time step
  subroutine micm_run(time_step, temperature, pressure, dry_air_density, &
                      rate_parameters, mixing_ratios, log_output_unit, errmsg, errcode)
    use musica_ccpp_micm_util, only: update_micm_state, extract_mixing_ratios_from_state
    use musica_micm,           only: solver_stats_t
    use musica_util,           only: string_t, error_t
    use iso_c_binding,         only: c_double, c_loc

    real(kind_phys),         intent(in)    :: time_step              ! s
    real(kind_phys), target, intent(in)    :: temperature(:,:)       ! K
    real(kind_phys), target, intent(in)    :: pressure(:,:)          ! Pa
    real(kind_phys), target, intent(in)    :: dry_air_density(:,:)   ! kg m-3
    real(kind_phys), target, intent(in)    :: rate_parameters(:,:,:) ! various units
    real(kind_phys), target, intent(inout) :: mixing_ratios(:,:,:)   ! kg kg-1
    integer,                 intent(in)    :: log_output_unit        ! file unit for logging output
    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errcode

    ! local variables
    integer                :: max_cells, i_state, state_size, state_1_size, offset
    type(state_t), pointer :: state
    type(string_t)         :: solver_state
    type(solver_stats_t)   :: solver_stats
    type(error_t)          :: error

    state_1_size = state_1%number_of_grid_cells
    do i_state = 1, ceiling( real( number_of_grid_cells ) / state_1_size )

      ! Determine which state to use for the current iteration
      state_size = min( number_of_grid_cells - ( i_state - 1 ) * state_1_size, &
                        state_1_size )
      if ( state_size == state_1_size ) then
        state => state_1 ! use the main state for the optimal number of grid cells
      else
        state => state_2 ! use the residual state for the remaining grid cells
        if (.not. associated( state )) then
          errmsg = "[MUSICA Error] Internal error. MICM residual state not initialized."
          errcode = 1
          return
        end if
        if (state%number_of_grid_cells /= state_size) then
          errmsg = "[MUSICA Error] Internal error. MICM residual state size mismatch."
          errcode = 1
          return
        end if
      end if
      offset = ( i_state - 1 ) * state_1_size ! number of grid cells already updated

      ! Update MICM state with the current conditions and mixing ratios
      call update_micm_state( state, offset, temperature, pressure, dry_air_density, &
          mixing_ratios, rate_parameters )
      
      ! Solve the system
      call micm%solve( time_step, state, solver_state, solver_stats, error )
      if (has_error_occurred(error, errmsg, errcode)) return
      if (solver_state%get_char_array() /= "Converged") then
        write(log_output_unit,*) &
          "[MUSICA Warning] MICM solver failure: '" // &
          trim(solver_state%get_char_array()) // "'. For grid cells ", &
          (i_state - 1) * state_1_size + 1, " to ", i_state * state_1_size, &
          " of ", number_of_grid_cells
      end if

      ! Update the mixing ratios with the results
      call extract_mixing_ratios_from_state( state, offset, mixing_ratios)

    end do

  end subroutine micm_run

  !> Finalizes MICM
  subroutine micm_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

    if (associated( state_1 )) then
      deallocate( state_1 )
      state_1 => null()
    end if
    if (associated( state_2 )) then
      deallocate( state_2 )
      state_2 => null()
    end if
    if (associated( micm )) then
      deallocate( micm )
      micm => null()
    end if

  end subroutine micm_final

end module musica_ccpp_micm
