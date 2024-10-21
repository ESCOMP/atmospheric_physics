module musica_ccpp_micm

  ! Note: "micm_t" is included in an external pre-built MICM library that the host
  ! model is responsible for linking to during compilation
  use musica_micm,          only: micm_t
  use musica_ccpp_util,     only: has_error_occurred
  use musica_ccpp_namelist, only: filename_of_micm_configuration
  use ccpp_kinds,           only: kind_phys

  implicit none
  private

  public :: micm_register, micm_init, micm_run, micm_final

  type(micm_t), pointer  :: micm => null( )

contains

  !> Register MICM constituent properties with the CCPP
  subroutine micm_register(solver_type, num_grid_cells, constituent_props, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_micm,               only: Rosenbrock, RosenbrockStandardOrder
    use musica_util,               only: error_t
    use iso_c_binding,             only: c_int

    integer(c_int),                                   intent(in)  :: solver_type
    integer(c_int),                                   intent(in)  :: num_grid_cells
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
    character(len=512),                               intent(out) :: errmsg
    integer,                                          intent(out) :: errcode

    ! local variables
    type(error_t)                 :: error
    real(kind=kind_phys)          :: molar_mass
    character(len=:), allocatable :: species_name
    logical                       :: is_advected
    integer                       :: i, species_index

    micm => micm_t(filename_of_micm_configuration, solver_type, num_grid_cells, error)
    if (has_error_occurred(error, errmsg, errcode)) return

    allocate(constituent_props(micm%species_ordering%size()), stat=errcode)
    if (errcode /= 0) then
      errmsg = "[MUSICA Error] Failed to allocate memory for constituent properties."
      return
    end if

    do i = 1, micm%species_ordering%size()
    associate( map => micm%species_ordering )
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

      call constituent_props(species_index)%instantiate( &
        std_name = species_name, &
        long_name = species_name, &
        units = 'kg kg-1', &
        vertical_dim = 'vertical_layer_dimension', &
        default_value = 0.0_kind_phys, &
        min_value = 0.0_kind_phys, &
        molar_mass = molar_mass, &
        advected = is_advected, &
        errcode = errcode, &
        errmsg = errmsg)
      if (errcode /= 0) return
    end associate ! map
    end do

  end subroutine micm_register

  !> Intitialize MICM
  subroutine micm_init(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

  end subroutine micm_init

  !> Solve chemistry at the current time step
  subroutine micm_run(time_step, temperature, pressure, dry_air_density, &
              user_defined_rate_parameters, constituents, errmsg, errcode)
    use musica_micm,   only: solver_stats_t
    use musica_util,   only: string_t, error_t
    use iso_c_binding, only: c_double

    real(kind_phys),    intent(in)    :: time_step          ! s
    real(c_double),     intent(in)    :: temperature(:)     ! K
    real(c_double),     intent(in)    :: pressure(:)        ! Pa
    real(c_double),     intent(in)    :: dry_air_density(:) ! kg m-3
    real(c_double),     intent(in)    :: user_defined_rate_parameters(:) ! various units
    real(c_double),     intent(inout) :: constituents(:)    ! mol m-3
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errcode

    ! local variables
    type(string_t)       :: solver_state
    type(solver_stats_t) :: solver_stats
    type(error_t)        :: error
    real(c_double)       :: c_time_step
    integer              :: i_elem

    c_time_step = real(time_step, c_double) 

    call micm%solve(c_time_step,                  &
                    temperature,                  &
                    pressure,                     &
                    dry_air_density,              &
                    constituents,                 &
                    user_defined_rate_parameters, &
                    solver_state,                 &
                    solver_stats,                 &
                    error)
    if (has_error_occurred(error, errmsg, errcode)) return

  end subroutine micm_run

  !> Finalize MICM
  subroutine micm_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

    deallocate( micm )

  end subroutine micm_final

end module musica_ccpp_micm
