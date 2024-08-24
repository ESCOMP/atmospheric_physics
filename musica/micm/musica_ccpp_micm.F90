module musica_ccpp_micm
  use iso_c_binding

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

  !> Register MICM constituents with the CCPP
  subroutine micm_register(constituents, solver_type, num_grid_cells, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_micm, only: Rosenbrock, RosenbrockStandardOrder
    use musica_util, only: error_t

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituents(:)
    integer(c_int),                                   intent(in)  :: solver_type
    integer(c_int),                                   intent(in)  :: num_grid_cells
    character(len=512),                               intent(out) :: errmsg
    integer,                                          intent(out) :: errcode

    ! local variables
    type(error_t)        :: error
    real(kind=kind_phys) :: molar_mass
    logical              :: is_advected
    integer              :: i

    errcode = 0
    errmsg = ''

    micm => micm_t(filename_of_micm_configuration, solver_type, num_grid_cells, error)
    if (has_error_occurred(error, errmsg, errcode)) return

    allocate(constituents(size(micm%species_ordering)), stat=errcode)
    if (errcode /= 0) then
      errmsg = "[MUSICA Error] Failed to allocate memory for constituents."
      return
    end if

    do i = 1, size(micm%species_ordering)
    associate( map => micm%species_ordering(i) )
      molar_mass = micm%get_species_property_double(map%name(), &
                                                    "molecular weight [kg mol-1]", &
                                                    error)
      if (has_error_occurred(error, errmsg, errcode)) return
      is_advected = micm%get_species_property_bool(map%name(), &
                                                      "__is advected", &
                                                      error)
      if (has_error_occurred(error, errmsg, errcode)) return

      call constituents(map%index())%instantiate( &
        std_name = map%name(), &
        long_name = map%name(), &
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
    integer, intent(out)            :: errcode

    errcode = 0
    errmsg = ''

  end subroutine micm_init

  !> Solve chemistry at the current time step
  subroutine micm_run(time_step, temperature, pressure, dry_air_density, constituents, &
                      rate_params, errmsg, errcode)
    use musica_micm, only: solver_stats_t
    use musica_util, only: string_t, error_t

    real(kind_phys),        intent(in)    :: time_step          ! s
    real(c_double), target, intent(in)    :: temperature(:)     ! K
    real(c_double), target, intent(in)    :: pressure(:)        ! Pa
    real(c_double), target, intent(in)    :: dry_air_density(:) ! kg m-3
    real(c_double), target, intent(inout) :: constituents(:)    ! mol m-3
    real(c_double), target, intent(inout) :: rate_params(:)
    character(len=512),     intent(out)   :: errmsg
    integer,                intent(out)   :: errcode

    ! local variables
    type(string_t)       :: solver_state
    type(solver_stats_t) :: solver_stats
    type(error_t)        :: error
    real(c_double)       :: c_time_step
    integer              :: i_elem
    logical              :: debug = .true.

    errcode = 0
    errmsg = ''
    c_time_step = real(time_step, c_double) 

    call micm%solve(c_time_step,     &
                    temperature,     &
                    pressure,        &
                    dry_air_density, &
                    constituents,    &
                    rate_params,     &
                    solver_state,    &
                    solver_stats,    &
                    error)
    if (has_error_occurred(error, errmsg, errcode)) return

    if (debug) then
      write(*,*) "[MUSICA DEBUG] Solver state: ", solver_state%get_char_array()
      write(*,*) "[MUSICA DEBUG] Function calls: ", solver_stats%function_calls()
      write(*,*) "[MUSICA DEBUG] Jacobian updates: ", solver_stats%jacobian_updates()
      write(*,*) "[MUSICA DEBUG] Number of steps: ", solver_stats%number_of_steps()
      write(*,*) "[MUSICA DEBUG] Accepted: ", solver_stats%accepted()
      write(*,*) "[MUSICA DEBUG] Rejected: ", solver_stats%rejected()
      write(*,*) "[MUSICA DEBUG] Decompositions: ", solver_stats%decompositions()
      write(*,*) "[MUSICA DEBUG] Solves: ", solver_stats%solves()
      write(*,*) "[MUSICA DEBUG] Singular: ", solver_stats%singular()
      write(*,*) "[MUSICA DEBUG] Final time: ", solver_stats%final_time()
    end if

  end subroutine micm_run

  !> Finalize MICM
  subroutine micm_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errcode

    errcode = 0
    errmsg = ''

  end subroutine micm_final

end module musica_ccpp_micm
