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
  subroutine micm_run(time_step, temperature, pressure, dry_air_density, constituent_props, &
                      constituents, rate_params, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use musica_micm, only: solver_stats_t
    use musica_util, only: string_t, error_t

    real(kind_phys),                   intent(in)    :: time_step             ! s
    real(c_double), target,           intent(in)    :: temperature(:)      ! K
    real(c_double), target,           intent(in)    :: pressure(:)         ! Pa
    real(c_double), target,           intent(in)    :: dry_air_density(:)  ! kg m-3
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: constituent_props(:)
    real(c_double), target,           intent(inout) :: constituents(:)   ! kg kg-1
    real(c_double), target,           intent(inout) :: rate_params(:)    ! kg kg-1
    character(len=512),                intent(out)   :: errmsg
    integer,                           intent(out)   :: errcode

    real(kind_phys), dimension(4)  :: molar_mass_arr ! kg mol-1
    type(string_t)       :: solver_state
    type(solver_stats_t) :: solver_stats
    type(error_t)        :: error
    integer              :: num_columns, num_layers, num_grid_cells
    integer              :: num_constituents, num_rate_params
    integer              :: i_column, i_layer, i_elem, start
    real(c_double)      :: c_time_step 

    num_columns = 2
    num_layers = 2
    num_grid_cells = 4
    num_constituents = 4
    num_rate_params = 2
    errcode = 0
    errmsg = ''

    c_time_step = real(time_step, c_double) 
    write(*,*) "num_columns ", num_columns
    write(*,*) "num_layers ", num_layers 
    write(*,*) "num_grid_cells ", num_grid_cells
    write(*,*) "num_constituents ", num_constituents
    write(*,*) "num_rate_params ", num_rate_params

    ! Get the molar_mass that is set in the call to instantiate()
    do i_elem = 1, num_constituents
      call constituent_props(i_elem)%molar_mass(molar_mass_arr(i_elem), errcode, errmsg)
      if (errcode /= 0) then
        errmsg = "[MUSICA Error] Unable to get molar mass."
        return
      end if
    end do

    ! ! TODO(jiwon) Check molar mass is non zero as it becomes a denominator for unit converison
    ! ! this code needs to go when ccpp framework does the check
    do i_elem = 1, num_constituents
      if (molar_mass_arr(i_elem) == 0) then
        errcode = 1
        errmsg = "[MUSICA Error] Molar mass must be a non zero value."
        return
      end if
    end do

    ! Convert CAM-SIMA unit to MICM unit (kg kg-1  ->  mol m-3)
    ! call convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)
    ! call convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, rate_params)
  
    c_time_step = real(time_step, c_double) 

    call micm%solve(c_time_step,          &
                    temperature,   &
                    pressure,        &
                    dry_air_density, &
                    constituents,    &
                    rate_params,     &
                    solver_state,         &
                    solver_stats,         &
                    error)
    if (has_error_occurred(error, errmsg, errcode)) then
      write(*,*) "[solve error]: ", errmsg
    end if

    write(*,*) "Solving done: "


    ! Convert MICM unit back to CAM-SIMA unit (mol m-3  ->  kg kg-1)
    ! call convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, constituents)
    ! call convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, rate_params)
  end subroutine micm_run

  !> Finalize MICM
  subroutine micm_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errcode

    errcode = 0
    errmsg = ''

  end subroutine micm_final

  ! Convert CAM-SIMA unit to MICM unit (kg kg-1  ->  mol m-3)
  subroutine convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)
    real(kind_phys), intent(in)    :: dry_air_density(:,:) ! kg m-3
    real(kind_phys), intent(in)    :: molar_mass_arr(:)    ! kg mol-1
    real(kind_phys), intent(inout) :: constituents(:,:,:)  ! in: kg kg-1 | out: mol m-3

    integer :: num_columns, num_layers, num_constituents
    integer :: i_column, i_layer, i_elem

    real(kind_phys) :: val

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)

    do i_column = 1, num_columns
      do i_layer = 1, num_layers
        do i_elem = 1, num_constituents
          val = constituents(i_column, i_layer, i_elem) * dry_air_density(i_column, i_layer) &
                / molar_mass_arr(i_elem)
          constituents(i_column, i_layer, i_elem) = val
        end do
      end do
    end do

  end subroutine convert_to_mol_per_cubic_meter

  ! Convert MICM unit to CAM-SIMA unit (mol m-3  ->  kg kg-1)
  subroutine convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, constituents)
    real(kind_phys), intent(in)    :: dry_air_density(:,:) ! kg m-3
    real(kind_phys), intent(in)    :: molar_mass_arr(:)    ! kg mol-1
    real(kind_phys), intent(inout) :: constituents(:,:,:)  ! in: mol m-3 | out: kg kg-1

    integer :: num_columns, num_layers, num_constituents
    integer :: i_column, i_layer, i_elem

    real(kind_phys) :: val

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)

    do i_column = 1, num_columns
      do i_layer = 1, num_layers
        do i_elem = 1, num_constituents
          val = constituents(i_column, i_layer, i_elem) / dry_air_density(i_column, i_layer) &
                * molar_mass_arr(i_elem)
          constituents(i_column, i_layer, i_elem) = val
        end do
      end do
    end do

  end subroutine convert_to_mass_mixing_ratio

end module musica_ccpp_micm
