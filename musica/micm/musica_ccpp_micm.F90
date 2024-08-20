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
                      constituents, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use musica_micm, only: solver_stats_t
    use musica_util, only: string_t, error_t

    real(kind_phys),                   intent(in)    :: time_step            ! [s]
    real(kind_phys),                   intent(in)    :: temperature(:,:)     ! [K]
    real(kind_phys),                   intent(in)    :: pressure(:,:)        ! [Pa]
    real(kind_phys),                   intent(in)    :: dry_air_density(:,:) ! [kg m-3]
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: constituent_props(:)
    real(kind_phys),                   intent(inout) :: constituents(:,:,:)  ! [kg kg-1]
    character(len=512),                intent(out)   :: errmsg
    integer,                           intent(out)   :: errcode

    ! local variables
    real(c_double)                                         :: c_time_step
    real(c_double), dimension(size(temperature, dim=1),  &
                              size(temperature, dim=2))    :: c_temperature
    real(c_double), dimension(size(pressure, dim=1),     &
                              size(pressure, dim=2))       :: c_pressure
    real(c_double), dimension(size(pressure, dim=1),     &
                              size(pressure, dim=2))       :: c_dry_air_density
    real(c_double), dimension(size(constituents, dim=1), &
                              size(constituents, dim=2), &
                              size(constituents, dim=3))   :: c_constituents
    real(c_double), dimension(size(constituents, dim=1), &
                              size(constituents, dim=2), &
                              0)                           :: c_rate_params
    real(kind_phys), dimension(size(constituents, dim=3))  :: molar_mass_arr ! [kg mol-1]

    ! 1-D array
    real(c_double), dimension(size(temperature, dim=1)  &
                            * size(temperature, dim=2))  :: c_temperature_1d
    real(c_double), dimension(size(pressure, dim=1)     &
                            * size(pressure, dim=2))     :: c_pressure_1d
    real(c_double), dimension(size(pressure, dim=1)     &
                            * size(pressure, dim=2))     :: c_dry_air_density_1d
    real(c_double), dimension(size(constituents, dim=1) &
                            * size(constituents, dim=2) & 
                            * size(constituents, dim=3)) :: c_constituents_1d
    real(c_double), dimension(size(constituents, dim=1) &
                            * size(constituents, dim=2)) :: c_rate_params_1d

    type(string_t)       :: solver_state
    type(solver_stats_t) :: solver_stats
    type(error_t)        :: error
    integer              :: num_columns, num_layers, num_constituents, num_grid_cells
    integer              :: i_column, i_layer, i_elem, start

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)
    num_grid_cells = num_columns * num_layers
    errcode = 0
    errmsg = ''

    ! Get the molar_mass that is set in the call to instantiate()
    do i_elem = 1, num_constituents
      call constituent_props(i_elem)%molar_mass(molar_mass_arr(i_elem), errcode, errmsg)

      if (errcode /= 0) then
        errmsg = "[MUSICA Error] Unable to get molar mass."
        return
      end if
    end do

    ! TODO(jiwon) Check molar mass is non zero as it becomes a denominator for unit converison
    ! this code needs to go when ccpp framework does the check
    do i_elem = 1, num_constituents
      if (molar_mass_arr(i_elem) == 0) then
        errcode = 1
        errmsg = "[MUSICA Error] Molar mass must be a non zero value."
        return
      end if
    end do

    ! Convert CAM-SIMA unit to MICM unit (kg kg-1  ->  mol m-3)
    call convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)

    c_time_step = real(time_step, c_double)
    c_temperature = real(temperature, c_double)
    c_pressure = real(pressure, c_double)
    c_dry_air_density = real(dry_air_density, c_double)
    c_constituents = real(constituents, c_double)

    ! Reshape to 1-D array
    c_temperature_1d = [c_temperature]
    c_pressure_1d = [c_pressure]
    c_dry_air_density_1d = [dry_air_density]
    c_rate_params_1d = [c_rate_params]
    
    ! Reshape to 1-D arry in species-column first order
    ! refers to: state.variables_[i_cell][i_species] = concentrations[i_species_elem++]
    start = 1
    do i_layer = 1, num_layers
      do i_column = 1, num_columns
        c_constituents_1d(start : start + num_constituents - 1) = c_constituents(i_column, i_layer, :)
        start = start + num_constituents
      end do
    end do

    write(*,*) "1d concentrations :   ", c_constituents_1d
    ! call micm%solve(c_time_step,     &
    !                 c_temperature_1d,    &
    !                 c_pressure_1d,      &
    !                 c_dry_air_density_1d, & 
    !                 c_constituents_1d,  &
    !                 c_rate_params_1d,   &
    !                 solver_state,                         &
    !                 solver_stats,                         &
    !                 error)
    ! if (has_error_occurred(error, errmsg, errcode)) return

    !TODO(jiwon) add something like this for test
    ! do start = 1, 20
    !     c_constituents_1d(start) = c_constituents_1d(start) + 100
    ! end do

    ! Reshape back to the original dimension
    start = 1
    do i_layer = 1, num_layers
      do i_column = 1, num_columns
        c_constituents(i_column, i_layer, :) = c_constituents_1d(start : start + num_constituents - 1)
        start = start + num_constituents
      end do
    end do

    write(*,*) "[musica] Jiwon Concentrations"
    ! the first digit in the format (4) indicates the number of column * the number of vertical layers
    write(*,fmt="(4(3x,f13.6))") c_constituents

    constituents = real(c_constituents, kind_phys)

    ! Convert MICM unit back to CAM-SIMA unit (mol m-3  ->  kg kg-1)
    call convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, constituents)

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

  subroutine reshape_to_1d_array(array_1d)
    real(c_double) array_1d(*)
  end subroutine

end module musica_ccpp_micm
