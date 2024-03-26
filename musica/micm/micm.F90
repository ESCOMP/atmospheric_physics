module micm
  use iso_c_binding

  ! Note: "micm_core" is included in an external pre-built MICM library that the host
  ! model is responsible for linking to during compilation
  use micm_core, only: micm_t
  use ccpp_kinds, only: kind_phys

  implicit none

  public :: micm_init, micm_run, micm_final
  private :: convert_to_mol_per_cubic_meter, convert_to_mass_mixing_ratio

  type(micm_t), pointer         :: micm_obj

contains

  !> \section arg_table_micm_init Argument Table
  !! \htmlinclude micm_init.html
  subroutine micm_init(config_path, iulog, errcode, errmsg)
    character(len=*), intent(in)    :: config_path
    integer, intent(in)             :: iulog
    integer, intent(out)            :: errcode
    character(len=512), intent(out) :: errmsg

    errcode = 0
    errmsg = ''

    micm_obj => micm_t(config_path, errcode)

    if (errcode /= 0) then
      errmsg = "[fatal] [micm] Failed to create MICM solver. Parsing configuration failed. &
                Please look over at MICM log file for further information."
      return
    endif

    write(iulog,*) "[info] [micm] Created MICM solver."

  end subroutine micm_init

  !> \section arg_table_micm_run Argument Table
  !! \htmlinclude micm_run.html
  subroutine micm_run(time_step, temperature, pressure, dry_air_density, constituent_props, &
                      constituents, iulog, errcode, errmsg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    real(kind_phys),                   intent(in)    :: time_step            ! s
    real(kind_phys),                   intent(in)    :: temperature(:,:)     ! K
    real(kind_phys),                   intent(in)    :: pressure(:,:)        ! Pa
    real(kind_phys),                   intent(in)    :: dry_air_density(:,:) ! kg m-3
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: constituent_props(:)
    real(kind_phys),                   intent(inout) :: constituents(:,:,:)  ! kg kg-1
    integer,                           intent(in)    :: iulog
    integer,                           intent(out)   :: errcode
    character(len=512),                intent(out)   :: errmsg

    ! local variables
    real(c_double)                                         :: c_time_step
    real(c_double), dimension(size(temperature, dim=1), &
                              size(temperature, dim=2))    :: c_temperature
    real(c_double), dimension(size(pressure, dim=1), &
                              size(pressure, dim=2))       :: c_pressure
    real(c_double), dimension(size(constituents, dim=1), &
                              size(constituents, dim=2), &
                              size(constituents, dim=3))   :: c_constituents

    real(kind_phys), dimension(size(constituents, dim=3))  :: molar_mass_arr ! kg mol-1

    integer :: num_columns, num_layers, num_constituents
    integer :: i_column, i_layer, i_elem

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)

    errcode = 0
    errmsg = ''

    ! Get the molar_mass that is set in the call to instantiate()
    do i_elem = 1, num_constituents
      call constituent_props(i_elem)%molar_mass(molar_mass_arr(i_elem), errcode, errmsg)

      if (errcode /= 0) then
        errmsg = "[error] [micm] Unable to get molar mass."
        return
      end if
    end do

    ! TODO(jiwon) Check molar mass is non zero as it becomes a denominator for unit converison
    ! this code needs to go when ccpp framework does the check
    do i_elem = 1, num_constituents
      if (molar_mass_arr(i_elem) == 0) then
        errcode = 1
        errmsg = "[error] [micm] Molar mass must be a non zero value."
        return
      end if
    end do

    ! Convert CAM-SIMA unit to MICM unit (kg kg-1  ->  mol m-3)
    call convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)

    c_time_step = real(time_step, c_double)
    c_temperature = real(temperature, c_double)
    c_pressure = real(pressure, c_double)
    c_constituents = real(constituents, c_double)

    write(iulog,*) "[info] [micm] Running MICM solver..."

    do i_column = 1, num_columns
      do i_layer = 1, num_layers
        call micm_obj%solve(c_temperature(i_column, i_layer), c_pressure(i_column, i_layer), &
                        c_time_step, num_constituents, c_constituents(i_column, i_layer, :))
      end do
    end do

    write(iulog,*) "[info] [micm] MICM solver has finished."

    constituents = real(c_constituents, kind_phys)

    ! Convert MICM unit back to CAM-SIMA unit (mol m-3  ->  kg kg-1)
    call convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, constituents)

  end subroutine micm_run

  !> \section arg_table_micm_final Argument Table
  !! \htmlinclude micm_final.html
  subroutine micm_final(iulog, errcode, errmsg)
    integer, intent(in)             :: iulog
    integer, intent(out)            :: errcode
    character(len=512), intent(out) :: errmsg

    errcode = 0
    errmsg = ''

    write(iulog,*) "[debug] [micm] Deallocating MICM object..."

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

end module micm