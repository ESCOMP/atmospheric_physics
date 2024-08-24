module micm_util
  implicit none

  private
  public :: reshape_into_micm_arr, reshape_into_ccpp_arr, convert_to_mol_per_cubic_meter, &
            convert_to_mass_mixing_ratio

contains

  subroutine reshape_into_micm_arr(temperature, pressure, dry_air_density, constituents, &
      rate_params, m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)
    use iso_c_binding, only: c_double
    use ccpp_kinds,    only: kind_phys

    real(kind_phys), target, intent(in)  :: temperature(:,:)     ! K
    real(kind_phys), target, intent(in)  :: pressure(:,:)        ! Pa
    real(kind_phys), target, intent(in)  :: dry_air_density(:,:) ! kg m-3
    real(kind_phys), target, intent(in)  :: constituents(:,:,:)  ! kg kg-1
    real(kind_phys), target, intent(in)  :: rate_params(:,:,:)
    real(c_double),  target, intent(out) :: m_temperature(:)     ! K
    real(c_double),  target, intent(out) :: m_pressure(:)        ! Pa
    real(c_double),  target, intent(out) :: m_dry_air_density(:) ! kg m-3
    real(c_double),  target, intent(out) :: m_constituents(:)    ! kg kg-1
    real(c_double),  target, intent(out) :: m_rate_params(:)

    ! local variables
    integer :: num_columns, num_layers
    integer :: num_constituents, num_rate_params
    integer :: i_column, i_layer, i_elem, i_constituents, i_rate_params

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)
    num_rate_params = size(rate_params, dim=3)
    
    ! Reshape into 1-D arry in species-column first order
    ! refers to: state.variables_[i_cell][i_species] = concentrations[i_species_elem++]
    i_elem = 1
    i_constituents = 1
    i_rate_params = 1
    do i_layer = 1, num_layers
      do i_column = 1, num_columns
        m_temperature(i_elem) = real(temperature(i_column, i_layer), c_double)
        m_pressure(i_elem) = real(pressure(i_column, i_layer), c_double)
        m_dry_air_density(i_elem) = real(dry_air_density(i_column, i_layer), c_double)
        m_constituents(i_constituents : i_constituents + num_constituents - 1) &
                            = real(constituents(i_column, i_layer, :), c_double)
        m_rate_params(i_rate_params : i_rate_params + num_constituents - 1) &
                            = real(rate_params(i_column, i_layer, :), c_double)
        i_elem = i_elem + 1
        i_constituents = i_constituents + num_constituents
        i_rate_params = i_rate_params + num_rate_params
      end do
    end do

  end subroutine reshape_into_micm_arr

  subroutine reshape_into_ccpp_arr(temperature, pressure, dry_air_density, constituents, &
      rate_params, m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)
    use iso_c_binding, only: c_double
    use ccpp_kinds,    only: kind_phys
    real(kind_phys), intent(out) :: temperature(:,:)     ! K
    real(kind_phys), intent(out) :: pressure(:,:)        ! Pa
    real(kind_phys), intent(out) :: dry_air_density(:,:) ! kg m-3
    real(kind_phys), intent(out) :: constituents(:,:,:)  ! kg kg-1
    real(kind_phys), intent(out) :: rate_params(:,:,:)
    real(c_double),  intent(in)  :: m_temperature(:)     ! K
    real(c_double),  intent(in)  :: m_pressure(:)        ! Pa
    real(c_double),  intent(in)  :: m_dry_air_density(:) ! kg m-3
    real(c_double),  intent(in)  :: m_constituents(:)    ! kg kg-1
    real(c_double),  intent(in)  :: m_rate_params(:)

    ! local variables
    integer :: num_columns, num_layers
    integer :: num_constituents, num_rate_params
    integer :: i_column, i_layer, i_elem, i_constituents, i_rate_params

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)
    num_rate_params = size(rate_params, dim=3)
    
    i_elem = 1
    i_constituents = 1
    i_rate_params = 1
    do i_layer = 1, num_layers
      do i_column = 1, num_columns
        temperature(i_column, i_layer) = real(m_temperature(i_elem), kind_phys)
        pressure(i_column, i_layer) = real(m_pressure(i_elem), kind_phys)
        dry_air_density(i_column, i_layer) = real(m_dry_air_density(i_elem), kind_phys)
        constituents(i_column, i_layer, :) &
            = real(m_constituents(i_constituents : i_constituents + num_constituents - 1), kind_phys)
        rate_params(i_column, i_layer, :) &
            = real(m_rate_params(i_rate_params : i_rate_params + num_constituents - 1), kind_phys)
        i_elem = i_elem + 1
        i_constituents = i_constituents + num_constituents
        i_rate_params = i_rate_params + num_rate_params
      end do
    end do

  end subroutine reshape_into_ccpp_arr

  ! Convert CAM-SIMA unit to MICM unit (kg kg-1  ->  mol m-3)
  subroutine convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)
    use ccpp_kinds,    only: kind_phys

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
    use ccpp_kinds,    only: kind_phys

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

end module micm_util