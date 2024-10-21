module musica_ccpp_micm_util
  implicit none

  private
  public :: reshape_into_micm_arr, reshape_into_ccpp_arr, convert_to_mol_per_cubic_meter, &
            convert_to_mass_mixing_ratio

contains

  !> Reshape array (2D/3D -> 1D) and convert type (kind_phys -> c_double)
  subroutine reshape_into_micm_arr(temperature, pressure, dry_air_density, constituents, &
      micm_temperature, micm_pressure, micm_dry_air_density, micm_constituents)
    use iso_c_binding, only: c_double
    use ccpp_kinds,    only: kind_phys

    real(kind_phys), intent(in)  :: temperature(:,:)        ! K
    real(kind_phys), intent(in)  :: pressure(:,:)           ! Pa
    real(kind_phys), intent(in)  :: dry_air_density(:,:)    ! kg m-3
    real(kind_phys), intent(in)  :: constituents(:,:,:)     ! kg kg-1
    real(c_double),  intent(out) :: micm_temperature(:)     ! K
    real(c_double),  intent(out) :: micm_pressure(:)        ! Pa
    real(c_double),  intent(out) :: micm_dry_air_density(:) ! kg m-3
    real(c_double),  intent(out) :: micm_constituents(:)    ! kg kg-1

    ! local variables
    integer :: num_columns, num_layers, num_constituents
    integer :: i_column, i_layer, i_elem, i_constituents

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)
    
    ! Reshape into 1-D arry in species-column first order, referring to
    ! state.variables_[i_cell][i_species] = concentrations[i_species_elem++]
    i_elem = 1
    i_constituents = 1
    do i_layer = 1, num_layers
      do i_column = 1, num_columns
        micm_temperature(i_elem) = real(temperature(i_column, i_layer), c_double)
        micm_pressure(i_elem) = real(pressure(i_column, i_layer), c_double)
        micm_dry_air_density(i_elem) = real(dry_air_density(i_column, i_layer), c_double)
        micm_constituents(i_constituents : i_constituents + num_constituents - 1) &
                            = real(constituents(i_column, i_layer, :), c_double)
        i_elem = i_elem + 1
        i_constituents = i_constituents + num_constituents
      end do
    end do

  end subroutine reshape_into_micm_arr

  !> Reshape array (1D -> 3D) and convert type (c_double -> kind_phys)
  subroutine reshape_into_ccpp_arr(micm_constituents, constituents)
    use iso_c_binding, only: c_double
    use ccpp_kinds,    only: kind_phys

    real(c_double),  intent(in)    :: micm_constituents(:) ! kg kg-1
    real(kind_phys), intent(inout) :: constituents(:,:,:)  ! kg kg-1

    ! local variables
    integer :: num_columns, num_layers, num_constituents
    integer :: i_column, i_layer, i_constituents

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)
    
    i_constituents = 1
    do i_layer = 1, num_layers
      do i_column = 1, num_columns
        constituents(i_column, i_layer, :) &
            = real(micm_constituents(i_constituents : i_constituents + num_constituents - 1), kind_phys)
        i_constituents = i_constituents + num_constituents
      end do
    end do

  end subroutine reshape_into_ccpp_arr

  ! Convert CAM-SIMA unit to MICM unit (kg kg-1  ->  mol m-3)
  subroutine convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)
    use ccpp_kinds,    only: kind_phys

    real(kind_phys), intent(in)    :: dry_air_density(:,:) ! kg m-3
    real(kind_phys), intent(in)    :: molar_mass_arr(:)    ! kg mol-1
    real(kind_phys), intent(inout) :: constituents(:,:,:)  ! in: kg kg-1 | out: mol m-3

    ! local variables
    integer         :: num_columns, num_layers, num_constituents
    integer         :: i_column, i_layer, i_elem
    real(kind_phys) :: val

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)

    do i_elem = 1, num_constituents
      do i_layer = 1, num_layers
        do i_column = 1, num_columns
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

    integer         :: num_columns, num_layers, num_constituents
    integer         :: i_column, i_layer, i_elem
    real(kind_phys) :: val

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_constituents = size(constituents, dim=3)

    do i_elem = 1, num_constituents
      do i_layer = 1, num_layers
        do i_column = 1, num_columns
          val = constituents(i_column, i_layer, i_elem) / dry_air_density(i_column, i_layer) &
                * molar_mass_arr(i_elem)
          constituents(i_column, i_layer, i_elem) = val
        end do
      end do
    end do

  end subroutine convert_to_mass_mixing_ratio

end module musica_ccpp_micm_util