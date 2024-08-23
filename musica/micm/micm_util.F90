module micm_util
  use iso_c_binding

  use ISO_FORTRAN_ENV, only: kind_phys => REAL64
  ! use ccpp_kinds,           only: kind_phys
  implicit none

  private

  public :: reshape_into_micm_arr, reshape_into_ccpp_arr
  ! kind_phyiscs -> c_double
  ! 3D -> 1D
contains

  subroutine reshape_into_micm_arr(temperature, pressure, dry_air_density, constituents, rate_params, &
    m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)
    use iso_c_binding

    use ISO_FORTRAN_ENV, only: kind_phys => REAL64
    implicit none
    real(kind_phys), target, intent(in) :: temperature(:,:)     ! K
    real(kind_phys), target, intent(in) :: pressure(:,:)        ! Pa
    real(kind_phys), target, intent(in) :: dry_air_density(:,:) ! kg m-3
    real(kind_phys), target, intent(in) :: constituents(:,:,:)  ! kg kg-1
    real(kind_phys), target, intent(in) :: rate_params(:,:,:)
    real(c_double), target, intent(out) :: m_temperature(:)     ! K
    real(c_double), target, intent(out) :: m_pressure(:)        ! Pa
    real(c_double), target, intent(out) :: m_dry_air_density(:) ! kg m-3
    real(c_double), target, intent(out) :: m_constituents(:)  ! kg kg-1
    real(c_double), target, intent(out) :: m_rate_params(:)

    ! local variables
    integer :: num_columns, num_layers, num_grid_cells
    integer :: num_constituents, num_rate_params
    integer :: i_column, i_layer, i_elem, i_constituents, i_rate_params

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_grid_cells = num_columns * num_layers
    num_constituents = size(constituents, dim=3)
    num_rate_params = size(rate_params, dim=3)

    write(*,*) "temperature2",num_columns
    write(*,*) "temperature3",num_layers
    write(*,*) "temperature4",num_grid_cells
    write(*,*) "temperature5",num_constituents
    write(*,*) "temperature6",num_rate_params
    
    i_elem = 1
    i_constituents = 1
    i_rate_params = 1
    do i_layer = 1, num_layers
      do i_column = 1, num_columns
        m_temperature(i_elem) = real(temperature(i_column, i_layer), c_double)
        m_pressure(i_elem) = real(pressure(i_column, i_layer), c_double)
        m_dry_air_density(i_elem) = real(dry_air_density(i_column, i_layer), c_double)
        m_constituents(i_constituents : i_constituents + num_constituents - 1) = real(constituents(i_column, i_layer, :), c_double)
        m_rate_params(i_rate_params : i_rate_params + num_constituents - 1) = real(rate_params(i_column, i_layer, :), c_double)
        i_elem = i_elem + 1
        i_constituents = i_constituents + num_constituents
        i_rate_params = i_rate_params + num_rate_params
      end do
    end do

  end subroutine reshape_into_micm_arr

  subroutine reshape_into_ccpp_arr(temperature, pressure, dry_air_density, constituents, rate_params, &
    m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)
    use iso_c_binding

    use ISO_FORTRAN_ENV, only: kind_phys => REAL64
    implicit none
    real(kind_phys), intent(out) :: temperature(:,:)     ! K
    real(kind_phys), intent(out) :: pressure(:,:)        ! Pa
    real(kind_phys), intent(out) :: dry_air_density(:,:) ! kg m-3
    real(kind_phys), intent(out) :: constituents(:,:,:)  ! kg kg-1
    real(kind_phys), intent(out) :: rate_params(:,:,:)
    real(c_double), intent(in) :: m_temperature(:)     ! K
    real(c_double), intent(in) :: m_pressure(:)        ! Pa
    real(c_double), intent(in) :: m_dry_air_density(:) ! kg m-3
    real(c_double), intent(in) :: m_constituents(:)  ! kg kg-1
    real(c_double), intent(in) :: m_rate_params(:)

        ! local variables
    integer :: num_columns, num_layers, num_grid_cells
    integer :: num_constituents, num_rate_params
    integer :: i_column, i_layer, i_elem, i_constituents, i_rate_params

    num_columns = size(constituents, dim=1)
    num_layers = size(constituents, dim=2)
    num_grid_cells = num_columns * num_layers
    num_constituents = size(constituents, dim=3)
    num_rate_params = size(rate_params, dim=3)
    
    i_elem = 1
    i_constituents = 1
    i_rate_params = 1
    do i_layer = 1, num_layers
      do i_column = 1, num_columns
        temperature(i_column, i_layer) = real(m_temperature(i_elem), c_double)
        pressure(i_column, i_layer) = real(m_pressure(i_elem), c_double)
        dry_air_density(i_column, i_layer) = real(m_dry_air_density(i_elem), c_double)
        constituents(i_column, i_layer, :) = real(m_constituents(i_constituents : i_constituents + num_constituents - 1), kind_phys)
        rate_params(i_column, i_layer, :) = real(m_rate_params(i_rate_params : i_rate_params + num_constituents - 1), c_double)
        i_elem = i_elem + 1
        i_constituents = i_constituents + num_constituents
        i_rate_params = i_rate_params + num_rate_params
      end do
    end do

  end subroutine reshape_into_ccpp_arr

end module micm_util

!   program ttt
!     use iso_c_binding

!     use ISO_FORTRAN_ENV, only: kind_phys => REAL64
!     use micm_util
!     integer, parameter                                                     :: NUM_SPECIES = 4
!     integer, parameter                                                     :: NUM_RATES = 2
!     integer, parameter                                                     :: NUM_COLUMNS = 2
!     integer, parameter                                                     :: NUM_LAYERS = 2
!     integer, parameter                                                     :: NUM_GRID_CELLS = 4
!     real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: temperature     ! K
!     real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: pressure        ! Pa
!     real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)             :: dry_air_density ! kg m-3
!     real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) :: constituents    ! kg kg-1
!     real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,NUM_RATES)   :: rate_params


!     real(c_double), target :: m_temperature(NUM_GRID_CELLS)     ! K
!     real(c_double), target :: m_pressure(NUM_GRID_CELLS)        ! Pa
!     real(c_double), target :: m_dry_air_density(NUM_GRID_CELLS) ! kg m-3
!     real(c_double), target :: m_constituents(NUM_GRID_CELLS*NUM_SPECIES)  ! kg kg-1
!     real(c_double), target :: m_rate_params(NUM_GRID_CELLS*NUM_RATES)

!     temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
!     temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)
!     pressure(:,1) = (/ 6000.04_kind_phys, 7000.04_kind_phys /)
!     pressure(:,2) = (/ 8000.04_kind_phys, 9000.04_kind_phys /)
!     dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
!     dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
!     constituents(1,1,:) = (/ 0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys, 0.4_kind_phys /)
!     constituents(1,2,:) = (/ 0.41_kind_phys, 0.42_kind_phys, 0.43_kind_phys, 0.44_kind_phys /)
!     constituents(2,1,:) = (/ 0.21_kind_phys, 0.22_kind_phys, 0.23_kind_phys, 0.24_kind_phys /)
!     constituents(2,2,:) = (/ 0.31_kind_phys, 0.32_kind_phys, 0.33_kind_phys, 0.34_kind_phys /)
!     rate_params(1,1,:) = (/2.7e-19_kind_phys, 1.13e-9_kind_phys/)
!     rate_params(1,2,:) = (/2.7e-19_kind_phys, 1.13e-9_kind_phys/)
!     rate_params(2,1,:) = (/2.7e-19_kind_phys, 1.13e-9_kind_phys/)
!     rate_params(2,2,:) = (/2.7e-19_kind_phys, 1.13e-9_kind_phys/)

!     write(*,*) "temperature: ", temperature

!     call reshape_into_micm_arr(temperature, pressure, dry_air_density, constituents, rate_params, &
!     m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)

!     write(*,*) "m_temperature: ", m_temperature

!     write(*,*) "-----------: "

!     call reshape_into_micm_arr(temperature, pressure, dry_air_density, constituents, rate_params, &
!     m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)

!     write(*,*) "ffftemperature: ", temperature

!   end program ttt



! ! end module micm_util