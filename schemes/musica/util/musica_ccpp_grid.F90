! Copyright (C) 2025 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_grid

  use ccpp_kinds, only: rk => kind_phys

  implicit none
  private

  public :: grid_t, GRID_INVALID

  !> grid_t defines the dimensions for gridded data used in the model.
  type :: grid_t
    real(rk), allocatable :: interaces_(:) !< Interfaces between grid sections
    real(rk), allocatable :: centers_(:)   !< Centers of grid sections
  contains
    procedure :: number_of_sections => grid_size
  end type grid_t

  interface grid_t
    module procedure grid_constructor_interfaces
    module procedure grid_constructor_interfaces_centers
    module procedure grid_constructor_evenly_spaced
  end interface grid_t

  integer, parameter :: GRID_INVALID = 1

contains

  !> @brief Constructor for grid_t based on interfaces only
  !> @param interfaces The interfaces between grid sections
  !> @param error_message The error message if an error occurs
  !> @param error_code The error code if an error occurs
  !> @return The grid_t instance
  function grid_constructor_interfaces(interfaces, error_message, error_code) &
      result(grid)
    type(grid_t) :: grid
    real(rk),           intent(in)  :: interfaces(:)
    character(len=512), intent(out) :: error_message
    integer,            intent(out) :: error_code
    integer :: i
    error_code = 0
    error_message = ''
    do i = 1, size(interfaces)-1
      if (interfaces(i) >= interfaces(i+1)) then
        error_code = GRID_INVALID
        error_message = 'Interfaces must be in increasing order'
        return
      end if
    end do
    grid%interaces_ = interfaces
    grid%centers_ = 0.5_rk * (grid%interaces_(1:size(interfaces)-1) &
                              + grid%interaces_(2:size(interfaces)))
  end function grid_constructor_interfaces

  !> @brief Constructor for grid_t based on interfaces and centers
  !> @param interfaces The interfaces between grid sections
  !> @param centers The centers of grid sections
  !> @param error_message The error message if an error occurs
  !> @param error_code The error code if an error occurs
  !> @return The grid_t instance
  function grid_constructor_interfaces_centers(interfaces, centers, &
      error_message, error_code) result(grid)
    type(grid_t) :: grid
    real(rk),           intent(in)  :: interfaces(:)
    real(rk),           intent(in)  :: centers(:)
    character(len=512), intent(out) :: error_message
    integer,            intent(out) :: error_code
    integer :: i
    error_code = 0
    error_message = ''
    if (size(interfaces) /= size(centers)+1) then
      error_code = GRID_INVALID
      error_message = 'Invalid dimensions for grid_t interfaces/centers'
    end if
    do i = 1, size(interfaces)-1
      if (interfaces(i) >= interfaces(i+1)) then
        error_code = GRID_INVALID
        error_message = 'Interfaces must be in increasing order'
        return
      end if
      if (centers(i) < interfaces(i) .or. &
          centers(i) > interfaces(i+1)) then
        error_code = GRID_INVALID
        error_message = 'Centers must be within grid interfaces'
        return
      end if
    end do
    grid%interaces_ = interfaces
    grid%centers_ = centers
  end function grid_constructor_interfaces_centers

  !> @brief Constructor for grid_t based on evenly spaced centers
  !> @param start The start of the grid
  !> @param end The end of the grid
  !> @param number_of_sections The number of sections in the grid
  !> @param error_message The error message if an error occurs
  !> @param error_code The error code if an error occurs
  !> @return The grid_t instance
  function grid_constructor_evenly_spaced(start, end, number_of_sections, &
      error_message, error_code) result(grid)
    type(grid_t) :: grid
    real(rk),           intent(in)  :: start
    real(rk),           intent(in)  :: end
    integer,            intent(in)  :: number_of_sections
    character(len=512), intent(out) :: error_message
    integer,            intent(out) :: error_code
    real(rk) :: delta
    integer :: i
    error_code = 0
    error_message = ''
    if (number_of_sections < 1) then
      error_code = GRID_INVALID
      error_message = 'Number of sections must be at least 1'
      return
    end if
    if (start >= end) then
      error_code = GRID_INVALID
      error_message = 'Start must be less than end'
      return
    end if
    delta = (end - start) / real(number_of_sections)
    allocate(grid%interaces_(number_of_sections+1), stat=error_code)
    if (error_code /= 0) then
      error_message = 'Failed to allocate memory for grid interfaces'
      return
    end if
    allocate(grid%centers_(number_of_sections), stat=error_code)
    if (error_code /= 0) then
      error_message = 'Failed to allocate memory for grid centers'
      return
    end if
    grid%interaces_ = (/ (start + real(i-1) * delta, &
                          i=1, number_of_sections+1) /)
    grid%centers_ = (/ (start + real(i-1) * delta + 0.5_rk * delta, &
                        i=1, number_of_sections) /)
  end function grid_constructor_evenly_spaced

  !> @brief Get the number of sections in the grid
  !> @param this The grid_t instance
  !> @return The number of sections
  function grid_size(this) result(number_of_sections)
    class(grid_t), intent(in) :: this
    integer :: number_of_sections
    number_of_sections = size(this%interaces_) - 1
  end function grid_size

end module musica_ccpp_grid