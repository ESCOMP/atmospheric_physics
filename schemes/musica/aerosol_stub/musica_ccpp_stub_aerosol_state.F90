! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_stub_aerosol_state

  use musica_ccpp_aerosol_state, only: aerosol_state_t

  implicit none
  private

  public :: stub_aerosol_state_t

  !> stub_aerosol_state_t defines the state of an aerosol system according to
  !! the aerosol representation of the stub aerosol package.
  type, extends(aerosol_state_t) :: stub_aerosol_state_t
    integer :: number_of_columns_ = 0 !< The number of columns in the model grid
    integer :: number_of_levels_ = 0  !< The number of levels in the model grid
  contains
    procedure :: number_of_columns => stub_aerosol_state_number_of_columns
    procedure :: number_of_levels => stub_aerosol_state_number_of_levels
  end type stub_aerosol_state_t

  interface stub_aerosol_state_t
    module procedure stub_aerosol_state_constructor
  end interface stub_aerosol_state_t

contains

  !> @brief Constructor for stub_aerosol_state_t
  !! @param number_of_columns The number of columns in the model grid.
  !! @param number_of_levels The number of levels in the model grid.
  !! @return The stub aerosol state instance.
  function stub_aerosol_state_constructor(number_of_columns, number_of_levels, &
      error_message, error_code) result(state)
    type(stub_aerosol_state_t), pointer :: state
    integer,            intent(in) :: number_of_columns
    integer,            intent(in) :: number_of_levels
    character(len=512), intent(out) :: error_message
    integer,            intent(out) :: error_code
    error_message = ''
    error_code = 0
    allocate(state, stat=error_code, errmsg=error_message)
    if (error_code /= 0) return
    error_message = ''
    state%number_of_columns_ = number_of_columns
    state%number_of_levels_ = number_of_levels
  end function stub_aerosol_state_constructor

  !> @brief Returns the number of columns in the model grid.
  !! @param this The stub aerosol state instance.
  !! @return The number of columns in the model grid.
  function stub_aerosol_state_number_of_columns(this) result(number_of_columns)
    class(stub_aerosol_state_t), intent(in) :: this
    integer :: number_of_columns
    number_of_columns = this%number_of_columns_
  end function stub_aerosol_state_number_of_columns

  !> @brief Returns the number of levels in the model grid.
  !! @param this The stub aerosol state instance.
  !! @return The number of levels in the model grid.
  function stub_aerosol_state_number_of_levels(this) result(number_of_levels)
    class(stub_aerosol_state_t), intent(in) :: this
    integer :: number_of_levels
    number_of_levels = this%number_of_levels_
  end function stub_aerosol_state_number_of_levels

end module musica_ccpp_stub_aerosol_state