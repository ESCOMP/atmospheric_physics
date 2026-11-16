! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!> Test-only stand-in for the host model's constituent registration.
!! Builds a CCPP model constituents object from a set of constituent properties
!! and makes it available to ccpp_scheme_utils, so that schemes can look up
!! constituent indices with ccpp_constituent_index.
module musica_test_constituents
  use ccpp_constituent_prop_mod, only: ccpp_model_constituents_t

  implicit none
  private

  public :: register_test_constituents, cleanup_test_constituents

  ! ccpp_scheme_utils can only be pointed at one constituents object,
  ! so every test in an executable shares this one
  type(ccpp_model_constituents_t), target, save :: model_constituents
  logical,                                 save :: is_scheme_utils_initialized = .false.

contains

  !> Registers constituent properties with the CCPP framework and returns
  !! the framework's constituent properties array.
  !! Constituent indices are assigned by the framework (advected constituents first),
  !! so they may not match the order of <constituent_props>.
  subroutine register_test_constituents(constituent_props, constituent_props_ptr, &
                                        errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t, &
                                         ccpp_constituent_prop_ptr_t
    use ccpp_scheme_utils,         only: ccpp_initialize_constituent_ptr

    type(ccpp_constituent_properties_t),            intent(in)  :: constituent_props(:)
    type(ccpp_constituent_prop_ptr_t), allocatable, intent(out) :: constituent_props_ptr(:)
    character(len=512),                             intent(out) :: errmsg
    integer,                                        intent(out) :: errcode

    ! local variables
    type(ccpp_model_constituents_t),     pointer :: model_constituents_ptr
    type(ccpp_constituent_properties_t), pointer :: const_prop
    integer :: i

    errmsg = ''
    errcode = 0

    ! Clears (and deallocates) anything registered by a previous test
    call model_constituents%initialize_table(size(constituent_props))

    do i = 1, size(constituent_props)
      ! The framework takes ownership of (and later deallocates) each property,
      ! so it gets its own copy
      allocate(const_prop)
      const_prop = constituent_props(i)
      call model_constituents%new_field(const_prop, errcode, errmsg)
      if (errcode /= 0) return
    end do

    call model_constituents%lock_table(errcode, errmsg)
    if (errcode /= 0) return

    if (.not. is_scheme_utils_initialized) then
      model_constituents_ptr => model_constituents
      call ccpp_initialize_constituent_ptr(model_constituents_ptr)
      is_scheme_utils_initialized = .true.
    end if

    constituent_props_ptr = model_constituents%const_metadata

  end subroutine register_test_constituents

  !> Deallocates all constituents registered with the CCPP framework
  subroutine cleanup_test_constituents()

    call model_constituents%reset()

  end subroutine cleanup_test_constituents

end module musica_test_constituents
