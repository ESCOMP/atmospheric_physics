! ccpp_const_utils contains utility functions that use
! the ccpp constituent properties pointer.
! this code was separated out to remove circular dependencies.
module ccpp_const_utils
  implicit none
  private

  public :: ccpp_const_get_idx

contains

  subroutine ccpp_const_get_idx(constituent_props, name, cindex, errmsg, errflg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! Input arguments
    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props(:)
    character(len=*),                  intent(in)  :: name   ! constituent name

    ! Output arguments
    integer,                           intent(out) :: cindex ! global constituent index
    character(len=512),                intent(out) :: errmsg ! error message
    integer,                           intent(out) :: errflg ! error flag

    ! Local variables
    integer            :: t_cindex
    character(len=256) :: t_const_name

    errmsg = ''
    errflg = 0

    cindex = -1

    ! This convoluted loop is brought to you in exchange for avoiding a
    ! circular dependency on cam_ccpp_cap::cam_const_get_index.
    const_props_loop: do t_cindex = lbound(constituent_props, 1), ubound(constituent_props, 1)
       call constituent_props(t_cindex)%standard_name(t_const_name, errflg, errmsg)
       if (errflg /= 0) then
          ! Abort subroutine and return with error.
          return
       end if

       ! Constituent standard names are case-insensitive.
       if (to_lower(trim(t_const_name)) == to_lower(trim(name))) then
          cindex = t_cindex
          exit const_props_loop
       end if
    enddo const_props_loop

  end subroutine ccpp_const_get_idx

  ! Return a copy of <str> with ASCII upper-case letters folded to lower case.
  ! This is inlined here to avoid host dependencies and should be retired together with the rest
  ! of this module once it is replaced by ccpp_scheme_utils.
  pure function to_lower(str) result(lower_str)
    character(len=*), intent(in) :: str
    character(len=len(str))      :: lower_str

    integer :: i, ascii_code
    integer, parameter :: upper_to_lower = iachar('a') - iachar('A')

    do i = 1, len(str)
      ascii_code = iachar(str(i:i))
      if (ascii_code >= iachar('A') .and. ascii_code <= iachar('Z')) then
        lower_str(i:i) = achar(ascii_code + upper_to_lower)
      else
        lower_str(i:i) = str(i:i)
      end if
    end do

  end function to_lower

end module ccpp_const_utils
