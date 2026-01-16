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

       if (trim(t_const_name) == trim(name)) then
          cindex = t_cindex
          exit const_props_loop
       end if
    enddo const_props_loop

  end subroutine ccpp_const_get_idx

end module ccpp_const_utils
