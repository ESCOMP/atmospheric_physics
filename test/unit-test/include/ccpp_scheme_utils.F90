! Minimal CCPP framework stub to satisfy dependencies.
module ccpp_scheme_utils
    implicit none

    private
    public :: ccpp_constituent_index
contains
    subroutine ccpp_constituent_index(standard_name, const_index, errcode, errmsg)
        use ccpp_constituent_prop_mod, only: int_unassigned

        character(*), intent(in) :: standard_name
        integer, intent(out) :: const_index
        integer, optional, intent(out) :: errcode
        character(*), optional, intent(out) :: errmsg

        if (standard_name /= '') then
            const_index = 1
        else
            const_index = int_unassigned
        end if

        if (present(errcode)) then
            errcode = 0
        end if

        if (present(errmsg)) then
            errmsg = ''
        end if
    end subroutine ccpp_constituent_index
end module ccpp_scheme_utils
