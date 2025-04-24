module ccpp_io_reader
    implicit none
    private

    public :: abstract_netcdf_reader_t
    public :: create_netcdf_reader_t

    type, abstract :: abstract_netcdf_reader_t
    contains
        procedure(open_file),    deferred :: open_file
        procedure(close_file),   deferred :: close_file
        procedure(get_var_int),  deferred :: get_var_int
        procedure(get_var_real), deferred :: get_var_real
        procedure(get_var_char), deferred :: get_var_char

        generic :: get_var => get_var_int, get_var_real, get_var_char
    end type abstract_netcdf_reader_t

    interface
        module function create_netcdf_reader_t() result(r)
            class(abstract_netcdf_reader_t), allocatable :: r
        end function create_netcdf_reader_t

        subroutine open_file(this, file_path, errmsg, errcode)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(inout) :: this
            character(len=*),                intent(in)    :: file_path
            character(len=*),                intent(out)   :: errmsg
            integer,                         intent(out)   :: errcode   !Error code
        end subroutine open_file

        subroutine close_file(this, errmsg, errcode)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(inout) :: this
            character(len=*),                intent(out)   :: errmsg
            integer,                         intent(out)   :: errcode !Error code
        end subroutine close_file

        subroutine get_var_int(this, varname, var, errmsg, errcode)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            integer, pointer,                intent(out) :: var(..) !Integer variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code
        end subroutine get_var_int

        subroutine get_var_real(this, varname, var, errmsg, errcode)
            use ccpp_kinds, only: kind_phys
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            real(kind_phys), pointer,        intent(out) :: var(..) !Floating-point variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code
        end subroutine get_var_real

        subroutine get_var_char(this, varname, var, errmsg, errcode)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            character(len=:), pointer,       intent(out) :: var(..) !Character variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code
        end subroutine get_var_char
        
    end interface

end module ccpp_io_reader
