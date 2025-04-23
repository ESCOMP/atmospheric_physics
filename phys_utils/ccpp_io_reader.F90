module ccpp_io_reader
    implicit none
    private

    public :: ccpp_io_reader_t
    public :: create_io_reader_t

    type, abstract :: ccpp_io_reader_t
    contains
        procedure(open_netcdf_file),    deferred :: open_netcdf_file
        procedure(close_netcdf_file),   deferred :: close_netcdf_file
        procedure(get_netcdf_var_int),  deferred :: get_netcdf_var_int
        procedure(get_netcdf_var_real), deferred :: get_netcdf_var_real
        procedure(get_netcdf_var_char), deferred :: get_netcdf_var_char

        generic :: get_netcdf_var => get_netcdf_var_int, get_netcdf_var_real, get_netcdf_var_char
    end type ccpp_io_reader_t

    interface
        module function create_io_reader_t() result(r)
            class(ccpp_io_reader_t), allocatable :: r
        end function create_io_reader_t

        subroutine open_netcdf_file(this, file_path, errcode, errmsg)
            import ccpp_io_reader_t

            class(ccpp_io_reader_t), intent(inout)  :: this
            character(len=*),        intent(in)  :: file_path
            integer,                 intent(out) :: errcode   !Error code
            character(len=*),        intent(out) :: errmsg
        end subroutine open_netcdf_file

        subroutine close_netcdf_file(this, errcode, errmsg)
            import ccpp_io_reader_t
            class(ccpp_io_reader_t), intent(inout)  :: this
            integer,                 intent(out) :: errcode !Error code
            character(len=*),        intent(out) :: errmsg
        end subroutine close_netcdf_file

        subroutine get_netcdf_var_int( this, varname, var, errcode, errmsg)
            import ccpp_io_reader_t
            class(ccpp_io_reader_t), intent(in)  :: this
            character(len=*),        intent(in)  :: varname
            integer, pointer,        intent(out) :: var(..) !Integer variable that file data will be read to.
            integer,                 intent(out) :: errcode !Error code
            character(len=*),        intent(out) :: errmsg  !Error message
        end subroutine get_netcdf_var_int

        subroutine get_netcdf_var_real( this, varname, var, errcode, errmsg)
            use ccpp_kinds, only: kind_phys
            import ccpp_io_reader_t
            class(ccpp_io_reader_t),       intent(in)  :: this
            character(len=*),              intent(in)  :: varname
            real(kind_phys), pointer,      intent(out) :: var(..) !Floating-point variable that file data will be read to.
            integer,                       intent(out) :: errcode !Error code
            character(len=*),              intent(out) :: errmsg  !Error message
        end subroutine get_netcdf_var_real

        subroutine get_netcdf_var_char( this, varname, var, errcode, errmsg)
            import ccpp_io_reader_t
            class(ccpp_io_reader_t),   intent(in)  :: this
            character(len=*),          intent(in)  :: varname
            character(len=:), pointer, intent(out) :: var(..) !Character variable that file data will be read to.
            integer,                   intent(out) :: errcode !Error code
            character(len=*),          intent(out) :: errmsg  !Error message
        end subroutine get_netcdf_var_char
        
    end interface

end module ccpp_io_reader
