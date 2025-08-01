module ccpp_io_reader
    implicit none
    private

    public :: abstract_netcdf_reader_t
    public :: create_netcdf_reader_t

    type, abstract :: abstract_netcdf_reader_t
    contains
        procedure(open_file),    deferred :: open_file
        procedure(close_file),   deferred :: close_file

        !Integer interfaces
        procedure(get_var_int_0d),   deferred :: get_var_int_0d
        procedure(get_var_int_1d),   deferred :: get_var_int_1d
        procedure(get_var_int_2d),   deferred :: get_var_int_2d
        procedure(get_var_int_3d),   deferred :: get_var_int_3d
        procedure(get_var_int_4d),   deferred :: get_var_int_4d
        procedure(get_var_int_5d),   deferred :: get_var_int_5d

        !Real interfaces
        procedure(get_var_real_0d),   deferred :: get_var_real_0d
        procedure(get_var_real_1d),   deferred :: get_var_real_1d
        procedure(get_var_real_2d),   deferred :: get_var_real_2d
        procedure(get_var_real_3d),   deferred :: get_var_real_3d
        procedure(get_var_real_4d),   deferred :: get_var_real_4d
        procedure(get_var_real_5d),   deferred :: get_var_real_5d

        !Character interfaces
        procedure(get_var_char_0d),   deferred :: get_var_char_0d
        procedure(get_var_char_1d),   deferred :: get_var_char_1d
        procedure(get_var_char_2d),   deferred :: get_var_char_2d
        procedure(get_var_char_3d),   deferred :: get_var_char_3d
        procedure(get_var_char_4d),   deferred :: get_var_char_4d
        procedure(get_var_char_5d),   deferred :: get_var_char_5d

        !Generic interface to routines that allocates a variable with data from
        !an opened NetCDF file variable.
        generic :: get_var => get_var_int_0d, get_var_int_1d, get_var_int_2d, get_var_int_3d, get_var_int_4d, get_var_int_5d, &
                              get_var_real_0d, get_var_real_1d, get_var_real_2d, get_var_real_3d, get_var_real_4d, get_var_real_5d,         &
                              get_var_char_0d, get_var_char_1d, get_var_char_2d, get_var_char_3d, get_var_char_4d, get_var_char_5d
    end type abstract_netcdf_reader_t

    interface
        module function create_netcdf_reader_t() result(r)
            class(abstract_netcdf_reader_t), pointer :: r
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

        ! ------------------------------------------------------------------
        ! Integer interfaces
        ! ------------------------------------------------------------------

        subroutine get_var_int_0d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            integer, allocatable,            intent(out) :: var     !Integer variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_int_0d

        subroutine get_var_int_1d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            integer, allocatable,            intent(out) :: var(:)  !Integer variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_int_1d

        subroutine get_var_int_2d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            integer, allocatable,            intent(out) :: var(:,:)!Integer variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_int_2d

        subroutine get_var_int_3d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            integer, allocatable,            intent(out) :: var(:,:,:) !Integer variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg     !Error message
            integer,                         intent(out) :: errcode    !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_int_3d

        subroutine get_var_int_4d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            integer, allocatable,            intent(out) :: var(:,:,:,:) !Integer variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg       !Error message
            integer,                         intent(out) :: errcode      !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_int_4d

        subroutine get_var_int_5d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            integer, allocatable,            intent(out) :: var(:,:,:,:,:) !Integer variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg         !Error message
            integer,                         intent(out) :: errcode        !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_int_5d

        ! ------------------------------------------------------------------
        ! Real interfaces
        ! ------------------------------------------------------------------

        subroutine get_var_real_0d(this, varname, var, errmsg, errcode, start, count)
            use ccpp_kinds, only: kind_phys
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            real(kind_phys), allocatable,    intent(out) :: var     !Floating-point variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_real_0d

        subroutine get_var_real_1d(this, varname, var, errmsg, errcode, start, count)
            use ccpp_kinds, only: kind_phys
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            real(kind_phys), allocatable,    intent(out) :: var(:)  !Floating-point variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_real_1d

        subroutine get_var_real_2d(this, varname, var, errmsg, errcode, start, count)
            use ccpp_kinds, only: kind_phys
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            real(kind_phys), allocatable,    intent(out) :: var(:,:)!Floating-point variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_real_2d

        subroutine get_var_real_3d(this, varname, var, errmsg, errcode, start, count)
            use ccpp_kinds, only: kind_phys
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            real(kind_phys), allocatable,    intent(out) :: var(:,:,:) !Floating-point variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg     !Error message
            integer,                         intent(out) :: errcode    !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_real_3d

        subroutine get_var_real_4d(this, varname, var, errmsg, errcode, start, count)
            use ccpp_kinds, only: kind_phys
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            real(kind_phys), allocatable,    intent(out) :: var(:,:,:,:) !Floating-point variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg       !Error message
            integer,                         intent(out) :: errcode      !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_real_4d

        subroutine get_var_real_5d(this, varname, var, errmsg, errcode, start, count)
            use ccpp_kinds, only: kind_phys
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            real(kind_phys), allocatable,    intent(out) :: var(:,:,:,:,:) !Floating-point variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg         !Error message
            integer,                         intent(out) :: errcode        !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_real_5d

        ! ------------------------------------------------------------------
        ! Character interfaces
        ! ------------------------------------------------------------------

        subroutine get_var_char_0d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            character(len=:), allocatable,   intent(out) :: var     !Character variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_char_0d

        subroutine get_var_char_1d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            character(len=:), allocatable,   intent(out) :: var(:)  !Character variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_char_1d

        subroutine get_var_char_2d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            character(len=:), allocatable,   intent(out) :: var(:,:)!Character variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg  !Error message
            integer,                         intent(out) :: errcode !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_char_2d

        subroutine get_var_char_3d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            character(len=:), allocatable,   intent(out) :: var(:,:,:) !Character variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg     !Error message
            integer,                         intent(out) :: errcode    !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_char_3d

        subroutine get_var_char_4d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            character(len=:), allocatable,   intent(out) :: var(:,:,:,:) !Character variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg       !Error message
            integer,                         intent(out) :: errcode      !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_char_4d

        subroutine get_var_char_5d(this, varname, var, errmsg, errcode, start, count)
            import abstract_netcdf_reader_t

            class(abstract_netcdf_reader_t), intent(in)  :: this
            character(len=*),                intent(in)  :: varname
            character(len=:), allocatable,   intent(out) :: var(:,:,:,:,:) !Character variable that file data will be read to.
            character(len=*),                intent(out) :: errmsg         !Error message
            integer,                         intent(out) :: errcode        !Error code

            !Optional arguments for reading a subset of the variable
            integer, optional,               intent(in)  :: start(:) !Start indices for each dimension
            integer, optional,               intent(in)  :: count(:) !Number of elements to read for each dimension
        end subroutine get_var_char_5d

    end interface

end module ccpp_io_reader
