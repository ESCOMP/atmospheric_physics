module cldfra_xu_randall_diag
    implicit none

    private
    public :: cldfra_xu_randall_diagnostics_init
    public :: cldfra_xu_randall_diagnostics_run
contains
    !> \section arg_table_cldfra_xu_randall_diagnostics_init Argument Table
    !! \htmlinclude cldfra_xu_randall_diagnostics_init.html
    subroutine cldfra_xu_randall_diagnostics_init( &
            errmsg, errflg)
        use cam_history, only: history_add_field

        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_add_field('cldfra_xu_randall', 'cloud_area_fraction', 'lev', 'avg', 'fraction')

        errmsg = ''
        errflg = 0
    end subroutine cldfra_xu_randall_diagnostics_init

    !> \section arg_table_cldfra_xu_randall_diagnostics_run Argument Table
    !! \htmlinclude cldfra_xu_randall_diagnostics_run.html
    subroutine cldfra_xu_randall_diagnostics_run( &
            cldfra, &
            errmsg, errflg)
        use cam_history, only: history_out_field
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: cldfra(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_out_field('cldfra_xu_randall', cldfra)

        errmsg = ''
        errflg = 0
    end subroutine cldfra_xu_randall_diagnostics_run
end module cldfra_xu_randall_diag
