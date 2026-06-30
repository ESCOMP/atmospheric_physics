!> This module contains diagnostic schemes that are specific to the YSU orographic gravity wave drag scheme,
!> which is part of the MMM physics.
module bl_gwdo_diag
    implicit none

    private
    public :: bl_gwdo_diagnostics_init
    public :: bl_gwdo_diagnostics_run
contains
    !> \section arg_table_bl_gwdo_diagnostics_init Argument Table
    !! \htmlinclude bl_gwdo_diagnostics_init.html
    subroutine bl_gwdo_diagnostics_init( &
            errmsg, errflg)
        use cam_history, only: history_add_field
        use cam_history_support, only: horiz_only

        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        ! The "bl_gwdo" physics scheme makes a distinction between X/Y winds and eastward/northward winds. See
        ! the "bl_gwdo_compat_pre" interstitial scheme for details. However, here we just refer to its diagnostics as
        ! eastward/northward to make them more familiar to CAM-SIMA users.
        call history_add_field('bl_gwdo_dtaux3d', 'tendency_of_eastward_wind_due_to_orographic_gwd', 'lev', 'avg', 'm s-2')
        call history_add_field('bl_gwdo_dtauy3d', 'tendency_of_northward_wind_due_to_orographic_gwd', 'lev', 'avg', 'm s-2')
        call history_add_field('bl_gwdo_dusfcg', 'atmosphere_eastward_stress_due_to_orographic_gwd', horiz_only, 'avg', 'Pa')
        call history_add_field('bl_gwdo_dvsfcg', 'atmosphere_northward_stress_due_to_orographic_gwd', horiz_only, 'avg', 'Pa')

        errmsg = ''
        errflg = 0
    end subroutine bl_gwdo_diagnostics_init

    !> \section arg_table_bl_gwdo_diagnostics_run Argument Table
    !! \htmlinclude bl_gwdo_diagnostics_run.html
    subroutine bl_gwdo_diagnostics_run( &
            dtaux3d, dtauy3d, dusfcg, dvsfcg, &
            errmsg, errflg)
        use cam_history, only: history_out_field
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: dtaux3d(:, :), dtauy3d(:, :), dusfcg(:), dvsfcg(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_out_field('bl_gwdo_dtaux3d', dtaux3d)
        call history_out_field('bl_gwdo_dtauy3d', dtauy3d)
        call history_out_field('bl_gwdo_dusfcg', dusfcg)
        call history_out_field('bl_gwdo_dvsfcg', dvsfcg)

        errmsg = ''
        errflg = 0
    end subroutine bl_gwdo_diagnostics_run
end module bl_gwdo_diag
