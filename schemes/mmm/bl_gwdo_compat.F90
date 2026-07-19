!> This module contains interstitial schemes that are specific to the YSU orographic gravity wave drag scheme,
!> which is part of the MMM physics.
module bl_gwdo_compat
    implicit none

    private
    public :: bl_gwdo_compat_pre_run
    public :: bl_gwdo_compat_init
    public :: bl_gwdo_compat_run
contains
    !> \section arg_table_bl_gwdo_compat_pre_run Argument Table
    !! \htmlinclude bl_gwdo_compat_pre_run.html
    pure subroutine bl_gwdo_compat_pre_run( &
            u, v, &
            uproj, vproj, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: u(:, :), v(:, :)
        real(kind_phys), intent(out) :: uproj(:, :), vproj(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! The "bl_gwdo" physics scheme was originally designed to be used with regional models like WRF, where the positive X and
        ! Y directions may not always point to the east and north, respectively. This is no longer the case for global models like
        ! CAM-SIMA.

        ! X and Y winds are just eastward and northward winds, respectively.
        uproj(:, :) = u(:, :)
        vproj(:, :) = v(:, :)
    end subroutine bl_gwdo_compat_pre_run

    !> \section arg_table_bl_gwdo_compat_init Argument Table
    !! \htmlinclude bl_gwdo_compat_init.html
    pure subroutine bl_gwdo_compat_init( &
            sina, cosa, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(out) :: sina(:), cosa(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! These variables do not change with time. Set them just once at model initialization for better performance.

        ! The "bl_gwdo" physics scheme was originally designed to be used with regional models like WRF, where the positive X and
        ! Y directions may not always point to the east and north, respectively. This is no longer the case for global models like
        ! CAM-SIMA.

        ! The angle of rotation from east to X is zero.
        sina(:) = 0.0_kind_phys
        cosa(:) = 1.0_kind_phys
    end subroutine bl_gwdo_compat_init

    !> \section arg_table_bl_gwdo_compat_run Argument Table
    !! \htmlinclude bl_gwdo_compat_run.html
    subroutine bl_gwdo_compat_run( &
            sina, cosa, &
            rublten, rvblten, &
            dtaux3d, dtauy3d, &
            dusfcg, dvsfcg, &
            uproj, vproj, &
            t1, q1, &
            prsi, prsl, prslk, zl, &
            var, oc1, &
            oa2d1, oa2d2, &
            oa2d3, oa2d4, &
            ol2d1, ol2d2, &
            ol2d3, ol2d4, &
            g_, cp_, rd_, rv_, fv_, pi_, &
            dxmeter, deltim, &
            its, ite, kte, kme, &
            errmsg, errflg)
        use bl_gwdo, only: bl_gwdo_run
        use ccpp_kinds, only: kind_phys

        integer, intent(in) :: its, ite, kte, kme
        real(kind_phys), intent(in) :: sina(:), cosa(:), &
                                       uproj(:, :), vproj(:, :), &
                                       t1(:, :), q1(:, :), &
                                       prsi(:, :), prsl(:, :), prslk(:, :), zl(:, :), &
                                       var(:), oc1(:), &
                                       oa2d1(:), oa2d2(:), &
                                       oa2d3(:), oa2d4(:), &
                                       ol2d1(:), ol2d2(:), &
                                       ol2d3(:), ol2d4(:), &
                                       g_, cp_, rd_, rv_, fv_, pi_, &
                                       dxmeter(:), deltim
        real(kind_phys), intent(inout) :: rublten(:, :), rvblten(:, :)
        real(kind_phys), intent(out) :: dtaux3d(:, :), dtauy3d(:, :), &
                                        dusfcg(:), dvsfcg(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! Some schemes of MMM physics expect vertical indexes to be in ascending order from bottom to top of atmosphere,
        ! which is the exact opposite to CAM-SIMA.
        !
        ! For all variables with a vertical dimension, they must be flipped upside down.
        ! This can be achieved by the `associate` construct with array bounds remapping so that the actual array bounds
        ! stays intact elsewhere.
        associate ( &
                rublten_r => rublten(:, size(rublten, 2):1:-1), &
                rvblten_r => rvblten(:, size(rvblten, 2):1:-1), &
                dtaux3d_r => dtaux3d(:, size(dtaux3d, 2):1:-1), &
                dtauy3d_r => dtauy3d(:, size(dtauy3d, 2):1:-1), &
                uproj_r => uproj(:, size(uproj, 2):1:-1), &
                vproj_r => vproj(:, size(vproj, 2):1:-1), &
                t1_r => t1(:, size(t1, 2):1:-1), &
                q1_r => q1(:, size(q1, 2):1:-1), &
                prsi_r => prsi(:, size(prsi, 2):1:-1), &
                prsl_r => prsl(:, size(prsl, 2):1:-1), &
                prslk_r => prslk(:, size(prslk, 2):1:-1), &
                zl_r => zl(:, size(zl, 2):1:-1))
            call bl_gwdo_run( &
                sina, cosa, &
                rublten_r, rvblten_r, &
                dtaux3d_r, dtauy3d_r, &
                dusfcg, dvsfcg, &
                uproj_r, vproj_r, &
                t1_r, q1_r, &
                prsi_r, prsl_r, prslk_r, zl_r, &
                var, oc1, &
                oa2d1, oa2d2, &
                oa2d3, oa2d4, &
                ol2d1, ol2d2, &
                ol2d3, ol2d4, &
                g_, cp_, rd_, rv_, fv_, pi_, &
                dxmeter, deltim, &
                its, ite, kte, kme, &
                errmsg, errflg)
        end associate
    end subroutine bl_gwdo_compat_run
end module bl_gwdo_compat
