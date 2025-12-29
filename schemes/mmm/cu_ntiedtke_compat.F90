!> This module contains interstitial schemes that are specific to new Tiedtke cumulus scheme,
!> which is part of MMM physics.
module cu_ntiedtke_compat
    implicit none

    private
    public :: cu_ntiedtke_compat_pre_run
    public :: cu_ntiedtke_compat_init
    public :: cu_ntiedtke_compat_run
    public :: cu_ntiedtke_diagnostics_init
    public :: cu_ntiedtke_diagnostics_run
contains
    !> \section arg_table_cu_ntiedtke_compat_pre_run Argument Table
    !! \htmlinclude cu_ntiedtke_compat_pre_run.html
    subroutine cu_ntiedtke_compat_pre_run( &
            cflx, exner, landfrac, &
            rthdynten, rthblten, rthratenlw, rthratensw, &
            rqvdynten, rqvblten, &
            lndj, &
            ptf, pqvf, evap, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys
        use ccpp_scheme_utils, only: ccpp_constituent_index

        real(kind_phys), intent(in) :: cflx(:, :), exner(:, :), landfrac(:), &
                                       rthdynten(:, :), rthblten(:, :), rthratenlw(:, :), rthratensw(:, :), &
                                       rqvdynten(:, :), rqvblten(:, :)
        integer, intent(out) :: lndj(:)
        real(kind_phys), intent(out) :: ptf(:, :), pqvf(:, :), evap(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        integer :: water_vapor_mixing_ratio_index

        where (landfrac >= 0.5_kind_phys)
            lndj = 1
        elsewhere
            lndj = 0
        end where

        ptf(:, :) = (rthdynten(:, :) + rthblten(:, :) + rthratenlw(:, :) + rthratensw(:, :)) * exner(:, :)
        pqvf(:, :) = rqvdynten(:, :) + rqvblten(:, :)
        evap(:) = 0.0_kind_phys

        call ccpp_constituent_index( &
            'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', water_vapor_mixing_ratio_index, errflg, errmsg)

        if (errflg /= 0 .or. &
            water_vapor_mixing_ratio_index < lbound(cflx, 2) .or. water_vapor_mixing_ratio_index > ubound(cflx, 2)) then
            errmsg = 'Failed to find desired constituent flux from cflx'
            errflg = 1

            return
        end if

        evap(:) = cflx(:, water_vapor_mixing_ratio_index)

        errmsg = ''
        errflg = 0
    end subroutine cu_ntiedtke_compat_pre_run

    !> \section arg_table_cu_ntiedtke_compat_init Argument Table
    !! \htmlinclude cu_ntiedtke_compat_init.html
    subroutine cu_ntiedtke_compat_init( &
            con_cp, con_rd, con_rv, con_xlv, con_xls, con_xlf, con_grav, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys
        use cu_ntiedtke, only: cu_ntiedtke_init

        real(kind_phys), intent(in) :: con_cp, con_rd, con_rv, con_xlv, con_xls, con_xlf, con_grav
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call cu_ntiedtke_init( &
            con_cp, con_rd, con_rv, con_xlv, con_xls, con_xlf, con_grav, &
            errmsg, errflg)

        errmsg = ''
        errflg = 0
    end subroutine cu_ntiedtke_compat_init

    !> \section arg_table_cu_ntiedtke_compat_run Argument Table
    !! \htmlinclude cu_ntiedtke_compat_run.html
    subroutine cu_ntiedtke_compat_run( &
            pu, pv, pt, pqv, pqc, pqi, &
            pqvf, ptf, poz, pzz, pomg, &
            pap, paph, evap, hfx, zprecc, lndj, lq, km, km1, dt, dx, &
            exner, &
            rucuten, rvcuten, rthcuten, rqvcuten, rqccuten, rqicuten, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys
        use cu_ntiedtke, only: cu_ntiedtke_run

        integer, intent(in) :: lndj(:), lq, km, km1
        real(kind_phys), intent(in) :: pu(:, :), pv(:, :), pt(:, :), pqv(:, :), pqc(:, :), pqi(:, :), &
                                       pqvf(:, :), ptf(:, :), poz(:, :), pzz(:, :), pomg(:, :), &
                                       pap(:, :), paph(:, :), evap(:), hfx(:), &
                                       dt, dx(:), &
                                       exner(:, :)
        real(kind_phys), intent(out) :: zprecc(:), &
                                        rucuten(:, :), rvcuten(:, :), &
                                        rthcuten(:, :), &
                                        rqvcuten(:, :), rqccuten(:, :), rqicuten(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        real(kind_phys), allocatable :: pu_local(:, :), pv_local(:, :), &
                                        pt_local(:, :), &
                                        pqv_local(:, :), pqc_local(:, :), pqi_local(:, :)

        zprecc(:) = 0.0_kind_phys

        rucuten(:, :) = 0.0_kind_phys
        rvcuten(:, :) = 0.0_kind_phys
        rthcuten(:, :) = 0.0_kind_phys
        rqvcuten(:, :) = 0.0_kind_phys
        rqccuten(:, :) = 0.0_kind_phys
        rqicuten(:, :) = 0.0_kind_phys

        ! The "cu_ntiedtke" physics scheme modifies model states directly, which is not ideal.
        ! Make local copies of the model states, pass them to the physics scheme, and compute the tendencies instead.

        allocate(pu_local, source=pu, errmsg=errmsg, stat=errflg)

        if (errflg /= 0) then
            return
        end if

        allocate(pv_local, source=pv, errmsg=errmsg, stat=errflg)

        if (errflg /= 0) then
            return
        end if

        allocate(pt_local, source=pt, errmsg=errmsg, stat=errflg)

        if (errflg /= 0) then
            return
        end if

        allocate(pqv_local, source=pqv, errmsg=errmsg, stat=errflg)

        if (errflg /= 0) then
            return
        end if

        allocate(pqc_local, source=pqc, errmsg=errmsg, stat=errflg)

        if (errflg /= 0) then
            return
        end if

        allocate(pqi_local, source=pqi, errmsg=errmsg, stat=errflg)

        if (errflg /= 0) then
            return
        end if

        call cu_ntiedtke_run( &
            pu_local, pv_local, pt_local, pqv_local, pqc_local, pqi_local, &
            pqvf, ptf, poz, pzz, pomg, &
            pap, paph, evap, hfx, zprecc, lndj, lq, km, km1, dt, dx, &
            errmsg, errflg)

        zprecc(:) = zprecc(:) * 0.001_kind_phys ! Convert from mm to m.

        rucuten(:, :) = (pu_local(:, :) - pu(:, :)) / dt
        rvcuten(:, :) = (pv_local(:, :) - pv(:, :)) / dt
        rthcuten(:, :) = (pt_local(:, :) - pt(:, :)) / exner(:, :) / dt
        rqvcuten(:, :) = (pqv_local(:, :) - pqv(:, :)) / dt
        rqccuten(:, :) = (pqc_local(:, :) - pqc(:, :)) / dt
        rqicuten(:, :) = (pqi_local(:, :) - pqi(:, :)) / dt

        errmsg = ''
        errflg = 0
    end subroutine cu_ntiedtke_compat_run

    !> \section arg_table_cu_ntiedtke_diagnostics_init Argument Table
    !! \htmlinclude cu_ntiedtke_diagnostics_init.html
    subroutine cu_ntiedtke_diagnostics_init( &
            errmsg, errflg)
        use cam_history, only: history_add_field
        use cam_history_support, only: horiz_only

        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_add_field('cu_ntiedtke_zprecc', &
            'lwe_thickness_of_convective_precipitation_amount', horiz_only, 'avg', 'm')
        call history_add_field('cu_ntiedtke_rucuten', &
            'tendency_of_eastward_wind_due_to_convection', 'lev', 'avg', 'm s-2')
        call history_add_field('cu_ntiedtke_rvcuten', &
            'tendency_of_northward_wind_due_to_convection', 'lev', 'avg', 'm s-2')
        call history_add_field('cu_ntiedtke_rthcuten', &
            'tendency_of_air_potential_temperature_due_to_convection', 'lev', 'avg', 'K s-1')
        call history_add_field('cu_ntiedtke_rqvcuten', &
            'tendency_of_water_vapor_mixing_ratio_wrt_dry_air_due_to_convection', 'lev', 'avg', 'kg kg-1 s-1')
        call history_add_field('cu_ntiedtke_rqccuten', &
            'tendency_of_cloud_liquid_water_mixing_ratio_wrt_dry_air_due_to_convection', 'lev', 'avg', 'kg kg-1 s-1')
        call history_add_field('cu_ntiedtke_rqicuten', &
            'tendency_of_cloud_ice_mixing_ratio_wrt_dry_air_due_to_convection', 'lev', 'avg', 'kg kg-1 s-1')

        errmsg = ''
        errflg = 0
    end subroutine cu_ntiedtke_diagnostics_init

    !> \section arg_table_cu_ntiedtke_diagnostics_run Argument Table
    !! \htmlinclude cu_ntiedtke_diagnostics_run.html
    subroutine cu_ntiedtke_diagnostics_run( &
            zprecc, &
            rucuten, rvcuten, rthcuten, rqvcuten, rqccuten, rqicuten, &
            errmsg, errflg)
        use cam_history, only: history_out_field
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: zprecc(:), &
                                       rucuten(:, :), rvcuten(:, :), &
                                       rthcuten(:, :), &
                                       rqvcuten(:, :), rqccuten(:, :), rqicuten(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_out_field('cu_ntiedtke_zprecc', zprecc)
        call history_out_field('cu_ntiedtke_rucuten', rucuten)
        call history_out_field('cu_ntiedtke_rvcuten', rvcuten)
        call history_out_field('cu_ntiedtke_rthcuten', rthcuten)
        call history_out_field('cu_ntiedtke_rqvcuten', rqvcuten)
        call history_out_field('cu_ntiedtke_rqccuten', rqccuten)
        call history_out_field('cu_ntiedtke_rqicuten', rqicuten)

        errmsg = ''
        errflg = 0
    end subroutine cu_ntiedtke_diagnostics_run
end module cu_ntiedtke_compat
