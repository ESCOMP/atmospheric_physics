!> This module contains interstitial schemes that are specific to the MMM physics.
module mmm_physics_compat
    implicit none

    private
    public :: mmm_physics_compat_init
    public :: mmm_physics_compat_run
    public :: mmm_physics_accumulate_tendencies_timestep_init
    public :: mmm_physics_accumulate_tendencies_run
    public :: mmm_physics_save_tendencies_timestep_init
    public :: mmm_physics_save_tendencies_run
    public :: mmm_physics_persist_states_init
    public :: mmm_physics_persist_states_timestep_final
    public :: compute_characteristic_grid_length_scale_init
    public :: compute_hydrostatic_upward_air_velocity_at_interface_run
    public :: compute_hydrostatic_upward_air_velocity_run
    public :: geopotential_height_wrt_sfc_at_interface_to_msl_run
    public :: geopotential_height_wrt_sfc_to_msl_run
    public :: lw_heating_rate_to_air_potential_temperature_tendency_run
    public :: sw_heating_rate_to_air_potential_temperature_tendency_run
contains
    !> \section arg_table_mmm_physics_compat_init Argument Table
    !! \htmlinclude mmm_physics_compat_init.html
    pure subroutine mmm_physics_compat_init( &
            initial_run, &
            icloud_bl, isfflx, isftcflx, iz0tlnd, spp_pbl, &
            cycling, restart, &
            xice_threshold, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        logical, intent(in) :: initial_run
        integer, intent(out) :: icloud_bl, isfflx, isftcflx, iz0tlnd, spp_pbl
        logical, intent(out) :: cycling, restart
        real(kind_phys), intent(out) :: xice_threshold
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! TODO:
        ! Should convert some of the following to namelist options after the convection-permitting
        ! physics suite is completed.

        icloud_bl = 0
        isfflx = 1
        isftcflx = 0
        iz0tlnd = 0
        spp_pbl = 0
        cycling = .false. ! There is no such thing as DA cycling in CAM-SIMA. Always false.
        restart = (.not. initial_run) ! Branch and restart runs are translated to restart runs for MMM physics.
        xice_threshold = 0.02_kind_phys
    end subroutine mmm_physics_compat_init

    !> \section arg_table_mmm_physics_compat_run Argument Table
    !! \htmlinclude mmm_physics_compat_run.html
    pure subroutine mmm_physics_compat_run( &
            nstep, &
            dt, &
            theta_curr, theta_prev, qv_curr, qv_prev, &
            icefrac, xice_threshold, landfrac, &
            scheme_name, &
            rthdynten, rqvdynten, &
            xland, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        integer, intent(in) :: nstep
        real(kind_phys), intent(in) :: dt, &
                                       theta_curr(:, :), theta_prev(:, :), qv_curr(:, :), qv_prev(:, :), &
                                       icefrac(:), xice_threshold, landfrac(:)
        character(*), intent(out) :: scheme_name
        real(kind_phys), intent(out) :: rthdynten(:, :), rqvdynten(:, :), &
                                        xland(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        scheme_name = 'mmm_physics_compat'

        if (nstep == 0) then
            rthdynten(:, :) = 0.0_kind_phys
            rqvdynten(:, :) = 0.0_kind_phys
        else
            rthdynten(:, :) = (theta_curr(:, :) - theta_prev(:, :)) / dt
            rqvdynten(:, :) = (qv_curr(:, :) - qv_prev(:, :)) / dt
        end if

        ! For MMM physics, land mask (`xland`) is defined as
        ! * xland = 1.0 for land cells, including sea ice cells.
        ! * xland = 2.0 for water cells.
        where (landfrac >= 0.5_kind_phys .or. &
               icefrac >= xice_threshold)
            xland = 1.0_kind_phys
        elsewhere
            xland = 2.0_kind_phys
        end where
    end subroutine mmm_physics_compat_run

    !> \section arg_table_mmm_physics_accumulate_tendencies_timestep_init Argument Table
    !! \htmlinclude mmm_physics_accumulate_tendencies_timestep_init.html
    pure subroutine mmm_physics_accumulate_tendencies_timestep_init( &
            dudt, dvdt, dtdt, &
            rublten, rucuten, rvblten, rvcuten, &
            rthblten, rthcuten, rthmpten, rthratenlw, rthratensw, &
            rqvblten, rqvcuten, rqvmpten, &
            rqcblten, rncblten, rqccuten, rqcmpten, rncmpten, &
            rqrmpten, rnrmpten, &
            rqiblten, rniblten, rqicuten, rqimpten, rnimpten, &
            rqsblten, rqsmpten, &
            rqgmpten, rngmpten, rvolgmpten, &
            rozblten, &
            rnwfablten, rnwfampten, rnifablten, rnifampten, rnbcablten, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(out) :: dudt(:, :), dvdt(:, :), dtdt(:, :), &
                                        rublten(:, :), rucuten(:, :), rvblten(:, :), rvcuten(:, :), &
                                        rthblten(:, :), rthcuten(:, :), rthmpten(:, :), rthratenlw(:, :), rthratensw(:, :), &
                                        rqvblten(:, :), rqvcuten(:, :), rqvmpten(:, :), &
                                        rqcblten(:, :), rncblten(:, :), rqccuten(:, :), rqcmpten(:, :), rncmpten(:, :), &
                                        rqrmpten(:, :), rnrmpten(:, :), &
                                        rqiblten(:, :), rniblten(:, :), rqicuten(:, :), rqimpten(:, :), rnimpten(:, :), &
                                        rqsblten(:, :), rqsmpten(:, :), &
                                        rqgmpten(:, :), rngmpten(:, :), rvolgmpten(:, :), &
                                        rozblten(:, :), &
                                        rnwfablten(:, :), rnwfampten(:, :), rnifablten(:, :), rnifampten(:, :), rnbcablten(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! Zero out tendencies at the beginning of each time step.

        ! Tendencies for feeding back to CAM-SIMA.
        dudt(:, :) = 0.0_kind_phys
        dvdt(:, :) = 0.0_kind_phys
        dtdt(:, :) = 0.0_kind_phys

        ! Tendencies collected from MMM physics schemes.
        rublten(:, :) = 0.0_kind_phys
        rucuten(:, :) = 0.0_kind_phys
        rvblten(:, :) = 0.0_kind_phys
        rvcuten(:, :) = 0.0_kind_phys

        rthblten(:, :) = 0.0_kind_phys
        rthcuten(:, :) = 0.0_kind_phys
        rthmpten(:, :) = 0.0_kind_phys
        rthratenlw(:, :) = 0.0_kind_phys
        rthratensw(:, :) = 0.0_kind_phys

        rqvblten(:, :) = 0.0_kind_phys
        rqvcuten(:, :) = 0.0_kind_phys
        rqvmpten(:, :) = 0.0_kind_phys

        rqcblten(:, :) = 0.0_kind_phys
        rncblten(:, :) = 0.0_kind_phys
        rqccuten(:, :) = 0.0_kind_phys
        rqcmpten(:, :) = 0.0_kind_phys
        rncmpten(:, :) = 0.0_kind_phys

        rqrmpten(:, :) = 0.0_kind_phys
        rnrmpten(:, :) = 0.0_kind_phys

        rqiblten(:, :) = 0.0_kind_phys
        rniblten(:, :) = 0.0_kind_phys
        rqicuten(:, :) = 0.0_kind_phys
        rqimpten(:, :) = 0.0_kind_phys
        rnimpten(:, :) = 0.0_kind_phys

        rqsblten(:, :) = 0.0_kind_phys
        rqsmpten(:, :) = 0.0_kind_phys

        rqgmpten(:, :) = 0.0_kind_phys
        rngmpten(:, :) = 0.0_kind_phys
        rvolgmpten(:, :) = 0.0_kind_phys

        rozblten(:, :) = 0.0_kind_phys

        rnwfablten(:, :) = 0.0_kind_phys
        rnwfampten(:, :) = 0.0_kind_phys
        rnifablten(:, :) = 0.0_kind_phys
        rnifampten(:, :) = 0.0_kind_phys
        rnbcablten(:, :) = 0.0_kind_phys
    end subroutine mmm_physics_accumulate_tendencies_timestep_init

    !> \section arg_table_mmm_physics_accumulate_tendencies_run Argument Table
    !! \htmlinclude mmm_physics_accumulate_tendencies_run.html
    pure subroutine mmm_physics_accumulate_tendencies_run( &
            dt, exner, &
            dudt, dvdt, dtdt, &
            theta, qv, qc, qr, qi, qs, qg, ozone, &
            nc, nr, ni, ng, nwfa, nifa, nbca, &
            volg, &
            rublten, rucuten, rvblten, rvcuten, &
            rthblten, rthcuten, rthmpten, rthratenlw, rthratensw, &
            rqvblten, rqvcuten, rqvmpten, &
            rqcblten, rncblten, rqccuten, rqcmpten, rncmpten, &
            rqrmpten, rnrmpten, &
            rqiblten, rniblten, rqicuten, rqimpten, rnimpten, &
            rqsblten, rqsmpten, &
            rqgmpten, rngmpten, rvolgmpten, &
            rozblten, &
            rnwfablten, rnwfampten, rnifablten, rnifampten, rnbcablten, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: dt, exner(:, :)
        real(kind_phys), intent(inout) :: dudt(:, :), dvdt(:, :), dtdt(:, :), &
                                          theta(:, :), qv(:, :), qc(:, :), qr(:, :), qi(:, :), qs(:, :), qg(:, :), ozone(:, :), &
                                          nc(:, :), nr(:, :), ni(:, :), ng(:, :), nwfa(:, :), nifa(:, :), nbca(:, :), &
                                          volg(:, :), &
                                          rublten(:, :), rucuten(:, :), rvblten(:, :), rvcuten(:, :), &
                                          rthblten(:, :), rthcuten(:, :), rthmpten(:, :), rthratenlw(:, :), rthratensw(:, :), &
                                          rqvblten(:, :), rqvcuten(:, :), rqvmpten(:, :), &
                                          rqcblten(:, :), rncblten(:, :), rqccuten(:, :), rqcmpten(:, :), rncmpten(:, :), &
                                          rqrmpten(:, :), rnrmpten(:, :), &
                                          rqiblten(:, :), rniblten(:, :), rqicuten(:, :), rqimpten(:, :), rnimpten(:, :), &
                                          rqsblten(:, :), rqsmpten(:, :), &
                                          rqgmpten(:, :), rngmpten(:, :), rvolgmpten(:, :), &
                                          rozblten(:, :), &
                                          rnwfablten(:, :), rnwfampten(:, :), rnifablten(:, :), rnifampten(:, :), rnbcablten(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! Accumulate the states and tendencies for feeding back to CAM-SIMA.
        dudt(:, :) = dudt(:, :) + (rublten(:, :) + rucuten(:, :))
        dvdt(:, :) = dvdt(:, :) + (rvblten(:, :) + rvcuten(:, :))
        dtdt(:, :) = dtdt(:, :) + (rthblten(:, :) + rthcuten(:, :) + rthmpten(:, :) + rthratenlw(:, :) + rthratensw(:, :)) * &
            exner(:, :)

        theta(:, :) = theta(:, :) + (rthblten(:, :) + rthcuten(:, :) + rthmpten(:, :) + rthratenlw(:, :) + rthratensw(:, :)) * dt
        qv(:, :) = qv(:, :) + (rqvblten(:, :) + rqvcuten(:, :) + rqvmpten(:, :)) * dt
        qc(:, :) = qc(:, :) + (rqcblten(:, :) + rqccuten(:, :) + rqcmpten(:, :)) * dt
        qr(:, :) = qr(:, :) + rqrmpten(:, :) * dt
        qi(:, :) = qi(:, :) + (rqiblten(:, :) + rqicuten(:, :) + rqimpten(:, :)) * dt
        qs(:, :) = qs(:, :) + (rqsblten(:, :) + rqsmpten(:, :)) * dt
        qg(:, :) = qg(:, :) + rqgmpten(:, :) * dt
        ozone(:, :) = ozone(:, :) + rozblten(:, :) * dt

        nc(:, :) = nc(:, :) + (rncblten(:, :) + rncmpten(:, :)) * dt
        nr(:, :) = nr(:, :) + rnrmpten(:, :) * dt
        ni(:, :) = ni(:, :) + (rniblten(:, :) + rnimpten(:, :)) * dt
        ng(:, :) = ng(:, :) + rngmpten(:, :) * dt
        nwfa(:, :) = nwfa(:, :) + (rnwfablten(:, :) + rnwfampten(:, :)) * dt
        nifa(:, :) = nifa(:, :) + (rnifablten(:, :) + rnifampten(:, :)) * dt
        nbca(:, :) = nbca(:, :) + rnbcablten(:, :) * dt

        volg(:, :) = volg(:, :) + rvolgmpten(:, :) * dt

        ! After accumulating, zero out the tendencies collected from MMM physics schemes to make this subroutine idempotent,
        ! preventing repeated application of the same tendencies.
        rublten(:, :) = 0.0_kind_phys
        rucuten(:, :) = 0.0_kind_phys
        rvblten(:, :) = 0.0_kind_phys
        rvcuten(:, :) = 0.0_kind_phys

        rthblten(:, :) = 0.0_kind_phys
        rthcuten(:, :) = 0.0_kind_phys
        rthmpten(:, :) = 0.0_kind_phys
        rthratenlw(:, :) = 0.0_kind_phys
        rthratensw(:, :) = 0.0_kind_phys

        rqvblten(:, :) = 0.0_kind_phys
        rqvcuten(:, :) = 0.0_kind_phys
        rqvmpten(:, :) = 0.0_kind_phys

        rqcblten(:, :) = 0.0_kind_phys
        rncblten(:, :) = 0.0_kind_phys
        rqccuten(:, :) = 0.0_kind_phys
        rqcmpten(:, :) = 0.0_kind_phys
        rncmpten(:, :) = 0.0_kind_phys

        rqrmpten(:, :) = 0.0_kind_phys
        rnrmpten(:, :) = 0.0_kind_phys

        rqiblten(:, :) = 0.0_kind_phys
        rniblten(:, :) = 0.0_kind_phys
        rqicuten(:, :) = 0.0_kind_phys
        rqimpten(:, :) = 0.0_kind_phys
        rnimpten(:, :) = 0.0_kind_phys

        rqsblten(:, :) = 0.0_kind_phys
        rqsmpten(:, :) = 0.0_kind_phys

        rqgmpten(:, :) = 0.0_kind_phys
        rngmpten(:, :) = 0.0_kind_phys
        rvolgmpten(:, :) = 0.0_kind_phys

        rozblten(:, :) = 0.0_kind_phys

        rnwfablten(:, :) = 0.0_kind_phys
        rnwfampten(:, :) = 0.0_kind_phys
        rnifablten(:, :) = 0.0_kind_phys
        rnifampten(:, :) = 0.0_kind_phys
        rnbcablten(:, :) = 0.0_kind_phys
    end subroutine mmm_physics_accumulate_tendencies_run

    !> \section arg_table_mmm_physics_save_tendencies_timestep_init Argument Table
    !! \htmlinclude mmm_physics_save_tendencies_timestep_init.html
    pure subroutine mmm_physics_save_tendencies_timestep_init( &
            rublten_p, rucuten_p, rvblten_p, rvcuten_p, &
            rthblten_p, rthcuten_p, rthratenlw_p, rthratensw_p, &
            rqvblten_p, rqvcuten_p, &
            rqcblten_p, rncblten_p, rqccuten_p, &
            rqiblten_p, rniblten_p, rqicuten_p, &
            rqsblten_p, &
            rozblten_p, &
            rnwfablten_p, rnifablten_p, rnbcablten_p, &
            rublten, rucuten, rvblten, rvcuten, &
            rthblten, rthcuten, rthratenlw, rthratensw, &
            rqvblten, rqvcuten, &
            rqcblten, rncblten, rqccuten, &
            rqiblten, rniblten, rqicuten, &
            rqsblten, &
            rozblten, &
            rnwfablten, rnifablten, rnbcablten, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(out) :: rublten_p(:, :), rucuten_p(:, :), rvblten_p(:, :), rvcuten_p(:, :), &
                                        rthblten_p(:, :), rthcuten_p(:, :), rthratenlw_p(:, :), rthratensw_p(:, :), &
                                        rqvblten_p(:, :), rqvcuten_p(:, :), &
                                        rqcblten_p(:, :), rncblten_p(:, :), rqccuten_p(:, :), &
                                        rqiblten_p(:, :), rniblten_p(:, :), rqicuten_p(:, :), &
                                        rqsblten_p(:, :), &
                                        rozblten_p(:, :), &
                                        rnwfablten_p(:, :), rnifablten_p(:, :), rnbcablten_p(:, :)
        real(kind_phys), intent(out) :: rublten(:, :), rucuten(:, :), rvblten(:, :), rvcuten(:, :), &
                                        rthblten(:, :), rthcuten(:, :), rthratenlw(:, :), rthratensw(:, :), &
                                        rqvblten(:, :), rqvcuten(:, :), &
                                        rqcblten(:, :), rncblten(:, :), rqccuten(:, :), &
                                        rqiblten(:, :), rniblten(:, :), rqicuten(:, :), &
                                        rqsblten(:, :), &
                                        rozblten(:, :), &
                                        rnwfablten(:, :), rnifablten(:, :), rnbcablten(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! Zero out tendencies at the beginning of each time step.

        ! "Pending" tendencies that serve as temporary working areas for each MMM physics scheme.
        rublten_p(:, :) = 0.0_kind_phys
        rucuten_p(:, :) = 0.0_kind_phys
        rvblten_p(:, :) = 0.0_kind_phys
        rvcuten_p(:, :) = 0.0_kind_phys

        rthblten_p(:, :) = 0.0_kind_phys
        rthcuten_p(:, :) = 0.0_kind_phys
        rthratenlw_p(:, :) = 0.0_kind_phys
        rthratensw_p(:, :) = 0.0_kind_phys

        rqvblten_p(:, :) = 0.0_kind_phys
        rqvcuten_p(:, :) = 0.0_kind_phys

        rqcblten_p(:, :) = 0.0_kind_phys
        rncblten_p(:, :) = 0.0_kind_phys
        rqccuten_p(:, :) = 0.0_kind_phys

        rqiblten_p(:, :) = 0.0_kind_phys
        rniblten_p(:, :) = 0.0_kind_phys
        rqicuten_p(:, :) = 0.0_kind_phys

        rqsblten_p(:, :) = 0.0_kind_phys

        rozblten_p(:, :) = 0.0_kind_phys

        rnwfablten_p(:, :) = 0.0_kind_phys
        rnifablten_p(:, :) = 0.0_kind_phys
        rnbcablten_p(:, :) = 0.0_kind_phys

        ! "Final" tendencies collected from MMM physics schemes.
        rublten(:, :) = 0.0_kind_phys
        rucuten(:, :) = 0.0_kind_phys
        rvblten(:, :) = 0.0_kind_phys
        rvcuten(:, :) = 0.0_kind_phys

        rthblten(:, :) = 0.0_kind_phys
        rthcuten(:, :) = 0.0_kind_phys
        rthratenlw(:, :) = 0.0_kind_phys
        rthratensw(:, :) = 0.0_kind_phys

        rqvblten(:, :) = 0.0_kind_phys
        rqvcuten(:, :) = 0.0_kind_phys

        rqcblten(:, :) = 0.0_kind_phys
        rncblten(:, :) = 0.0_kind_phys
        rqccuten(:, :) = 0.0_kind_phys

        rqiblten(:, :) = 0.0_kind_phys
        rniblten(:, :) = 0.0_kind_phys
        rqicuten(:, :) = 0.0_kind_phys

        rqsblten(:, :) = 0.0_kind_phys

        rozblten(:, :) = 0.0_kind_phys

        rnwfablten(:, :) = 0.0_kind_phys
        rnifablten(:, :) = 0.0_kind_phys
        rnbcablten(:, :) = 0.0_kind_phys
    end subroutine mmm_physics_save_tendencies_timestep_init

    !> \section arg_table_mmm_physics_save_tendencies_run Argument Table
    !! \htmlinclude mmm_physics_save_tendencies_run.html
    pure subroutine mmm_physics_save_tendencies_run( &
            rublten_p, rucuten_p, rvblten_p, rvcuten_p, &
            rthblten_p, rthcuten_p, rthratenlw_p, rthratensw_p, &
            rqvblten_p, rqvcuten_p, &
            rqcblten_p, rncblten_p, rqccuten_p, &
            rqiblten_p, rniblten_p, rqicuten_p, &
            rqsblten_p, &
            rozblten_p, &
            rnwfablten_p, rnifablten_p, rnbcablten_p, &
            rublten, rucuten, rvblten, rvcuten, &
            rthblten, rthcuten, rthratenlw, rthratensw, &
            rqvblten, rqvcuten, &
            rqcblten, rncblten, rqccuten, &
            rqiblten, rniblten, rqicuten, &
            rqsblten, &
            rozblten, &
            rnwfablten, rnifablten, rnbcablten, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(inout) :: rublten_p(:, :), rucuten_p(:, :), rvblten_p(:, :), rvcuten_p(:, :), &
                                          rthblten_p(:, :), rthcuten_p(:, :), rthratenlw_p(:, :), rthratensw_p(:, :), &
                                          rqvblten_p(:, :), rqvcuten_p(:, :), &
                                          rqcblten_p(:, :), rncblten_p(:, :), rqccuten_p(:, :), &
                                          rqiblten_p(:, :), rniblten_p(:, :), rqicuten_p(:, :), &
                                          rqsblten_p(:, :), &
                                          rozblten_p(:, :), &
                                          rnwfablten_p(:, :), rnifablten_p(:, :), rnbcablten_p(:, :)
        real(kind_phys), intent(inout) :: rublten(:, :), rucuten(:, :), rvblten(:, :), rvcuten(:, :), &
                                          rthblten(:, :), rthcuten(:, :), rthratenlw(:, :), rthratensw(:, :), &
                                          rqvblten(:, :), rqvcuten(:, :), &
                                          rqcblten(:, :), rncblten(:, :), rqccuten(:, :), &
                                          rqiblten(:, :), rniblten(:, :), rqicuten(:, :), &
                                          rqsblten(:, :), &
                                          rozblten(:, :), &
                                          rnwfablten(:, :), rnifablten(:, :), rnbcablten(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! Save pending tendencies to final tendencies.
        !
        ! The reason for this mechanism is to mitigate the inconsistent tendency behaviors across various MMM physics schemes.
        ! Although some schemes correctly accumulate their tendencies with others, there are certain schemes that outright
        ! overwrite the tendencies from others. Instead of wrestling with the problematic schemes individually like MPAS does,
        ! which is very error-prone, a suite-level scheme is introduced to blanketly enforce the correct behavior.
        rublten(:, :) = rublten(:, :) + rublten_p(:, :)
        rucuten(:, :) = rucuten(:, :) + rucuten_p(:, :)
        rvblten(:, :) = rvblten(:, :) + rvblten_p(:, :)
        rvcuten(:, :) = rvcuten(:, :) + rvcuten_p(:, :)

        rthblten(:, :) = rthblten(:, :) + rthblten_p(:, :)
        rthcuten(:, :) = rthcuten(:, :) + rthcuten_p(:, :)
        rthratenlw(:, :) = rthratenlw(:, :) + rthratenlw_p(:, :)
        rthratensw(:, :) = rthratensw(:, :) + rthratensw_p(:, :)

        rqvblten(:, :) = rqvblten(:, :) + rqvblten_p(:, :)
        rqvcuten(:, :) = rqvcuten(:, :) + rqvcuten_p(:, :)

        rqcblten(:, :) = rqcblten(:, :) + rqcblten_p(:, :)
        rncblten(:, :) = rncblten(:, :) + rncblten_p(:, :)
        rqccuten(:, :) = rqccuten(:, :) + rqccuten_p(:, :)

        rqiblten(:, :) = rqiblten(:, :) + rqiblten_p(:, :)
        rniblten(:, :) = rniblten(:, :) + rniblten_p(:, :)
        rqicuten(:, :) = rqicuten(:, :) + rqicuten_p(:, :)

        rqsblten(:, :) = rqsblten(:, :) + rqsblten_p(:, :)

        rozblten(:, :) = rozblten(:, :) + rozblten_p(:, :)

        rnwfablten(:, :) = rnwfablten(:, :) + rnwfablten_p(:, :)
        rnifablten(:, :) = rnifablten(:, :) + rnifablten_p(:, :)
        rnbcablten(:, :) = rnbcablten(:, :) + rnbcablten_p(:, :)

        ! After saving, zero out pending tendencies to make this subroutine idempotent,
        ! preventing repeated application of the same tendencies.
        rublten_p(:, :) = 0.0_kind_phys
        rucuten_p(:, :) = 0.0_kind_phys
        rvblten_p(:, :) = 0.0_kind_phys
        rvcuten_p(:, :) = 0.0_kind_phys

        rthblten_p(:, :) = 0.0_kind_phys
        rthcuten_p(:, :) = 0.0_kind_phys
        rthratenlw_p(:, :) = 0.0_kind_phys
        rthratensw_p(:, :) = 0.0_kind_phys

        rqvblten_p(:, :) = 0.0_kind_phys
        rqvcuten_p(:, :) = 0.0_kind_phys

        rqcblten_p(:, :) = 0.0_kind_phys
        rncblten_p(:, :) = 0.0_kind_phys
        rqccuten_p(:, :) = 0.0_kind_phys

        rqiblten_p(:, :) = 0.0_kind_phys
        rniblten_p(:, :) = 0.0_kind_phys
        rqicuten_p(:, :) = 0.0_kind_phys

        rqsblten_p(:, :) = 0.0_kind_phys

        rozblten_p(:, :) = 0.0_kind_phys

        rnwfablten_p(:, :) = 0.0_kind_phys
        rnifablten_p(:, :) = 0.0_kind_phys
        rnbcablten_p(:, :) = 0.0_kind_phys
    end subroutine mmm_physics_save_tendencies_run

    !> \section arg_table_mmm_physics_persist_states_init Argument Table
    !! \htmlinclude mmm_physics_persist_states_init.html
    pure subroutine mmm_physics_persist_states_init( &
            theta_prev, qv_prev, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(out) :: theta_prev(:, :), qv_prev(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! For remembering the model states from the previous time step. They must be allocated at model initialization
        ! because they need to persist across time steps.
        theta_prev(:, :) = 0.0_kind_phys
        qv_prev(:, :) = 0.0_kind_phys
    end subroutine mmm_physics_persist_states_init

    !> \section arg_table_mmm_physics_persist_states_timestep_final Argument Table
    !! \htmlinclude mmm_physics_persist_states_timestep_final.html
    pure subroutine mmm_physics_persist_states_timestep_final( &
            theta_curr, theta_prev, qv_curr, qv_prev, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: theta_curr(:, :), qv_curr(:, :)
        real(kind_phys), intent(out) :: theta_prev(:, :), qv_prev(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! Remember the model states at this time step. When the next time step comes, they will become
        ! the model states from the previous one.
        theta_prev(:, :) = theta_curr(:, :)
        qv_prev(:, :) = qv_curr(:, :)
    end subroutine mmm_physics_persist_states_timestep_final

    !> \section arg_table_compute_characteristic_grid_length_scale_init Argument Table
    !! \htmlinclude compute_characteristic_grid_length_scale_init.html
    pure subroutine compute_characteristic_grid_length_scale_init( &
            omega, rearth, &
            dx, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: omega(:), rearth
        real(kind_phys), intent(out) :: dx(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! Grid sizes do not change with time. Set them just once at model initialization for better performance.

        ! Compute grid sizes in meters. This is trivial for models with regular grids like WRF,
        ! but not so straightforward for models with unstructured grids like CAM-SIMA. Here, the square root of cell area is used.
        dx(:) = sqrt(omega(:) * (rearth ** 2))
    end subroutine compute_characteristic_grid_length_scale_init

    !> \section arg_table_compute_hydrostatic_upward_air_velocity_at_interface_run Argument Table
    !! \htmlinclude compute_hydrostatic_upward_air_velocity_at_interface_run.html
    pure subroutine compute_hydrostatic_upward_air_velocity_at_interface_run( &
            gravit, omega, rho, &
            wint, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: gravit, omega(:, :), rho(:, :)
        real(kind_phys), intent(out) :: wint(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        integer :: k
        real(kind_phys) :: wmid(size(omega, 1), size(omega, 2))

        ! Compute upward air velocity by hydrostatic equation.
        wmid(:, :) = -omega(:, :) / (rho(:, :) * gravit)

        ! Values at upper and lower boundaries are assumed to be zero.
        wint(:, :) = 0.0_kind_phys

        ! Compute upward air velocity at interface.
        do k = 2, size(omega, 2)
            wint(:, k) = 0.5_kind_phys * (wmid(:, k - 1) + wmid(:, k))
        end do

        errmsg = ''
        errflg = 0
    end subroutine compute_hydrostatic_upward_air_velocity_at_interface_run

    !> \section arg_table_compute_hydrostatic_upward_air_velocity_run Argument Table
    !! \htmlinclude compute_hydrostatic_upward_air_velocity_run.html
    pure subroutine compute_hydrostatic_upward_air_velocity_run( &
            gravit, omega, rho, &
            wmid, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: gravit, omega(:, :), rho(:, :)
        real(kind_phys), intent(out) :: wmid(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        ! Compute upward air velocity by hydrostatic equation.
        wmid(:, :) = -omega(:, :) / (rho(:, :) * gravit)

        errmsg = ''
        errflg = 0
    end subroutine compute_hydrostatic_upward_air_velocity_run

    !> \section arg_table_geopotential_height_wrt_sfc_at_interface_to_msl_run Argument Table
    !! \htmlinclude geopotential_height_wrt_sfc_at_interface_to_msl_run.html
    pure subroutine geopotential_height_wrt_sfc_at_interface_to_msl_run( &
            ncol, &
            gravit, phis, zisfc, &
            zimsl, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        integer, intent(in) :: ncol
        real(kind_phys), intent(in) :: gravit, phis(:), zisfc(:, :)
        real(kind_phys), intent(out) :: zimsl(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        integer :: i

        errmsg = ''
        errflg = 0

        ! Convert geopotential height wrt surface at interface to geopotential height wrt mean sea level at interface,
        ! in accordance with its normal definition.
        do i = 1, ncol
            zimsl(i, :) = phis(i) / gravit + zisfc(i, :)
        end do
    end subroutine geopotential_height_wrt_sfc_at_interface_to_msl_run

    !> \section arg_table_geopotential_height_wrt_sfc_to_msl_run Argument Table
    !! \htmlinclude geopotential_height_wrt_sfc_to_msl_run.html
    pure subroutine geopotential_height_wrt_sfc_to_msl_run( &
            ncol, &
            gravit, phis, zmsfc, &
            zmmsl, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        integer, intent(in) :: ncol
        real(kind_phys), intent(in) :: gravit, phis(:), zmsfc(:, :)
        real(kind_phys), intent(out) :: zmmsl(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        integer :: i

        errmsg = ''
        errflg = 0

        ! Convert geopotential height wrt surface to geopotential height wrt mean sea level, in accordance with
        ! its normal definition.
        do i = 1, ncol
            zmmsl(i, :) = phis(i) / gravit + zmsfc(i, :)
        end do
    end subroutine geopotential_height_wrt_sfc_to_msl_run

    !> \section arg_table_lw_heating_rate_to_air_potential_temperature_tendency_run Argument Table
    !! \htmlinclude lw_heating_rate_to_air_potential_temperature_tendency_run.html
    pure subroutine lw_heating_rate_to_air_potential_temperature_tendency_run( &
            cpairv, exner, lw_heating_rate, &
            rthratenlw, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: cpairv(:, :), exner(:, :), lw_heating_rate(:, :)
        real(kind_phys), intent(out) :: rthratenlw(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        rthratenlw(:, :) = lw_heating_rate(:, :) / cpairv(:, :) / exner(:, :)
    end subroutine lw_heating_rate_to_air_potential_temperature_tendency_run

    !> \section arg_table_sw_heating_rate_to_air_potential_temperature_tendency_run Argument Table
    !! \htmlinclude sw_heating_rate_to_air_potential_temperature_tendency_run.html
    pure subroutine sw_heating_rate_to_air_potential_temperature_tendency_run( &
            cpairv, exner, sw_heating_rate, &
            rthratensw, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: cpairv(:, :), exner(:, :), sw_heating_rate(:, :)
        real(kind_phys), intent(out) :: rthratensw(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        rthratensw(:, :) = sw_heating_rate(:, :) / cpairv(:, :) / exner(:, :)
    end subroutine sw_heating_rate_to_air_potential_temperature_tendency_run
end module mmm_physics_compat
