!> This module contains interstitial schemes that are specific to MYNN PBL scheme,
!> which is part of MMM physics.
module bl_mynn_compat
    implicit none

    private
    public :: bl_mynn_compat_pre_run
    public :: bl_mynn_compat_init
    public :: bl_mynn_compat_run
    public :: bl_mynn_diagnostics_init
    public :: bl_mynn_diagnostics_run
contains
    !> \section arg_table_bl_mynn_compat_pre_run Argument Table
    !! \htmlinclude bl_mynn_compat_pre_run.html
    pure subroutine bl_mynn_compat_pre_run( &
            itimestep, spp_pbl, &
            restart, &
            rthratenlw, rthratensw, &
            initflag, &
            kpbl, ktop_plume, &
            edmf_a, edmf_w, &
            edmf_qt, edmf_thl, edmf_ent, edmf_qc, &
            sub_thl, sub_sqv, &
            det_thl, det_sqv, &
            edmf_a_dd, edmf_w_dd, &
            edmf_qt_dd, edmf_thl_dd, edmf_ent_dd, edmf_qc_dd, &
            pattern_spp_pbl, rthraten, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        integer, intent(in) :: itimestep, spp_pbl
        logical, intent(in) :: restart
        real(kind_phys), intent(in) :: rthratenlw(:, :), rthratensw(:, :)
        integer, intent(out) :: initflag, &
                                kpbl(:), ktop_plume(:)
        real(kind_phys), intent(out) :: edmf_a(:, :), edmf_w(:, :), &
                                        edmf_qt(:, :), edmf_thl(:, :), edmf_ent(:, :), edmf_qc(:, :), &
                                        sub_thl(:, :), sub_sqv(:, :), &
                                        det_thl(:, :), det_sqv(:, :), &
                                        edmf_a_dd(:, :), edmf_w_dd(:, :), &
                                        edmf_qt_dd(:, :), edmf_thl_dd(:, :), edmf_ent_dd(:, :), edmf_qc_dd(:, :), &
                                        pattern_spp_pbl(:, :), rthraten(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        if (itimestep > 1 .or. restart) then
            initflag = 0
        else
            initflag = 1
        end if

        kpbl(:) = 0
        ktop_plume(:) = 0

        edmf_a(:, :) = 0.0_kind_phys
        edmf_w(:, :) = 0.0_kind_phys
        edmf_qt(:, :) = 0.0_kind_phys
        edmf_thl(:, :) = 0.0_kind_phys
        edmf_ent(:, :) = 0.0_kind_phys
        edmf_qc(:, :) = 0.0_kind_phys
        sub_thl(:, :) = 0.0_kind_phys
        sub_sqv(:, :) = 0.0_kind_phys
        det_thl(:, :) = 0.0_kind_phys
        det_sqv(:, :) = 0.0_kind_phys
        edmf_a_dd(:, :) = 0.0_kind_phys
        edmf_w_dd(:, :) = 0.0_kind_phys
        edmf_qt_dd(:, :) = 0.0_kind_phys
        edmf_thl_dd(:, :) = 0.0_kind_phys
        edmf_ent_dd(:, :) = 0.0_kind_phys
        edmf_qc_dd(:, :) = 0.0_kind_phys

        if (spp_pbl /= 0) then
            errmsg = 'bl_mynn_compat_pre_run: Stochastically perturbed parameterization is not supported'
            errflg = 1

            return
        else
            pattern_spp_pbl(:, :) = 0.0_kind_phys
        end if

        rthraten(:, :) = rthratenlw(:, :) + rthratensw(:, :)

        errmsg = ''
        errflg = 0
    end subroutine bl_mynn_compat_pre_run

    !> \section arg_table_bl_mynn_compat_init Argument Table
    !! \htmlinclude bl_mynn_compat_init.html
    subroutine bl_mynn_compat_init( &
            con_cp, con_cpv, con_cice, con_cliq, con_ep1, con_ep2, con_grav, con_karman, con_p0, &
            con_rd, con_rv, con_svp1, con_svp2, con_svp3, con_svpt0, con_xlf, con_xls, con_xlv, &
            bl_mynn_cloudpdf, bl_mynn_mixlength, bl_mynn_stfunc, &
            bl_mynn_tkeadvect, &
            bl_mynn_tkebudget, &
            bl_mynn_topdown, &
            bl_mynn_edmf, bl_mynn_edmf_dd, bl_mynn_edmf_mom, &
            bl_mynn_edmf_tke, bl_mynn_mixscalars, bl_mynn_output, &
            bl_mynn_cloudmix, bl_mynn_mixqt, bl_mynn_scaleaware, &
            bl_mynn_dheatopt, &
            bl_mynn_closure, &
            flag_qc, flag_qi, flag_qs, flag_qoz, &
            flag_qnc, flag_qni, flag_qnwfa, flag_qnifa, flag_qnbca, &
            qke, qke_adv, &
            tsq, qsq, cov, &
            el_pbl, &
            sh, sm, &
            qc_bl, qi_bl, cldfra_bl, &
            errmsg, errflg)
        use bl_mynn, only: bl_mynn_init
        use ccpp_kinds, only: kind_phys
        use ccpp_scheme_utils, only: ccpp_constituent_index

        real(kind_phys), intent(in) :: con_cp, con_cpv, con_cice, con_cliq, con_ep1, con_ep2, con_grav, con_karman, con_p0, &
                                       con_rd, con_rv, con_svp1, con_svp2, con_svp3, con_svpt0, con_xlf, con_xls, con_xlv
        integer, intent(out) :: bl_mynn_cloudpdf, bl_mynn_mixlength, bl_mynn_stfunc
        logical, intent(out) :: bl_mynn_tkeadvect, &
                                bl_mynn_tkebudget, &
                                bl_mynn_topdown, &
                                bl_mynn_edmf, bl_mynn_edmf_dd, bl_mynn_edmf_mom, &
                                bl_mynn_edmf_tke, bl_mynn_mixscalars, bl_mynn_output, &
                                bl_mynn_cloudmix, bl_mynn_mixqt, bl_mynn_scaleaware, &
                                bl_mynn_dheatopt, &
                                flag_qc, flag_qi, flag_qs, flag_qoz, &
                                flag_qnc, flag_qni, flag_qnwfa, flag_qnifa, flag_qnbca
        real(kind_phys), intent(out) :: bl_mynn_closure, &
                                        qke(:, :), qke_adv(:, :), &
                                        tsq(:, :), qsq(:, :), cov(:, :), &
                                        el_pbl(:, :), &
                                        sh(:, :), sm(:, :), &
                                        qc_bl(:, :), qi_bl(:, :), cldfra_bl(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        integer :: i

        errmsg = ''
        errflg = 0

        ! Initialize constants in MYNN PBL scheme.
        call bl_mynn_init( &
            con_cp, con_cpv, con_cice, con_cliq, con_ep1, con_ep2, con_grav, con_karman, con_p0, &
            con_rd, con_rv, con_svp1, con_svp2, con_svp3, con_svpt0, con_xlf, con_xls, con_xlv, &
            errmsg, errflg)

        if (errflg /= 0) then
            errflg = 1
            errmsg = 'bl_mynn_compat_init: Failed to call "bl_mynn_init"' // new_line('') // &
                'External procedure returned with error: ' // trim(adjustl(errmsg))

            return
        end if

        bl_mynn_cloudpdf = 2
        bl_mynn_mixlength = 2
        bl_mynn_stfunc = 1
        bl_mynn_tkeadvect = .false.
        bl_mynn_tkebudget = .false.
        bl_mynn_topdown = .false.
        bl_mynn_edmf = .true.
        bl_mynn_edmf_dd = .false.
        bl_mynn_edmf_mom = .false.
        bl_mynn_edmf_tke = .false.
        bl_mynn_mixscalars = .true.
        bl_mynn_output = .false.
        bl_mynn_cloudmix = .true.
        bl_mynn_mixqt = .false.
        bl_mynn_scaleaware = .true.
        bl_mynn_dheatopt = .true.
        bl_mynn_closure = 2.5_kind_phys

        flag_qc = .false.
        flag_qi = .false.
        flag_qs = .false.
        flag_qoz = .false.

        call ccpp_constituent_index( &
            'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', i, errflg, errmsg)

        if (errflg == 0) then
            flag_qc = .true.
        end if

        call ccpp_constituent_index( &
            'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', i, errflg, errmsg)

        if (errflg == 0) then
            flag_qi = .true.
        end if

        call ccpp_constituent_index( &
            'snow_mixing_ratio_wrt_moist_air_and_condensed_water', i, errflg, errmsg)

        if (errflg == 0) then
            flag_qs = .true.
        end if

        call ccpp_constituent_index( &
            'ozone_mixing_ratio_wrt_moist_air_and_condensed_water', i, errflg, errmsg)

        if (errflg == 0) then
            flag_qoz = .true.
        end if

        flag_qnc = .false.
        flag_qni = .false.
        flag_qnwfa = .false.
        flag_qnifa = .false.
        flag_qnbca = .false.

        call ccpp_constituent_index( &
            'mass_number_concentration_of_cloud_liquid_water_particles_in_air', i, errflg, errmsg)

        if (errflg == 0) then
            flag_qnc = .true.
        end if

        call ccpp_constituent_index( &
            'mass_number_concentration_of_cloud_ice_water_crystals_in_air', i, errflg, errmsg)

        if (errflg == 0) then
            flag_qni = .true.
        end if

        call ccpp_constituent_index( &
            'mass_number_concentration_of_hygroscopic_aerosols_in_air', i, errflg, errmsg)

        if (errflg == 0) then
            flag_qnwfa = .true.
        end if

        call ccpp_constituent_index( &
            'mass_number_concentration_of_nonhygroscopic_ice_nucleating_aerosols_in_air', i, errflg, errmsg)

        if (errflg == 0) then
            flag_qnifa = .true.
        end if

        call ccpp_constituent_index( &
            'mass_number_concentration_of_hydrophobic_black_carbon_in_air', i, errflg, errmsg)

        if (errflg == 0) then
            flag_qnbca = .true.
        end if

        ! MYNN PBL scheme reuses these variables internally.
        ! As a result, they must be able to persist across time steps.
        qke(:, :) = 0.0_kind_phys
        qke_adv(:, :) = 0.0_kind_phys
        tsq(:, :) = 0.0_kind_phys
        qsq(:, :) = 0.0_kind_phys
        cov(:, :) = 0.0_kind_phys
        el_pbl(:, :) = 0.0_kind_phys
        sh(:, :) = 0.0_kind_phys
        sm(:, :) = 0.0_kind_phys
        qc_bl(:, :) = 0.0_kind_phys
        qi_bl(:, :) = 0.0_kind_phys
        cldfra_bl(:, :) = 0.0_kind_phys

        errmsg = ''
        errflg = 0
    end subroutine bl_mynn_compat_init

    !> \section arg_table_bl_mynn_compat_run Argument Table
    !! \htmlinclude bl_mynn_compat_run.html
    subroutine bl_mynn_compat_run( &
            initflag, restart, cycling, &
            delt, dz, dx, &
            znt, u, v, &
            w, th, sqv_dry, &
            sqc_dry, sqi_dry, sqs_dry, &
            qnc, qni, qnwfa, &
            qnifa, qnbca, qozone, &
            p, exner, rho, &
            tt, xland, ts, &
            qsfc, ps, ust, &
            ch, hfx, qfx, &
            rmol, wspd, uoce, &
            voce, qke, qke_adv, &
            tsq, qsq, cov, &
            rublten, rvblten, rthblten, &
            rqvblten, rqcblten, rqiblten, &
            rqsblten, rqncblten, rqniblten, &
            rqnwfablten, rqnifablten, rqnbcablten, &
            rqozblten, exch_h, exch_m, &
            pblh, kpbl, el_pbl, &
            dqke, qwt, qshear, &
            qbuoy, qdiss, sh, &
            sm, qc_bl, qi_bl, &
            cldfra_bl, icloud_bl, bl_mynn_tkeadvect, &
            bl_mynn_tkebudget, bl_mynn_cloudpdf, bl_mynn_mixlength, &
            bl_mynn_closure, bl_mynn_stfunc, bl_mynn_topdown, &
            bl_mynn_edmf, bl_mynn_edmf_dd, bl_mynn_edmf_mom, &
            bl_mynn_edmf_tke, bl_mynn_mixscalars, bl_mynn_output, &
            bl_mynn_cloudmix, bl_mynn_mixqt, bl_mynn_scaleaware, &
            bl_mynn_dheatopt, edmf_a, edmf_w, &
            edmf_qt, edmf_thl, edmf_ent, &
            edmf_qc, sub_thl, sub_sqv, &
            det_thl, det_sqv, edmf_a_dd, &
            edmf_w_dd, edmf_qt_dd, edmf_thl_dd, &
            edmf_ent_dd, edmf_qc_dd, maxwidth, &
            maxmf, ztop_plume, ktop_plume, &
            spp_pbl, pattern_spp_pbl, rthraten, &
            flag_qc, flag_qi, flag_qs, &
            flag_qnc, flag_qni, flag_qnwfa, &
            flag_qnifa, flag_qnbca, flag_qoz, &
            its, ite, kte, kme, errmsg, errflg)
        use bl_mynn, only: bl_mynn_run
        use ccpp_kinds, only: kind_phys

        integer, intent(in) :: initflag, &
                               icloud_bl, &
                               bl_mynn_cloudpdf, bl_mynn_mixlength, &
                               bl_mynn_stfunc, &
                               spp_pbl, &
                               its, ite, kte, kme
        logical, intent(in) :: restart, cycling, &
                               bl_mynn_tkeadvect, &
                               bl_mynn_tkebudget, &
                               bl_mynn_topdown, &
                               bl_mynn_edmf, bl_mynn_edmf_dd, bl_mynn_edmf_mom, &
                               bl_mynn_edmf_tke, bl_mynn_mixscalars, bl_mynn_output, &
                               bl_mynn_cloudmix, bl_mynn_mixqt, bl_mynn_scaleaware, &
                               bl_mynn_dheatopt, &
                               flag_qc, flag_qi, flag_qs, &
                               flag_qnc, flag_qni, flag_qnwfa, &
                               flag_qnifa, flag_qnbca, flag_qoz
        real(kind_phys), intent(in) :: delt, dz(:, :), dx(:), &
                                       znt(:), u(:, :), v(:, :), &
                                       w(:, :), th(:, :), sqv_dry(:, :), &
                                       sqc_dry(:, :), sqi_dry(:, :), sqs_dry(:, :), &
                                       qnc(:, :), qni(:, :), qnwfa(:, :), &
                                       qnifa(:, :), qnbca(:, :), qozone(:, :), &
                                       p(:, :), exner(:, :), rho(:, :), &
                                       tt(:, :), xland(:), ts(:), &
                                       qsfc(:), ps(:), ust(:), &
                                       ch(:), hfx(:), qfx(:), &
                                       rmol(:), wspd(:), uoce(:), &
                                       voce(:), &
                                       bl_mynn_closure, &
                                       pattern_spp_pbl(:, :), &
                                       rthraten(:, :)
        integer, intent(inout) :: kpbl(:), &
                                  ktop_plume(:)
        real(kind_phys), intent(inout) :: qke(:, :), qke_adv(:, :), &
                                          tsq(:, :), qsq(:, :), cov(:, :), &
                                          rublten(:, :), rvblten(:, :), rthblten(:, :), &
                                          rqvblten(:, :), rqcblten(:, :), rqiblten(:, :), &
                                          rqsblten(:, :), rqncblten(:, :), rqniblten(:, :), &
                                          rqnwfablten(:, :), rqnifablten(:, :), rqnbcablten(:, :), &
                                          rqozblten(:, :), &
                                          pblh(:), el_pbl(:, :), &
                                          sh(:, :), &
                                          sm(:, :), qc_bl(:, :), qi_bl(:, :), &
                                          cldfra_bl(:, :), &
                                          edmf_a(:, :), edmf_w(:, :), &
                                          edmf_qt(:, :), edmf_thl(:, :), edmf_ent(:, :), &
                                          edmf_qc(:, :), sub_thl(:, :), sub_sqv(:, :), &
                                          det_thl(:, :), det_sqv(:, :), edmf_a_dd(:, :), &
                                          edmf_w_dd(:, :), edmf_qt_dd(:, :), edmf_thl_dd(:, :), &
                                          edmf_ent_dd(:, :), edmf_qc_dd(:, :)
        real(kind_phys), intent(out) :: exch_h(:, :), exch_m(:, :), &
                                        dqke(:, :), qwt(:, :), qshear(:, :), &
                                        qbuoy(:, :), qdiss(:, :), &
                                        maxwidth(:), &
                                        maxmf(:), ztop_plume(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        integer, parameter :: kts = 1
        real(kind_phys), allocatable :: sqv(:, :), sqc(:, :), sqi(:, :), sqs(:, :)

        errmsg = ''
        errflg = 0

        allocate( &
            sqv(size(sqv_dry, 1), size(sqv_dry, 2)), &
            sqc(size(sqc_dry, 1), size(sqc_dry, 2)), &
            sqi(size(sqi_dry, 1), size(sqi_dry, 2)), &
            sqs(size(sqs_dry, 1), size(sqs_dry, 2)), &
            errmsg=errmsg, stat=errflg)

        if (errflg /= 0) then
            errmsg = 'bl_mynn_compat_run: Failed to allocate "sqv", "sqc", "sqi", "sqs"' // new_line('') // &
                'Allocation returned with error: ' // trim(adjustl(errmsg))

            return
        end if

        ! Convert constituents from dry to moist basis. These are what MYNN PBL scheme wants.
        sqv(:, :) = sqv_dry(:, :) / (1.0_kind_phys + sqv_dry(:, :))
        sqc(:, :) = sqc_dry(:, :) / (1.0_kind_phys + sqv_dry(:, :))
        sqi(:, :) = sqi_dry(:, :) / (1.0_kind_phys + sqv_dry(:, :))
        sqs(:, :) = sqs_dry(:, :) / (1.0_kind_phys + sqv_dry(:, :))

        ! Some schemes of MMM physics expect vertical indexes to be in ascending order from bottom to top of atmosphere,
        ! which is the exact opposite to CAM-SIMA.
        !
        ! For all variables with a vertical dimension, they must be flipped upside down.
        ! This can be achieved by the `associate` construct with array bounds remapping so that the actual array bounds
        ! stays intact elsewhere.
        associate ( &
            dz_r => dz(:, size(dz, 2):1:-1), &
            u_r => u(:, size(u, 2):1:-1), &
            v_r => v(:, size(v, 2):1:-1), &
            w_r => w(:, size(w, 2):1:-1), &
            th_r => th(:, size(th, 2):1:-1), &
            sqv_r => sqv(:, size(sqv, 2):1:-1), &
            sqc_r => sqc(:, size(sqc, 2):1:-1), &
            sqi_r => sqi(:, size(sqi, 2):1:-1), &
            sqs_r => sqs(:, size(sqs, 2):1:-1), &
            qnc_r => qnc(:, size(qnc, 2):1:-1), &
            qni_r => qni(:, size(qni, 2):1:-1), &
            qnwfa_r => qnwfa(:, size(qnwfa, 2):1:-1), &
            qnifa_r => qnifa(:, size(qnifa, 2):1:-1), &
            qnbca_r => qnbca(:, size(qnbca, 2):1:-1), &
            qozone_r => qozone(:, size(qozone, 2):1:-1), &
            p_r => p(:, size(p, 2):1:-1), &
            exner_r => exner(:, size(exner, 2):1:-1), &
            rho_r => rho(:, size(rho, 2):1:-1), &
            tt_r => tt(:, size(tt, 2):1:-1), &
            qke_r => qke(:, size(qke, 2):1:-1), &
            qke_adv_r => qke_adv(:, size(qke_adv, 2):1:-1), &
            tsq_r => tsq(:, size(tsq, 2):1:-1), &
            qsq_r => qsq(:, size(qsq, 2):1:-1), &
            cov_r => cov(:, size(cov, 2):1:-1), &
            rublten_r => rublten(:, size(rublten, 2):1:-1), &
            rvblten_r => rvblten(:, size(rvblten, 2):1:-1), &
            rthblten_r => rthblten(:, size(rthblten, 2):1:-1), &
            rqvblten_r => rqvblten(:, size(rqvblten, 2):1:-1), &
            rqcblten_r => rqcblten(:, size(rqcblten, 2):1:-1), &
            rqiblten_r => rqiblten(:, size(rqiblten, 2):1:-1), &
            rqsblten_r => rqsblten(:, size(rqsblten, 2):1:-1), &
            rqncblten_r => rqncblten(:, size(rqncblten, 2):1:-1), &
            rqniblten_r => rqniblten(:, size(rqniblten, 2):1:-1), &
            rqnwfablten_r => rqnwfablten(:, size(rqnwfablten, 2):1:-1), &
            rqnifablten_r => rqnifablten(:, size(rqnifablten, 2):1:-1), &
            rqnbcablten_r => rqnbcablten(:, size(rqnbcablten, 2):1:-1), &
            rqozblten_r => rqozblten(:, size(rqozblten, 2):1:-1), &
            exch_h_r => exch_h(:, size(exch_h, 2):1:-1), &
            exch_m_r => exch_m(:, size(exch_m, 2):1:-1), &
            el_pbl_r => el_pbl(:, size(el_pbl, 2):1:-1), &
            dqke_r => dqke(:, size(dqke, 2):1:-1), &
            qwt_r => qwt(:, size(qwt, 2):1:-1), &
            qshear_r => qshear(:, size(qshear, 2):1:-1), &
            qbuoy_r => qbuoy(:, size(qbuoy, 2):1:-1), &
            qdiss_r => qdiss(:, size(qdiss, 2):1:-1), &
            sh_r => sh(:, size(sh, 2):1:-1), &
            sm_r => sm(:, size(sm, 2):1:-1), &
            qc_bl_r => qc_bl(:, size(qc_bl, 2):1:-1), &
            qi_bl_r => qi_bl(:, size(qi_bl, 2):1:-1), &
            cldfra_bl_r => cldfra_bl(:, size(cldfra_bl, 2):1:-1), &
            edmf_a_r => edmf_a(:, size(edmf_a, 2):1:-1), &
            edmf_w_r => edmf_w(:, size(edmf_w, 2):1:-1), &
            edmf_qt_r => edmf_qt(:, size(edmf_qt, 2):1:-1), &
            edmf_thl_r => edmf_thl(:, size(edmf_thl, 2):1:-1), &
            edmf_ent_r => edmf_ent(:, size(edmf_ent, 2):1:-1), &
            edmf_qc_r => edmf_qc(:, size(edmf_qc, 2):1:-1), &
            sub_thl_r => sub_thl(:, size(sub_thl, 2):1:-1), &
            sub_sqv_r => sub_sqv(:, size(sub_sqv, 2):1:-1), &
            det_thl_r => det_thl(:, size(det_thl, 2):1:-1), &
            det_sqv_r => det_sqv(:, size(det_sqv, 2):1:-1), &
            edmf_a_dd_r => edmf_a_dd(:, size(edmf_a_dd, 2):1:-1), &
            edmf_w_dd_r => edmf_w_dd(:, size(edmf_w_dd, 2):1:-1), &
            edmf_qt_dd_r => edmf_qt_dd(:, size(edmf_qt_dd, 2):1:-1), &
            edmf_thl_dd_r => edmf_thl_dd(:, size(edmf_thl_dd, 2):1:-1), &
            edmf_ent_dd_r => edmf_ent_dd(:, size(edmf_ent_dd, 2):1:-1), &
            edmf_qc_dd_r => edmf_qc_dd(:, size(edmf_qc_dd, 2):1:-1), &
            pattern_spp_pbl_r => pattern_spp_pbl(:, size(pattern_spp_pbl, 2):1:-1), &
            rthraten_r => rthraten(:, size(rthraten, 2):1:-1))
            call bl_mynn_run( &
                initflag, restart, cycling, &
                delt, dz_r, dx, &
                znt, u_r, v_r, &
                w_r, th_r, sqv_r, &
                sqc_r, sqi_r, sqs_r, &
                qnc_r, qni_r, qnwfa_r, &
                qnifa_r, qnbca_r, qozone_r, &
                p_r, exner_r, rho_r, &
                tt_r, xland, ts, &
                qsfc, ps, ust, &
                ch, hfx, qfx, &
                rmol, wspd, uoce, &
                voce, qke_r, qke_adv_r, &
                tsq_r, qsq_r, cov_r, &
                rublten_r, rvblten_r, rthblten_r, &
                rqvblten_r, rqcblten_r, rqiblten_r, &
                rqsblten_r, rqncblten_r, rqniblten_r, &
                rqnwfablten_r, rqnifablten_r, rqnbcablten_r, &
                rqozblten_r, exch_h_r, exch_m_r, &
                pblh, kpbl, el_pbl_r, &
                dqke_r, qwt_r, qshear_r, &
                qbuoy_r, qdiss_r, sh_r, &
                sm_r, qc_bl_r, qi_bl_r, &
                cldfra_bl_r, icloud_bl, bl_mynn_tkeadvect, &
                bl_mynn_tkebudget, bl_mynn_cloudpdf, bl_mynn_mixlength, &
                bl_mynn_closure, bl_mynn_stfunc, bl_mynn_topdown, &
                bl_mynn_edmf, bl_mynn_edmf_dd, bl_mynn_edmf_mom, &
                bl_mynn_edmf_tke, bl_mynn_mixscalars, bl_mynn_output, &
                bl_mynn_cloudmix, bl_mynn_mixqt, bl_mynn_scaleaware, &
                bl_mynn_dheatopt, edmf_a_r, edmf_w_r, &
                edmf_qt_r, edmf_thl_r, edmf_ent_r, &
                edmf_qc_r, sub_thl_r, sub_sqv_r, &
                det_thl_r, det_sqv_r, edmf_a_dd_r, &
                edmf_w_dd_r, edmf_qt_dd_r, edmf_thl_dd_r, &
                edmf_ent_dd_r, edmf_qc_dd_r, maxwidth, &
                maxmf, ztop_plume, ktop_plume, &
                spp_pbl, pattern_spp_pbl_r, rthraten_r, &
                flag_qc, flag_qi, flag_qs, &
                flag_qnc, flag_qni, flag_qnwfa, &
                flag_qnifa, flag_qnbca, flag_qoz, &
                its, ite, kts, kte, kme, errmsg, errflg)
        end associate

        if (errflg /= 0) then
            errflg = 1
            errmsg = 'bl_mynn_compat_run: Failed to call "bl_mynn_run"' // new_line('') // &
                'External procedure returned with error: ' // trim(adjustl(errmsg))

            return
        end if

        sqv(:, :) = sqv(:, :) + rqvblten(:, :) * delt
        sqc(:, :) = sqc(:, :) + rqcblten(:, :) * delt
        sqi(:, :) = sqi(:, :) + rqiblten(:, :) * delt
        sqs(:, :) = sqs(:, :) + rqsblten(:, :) * delt

        ! Convert tendencies from moist to dry basis.
        rqvblten(:, :) = (sqv(:, :) / (1.0_kind_phys - sqv(:, :)) - sqv_dry(:, :)) / delt
        rqcblten(:, :) = (sqc(:, :) / (1.0_kind_phys - sqv(:, :)) - sqc_dry(:, :)) / delt
        rqiblten(:, :) = (sqi(:, :) / (1.0_kind_phys - sqv(:, :)) - sqi_dry(:, :)) / delt
        rqsblten(:, :) = (sqs(:, :) / (1.0_kind_phys - sqv(:, :)) - sqs_dry(:, :)) / delt

        errmsg = ''
        errflg = 0
    end subroutine bl_mynn_compat_run

    !> \section arg_table_bl_mynn_diagnostics_init Argument Table
    !! \htmlinclude bl_mynn_diagnostics_init.html
    subroutine bl_mynn_diagnostics_init( &
            errmsg, errflg)
        use cam_history, only: history_add_field
        use cam_history_support, only: horiz_only

        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_add_field( &
            'bl_mynn_qke', &
            'tke_multiplied_by_two', &
            'lev', 'avg', 'm2 s-2')
        call history_add_field( &
            'bl_mynn_qke_adv', &
            'tke_multiplied_by_two_due_to_advection', &
            'lev', 'avg', 'm2 s-2')
        call history_add_field( &
            'bl_mynn_tsq', &
            'variance_of_liquid_water_potential_temperature', &
            'lev', 'avg', 'K2')
        call history_add_field( &
            'bl_mynn_qsq', &
            'variance_of_water_vapor_mixing_ratio_wrt_moist_air', &
            'lev', 'avg', 'kg2 kg-2')
        call history_add_field( &
            'bl_mynn_cov', &
            'covariance_of_liquid_water_potential_temperature_and_water_vapor_mixing_ratio_wrt_moist_air', &
            'lev', 'avg', 'K kg kg-1')
        call history_add_field( &
            'bl_mynn_pblh', &
            'atmosphere_boundary_layer_thickness', &
            horiz_only, 'avg', 'm')
        call history_add_field( &
            'bl_mynn_el_pbl', &
            'turbulent_mixing_length', &
            'lev', 'avg', 'm')
        call history_add_field( &
            'bl_mynn_sh', &
            'stability_function_for_heat', &
            'lev', 'avg', '1')
        call history_add_field( &
            'bl_mynn_sm', &
            'stability_function_for_momentum', &
            'lev', 'avg', '1')
        call history_add_field( &
            'bl_mynn_qc_bl', &
            'subgrid_scale_cloud_liquid_water_mixing_ratio_wrt_moist_air', &
            'lev', 'avg', 'kg kg-1')
        call history_add_field( &
            'bl_mynn_qi_bl', &
            'subgrid_scale_cloud_ice_mixing_ratio_wrt_moist_air', &
            'lev', 'avg', 'kg kg-1')
        call history_add_field( &
            'bl_mynn_cldfra_bl', &
            'subgrid_scale_cloud_area_fraction_in_atmosphere_layer', &
            'lev', 'avg', 'fraction')
        call history_add_field( &
            'bl_mynn_edmf_a', &
            'updraft_area_fraction_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'fraction')
        call history_add_field( &
            'bl_mynn_edmf_w', &
            'updraft_vertical_velocity_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'm s-1')
        call history_add_field( &
            'bl_mynn_edmf_qt', &
            'updraft_total_water_mixing_ratio_wrt_moist_air_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'kg kg-1')
        call history_add_field( &
            'bl_mynn_edmf_thl', &
            'updraft_liquid_water_potential_temperature_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'K')
        call history_add_field( &
            'bl_mynn_edmf_ent', &
            'updraft_entrainment_rate_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'm-1')
        call history_add_field( &
            'bl_mynn_edmf_qc', &
            'updraft_cloud_liquid_water_mixing_ratio_wrt_moist_air_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'kg kg-1')
        call history_add_field( &
            'bl_mynn_sub_thl', &
            'tendency_of_liquid_water_potential_temperature_due_to_subsidence_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'K s-1')
        call history_add_field( &
            'bl_mynn_sub_sqv', &
            'tendency_of_water_vapor_mixing_ratio_wrt_moist_air_due_to_subsidence_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'kg kg-1 s-1')
        call history_add_field( &
            'bl_mynn_det_thl', &
            'tendency_of_liquid_water_potential_temperature_due_to_detrainment_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'K s-1')
        call history_add_field( &
            'bl_mynn_det_sqv', &
            'tendency_of_water_vapor_mixing_ratio_wrt_moist_air_due_to_detrainment_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'kg kg-1 s-1')
        call history_add_field( &
            'bl_mynn_edmf_a_dd', &
            'downdraft_area_fraction_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'fraction')
        call history_add_field( &
            'bl_mynn_edmf_w_dd', &
            'downdraft_vertical_velocity_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'm s-1')
        call history_add_field( &
            'bl_mynn_edmf_qt_dd', &
            'downdraft_total_water_mixing_ratio_wrt_moist_air_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'kg kg-1')
        call history_add_field( &
            'bl_mynn_edmf_thl_dd', &
            'downdraft_liquid_water_potential_temperature_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'K')
        call history_add_field( &
            'bl_mynn_edmf_ent_dd', &
            'downdraft_entrainment_rate_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'm-1')
        call history_add_field( &
            'bl_mynn_edmf_qc_dd', &
            'downdraft_cloud_liquid_water_mixing_ratio_wrt_moist_air_in_scale_aware_tke_moist_edmf_pbl_scheme', &
            'lev', 'avg', 'kg kg-1')
        call history_add_field( &
            'bl_mynn_exch_h', &
            'atmosphere_heat_diffusivity', &
            'lev', 'avg', 'm2 s-1')
        call history_add_field( &
            'bl_mynn_exch_m', &
            'atmosphere_momentum_diffusivity', &
            'lev', 'avg', 'm2 s-1')
        call history_add_field( &
            'bl_mynn_dqke', &
            'total_tendency_of_tke', &
            'lev', 'avg', 'm2 s-3')
        call history_add_field( &
            'bl_mynn_qwt', &
            'tendency_of_tke_due_to_vertical_transport', &
            'lev', 'avg', 'm2 s-3')
        call history_add_field( &
            'bl_mynn_qshear', &
            'tendency_of_tke_due_to_shear', &
            'lev', 'avg', 'm2 s-3')
        call history_add_field( &
            'bl_mynn_qbuoy', &
            'tendency_of_tke_due_to_buoyancy', &
            'lev', 'avg', 'm2 s-3')
        call history_add_field( &
            'bl_mynn_qdiss', &
            'tendency_of_tke_due_to_dissipation', &
            'lev', 'avg', 'm2 s-3')
        call history_add_field( &
            'bl_mynn_maxwidth', &
            'maximum_plume_width_for_mynn_pbl_scheme', &
            horiz_only, 'avg', 'm')
        call history_add_field( &
            'bl_mynn_maxmf', &
            'maximum_mass_flux_for_mynn_pbl_scheme', &
            horiz_only, 'avg', 'm s-1')
        call history_add_field( &
            'bl_mynn_ztop_plume', &
            'height_of_highest_plume_for_mynn_pbl_scheme', &
            horiz_only, 'avg', 'm')

        errmsg = ''
        errflg = 0
    end subroutine bl_mynn_diagnostics_init

    !> \section arg_table_bl_mynn_diagnostics_run Argument Table
    !! \htmlinclude bl_mynn_diagnostics_run.html
    subroutine bl_mynn_diagnostics_run( &
            qke, qke_adv, &
            tsq, qsq, cov, &
            pblh, el_pbl, &
            sh, sm, &
            qc_bl, qi_bl, cldfra_bl, &
            edmf_a, edmf_w, &
            edmf_qt, edmf_thl, edmf_ent, edmf_qc, &
            sub_thl, sub_sqv, &
            det_thl, det_sqv, &
            edmf_a_dd, edmf_w_dd, &
            edmf_qt_dd, edmf_thl_dd, edmf_ent_dd, edmf_qc_dd, &
            exch_h, exch_m, &
            dqke, qwt, qshear, qbuoy, qdiss, &
            maxwidth, maxmf, ztop_plume, &
            errmsg, errflg)
        use cam_history, only: history_out_field
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: qke(:, :), qke_adv(:, :), &
                                       tsq(:, :), qsq(:, :), cov(:, :), &
                                       pblh(:), el_pbl(:, :), &
                                       sh(:, :), sm(:, :), &
                                       qc_bl(:, :), qi_bl(:, :), cldfra_bl(:, :), &
                                       edmf_a(:, :), edmf_w(:, :), &
                                       edmf_qt(:, :), edmf_thl(:, :), edmf_ent(:, :), edmf_qc(:, :), &
                                       sub_thl(:, :), sub_sqv(:, :), &
                                       det_thl(:, :), det_sqv(:, :), &
                                       edmf_a_dd(:, :), edmf_w_dd(:, :), &
                                       edmf_qt_dd(:, :), edmf_thl_dd(:, :), edmf_ent_dd(:, :), edmf_qc_dd(:, :), &
                                       exch_h(:, :), exch_m(:, :), &
                                       dqke(:, :), qwt(:, :), qshear(:, :), qbuoy(:, :), qdiss(:, :), &
                                       maxwidth(:), maxmf(:), ztop_plume(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_out_field('bl_mynn_qke', qke)
        call history_out_field('bl_mynn_qke_adv', qke_adv)
        call history_out_field('bl_mynn_tsq', tsq)
        call history_out_field('bl_mynn_qsq', qsq)
        call history_out_field('bl_mynn_cov', cov)
        call history_out_field('bl_mynn_pblh', pblh)
        call history_out_field('bl_mynn_el_pbl', el_pbl)
        call history_out_field('bl_mynn_sh', sh)
        call history_out_field('bl_mynn_sm', sm)
        call history_out_field('bl_mynn_qc_bl', qc_bl)
        call history_out_field('bl_mynn_qi_bl', qi_bl)
        call history_out_field('bl_mynn_cldfra_bl', cldfra_bl)
        call history_out_field('bl_mynn_edmf_a', edmf_a)
        call history_out_field('bl_mynn_edmf_w', edmf_w)
        call history_out_field('bl_mynn_edmf_qt', edmf_qt)
        call history_out_field('bl_mynn_edmf_thl', edmf_thl)
        call history_out_field('bl_mynn_edmf_ent', edmf_ent)
        call history_out_field('bl_mynn_edmf_qc', edmf_qc)
        call history_out_field('bl_mynn_sub_thl', sub_thl)
        call history_out_field('bl_mynn_sub_sqv', sub_sqv)
        call history_out_field('bl_mynn_det_thl', det_thl)
        call history_out_field('bl_mynn_det_sqv', det_sqv)
        call history_out_field('bl_mynn_edmf_a_dd', edmf_a_dd)
        call history_out_field('bl_mynn_edmf_w_dd', edmf_w_dd)
        call history_out_field('bl_mynn_edmf_qt_dd', edmf_qt_dd)
        call history_out_field('bl_mynn_edmf_thl_dd', edmf_thl_dd)
        call history_out_field('bl_mynn_edmf_ent_dd', edmf_ent_dd)
        call history_out_field('bl_mynn_edmf_qc_dd', edmf_qc_dd)
        call history_out_field('bl_mynn_exch_h', exch_h)
        call history_out_field('bl_mynn_exch_m', exch_m)
        call history_out_field('bl_mynn_dqke', dqke)
        call history_out_field('bl_mynn_qwt', qwt)
        call history_out_field('bl_mynn_qshear', qshear)
        call history_out_field('bl_mynn_qbuoy', qbuoy)
        call history_out_field('bl_mynn_qdiss', qdiss)
        call history_out_field('bl_mynn_maxwidth', maxwidth)
        call history_out_field('bl_mynn_maxmf', maxmf)
        call history_out_field('bl_mynn_ztop_plume', ztop_plume)

        errmsg = ''
        errflg = 0
    end subroutine bl_mynn_diagnostics_run
end module bl_mynn_compat
