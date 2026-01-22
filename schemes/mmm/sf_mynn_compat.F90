!> This module contains interstitial schemes that are specific to MYNN surface layer scheme,
!> which is part of MMM physics.
module sf_mynn_compat
    use ccpp_kinds, only: kind_phys

    implicit none

    private
    public :: sf_mynn_compat_pre_run
    public :: sf_mynn_compat_init
    public :: sf_mynn_compat_run
    public :: sf_mynn_diagnostics_init
    public :: sf_mynn_diagnostics_run

    ! This threshold is hardcoded to the same value as in MMM physics.
    ! It is named `xice_threshold` there.
    real(kind_phys), parameter :: sea_ice_area_fraction_threshold = 0.02_kind_phys
contains
    !> \section arg_table_sf_mynn_compat_pre_run Argument Table
    !! \htmlinclude sf_mynn_compat_pre_run.html
    subroutine sf_mynn_compat_pre_run( &
            itimestep, &
            u, v, t, qv, p, dz, rho, &
            icefrac, landfrac, snowhice, snowhland, &
            u1d, v1d, t1d, qv1d, p1d, dz8w1d, rho1d, &
            u1d2, v1d2, dz2w1d, &
            chs, chs2, cqs2, cpm, rmol, &
            znt, ust, zol, mol, regime, psim, &
            psih, qfx, &
            flhc, flqc, snowh, qgh, qsfc, &
            gz1oz0, wspd, br, svp1, svp2, &
            svp3, svpt0, qcg, &
            spp_pbl, rstoch1d, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        integer, intent(in) :: itimestep
        real(kind_phys), intent(in) :: u(:, :), v(:, :), t(:, :), qv(:, :), p(:, :), dz(:, :), rho(:, :), &
                                       icefrac(:), landfrac(:), snowhice(:), snowhland(:)
        logical, intent(out) :: spp_pbl
        real(kind_phys), intent(out) :: u1d(:), v1d(:), t1d(:), qv1d(:), p1d(:), dz8w1d(:), rho1d(:), &
                                        u1d2(:), v1d2(:), dz2w1d(:), &
                                        chs(:), chs2(:), cqs2(:), cpm(:), rmol(:), &
                                        znt(:), ust(:), zol(:), mol(:), regime(:), psim(:), &
                                        psih(:), qfx(:), &
                                        flhc(:), flqc(:), snowh(:), qgh(:), qsfc(:), &
                                        gz1oz0(:), wspd(:), br(:), svp1, svp2, &
                                        svp3, svpt0, qcg(:), &
                                        rstoch1d(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        ! Provide first guess at the first time step.
        if (itimestep == 1) then
            ust(:) = max(0.04_kind_phys * sqrt(u(:, 1) ** 2 + v(:, 1) ** 2), 0.001_kind_phys)
            mol(:) = 0.0_kind_phys
            qsfc(:) = qv(:, 1) ! Note that this `qv` is wet.
        end if

        ! In CAM-SIMA, the first vertical index is at top of atmosphere.
        ! The last one is at bottom of atmosphere, which is what we want here.

        u1d(:) = u(:, size(u, 2))
        v1d(:) = v(:, size(v, 2))
        t1d(:) = t(:, size(t, 2))
        qv1d(:) = qv(:, size(qv, 2))
        p1d(:) = p(:, size(p, 2))
        dz8w1d(:) = dz(:, size(dz, 2))
        rho1d(:) = rho(:, size(rho, 2))

        u1d2(:) = u(:, size(u, 2) - 1)
        v1d2(:) = v(:, size(v, 2) - 1)
        dz2w1d(:) = dz(:, size(dz, 2) - 1)

        chs(:) = 0.0_kind_phys
        chs2(:) = 0.0_kind_phys
        cqs2(:) = 0.0_kind_phys
        cpm(:) = 0.0_kind_phys
        rmol(:) = 0.0_kind_phys

        znt(:) = 0.0_kind_phys
        zol(:) = 0.0_kind_phys
        regime(:) = 0.0_kind_phys
        psim(:) = 0.0_kind_phys

        psih(:) = 0.0_kind_phys
        qfx(:) = 0.0_kind_phys

        flhc(:) = 0.0_kind_phys
        flqc(:) = 0.0_kind_phys
        snowh(:) = 0.0_kind_phys

        where (landfrac >= 0.5_kind_phys)
            snowh = snowhland
        end where

        where (icefrac >= sea_ice_area_fraction_threshold)
            snowh = snowhice
        end where

        qgh(:) = 0.0_kind_phys

        gz1oz0(:) = 0.0_kind_phys
        wspd(:) = 0.0_kind_phys
        br(:) = 0.0_kind_phys

        ! Constants in equation 10 from Bolton (1980). See
        ! doi:10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2.
        svp1 = 6.112_kind_phys
        svp2 = 17.67_kind_phys
        svp3 = 29.65_kind_phys
        svpt0 = 273.15_kind_phys

        qcg(:) = 0.0_kind_phys ! Not used but still appear in the argument list...

        spp_pbl = .false.
        rstoch1d(:) = 0.0_kind_phys

        errmsg = ''
        errflg = 0
    end subroutine sf_mynn_compat_pre_run

    !> \section arg_table_sf_mynn_compat_init Argument Table
    !! \htmlinclude sf_mynn_compat_init.html
    subroutine sf_mynn_compat_init( &
            ust, mol, qsfc, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys
        use sf_mynn, only: sf_mynn_init

        real(kind_phys), intent(out) :: ust(:), mol(:), qsfc(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        ! Precompute lookup tables in MYNN surface layer scheme.
        call sf_mynn_init(errmsg, errflg)

        if (errflg /= 0) then
            return
        end if

        ! MYNN surface layer scheme takes time averages of these variables internally.
        ! As a result, they must be able to persist across time steps.
        ust(:) = 0.0_kind_phys
        mol(:) = 0.0_kind_phys
        qsfc(:) = 0.0_kind_phys

        errmsg = ''
        errflg = 0
    end subroutine sf_mynn_compat_init

    !> \section arg_table_sf_mynn_compat_run Argument Table
    !! \htmlinclude sf_mynn_compat_run.html
    subroutine sf_mynn_compat_run( &
            ncol, cflx, icefrac, sst, &
            u1d, v1d, t1d, qv1d, p1d, dz8w1d, rho1d, &
            u1d2, v1d2, dz2w1d, cp, g, rovcp, r, xlv, &
            psfcpa, chs, chs2, cqs2, cpm, pblh, rmol, &
            znt, ust, mavail, zol, mol, regime, psim, &
            psih, xland, hfx, qfx, tsk, u10, v10, th2, &
            t2, q2, flhc, flqc, snowh, qgh, qsfc, lh, &
            gz1oz0, wspd, br, isfflx, dx, svp1, svp2, &
            svp3, svpt0, ep1, ep2, karman, qcg, &
            itimestep, wstar, qstar, ustm, ck, cka, &
            cd, cda, spp_pbl, rstoch1d, isftcflx, &
            iz0tlnd, its, ite, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys
        use ccpp_scheme_utils, only: ccpp_constituent_index
        use sf_mynn, only: sf_mynn_run

        ! Typical threshold between sea surface temperature and ice surface temperature. See
        ! Algorithm Theoretical Basis Document (ATBD) for the MODIS Snow and Sea Ice-Mapping Algorithms,
        ! Section 4.4.3 Ice Surface Temperature (IST) Algorithm.
        real(kind_phys), parameter :: sea_ice_temperature_threshold = 271.4_kind_phys

        integer, intent(in) :: ncol, &
                               isfflx, &
                               itimestep, &
                               isftcflx, &
                               iz0tlnd, its, ite
        logical, intent(in) :: spp_pbl
        real(kind_phys), intent(in) :: icefrac(:), sst(:), &
                                       u1d(:), v1d(:), t1d(:), qv1d(:), p1d(:), dz8w1d(:), rho1d(:), &
                                       u1d2(:), v1d2(:), dz2w1d(:), cp, g, rovcp, r, xlv, &
                                       psfcpa(:), pblh(:), &
                                       mavail(:), &
                                       xland(:), tsk(:), &
                                       snowh(:), &
                                       dx(:), svp1, svp2, &
                                       svp3, svpt0, ep1, ep2, karman, qcg(:), &
                                       rstoch1d(:)
        real(kind_phys), intent(inout) :: cflx(:, :), &
                                          chs(:), chs2(:), cqs2(:), cpm(:), rmol(:), &
                                          znt(:), ust(:), zol(:), mol(:), regime(:), psim(:), &
                                          psih(:), hfx(:), qfx(:), &
                                          flhc(:), flqc(:), qgh(:), qsfc(:), lh(:), &
                                          gz1oz0(:), wspd(:), br(:)
        real(kind_phys), intent(out) :: u10(:), v10(:), th2(:), &
                                        t2(:), q2(:), &
                                        wstar(:), qstar(:), ustm(:), ck(:), cka(:), &
                                        cd(:), cda(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        integer :: water_vapor_mixing_ratio_index
        real(kind_phys), allocatable :: ch(:)

        ! Special handling for sea ice cells like in MMM physics.
        logical, allocatable :: mask_sea_ice_cell(:)
        real(kind_phys), allocatable :: mavail_sea(:), &
                                        xland_sea(:), tsk_sea(:)
        real(kind_phys), allocatable :: chs_sea(:), chs2_sea(:), cqs2_sea(:), cpm_sea(:), rmol_sea(:), &
                                        znt_sea(:), ust_sea(:), zol_sea(:), mol_sea(:), regime_sea(:), psim_sea(:), &
                                        psih_sea(:), hfx_sea(:), qfx_sea(:), &
                                        flhc_sea(:), flqc_sea(:), qgh_sea(:), qsfc_sea(:), lh_sea(:), &
                                        gz1oz0_sea(:), wspd_sea(:), br_sea(:), &
                                        ch_sea(:)
        real(kind_phys), allocatable :: u10_sea(:), v10_sea(:), th2_sea(:), &
                                        t2_sea(:), q2_sea(:), &
                                        wstar_sea(:), qstar_sea(:), ustm_sea(:), ck_sea(:), cka_sea(:), &
                                        cd_sea(:), cda_sea(:)

        ! `ch` is duplicate of `chs`, but for unknown reasons it is passed separately in the argument list
        ! to MYNN surface layer scheme.
        allocate(ch(ncol), errmsg=errmsg, stat=errflg)

        if (errflg /= 0) then
            return
        end if

        ch(:) = chs(:)

        qfx(:) = 0.0_kind_phys

        call ccpp_constituent_index( &
            'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', water_vapor_mixing_ratio_index, errflg, errmsg)

        if (errflg /= 0 .or. &
            water_vapor_mixing_ratio_index < lbound(cflx, 2) .or. water_vapor_mixing_ratio_index > ubound(cflx, 2)) then
            errmsg = 'Failed to find desired constituent flux from cflx'
            errflg = 1

            return
        end if

        qfx(:) = cflx(:, water_vapor_mixing_ratio_index)

        allocate(mask_sea_ice_cell(ncol), errmsg=errmsg, stat=errflg)

        if (errflg /= 0) then
            return
        end if

        mask_sea_ice_cell(:) = (icefrac >= sea_ice_area_fraction_threshold)

        ! If there are sea ice cells, make local copies of the variables in preparation for the second pass.
        if (any(mask_sea_ice_cell)) then
            allocate(mavail_sea(ncol), xland_sea(ncol), tsk_sea(ncol), &
                errmsg=errmsg, stat=errflg)

            if (errflg /= 0) then
                return
            end if

            mavail_sea(:) = mavail(:)
            xland_sea(:) = xland(:)
            tsk_sea(:) = tsk(:)

            allocate(chs_sea(ncol), chs2_sea(ncol), cqs2_sea(ncol), cpm_sea(ncol), rmol_sea(ncol), &
                znt_sea(ncol), ust_sea(ncol), zol_sea(ncol), mol_sea(ncol), regime_sea(ncol), psim_sea(ncol), &
                psih_sea(ncol), hfx_sea(ncol), qfx_sea(ncol), &
                flhc_sea(ncol), flqc_sea(ncol), qgh_sea(ncol), qsfc_sea(ncol), lh_sea(ncol), &
                gz1oz0_sea(ncol), wspd_sea(ncol), br_sea(ncol), &
                ch_sea(ncol), &
                errmsg=errmsg, stat=errflg)

            if (errflg /= 0) then
                return
            end if

            chs_sea(:) = chs(:)
            chs2_sea(:) = chs2(:)
            cqs2_sea(:) = cqs2(:)
            cpm_sea(:) = cpm(:)
            rmol_sea(:) = rmol(:)
            znt_sea(:) = znt(:)
            ust_sea(:) = ust(:)
            zol_sea(:) = zol(:)
            mol_sea(:) = mol(:)
            regime_sea(:) = regime(:)
            psim_sea(:) = psim(:)
            psih_sea(:) = psih(:)
            hfx_sea(:) = hfx(:)
            qfx_sea(:) = qfx(:)
            flhc_sea(:) = flhc(:)
            flqc_sea(:) = flqc(:)
            qgh_sea(:) = qgh(:)
            qsfc_sea(:) = qsfc(:)
            lh_sea(:) = lh(:)
            gz1oz0_sea(:) = gz1oz0(:)
            wspd_sea(:) = wspd(:)
            br_sea(:) = br(:)
            ch_sea(:) = ch(:)

            allocate(u10_sea(ncol), v10_sea(ncol), th2_sea(ncol), &
                t2_sea(ncol), q2_sea(ncol), &
                wstar_sea(ncol), qstar_sea(ncol), ustm_sea(ncol), ck_sea(ncol), cka_sea(ncol), &
                cd_sea(ncol), cda_sea(ncol), &
                errmsg=errmsg, stat=errflg)

            if (errflg /= 0) then
                return
            end if

            u10_sea(:) = u10(:)
            v10_sea(:) = v10(:)
            th2_sea(:) = th2(:)
            t2_sea(:) = t2(:)
            q2_sea(:) = q2(:)
            wstar_sea(:) = wstar(:)
            qstar_sea(:) = qstar(:)
            ustm_sea(:) = ustm(:)
            ck_sea(:) = ck(:)
            cka_sea(:) = cka(:)
            cd_sea(:) = cd(:)
            cda_sea(:) = cda(:)

            where (mask_sea_ice_cell)
                ! Set surface moisture availability to maximum.
                mavail_sea = 1.0_kind_phys

                ! Impose minimum skin temperature.
                tsk_sea = max(sea_ice_temperature_threshold, sst)

                ! Treat as water cells.
                xland_sea = 2.0_kind_phys

                ! Set surface roughness length to calm water.
                znt_sea = 0.0001_kind_phys
            end where
        end if

        ! First pass for all cells.
        call sf_mynn_run( &
            u1d, v1d, t1d, qv1d, p1d, dz8w1d, rho1d, &
            u1d2, v1d2, dz2w1d, cp, g, rovcp, r, xlv, &
            psfcpa, chs, chs2, cqs2, cpm, pblh, rmol, &
            znt, ust, mavail, zol, mol, regime, psim, &
            psih, xland, hfx, qfx, tsk, u10, v10, th2, &
            t2, q2, flhc, flqc, snowh, qgh, qsfc, lh, &
            gz1oz0, wspd, br, isfflx, dx, svp1, svp2, &
            svp3, svpt0, ep1, ep2, karman, ch, qcg, &
            itimestep, wstar, qstar, ustm, ck, cka, &
            cd, cda, spp_pbl, rstoch1d, isftcflx, &
            iz0tlnd, its, ite, &
            errmsg, errflg)

        if (errflg /= 0) then
            return
        end if

        ! The first pass treated sea ice cells as land cells due to how `xland` is defined.
        ! The second pass treats sea ice cells as water cells instead.
        ! Then, for sea ice cells only, the final results are weighted averages according to their area fractions.
        if (any(mask_sea_ice_cell)) then
            call sf_mynn_run( &
                u1d, v1d, t1d, qv1d, p1d, dz8w1d, rho1d, &
                u1d2, v1d2, dz2w1d, cp, g, rovcp, r, xlv, &
                psfcpa, chs_sea, chs2_sea, cqs2_sea, cpm_sea, pblh, rmol_sea, &
                znt_sea, ust_sea, mavail_sea, zol_sea, mol_sea, regime_sea, psim_sea, &
                psih_sea, xland_sea, hfx_sea, qfx_sea, tsk_sea, u10_sea, v10_sea, th2_sea, &
                t2_sea, q2_sea, flhc_sea, flqc_sea, snowh, qgh_sea, qsfc_sea, lh_sea, &
                gz1oz0_sea, wspd_sea, br_sea, isfflx, dx, svp1, svp2, &
                svp3, svpt0, ep1, ep2, karman, ch_sea, qcg, &
                itimestep, wstar_sea, qstar_sea, ustm_sea, ck_sea, cka_sea, &
                cd_sea, cda_sea, spp_pbl, rstoch1d, isftcflx, &
                iz0tlnd, its, ite, &
                errmsg, errflg)

            if (errflg /= 0) then
                return
            end if

            ! Blend the final results between sea ice and sea water.
            where (mask_sea_ice_cell)
                chs = chs * icefrac + (1.0_kind_phys - icefrac) * chs_sea
                chs2 = chs2 * icefrac + (1.0_kind_phys - icefrac) * chs2_sea
                cqs2 = cqs2 * icefrac + (1.0_kind_phys - icefrac) * cqs2_sea
                cpm = cpm * icefrac + (1.0_kind_phys - icefrac) * cpm_sea
                rmol = rmol * icefrac + (1.0_kind_phys - icefrac) * rmol_sea
                znt = znt * icefrac + (1.0_kind_phys - icefrac) * znt_sea
                ust = ust * icefrac + (1.0_kind_phys - icefrac) * ust_sea
                zol = zol * icefrac + (1.0_kind_phys - icefrac) * zol_sea
                mol = mol * icefrac + (1.0_kind_phys - icefrac) * mol_sea
                psim = psim * icefrac + (1.0_kind_phys - icefrac) * psim_sea
                psih = psih * icefrac + (1.0_kind_phys - icefrac) * psih_sea
                hfx = hfx * icefrac + (1.0_kind_phys - icefrac) * hfx_sea
                qfx = qfx * icefrac + (1.0_kind_phys - icefrac) * qfx_sea
                flhc = flhc * icefrac + (1.0_kind_phys - icefrac) * flhc_sea
                flqc = flqc * icefrac + (1.0_kind_phys - icefrac) * flqc_sea
                qgh = qgh * icefrac + (1.0_kind_phys - icefrac) * qgh_sea
                qsfc = qsfc * icefrac + (1.0_kind_phys - icefrac) * qsfc_sea
                lh = lh * icefrac + (1.0_kind_phys - icefrac) * lh_sea
                gz1oz0 = gz1oz0 * icefrac + (1.0_kind_phys - icefrac) * gz1oz0_sea
                wspd = wspd * icefrac + (1.0_kind_phys - icefrac) * wspd_sea
                br = br * icefrac + (1.0_kind_phys - icefrac) * br_sea

                u10 = u10 * icefrac + (1.0_kind_phys - icefrac) * u10_sea
                v10 = v10 * icefrac + (1.0_kind_phys - icefrac) * v10_sea
                th2 = th2 * icefrac + (1.0_kind_phys - icefrac) * th2_sea
                t2 = t2 * icefrac + (1.0_kind_phys - icefrac) * t2_sea
                q2 = q2 * icefrac + (1.0_kind_phys - icefrac) * q2_sea
                wstar = wstar * icefrac + (1.0_kind_phys - icefrac) * wstar_sea
                qstar = qstar * icefrac + (1.0_kind_phys - icefrac) * qstar_sea
                ustm = ustm * icefrac + (1.0_kind_phys - icefrac) * ustm_sea
                ck = ck * icefrac + (1.0_kind_phys - icefrac) * ck_sea
                cka = cka * icefrac + (1.0_kind_phys - icefrac) * cka_sea
                cd = cd * icefrac + (1.0_kind_phys - icefrac) * cd_sea
                cda = cda * icefrac + (1.0_kind_phys - icefrac) * cda_sea
            end where

            ! `regime` is a categorical variable. Assign according to whether sea ice or sea water is dominant.
            where (mask_sea_ice_cell .and. icefrac < 0.5_kind_phys)
                regime = regime_sea
            end where
        end if

        cflx(:, water_vapor_mixing_ratio_index) = qfx(:)

        errmsg = ''
        errflg = 0
    end subroutine sf_mynn_compat_run

    !> \section arg_table_sf_mynn_diagnostics_init Argument Table
    !! \htmlinclude sf_mynn_diagnostics_init.html
    subroutine sf_mynn_diagnostics_init( &
            errmsg, errflg)
        use cam_history, only: history_add_field
        use cam_history_support, only: horiz_only

        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_add_field('sf_mynn_lh', &
            'surface_upward_latent_heat_flux_from_coupler', horiz_only, 'avg', 'W m-2')
        call history_add_field('sf_mynn_hfx', &
            'surface_upward_sensible_heat_flux_from_coupler', horiz_only, 'avg', 'W m-2')
        call history_add_field('sf_mynn_qfx', &
            'surface_upward_water_vapor_flux', horiz_only, 'avg', 'kg m-2 s-1')

        call history_add_field('sf_mynn_cd', &
            'drag_coefficient_for_momentum_at_10m', horiz_only, 'avg', '1')
        call history_add_field('sf_mynn_cda', &
            'drag_coefficient_for_momentum_at_surface_adjacent_layer', horiz_only, 'avg', '1')
        call history_add_field('sf_mynn_ck', &
            'bulk_exchange_coefficient_for_enthalpy_at_2m', horiz_only, 'avg', '1')
        call history_add_field('sf_mynn_cka', &
            'bulk_exchange_coefficient_for_enthalpy_at_surface_adjacent_layer', horiz_only, 'avg', '1')

        call history_add_field('sf_mynn_q2', &
            'water_vapor_mixing_ratio_wrt_dry_air_at_2m', horiz_only, 'avg', 'kg kg-1')
        call history_add_field('sf_mynn_t2', &
            'air_temperature_at_2m', horiz_only, 'avg', 'K')
        call history_add_field('sf_mynn_th2', &
            'air_potential_temperature_at_2m', horiz_only, 'avg', 'K')
        call history_add_field('sf_mynn_u10', &
            'eastward_wind_at_10m', horiz_only, 'avg', 'm s-1')
        call history_add_field('sf_mynn_v10', &
            'northward_wind_at_10m', horiz_only, 'avg', 'm s-1')

        call history_add_field('sf_mynn_ustm', &
            'surface_friction_velocity_assuming_no_correction_for_surface_convective_velocity_scale', horiz_only, 'avg', 'm s-1')
        call history_add_field('sf_mynn_wstar', &
            'surface_convective_velocity_scale', horiz_only, 'avg', 'm s-1')
        call history_add_field('sf_mynn_qstar', &
            'surface_moisture_scale', horiz_only, 'avg', 'g kg-1')

        errmsg = ''
        errflg = 0
    end subroutine sf_mynn_diagnostics_init

    !> \section arg_table_sf_mynn_diagnostics_run Argument Table
    !! \htmlinclude sf_mynn_diagnostics_run.html
    subroutine sf_mynn_diagnostics_run( &
            lh, hfx, qfx, &
            cd, cda, ck, cka, &
            q2, t2, th2, u10, v10, &
            ustm, wstar, qstar, &
            errmsg, errflg)
        use cam_history, only: history_out_field
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: lh(:), hfx(:), qfx(:), &
                                       cd(:), cda(:), ck(:), cka(:), &
                                       q2(:), t2(:), th2(:), u10(:), v10(:), &
                                       ustm(:), wstar(:), qstar(:)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_out_field('sf_mynn_lh', lh)
        call history_out_field('sf_mynn_hfx', hfx)
        call history_out_field('sf_mynn_qfx', qfx)

        call history_out_field('sf_mynn_cd', cd)
        call history_out_field('sf_mynn_cda', cda)
        call history_out_field('sf_mynn_ck', ck)
        call history_out_field('sf_mynn_cka', cka)

        call history_out_field('sf_mynn_q2', q2)
        call history_out_field('sf_mynn_t2', t2)
        call history_out_field('sf_mynn_th2', th2)
        call history_out_field('sf_mynn_u10', u10)
        call history_out_field('sf_mynn_v10', v10)

        call history_out_field('sf_mynn_ustm', ustm)
        call history_out_field('sf_mynn_wstar', wstar)
        call history_out_field('sf_mynn_qstar', qstar)

        errmsg = ''
        errflg = 0
    end subroutine sf_mynn_diagnostics_run
end module sf_mynn_compat
