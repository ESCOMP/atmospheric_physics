!> This module contains diagnostic schemes that are specific to the MYNN PBL scheme,
!> which is part of the MMM physics.
module bl_mynn_diag
    implicit none

    private
    public :: bl_mynn_diagnostics_init
    public :: bl_mynn_diagnostics_run
contains
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
end module bl_mynn_diag
