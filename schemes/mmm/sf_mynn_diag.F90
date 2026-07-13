!> This module contains diagnostic schemes that are specific to the MYNN surface layer scheme,
!> which is part of the MMM physics.
module sf_mynn_diag
    implicit none

    private
    public :: sf_mynn_diagnostics_init
    public :: sf_mynn_diagnostics_run
contains
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
end module sf_mynn_diag
