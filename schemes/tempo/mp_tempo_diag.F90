module mp_tempo_diag
    implicit none

    private
    public :: mp_tempo_diagnostics_init
    public :: mp_tempo_diagnostics_run
contains
    !> \section arg_table_mp_tempo_diagnostics_init Argument Table
    !! \htmlinclude mp_tempo_diagnostics_init.html
    subroutine mp_tempo_diagnostics_init( &
            errmsg, errflg)
        use cam_history, only: history_add_field
        use cam_history_support, only: horiz_only

        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_add_field( &
            'mp_tempo_rain', &
            'lwe_thickness_of_rainfall_amount', &
            horiz_only, 'avg', 'm')
        call history_add_field( &
            'mp_tempo_ice', &
            'lwe_thickness_of_ice_amount', &
            horiz_only, 'avg', 'm')
        call history_add_field( &
            'mp_tempo_snow', &
            'lwe_thickness_of_snowfall_amount', &
            horiz_only, 'avg', 'm')
        call history_add_field( &
            'mp_tempo_graupel', &
            'lwe_thickness_of_graupel_amount', &
            horiz_only, 'avg', 'm')
        call history_add_field( &
            'mp_tempo_prcp', &
            'lwe_thickness_of_precipitation_amount', &
            horiz_only, 'avg', 'm')
        call history_add_field( &
            'mp_tempo_frozen_fraction', &
            'ratio_of_snowfall_to_rainfall', &
            horiz_only, 'avg', 'fraction')

        call history_add_field( &
            'mp_tempo_refl10cm', &
            'radar_reflectivity_at_10cm_wavelength', &
            'lev', 'avg', 'dBZ')
        call history_add_field( &
            'mp_tempo_re_cloud', &
            'effective_radius_of_stratiform_cloud_liquid_water_droplet', &
            'lev', 'avg', 'um')
        call history_add_field( &
            'mp_tempo_re_ice', &
            'effective_radius_of_stratiform_cloud_ice_crystal', &
            'lev', 'avg', 'um')
        call history_add_field( &
            'mp_tempo_re_snow', &
            'effective_radius_of_stratiform_cloud_snow_crystal', &
            'lev', 'avg', 'um')

        errmsg = ''
        errflg = 0
    end subroutine mp_tempo_diagnostics_init

    !> \section arg_table_mp_tempo_diagnostics_run Argument Table
    !! \htmlinclude mp_tempo_diagnostics_run.html
    subroutine mp_tempo_diagnostics_run( &
            rain, ice, snow, graupel, prcp, &
            frozen_fraction, refl10cm, re_cloud, re_ice, re_snow, &
            errmsg, errflg)
        use cam_history, only: history_out_field
        use ccpp_kinds, only: kind_phys

        real(kind_phys), intent(in) :: rain(:), ice(:), snow(:), graupel(:), prcp(:), frozen_fraction(:), &
                                       refl10cm(:, :), re_cloud(:, :), re_ice(:, :), re_snow(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        call history_out_field('mp_tempo_rain', rain)
        call history_out_field('mp_tempo_ice', ice)
        call history_out_field('mp_tempo_snow', snow)
        call history_out_field('mp_tempo_graupel', graupel)
        call history_out_field('mp_tempo_prcp', prcp)
        call history_out_field('mp_tempo_frozen_fraction', frozen_fraction)

        call history_out_field('mp_tempo_refl10cm', refl10cm)
        call history_out_field('mp_tempo_re_cloud', re_cloud)
        call history_out_field('mp_tempo_re_ice', re_ice)
        call history_out_field('mp_tempo_re_snow', re_snow)

        errmsg = ''
        errflg = 0
    end subroutine mp_tempo_diagnostics_run
end module mp_tempo_diag
