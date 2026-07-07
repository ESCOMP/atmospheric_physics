!> This module contains diagnostic schemes that are specific to the new Tiedtke cumulus scheme,
!> which is part of the MMM physics.
module cu_ntiedtke_diag
    implicit none

    private
    public :: cu_ntiedtke_diagnostics_init
    public :: cu_ntiedtke_diagnostics_run
contains
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
            rucuten, rvcuten, &
            rthcuten, &
            rqvcuten, rqccuten, rqicuten, &
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
end module cu_ntiedtke_diag
