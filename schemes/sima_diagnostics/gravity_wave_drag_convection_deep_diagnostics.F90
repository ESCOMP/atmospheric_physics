! Diagnostics for gravity wave drag from deep convection (Beres scheme)
module gravity_wave_drag_convection_deep_diagnostics
  implicit none
  private

  public :: gravity_wave_drag_convection_deep_diagnostics_init
  public :: gravity_wave_drag_convection_deep_diagnostics_run

contains

!> \section arg_table_gravity_wave_drag_convection_deep_diagnostics_init  Argument Table
!! \htmlinclude gravity_wave_drag_convection_deep_diagnostics_init.html
  subroutine gravity_wave_drag_convection_deep_diagnostics_init(errmsg, errflg)
    use cam_history, only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_add_field('BTTGWSDF', 'tendency_of_air_temperature_due_to_diffusion_due_to_deep_convective_gravity_wave_drag', 'lev', 'avg', 'K s-1')
    call history_add_field('BTTGWSKE', 'tendency_of_air_temperature_due_to_kinetic_energy_dissipation_due_to_deep_convective_gravity_wave_drag', 'lev', 'avg', 'K s-1')

    ! Wind tendencies
    ! These are vector fields.
    ! TODO: register_vector_field equivalent in CAM-SIMA?
    call history_add_field('BUTGWSPEC', 'tendency_of_eastward_wind_due_to_deep_convective_gravity_wave_drag', 'lev', 'avg', 'm s-2')
    call history_add_field('BVTGWSPEC', 'tendency_of_northward_wind_due_to_deep_convective_gravity_wave_drag', 'lev', 'avg', 'm s-2')

    call history_add_field('BTTGWSPEC', 'tendency_of_air_temperature_due_to_deep_convective_gravity_wave_drag', 'lev', 'avg', 'K s-1')

    call history_add_field('BTAUE', 'eastward_reynolds_stress_due_to_deep_convective_gravity_wave_drag', 'ilev', 'avg', 'Pa')
    call history_add_field('BTAUW', 'westward_reynolds_stress_due_to_deep_convective_gravity_wave_drag', 'ilev', 'avg', 'Pa')
    call history_add_field('BTAUN', 'northward_reynolds_stress_due_to_deep_convective_gravity_wave_drag', 'ilev', 'avg', 'Pa')
    call history_add_field('BTAUS', 'southward_reynolds_stress_due_to_deep_convective_gravity_wave_drag', 'ilev', 'avg', 'Pa')
    call history_add_field('BTAUNET', 'net_eastward_reynolds_stress_due_to_deep_convective_gravity_wave_drag', 'ilev', 'avg', 'Pa')

    call history_add_field('NETDT', 'tendency_of_air_temperature_due_to_deep_convection', 'lev', 'avg', 'K s-1')
    call history_add_field('HDEPTH', 'convective_heating_depth_for_deep_convective_gravity_wave_drag', horiz_only, 'avg', 'km')
    call history_add_field('MAXQ0', 'maximum_column_heating_rate_due_to_deep_convective_gravity_wave_drag', horiz_only, 'avg', 'K day-1')

  end subroutine gravity_wave_drag_convection_deep_diagnostics_init

!> \section arg_table_gravity_wave_drag_convection_deep_diagnostics_run  Argument Table
!! \htmlinclude gravity_wave_drag_convection_deep_diagnostics_run.html
  subroutine gravity_wave_drag_convection_deep_diagnostics_run( &
    cpair, &
    dttdf, dttke, &
    utgw, vtgw, ttgw, &
    taucd_east, taucd_west, taucd_north, taucd_south, &
    ttend_dp, hdepth, maxq0, &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys
    use cam_history, only: history_out_field

    real(kind_phys), intent(in) :: cpair
    real(kind_phys), intent(in) :: dttdf(:,:)
    real(kind_phys), intent(in) :: dttke(:,:)
    real(kind_phys), intent(in) :: utgw(:,:)
    real(kind_phys), intent(in) :: vtgw(:,:)
    real(kind_phys), intent(in) :: ttgw(:,:)
    real(kind_phys), intent(in) :: taucd_east(:,:)
    real(kind_phys), intent(in) :: taucd_west(:,:)
    real(kind_phys), intent(in) :: taucd_north(:,:)
    real(kind_phys), intent(in) :: taucd_south(:,:)
    real(kind_phys), intent(in) :: ttend_dp(:,:)  ! [K s-1]
    real(kind_phys), intent(in) :: hdepth(:)      ! [m]
    real(kind_phys), intent(in) :: maxq0(:)       ! [K day-1]

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Output temperature tendencies (convert from J/kg/s to K/s)
    ! Note: this uses cpair instead of cpairv.
    ! Also see atmospheric_physics issue #320.
    call history_out_field('BTTGWSDF', dttdf / cpair)
    call history_out_field('BTTGWSKE', dttke / cpair)

    call history_out_field('BUTGWSPEC', utgw)
    call history_out_field('BVTGWSPEC', vtgw)
    call history_out_field('BTTGWSPEC', ttgw)

    call history_out_field('BTAUE', taucd_east)
    call history_out_field('BTAUW', taucd_west)
    call history_out_field('BTAUN', taucd_north)
    call history_out_field('BTAUS', taucd_south)
    call history_out_field('BTAUNET', taucd_east + taucd_west)

    call history_out_field('NETDT', ttend_dp)

    ! Output heating diagnostics (convert depth to km)
    call history_out_field('HDEPTH', hdepth / 1000.0_kind_phys)
    call history_out_field('MAXQ0', maxq0)

  end subroutine gravity_wave_drag_convection_deep_diagnostics_run

end module gravity_wave_drag_convection_deep_diagnostics
