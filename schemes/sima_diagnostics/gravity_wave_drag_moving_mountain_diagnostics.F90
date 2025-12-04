! Diagnostics for gravity wave drag from moving mountain parameterization
module gravity_wave_drag_moving_mountain_diagnostics
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gravity_wave_drag_moving_mountain_diagnostics_init
  public :: gravity_wave_drag_moving_mountain_diagnostics_run

contains

!> \section arg_table_gravity_wave_drag_moving_mountain_diagnostics_init  Argument Table
!! \htmlinclude gravity_wave_drag_moving_mountain_diagnostics_init.html
  subroutine gravity_wave_drag_moving_mountain_diagnostics_init(errmsg, errflg)
    use cam_history, only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_add_field('U_MOVMTN_IN', 'eastward_wind_before_moving_mountain_gravity_wave_drag', 'lev', 'avg', 'm s-1')
    call history_add_field('V_MOVMTN_IN', 'northward_wind_before_moving_mountain_gravity_wave_drag', 'lev', 'avg', 'm s-1')

    call history_add_field('SRC_LEVEL_MOVMTN', 'vertical_index_of_gravity_wave_source_level', horiz_only, 'avg', 'index')
    call history_add_field('TND_LEVEL_MOVMTN', 'vertical_index_of_lowest_gravity_wave_tendency_level', horiz_only, 'avg', 'index')
    call history_add_field('STEER_LEVEL_MOVMTN', 'index_of_steering_level_due_to_moving_mountain_gravity_wave_drag', horiz_only, 'avg', 'index')

    call history_add_field('UBI_MOVMTN', 'wind_speed_along_gravity_wave_source_direction_at_interfaces', 'ilev', 'inst', 'm s-1')
    call history_add_field('UBM_MOVMTN', 'wind_speed_along_gravity_wave_source_direction', 'lev', 'inst', 'm s-1')

    call history_add_field('TAU_MOVMTN', 'total_stress_at_interface_from_reference_wavenumber_due_to_moving_mountain_gravity_wave_drag', 'ilev', 'inst', 'N m-2')
    call history_add_field('GWUT_MOVMTN', 'wind_tendency_in_source_direction_due_to_moving_mountain_gravity_wave_drag', 'lev', 'inst', 'm s-2')
    call history_add_field('UTGW_MOVMTN', 'tendency_of_eastward_wind_due_to_moving_mountain_gravity_wave_drag', 'lev', 'inst', 'm s-2')
    call history_add_field('VTGW_MOVMTN', 'tendency_of_northward_wind_due_to_moving_mountain_gravity_wave_drag', 'lev', 'inst', 'm s-2')

    call history_add_field('HDEPTH_MOVMTN', 'convective_heating_depth_due_to_moving_mountain_gravity_wave_drag', horiz_only, 'inst', 'km')
    call history_add_field('NETDT_MOVMTN', 'tendency_of_air_temperature_due_to_deep_convection', 'lev', 'inst', 'K s-1')

    call history_add_field('TTEND_CLUBB', 'tendency_of_air_temperature_due_to_clubb', 'lev', 'avg', 'K s-1')
    call history_add_field('UPWP_CLUBB_GW', 'X-momflux from CLUBB to GW', 'ilev', 'avg', 'm2 s-2')
    call history_add_field('VPWP_CLUBB_GW', 'Y-momflux from CLUBB to GW', 'ilev', 'avg', 'm2 s-2')
    call history_add_field('VORT4GW', 'relative_vorticity', 'lev', 'avg', 's-1')

    call history_add_field('UCELL_MOVMTN', 'eastward_wind_at_steering_level_due_to_moving_mountain_gravity_wave_drag', horiz_only, 'inst', 'm s-1')
    call history_add_field('VCELL_MOVMTN', 'northward_wind_at_steering_level_due_to_moving_mountain_gravity_wave_drag', horiz_only, 'inst', 'm s-1')
    call history_add_field('CS_MOVMTN', 'gravity_wave_phase_speed_in_source_direction_due_to_moving_mountain_gravity_wave_drag', horiz_only, 'inst', 'm s-1')
    call history_add_field('XPWP_SRC_MOVMTN', 'momentum_flux_source_due_to_moving_mountain_gravity_wave_drag', horiz_only, 'inst', 'm2 s-2')

  end subroutine gravity_wave_drag_moving_mountain_diagnostics_init

!> \section arg_table_gravity_wave_drag_moving_mountain_diagnostics_run  Argument Table
!! \htmlinclude gravity_wave_drag_moving_mountain_diagnostics_run.html
  subroutine gravity_wave_drag_moving_mountain_diagnostics_run( &
    u, v, &
    src_level, tend_level, steer_level, &
    ubi, ubm, &
    tau0, gwut0, utgw, vtgw, &
    hdepth, ttend_dp, &
    ttend_clubb, &
    upwp_clubb, vpwp_clubb, vorticity, &
    usteer, vsteer, CS, xpwp_src, &
    errmsg, errflg)

    use cam_history, only: history_out_field

    ! Input winds [m s-1]
    real(kind_phys), intent(in) :: u(:,:)
    real(kind_phys), intent(in) :: v(:,:)

    ! Outputs from gravity wave scheme, including inputs from CLUBB:
    integer,         intent(in) :: src_level(:)
    integer,         intent(in) :: tend_level(:)
    real(kind_phys), intent(in) :: steer_level(:)
    real(kind_phys), intent(in) :: ubi(:,:)
    real(kind_phys), intent(in) :: ubm(:,:)
    real(kind_phys), intent(in) :: tau0(:,:)
    real(kind_phys), intent(in) :: gwut0(:,:)
    real(kind_phys), intent(in) :: utgw(:,:)
    real(kind_phys), intent(in) :: vtgw(:,:)
    real(kind_phys), intent(in) :: hdepth(:)
    real(kind_phys), intent(in) :: ttend_dp(:,:)
    real(kind_phys), intent(in) :: ttend_clubb(:,:)     ! [K s-1]
    real(kind_phys), intent(in) :: upwp_clubb(:,:)   ! [m2 s-2]
    real(kind_phys), intent(in) :: vpwp_clubb(:,:)   ! [m2 s-2]
    real(kind_phys), intent(in) :: vorticity(:,:)         ! [s-1]

    real(kind_phys), intent(in) :: usteer(:)   ! [m s-1]
    real(kind_phys), intent(in) :: vsteer(:)   ! [m s-1]
    real(kind_phys), intent(in) :: CS(:)       ! [m s-1]
    real(kind_phys), intent(in) :: xpwp_src(:) ! [m2 s-2]

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_out_field('U_MOVMTN_IN', u)
    call history_out_field('V_MOVMTN_IN', v)
    call history_out_field('SRC_LEVEL_MOVMTN', real(src_level, kind_phys))
    call history_out_field('TND_LEVEL_MOVMTN', real(tend_level, kind_phys))
    call history_out_field('STEER_LEVEL_MOVMTN', steer_level)
    call history_out_field('UBI_MOVMTN', ubi)
    call history_out_field('UBM_MOVMTN', ubm)
    call history_out_field('TAU_MOVMTN', tau0)
    call history_out_field('GWUT_MOVMTN', gwut0)
    call history_out_field('UTGW_MOVMTN', utgw)
    call history_out_field('VTGW_MOVMTN', vtgw)
    call history_out_field('HDEPTH_MOVMTN', hdepth / 1000.0_kind_phys)
    call history_out_field('NETDT_MOVMTN', ttend_dp)
    call history_out_field('TTEND_CLUBB', ttend_clubb)
    call history_out_field('UPWP_CLUBB_GW', upwp_clubb)
    call history_out_field('VPWP_CLUBB_GW', vpwp_clubb)
    call history_out_field('VORT4GW', vorticity)
    call history_out_field('UCELL_MOVMTN', usteer)
    call history_out_field('VCELL_MOVMTN', vsteer)
    call history_out_field('CS_MOVMTN', CS)
    call history_out_field('XPWP_SRC_MOVMTN', xpwp_src)

  end subroutine gravity_wave_drag_moving_mountain_diagnostics_run
end module gravity_wave_drag_moving_mountain_diagnostics
