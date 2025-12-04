! Diagnostics for ridge-based gravity wave drag (meso-Beta)
module gravity_wave_drag_ridge_beta_diagnostics
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gravity_wave_drag_ridge_beta_diagnostics_init
  public :: gravity_wave_drag_ridge_beta_diagnostics_run

contains

!> \section arg_table_gravity_wave_drag_ridge_beta_diagnostics_init  Argument Table
!! \htmlinclude gravity_wave_drag_ridge_beta_diagnostics_init.html
  subroutine gravity_wave_drag_ridge_beta_diagnostics_init(errmsg, errflg)
    use cam_history, only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! State at entry to GW scheme
    call history_add_field('UEGW', 'eastward_wind_before_beta_ridge_gravity_wave_drag', 'lev', 'avg', 'm s-1')
    call history_add_field('VEGW', 'northward_wind_before_beta_ridge_gravity_wave_drag', 'lev', 'avg', 'm s-1')
    call history_add_field('TEGW', 'air_temperature_before_beta_ridge_gravity_wave_drag', 'lev', 'avg', 'K')
    call history_add_field('ZEGW', 'geopotential_height_wrt_surface_at_interfaces_before_beta_ridge_gravity_wave_drag', 'ilev', 'avg', 'm')
    call history_add_field('ZMGW', 'geopotential_height_wrt_surface_before_beta_ridge_gravity_wave_drag', 'lev', 'avg', 'm')

    ! Momentum flux profile
    call history_add_field('TAUARDGBETAX', 'eastward_stress_at_interfaces_due_to_beta_ridge_gravity_wave_drag_before_residual', 'ilev', 'inst', 'N m-2')
    call history_add_field('TAUARDGBETAY', 'northward_stress_at_interfaces_due_to_beta_ridge_gravity_wave_drag_before_residual', 'ilev', 'inst', 'N m-2')

    ! Wind tendencies
    call history_add_field('UTGWORO', 'tendency_of_eastward_wind_due_to_beta_ridge_gravity_wave_drag', 'lev', 'avg', 'm s-2')
    call history_add_field('VTGWORO', 'tendency_of_northward_wind_due_to_beta_ridge_gravity_wave_drag', 'lev', 'avg', 'm s-2')
    call history_add_field('TTGWORO', 'tendency_of_air_temperature_due_to_beta_ridge_gravity_wave_drag', 'lev', 'avg', 'K s-1')

    ! Surface stress (retrieved from surface layer of stress at interfaces output from scheme.)
    call history_add_field('TAUGWX', 'eastward_stress_at_surface_due_to_beta_ridge_gravity_wave_drag', horiz_only, 'avg', 'N m-2')
    call history_add_field('TAUGWY', 'northward_stress_at_surface_due_to_beta_ridge_gravity_wave_drag', horiz_only, 'avg', 'N m-2')

  end subroutine gravity_wave_drag_ridge_beta_diagnostics_init

!> \section arg_table_gravity_wave_drag_ridge_beta_diagnostics_run  Argument Table
!! \htmlinclude gravity_wave_drag_ridge_beta_diagnostics_run.html
  subroutine gravity_wave_drag_ridge_beta_diagnostics_run( &
    pverp, &
    u, v, t, zi, zm, &
    taurx, taury, &
    tauardgx, tauardgy, &
    utgw, vtgw, ttgw, &
    errmsg, errflg)

    use cam_history, only: history_out_field

    integer,         intent(in) :: pverp

    ! State variables at entry [m s-1], [K], [m]
    real(kind_phys), intent(in) :: u(:,:)
    real(kind_phys), intent(in) :: v(:,:)
    real(kind_phys), intent(in) :: t(:,:)
    real(kind_phys), intent(in) :: zi(:,:)
    real(kind_phys), intent(in) :: zm(:,:)

    ! Total momentum flux [N m-2]
    real(kind_phys), intent(in) :: taurx(:,:)
    real(kind_phys), intent(in) :: taury(:,:)

    ! Ridge-based momentum flux [N m-2]
    real(kind_phys), intent(in) :: tauardgx(:,:)
    real(kind_phys), intent(in) :: tauardgy(:,:)

    real(kind_phys), intent(in) :: utgw(:,:)
    real(kind_phys), intent(in) :: vtgw(:,:)
    real(kind_phys), intent(in) :: ttgw(:,:)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Output state before gravity wave drag parameterization.
    call history_out_field('UEGW', u)
    call history_out_field('VEGW', v)
    call history_out_field('TEGW', t)
    call history_out_field('ZEGW', zi)
    call history_out_field('ZMGW', zm)

    ! Output momentum flux profiles
    call history_out_field('TAUARDGBETAX', tauardgx)
    call history_out_field('TAUARDGBETAY', tauardgy)

    ! Output wind and temperature tendencies
    call history_out_field('UTGWORO', utgw)
    call history_out_field('VTGWORO', vtgw)
    call history_out_field('TTGWORO', ttgw)

    ! Output surface stress (bottom level)
    call history_out_field('TAUGWX', taurx(:,pverp))
    call history_out_field('TAUGWY', taury(:,pverp))

  end subroutine gravity_wave_drag_ridge_beta_diagnostics_run

end module gravity_wave_drag_ridge_beta_diagnostics
