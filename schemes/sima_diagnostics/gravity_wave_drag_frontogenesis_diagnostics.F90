! Diagnostics for gravity wave drag from frontogenesis.
module gravity_wave_drag_frontogenesis_diagnostics
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gravity_wave_drag_frontogenesis_diagnostics_init
  public :: gravity_wave_drag_frontogenesis_diagnostics_run

contains

!> \section arg_table_gravity_wave_drag_frontogenesis_diagnostics_init  Argument Table
!! \htmlinclude gravity_wave_drag_frontogenesis_diagnostics_init.html
  subroutine gravity_wave_drag_frontogenesis_diagnostics_init(errmsg, errflg)
    use cam_history, only: history_add_field

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Frontogenesis function
    call history_add_field('FRONTGF', 'frontogenesis_function', 'lev', 'avg', 'K2 m-2 s-1')
    call history_add_field('FRONTGFA', 'frontogenesis_angle', 'lev', 'avg', 'radians')

    ! Wind tendency bins
    call history_add_field('UTEND1', 'wind_tendency_in_first_speed_bin_due_to_frontogenesis_gravity_wave_drag c<-40', 'lev', 'avg', 'm s-2')
    call history_add_field('UTEND2', 'wind_tendency_in_second_speed_bin_due_to_frontogenesis_gravity_wave_drag -40<c<-15', 'lev', 'avg', 'm s-2')
    call history_add_field('UTEND3', 'wind_tendency_in_third_speed_bin_due_to_frontogenesis_gravity_wave_drag -15<c<15', 'lev', 'avg', 'm s-2')
    call history_add_field('UTEND4', 'wind_tendency_in_fourth_speed_bin_due_to_frontogenesis_gravity_wave_drag 15<c<40', 'lev', 'avg', 'm s-2')
    call history_add_field('UTEND5', 'wind_tendency_in_fifth_speed_bin_due_to_frontogenesis_gravity_wave_drag c>40', 'lev', 'avg', 'm s-2')

    ! Temperature tendencies (converted from J kg-1 s-1 to K s-1 for output)
    call history_add_field('TTGWSDF', 'tendency_of_air_temperature_due_to_diffusion_due_to_frontogenesis_gravity_wave_drag', 'lev', 'avg', 'K s-1')
    call history_add_field('TTGWSKE', 'tendency_of_air_temperature_due_to_kinetic_energy_dissipation_due_to_frontogenesis_gravity_wave_drag', 'lev', 'avg', 'K s-1')

    ! Total wind tendencies
    call history_add_field('UTGWSPEC', 'tendency_of_eastward_wind_due_to_frontogenesis_gravity_wave_drag', 'lev', 'avg', 'm s-2')
    call history_add_field('VTGWSPEC', 'tendency_of_northward_wind_due_to_frontogenesis_gravity_wave_drag', 'lev', 'avg', 'm s-2')
    call history_add_field('TTGWSPEC', 'tendency_of_air_temperature_due_to_frontogenesis_gravity_wave_drag', 'lev', 'avg', 'K s-1')

    ! Directional momentum flux
    call history_add_field('TAUE', 'eastward_reynolds_stress_due_to_frontogenesis_gravity_wave_drag', 'ilev', 'avg', 'N m-2')
    call history_add_field('TAUW', 'westward_reynolds_stress_due_to_frontogenesis_gravity_wave_drag', 'ilev', 'avg', 'N m-2')
    call history_add_field('TAUN', 'northward_reynolds_stress_due_to_frontogenesis_gravity_wave_drag', 'ilev', 'avg', 'N m-2')
    call history_add_field('TAUS', 'southward_reynolds_stress_due_to_frontogenesis_gravity_wave_drag', 'ilev', 'avg', 'N m-2')
    call history_add_field('TAUNET', 'net_eastward_reynolds_stress_due_to_frontogenesis_gravity_wave_drag', 'ilev', 'avg', 'N m-2')

  end subroutine gravity_wave_drag_frontogenesis_diagnostics_init

!> \section arg_table_gravity_wave_drag_frontogenesis_diagnostics_run  Argument Table
!! \htmlinclude gravity_wave_drag_frontogenesis_diagnostics_run.html
  subroutine gravity_wave_drag_frontogenesis_diagnostics_run( &
    cpair, &
    frontgf, frontga, &
    utend1, utend2, utend3, utend4, utend5, &
    dttdf, dttke, &
    utgw, vtgw, ttgw, &
    taucd_east, taucd_west, taucd_north, taucd_south, &
    errmsg, errflg)

    use cam_history, only: history_out_field

    real(kind_phys), intent(in) :: cpair

    real(kind_phys), intent(in) :: frontgf(:,:)
    real(kind_phys), intent(in) :: frontga(:,:)

    real(kind_phys), intent(in) :: utend1(:,:)
    real(kind_phys), intent(in) :: utend2(:,:)
    real(kind_phys), intent(in) :: utend3(:,:)
    real(kind_phys), intent(in) :: utend4(:,:)
    real(kind_phys), intent(in) :: utend5(:,:)

    real(kind_phys), intent(in) :: dttdf(:,:)
    real(kind_phys), intent(in) :: dttke(:,:)

    real(kind_phys), intent(in) :: utgw(:,:)
    real(kind_phys), intent(in) :: vtgw(:,:)
    real(kind_phys), intent(in) :: ttgw(:,:)

    real(kind_phys), intent(in) :: taucd_east(:,:)
    real(kind_phys), intent(in) :: taucd_west(:,:)
    real(kind_phys), intent(in) :: taucd_north(:,:)
    real(kind_phys), intent(in) :: taucd_south(:,:)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Output frontogenesis function and angle.
    call history_out_field('FRONTGF', frontgf)
    call history_out_field('FRONTGFA', frontga)

    ! Output wind tendency bins
    call history_out_field('UTEND1', utend1)
    call history_out_field('UTEND2', utend2)
    call history_out_field('UTEND3', utend3)
    call history_out_field('UTEND4', utend4)
    call history_out_field('UTEND5', utend5)

    ! Output temperature tendencies (convert from J/kg/s to K/s)
    ! Note: this uses cpair instead of cpairv.
    ! Also see atmospheric_physics issue #320.
    call history_out_field('TTGWSDF', dttdf / cpair)
    call history_out_field('TTGWSKE', dttke / cpair)

    ! Output wind and temperature tendencies
    call history_out_field('UTGWSPEC', utgw)
    call history_out_field('VTGWSPEC', vtgw)
    call history_out_field('TTGWSPEC', ttgw)

    ! Output directional momentum flux
    call history_out_field('TAUE', taucd_east)
    call history_out_field('TAUW', taucd_west)
    call history_out_field('TAUN', taucd_north)
    call history_out_field('TAUS', taucd_south)
    call history_out_field('TAUNET', taucd_east + taucd_west)

  end subroutine gravity_wave_drag_frontogenesis_diagnostics_run

end module gravity_wave_drag_frontogenesis_diagnostics
