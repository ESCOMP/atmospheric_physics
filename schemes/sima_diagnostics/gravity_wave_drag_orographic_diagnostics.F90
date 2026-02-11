! Diagnostics for orographic gravity wave drag
module gravity_wave_drag_orographic_diagnostics
  implicit none
  private

  public :: gravity_wave_drag_orographic_diagnostics_init
  public :: gravity_wave_drag_orographic_diagnostics_run

contains

!> \section arg_table_gravity_wave_drag_orographic_diagnostics_init  Argument Table
!! \htmlinclude gravity_wave_drag_orographic_diagnostics_init.html
  subroutine gravity_wave_drag_orographic_diagnostics_init(errmsg, errflg)
    use cam_history, only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_add_field('TAUAORO', 'total_stress_at_interface_from_reference_wavenumber_due_to_orographic_gravity_wave_drag', 'ilev', 'inst', 'N m-2')

    call history_add_field('UTGWORO', 'tendency_of_eastward_wind_due_to_orographic_gravity_wave_drag', 'lev', 'avg', 'm s-2')
    call history_add_field('VTGWORO', 'tendency_of_northward_wind_due_to_orographic_gravity_wave_drag', 'lev', 'avg', 'm s-2')
    call history_add_field('TTGWORO', 'tendency_of_air_temperature_due_to_orographic_gravity_wave_drag', 'lev', 'avg', 'K s-1')

    call history_add_field('TTGWSDFORO', 'tendency_of_air_temperature_due_to_diffusion_due_to_orographic_gravity_wave_drag', 'lev', 'avg', 'K s-1')
    call history_add_field('TTGWSKEORO', 'tendency_of_air_temperature_due_to_kinetic_energy_dissipation_due_to_orographic_gravity_wave_drag', 'lev', 'avg', 'K s-1')

    call history_add_field('TAUGWX', 'eastward_stress_at_surface_due_to_orographic_gravity_wave_drag', horiz_only, 'avg', 'N m-2')
    call history_add_field('TAUGWY', 'northward_stress_at_surface_due_to_orographic_gravity_wave_drag', horiz_only, 'avg', 'N m-2')

  end subroutine gravity_wave_drag_orographic_diagnostics_init

!> \section arg_table_gravity_wave_drag_orographic_diagnostics_run  Argument Table
!! \htmlinclude gravity_wave_drag_orographic_diagnostics_run.html
  subroutine gravity_wave_drag_orographic_diagnostics_run( &
    cpair, &
    taua, &
    utgw, vtgw, ttgw, &
    dttdf, dttke, &
    tau0x, tau0y, &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    use cam_history, only: history_out_field

    ! Thermodynamic constant [J kg-1 K-1]
    real(kind_phys), intent(in) :: cpair

    ! Momentum flux profile [N m-2]
    real(kind_phys), intent(in) :: taua(:,:)

    real(kind_phys), intent(in) :: utgw(:,:)
    real(kind_phys), intent(in) :: vtgw(:,:)
    real(kind_phys), intent(in) :: ttgw(:,:) ! [K s-1] converted inside gw scheme using cpairv.

    ! Dry air enthalpy tendencies due to diffusion and KE dissipation [J kg-1 s-1] -- unit converted for output.
    real(kind_phys), intent(in) :: dttdf(:,:)
    real(kind_phys), intent(in) :: dttke(:,:)

    ! Surface stress components [N m-2]
    real(kind_phys), intent(in) :: tau0x(:)
    real(kind_phys), intent(in) :: tau0y(:)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_out_field('TAUAORO', taua)

    call history_out_field('UTGWORO', utgw)
    call history_out_field('VTGWORO', vtgw)
    call history_out_field('TTGWORO', ttgw)

    ! Output temperature tendencies from diffusion (convert from J/kg/s to K/s)
    ! Note: this uses cpair instead of cpairv (inconsistent with ttgw, but consistent with other gw schemes).
    ! Also see atmospheric_physics issue #320.
    call history_out_field('TTGWSDFORO', dttdf / cpair)
    call history_out_field('TTGWSKEORO', dttke / cpair)

    call history_out_field('TAUGWX', tau0x)
    call history_out_field('TAUGWY', tau0y)

  end subroutine gravity_wave_drag_orographic_diagnostics_run

end module gravity_wave_drag_orographic_diagnostics
