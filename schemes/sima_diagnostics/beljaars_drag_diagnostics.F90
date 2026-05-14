! Diagnostic scheme for Beljaars SGO form drag
module beljaars_drag_diagnostics
  implicit none
  private

  public :: beljaars_drag_diagnostics_init
  public :: beljaars_drag_diagnostics_run

contains

!> \section arg_table_beljaars_drag_diagnostics_init  Argument Table
!! \htmlinclude beljaars_drag_diagnostics_init.html
  subroutine beljaars_drag_diagnostics_init(errmsg, errflg)
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! 3D field (layer centers)
    call history_add_field('DRAGBLJ', 'turbulent_orographic_form_drag_coefficent', &
                           'lev', 'inst', 's-1')

    ! 2D fields (horizontal only)
    call history_add_field('TAUBLJX', 'eastward_beljaars_surface_stress_tbd', &
                           horiz_only, 'inst', 'N m-2')
    call history_add_field('TAUBLJY', 'northward_beljaars_surface_stress_tbd', &
                           horiz_only, 'inst', 'N m-2')

  end subroutine beljaars_drag_diagnostics_init

!> \section arg_table_beljaars_drag_diagnostics_run  Argument Table
!! \htmlinclude beljaars_drag_diagnostics_run.html
  subroutine beljaars_drag_diagnostics_run( &
    drag, taux, tauy, &
    errmsg, errflg)
    use ccpp_kinds,  only: kind_phys
    use cam_history, only: history_out_field

    real(kind_phys),    intent(in)  :: drag(:, :)   ! SGO drag profile [s-1]
    real(kind_phys),    intent(in)  :: taux(:)      ! surface zonal wind stress [N m-2]
    real(kind_phys),    intent(in)  :: tauy(:)      ! surface meridional wind stress [N m-2]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_out_field('DRAGBLJ', drag)
    call history_out_field('TAUBLJX', taux)
    call history_out_field('TAUBLJY', tauy)

  end subroutine beljaars_drag_diagnostics_run

end module beljaars_drag_diagnostics
