! Surface precipitation rate diagnostics: PRECC/PRECL/PRECT
!
! Uses the same per-process rates sent to the coupler by set_surface_coupling_vars
! to match what the surface components receives.
module precipitation_diagnostics
  implicit none
  private

  public :: precipitation_diagnostics_init
  public :: precipitation_diagnostics_run

contains

!> \section arg_table_precipitation_diagnostics_init  Argument Table
!! \htmlinclude precipitation_diagnostics_init.html
  subroutine precipitation_diagnostics_init(errmsg, errflg)
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=*),  intent(out) :: errmsg
    integer,           intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Convective precipitation rate (liq + ice)
    call history_add_field('PRECC',  'lwe_convective_precipitation_rate_at_surface', horiz_only, 'avg', 'm s-1')
    ! Large-scale (stable) precipitation rate (liq + ice)
    call history_add_field('PRECL',  'lwe_large_scale_precipitation_rate_at_surface', horiz_only, 'avg', 'm s-1')
    ! Total (convective and large-scale) precipitation rate (liq + ice)
    call history_add_field('PRECT',  'total_precipitation_rate_at_surface', horiz_only, 'avg', 'm s-1')
    ! Convective snow rate (water equivalent)
    call history_add_field('PRECSC', 'lwe_convective_snowfall_rate_at_surface', horiz_only, 'avg', 'm s-1')
    ! Large-scale (stable) snow rate (water equivalent)
    call history_add_field('PRECSL', 'lwe_snow_and_cloud_ice_precipitation_rate_at_surface_due_to_microphysics', horiz_only, 'avg', 'm s-1')

  end subroutine precipitation_diagnostics_init

!> \section arg_table_precipitation_diagnostics_run  Argument Table
!! \htmlinclude precipitation_diagnostics_run.html
  subroutine precipitation_diagnostics_run(ncol, prec_dp, snow_dp, prec_sh, snow_sh, &
      prec_str, snow_str, errmsg, errflg)
    use ccpp_kinds,  only: kind_phys
    use cam_history, only: history_out_field

    integer,           intent(in)  :: ncol
    real(kind_phys),   intent(in)  :: prec_dp(:)   ! Deep convective precipitation rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: snow_dp(:)   ! Deep convective snow rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: prec_sh(:)   ! Shallow convective precipitation rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: snow_sh(:)   ! Shallow convective snow rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: prec_str(:)  ! Stratiform precipitation rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: snow_str(:)  ! Stratiform snow rate at surface [m s-1]
    character(len=*),  intent(out) :: errmsg
    integer,           intent(out) :: errflg

    ! Local variables
    real(kind_phys) :: precc(ncol)   ! Total convective precipitation rate [m s-1]
    real(kind_phys) :: precsc(ncol)  ! Total convective snow rate [m s-1]

    errmsg = ''
    errflg = 0

    precc(:)  = prec_dp(:) + prec_sh(:)
    precsc(:) = snow_dp(:) + snow_sh(:)

    call history_out_field('PRECC',  precc)
    call history_out_field('PRECL',  prec_str)
    call history_out_field('PRECT',  precc(:) + prec_str(:))
    call history_out_field('PRECSC', precsc)
    call history_out_field('PRECSL', snow_str)

  end subroutine precipitation_diagnostics_run

end module precipitation_diagnostics
