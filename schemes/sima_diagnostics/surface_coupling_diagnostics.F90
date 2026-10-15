! Surface coupling diagnostics:
! precipitation rates (PRECC/PRECL/PRECT/PRECSC/PRECSL) and the surface fluxes received from the coupler (SHFLX/LHFLX/QFLX/TAUX/TAUY).
!
! Precipitation uses the same per-process rates sent to the coupler by set_surface_coupling_vars
! to match what the surface components receive.
module surface_coupling_diagnostics
  implicit none
  private

  public :: surface_coupling_diagnostics_init
  public :: surface_coupling_diagnostics_run

  ! Index of water vapor in the constituent array.
  ! QFLX is only written when it is present.
  integer :: wv_idx

contains

!> \section arg_table_surface_coupling_diagnostics_init  Argument Table
!! \htmlinclude arg_table_surface_coupling_diagnostics_init.html
  subroutine surface_coupling_diagnostics_init(errmsg, errflg)
    use cam_history,               only: history_add_field
    use cam_history_support,       only: horiz_only
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: int_unassigned

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

    ! Surface fluxes received from the coupler
    call history_add_field('SHFLX',  'surface_upward_sensible_heat_flux_from_coupler', horiz_only, 'avg', 'W m-2')
    call history_add_field('LHFLX',  'surface_upward_latent_heat_flux_from_coupler', horiz_only, 'avg', 'W m-2')
    call history_add_field('TAUX',   'surface_eastward_wind_stress_from_coupler', horiz_only, 'avg', 'N m-2')
    call history_add_field('TAUY',   'surface_northward_wind_stress_from_coupler', horiz_only, 'avg', 'N m-2')

    call ccpp_constituent_index('water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', wv_idx, errflg, errmsg)
    if (errflg /= 0) then
      return
    end if
    if (wv_idx /= int_unassigned) then
      call history_add_field('QFLX', 'surface_upward_water_vapor_flux_from_coupler', horiz_only, 'avg', 'kg m-2 s-1')
    end if

  end subroutine surface_coupling_diagnostics_init

!> \section arg_table_surface_coupling_diagnostics_run  Argument Table
!! \htmlinclude arg_table_surface_coupling_diagnostics_run.html
  subroutine surface_coupling_diagnostics_run(ncol, prec_dp, snow_dp, prec_sh, snow_sh, &
      prec_str, snow_str, shf, lhf, cflx, wsx, wsy, errmsg, errflg)
    use ccpp_kinds,  only: kind_phys
    use cam_history, only: history_out_field

    integer,           intent(in)  :: ncol
    real(kind_phys),   intent(in)  :: prec_dp(:)   ! Deep convective precipitation rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: snow_dp(:)   ! Deep convective snow rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: prec_sh(:)   ! Shallow convective precipitation rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: snow_sh(:)   ! Shallow convective snow rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: prec_str(:)  ! Stratiform precipitation rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: snow_str(:)  ! Stratiform snow rate at surface [m s-1]
    real(kind_phys),   intent(in)  :: shf(:)       ! Surface upward sensible heat flux [W m-2]
    real(kind_phys),   intent(in)  :: lhf(:)       ! Surface upward latent heat flux [W m-2]
    real(kind_phys),   intent(in)  :: cflx(:,:)    ! Surface upward constituent fluxes [kg m-2 s-1]
    real(kind_phys),   intent(in)  :: wsx(:)       ! Surface eastward wind stress [N m-2]
    real(kind_phys),   intent(in)  :: wsy(:)       ! Surface northward wind stress [N m-2]
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

    call history_out_field('SHFLX',  shf)
    call history_out_field('LHFLX',  lhf)
    call history_out_field('TAUX',   wsx)
    call history_out_field('TAUY',   wsy)
    if (wv_idx /= int_unassigned) then
      call history_out_field('QFLX', cflx(:, wv_idx))
    end if

  end subroutine surface_coupling_diagnostics_run

end module surface_coupling_diagnostics
