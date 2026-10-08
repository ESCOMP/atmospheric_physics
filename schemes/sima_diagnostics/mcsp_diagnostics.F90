! Diagnostic scheme for the mesoscale coherent structure parameterization (MCSP)
module mcsp_diagnostics
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mcsp_diagnostics_init
  public :: mcsp_diagnostics_run

contains

!> \section arg_table_mcsp_diagnostics_init Argument Table
!! \htmlinclude mcsp_diagnostics_init.html
  subroutine mcsp_diagnostics_init(errmsg, errflg)
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_add_field('MCSP_DT', 'tendency_of_air_temperature_due_to_mesoscale_coherent_structure_parameterization', &
                           'lev', 'avg', 'K s-1')
    call history_add_field('MCSP_DQ', &
                           'tendency_of_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_mesoscale_coherent_structure_parameterization', &
                           'lev', 'avg', 'kg kg-1 s-1')
    call history_add_field('MCSP_DU', 'tendency_of_eastward_wind_due_to_mesoscale_coherent_structure_parameterization', &
                           'lev', 'avg', 'm s-2')
    call history_add_field('MCSP_DV', 'tendency_of_northward_wind_due_to_mesoscale_coherent_structure_parameterization', &
                           'lev', 'avg', 'm s-2')
    call history_add_field('MCSP_DT_max', 'heating_amplitude_of_mesoscale_coherent_structure_parameterization', &
                           horiz_only, 'avg', 'K s-1')
    call history_add_field('MCSP_freq', 'frequency_of_activation_of_mesoscale_coherent_structure_parameterization', &
                           horiz_only, 'avg', '1')
    call history_add_field('MCSP_shear', 'low_level_eastward_wind_shear_for_mesoscale_coherent_structure_parameterization', &
                           horiz_only, 'avg', 'm s-1')
    call history_add_field('MCSP_conv_depth', 'pressure_depth_of_deep_convection', &
                           horiz_only, 'avg', 'Pa')

  end subroutine mcsp_diagnostics_init

!> \section arg_table_mcsp_diagnostics_run Argument Table
!! \htmlinclude mcsp_diagnostics_run.html
  subroutine mcsp_diagnostics_run(mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
                                  mcsp_dt_max, mcsp_freq, mcsp_shear, conv_depth, errmsg, errflg)
    use cam_history, only: history_out_field

    real(kind_phys),    intent(in)  :: mcsp_dt_out(:,:)  ! MCSP temperature tendency [K s-1]
    real(kind_phys),    intent(in)  :: mcsp_dq_out(:,:)  ! MCSP water vapor tendency [kg kg-1 s-1]
    real(kind_phys),    intent(in)  :: mcsp_du_out(:,:)  ! MCSP zonal wind tendency [m s-2]
    real(kind_phys),    intent(in)  :: mcsp_dv_out(:,:)  ! MCSP meridional wind tendency [m s-2]
    real(kind_phys),    intent(in)  :: mcsp_dt_max(:)    ! MCSP heating amplitude [K s-1]
    real(kind_phys),    intent(in)  :: mcsp_freq(:)      ! MCSP frequency of activation [1]
    real(kind_phys),    intent(in)  :: mcsp_shear(:)     ! low-level zonal wind shear [m s-1]
    real(kind_phys),    intent(in)  :: conv_depth(:)     ! pressure depth of deep convection [Pa]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_out_field('MCSP_DT',         mcsp_dt_out)
    call history_out_field('MCSP_DQ',         mcsp_dq_out)
    call history_out_field('MCSP_DU',         mcsp_du_out)
    call history_out_field('MCSP_DV',         mcsp_dv_out)
    call history_out_field('MCSP_DT_max',     mcsp_dt_max)
    call history_out_field('MCSP_freq',       mcsp_freq)
    call history_out_field('MCSP_shear',      mcsp_shear)
    call history_out_field('MCSP_conv_depth', conv_depth)

  end subroutine mcsp_diagnostics_run

end module mcsp_diagnostics
