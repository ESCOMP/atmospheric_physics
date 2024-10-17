module tropopause_diagnostics
  ! ... output tropopause diagnostics within CAM-SIMA
  use ccpp_kinds,           only: kind_phys

  implicit none
  private

  ! CCPP-compliant subroutines
  public :: tropopause_diagnostics_init
  public :: tropopause_diagnostics_run

  integer, parameter :: NOTFOUND = -1

contains
  ! Initialize the output history fields.
!> \section arg_table_tropopause_diagnostics_init Argument Table
!! \htmlinclude tropopause_diagnostics_init.html
  subroutine tropopause_diagnostics_init(errmsg, errflg)

    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! Define the output fields.

    ! Primary (Lapse rate) + backup (climatology) method - will always find a tropopause
    call history_add_field('TROP_P',   'tropopause_air_pressure',    horiz_only,  'avg', 'Pa')
    call history_add_field('TROP_T',   'tropopause_air_temperature', horiz_only,  'avg', 'K' )
    call history_add_field('TROP_Z',   'tropopause_geopotential_height_wrt_surface', horiz_only,  'avg', 'm' )
    call history_add_field('TROP_DZ',  'geopotential_height_difference_between_atmosphere_layer_and_tropopause', 'lev', 'avg', 'm')
    call history_add_field('TROP_PD',  'probability_distribution_of_tropopause_vertical_layer_index', 'lev', 'avg', 'probability')
    call history_add_field('TROP_FD',  'tropopause_found', horiz_only,  'avg', 'probability')

    ! Primary (Lapse rate) only
    call history_add_field('TROPP_P',  'tropopause_air_pressure_from_lapse_rate_method', horiz_only,  'avg', 'Pa')
    call history_add_field('TROPP_T',  'tropopause_air_temperature_from_lapse_rate_method', horiz_only,  'avg', 'K' )
    call history_add_field('TROPP_Z',  'tropopause_geopotential_height_wrt_surface_from_lapse_rate_method', horiz_only, 'avg', 'm' )
    call history_add_field('TROPP_DZ', 'geopotential_height_difference_between_atmosphere_layer_and_tropopause_from_lapse_rate_method', 'lev', 'avg', 'm')
    call history_add_field('TROPP_PD', 'probability_distribution_of_tropopause_vertical_layer_index_from_lapse_rate_method', 'lev', 'avg', 'probability')
    call history_add_field('TROPP_FD', 'tropopause_found_from_lapse_rate_method', horiz_only,  'avg', 'probability')

    ! Cold point (CPP) only
    call history_add_field('TROPF_P',  'tropopause_air_pressure_from_cold_point_method', horiz_only,  'avg', 'Pa')
    call history_add_field('TROPF_T',  'tropopause_air_temperature_from_cold_point_method', horiz_only,  'avg', 'K' )
    call history_add_field('TROPF_Z',  'tropopause_geopotential_height_wrt_surface_from_cold_point_method', horiz_only,  'avg', 'm' )
    call history_add_field('TROPF_DZ', 'geopotential_height_difference_between_atmosphere_layer_and_tropopause_from_cold_point_method', 'lev', 'avg', 'm')
    call history_add_field('TROPF_PD', 'probability_distribution_of_tropopause_vertical_layer_index_from_cold_point_method', 'lev', 'avg', 'probability')
    call history_add_field('TROPF_FD', 'tropopause_found_from_cold_point_method', horiz_only,  'avg', 'probability')

    ! Climatology only - will always find a tropopause
    call history_add_field('TROPC_P',  'tropopause_air_pressure_from_climatological_method', horiz_only,  'avg', 'Pa')
    call history_add_field('TROPC_T',  'tropopause_air_temperature_from_climatological_method', horiz_only,  'avg', 'K' )
    call history_add_field('TROPC_Z',  'tropopause_geopotential_height_wrt_surface_from_climatological_method', horiz_only, 'avg', 'm' )
    call history_add_field('TROPC_DZ', 'geopotential_height_difference_between_atmosphere_layer_and_tropopause_from_climatological_method', 'lev', 'avg', 'm')
    call history_add_field('TROPC_PD', 'probability_distribution_of_tropopause_vertical_layer_index_from_climatological_method', 'lev', 'avg', 'probability')
    call history_add_field('TROPC_FD', 'tropopause_found_from_climatological_method', horiz_only,  'avg', 'probability')

    ! Hybridstobie output fields
    call history_add_field('hstobie_trop', 'vertical_layer_index_lower_bound_from_hybrid_stobie_linoz_with_climatological_backup_method_for_stratospheric_chemistry', 'lev', 'inst', 'fraction of model time')
    call history_add_field('hstobie_linoz', 'vertical_layer_index_lower_bound_from_hybrid_stobie_linoz_with_climatological_backup_method_for_linearized_ozone_chemistry', 'lev', 'inst', 'fraction of model time')
    call history_add_field('hstobie_tropop', 'tropopause_vertical_layer_index_from_hybrid_stobie_linoz_with_climatological_backup_method_for_chemistry', 'lev', 'inst', 'fraction of model time')

  end subroutine tropopause_diagnostics_init

  ! Output the tropopause pressure and temperature to the history files. Two sets
  ! of output will be generated, one for the default algorithm and another one
  ! using the default routine, but backed by a climatology when the default
  ! algorithm fails.
!> \section arg_table_tropopause_diagnostics_run Argument Table
!! \htmlinclude tropopause_diagnostics_run.html
  subroutine tropopause_diagnostics_run(ncol, pver, &
                                        zm, &
                                        tropLev, tropP, tropT, tropZ, & ! Default primary+backup (twmo+climate)
                                        tropLev_twmo, tropP_twmo, tropT_twmo, tropZ_twmo, & ! Primary only (twmo)
                                        tropLev_clim, tropP_clim, tropT_clim, tropZ_clim, & ! Climate-only
                                        tropLev_hybstob, tropP_hybstob, tropT_hybstob, tropZ_hybstob, & !      Hybridstobie + climate backup
                                        tropLev_cpp, tropP_cpp, tropT_cpp, tropZ_cpp, & ! Cold point only
                                        hstobie_trop, hstobie_linoz, hstobie_tropop, & ! Hybridstobie only for chemistry diagnostics
                                        errmsg, errflg)
    use cam_history,          only: history_out_field
    use cam_history_support,  only: fillvalue

    integer,         intent(in)      :: ncol               ! Number of atmospheric columns
    integer,         intent(in)      :: pver               ! Number of vertical levels

    real(kind_phys), intent(in)      :: zm(:,:)            ! Geopotential height above surface at midpoints (m), pver

    integer,         intent(in)      :: tropLev(:)         ! tropopause level index
    real(kind_phys), intent(in)      :: tropP(:)           ! tropopause pressure (Pa)
    real(kind_phys), intent(in)      :: tropT(:)           ! tropopause temperature (K)
    real(kind_phys), intent(in)      :: tropZ(:)           ! tropopause height (m)

    integer,         intent(in)      :: tropLev_twmo(:)    ! lapse-rate tropopause level index
    real(kind_phys), intent(in)      :: tropP_twmo(:)      ! lapse-rate tropopause pressure (Pa)
    real(kind_phys), intent(in)      :: tropT_twmo(:)      ! lapse-rate tropopause temperature (K)
    real(kind_phys), intent(in)      :: tropZ_twmo(:)      ! lapse-rate tropopause height (m)

    integer,         intent(in)      :: tropLev_clim(:)    ! climatology-backed tropopause level index
    real(kind_phys), intent(in)      :: tropP_clim(:)      ! climatology-backed tropopause pressure (Pa)
    real(kind_phys), intent(in)      :: tropT_clim(:)      ! climatology-backed tropopause temperature (K)
    real(kind_phys), intent(in)      :: tropZ_clim(:)      ! climatology-backed tropopause height (m)

    integer,         intent(in)      :: tropLev_hybstob(:) ! hybridstobie climatology-backed tropopause level index
    real(kind_phys), intent(in)      :: tropP_hybstob(:)   ! hybridstobie climatology-backed tropopause pressure (Pa)
    real(kind_phys), intent(in)      :: tropT_hybstob(:)   ! hybridstobie climatology-backed tropopause temperature (K)
    real(kind_phys), intent(in)      :: tropZ_hybstob(:)   ! hybridstobie climatology-backed tropopause height (m)

    integer,         intent(in)      :: tropLev_cpp(:)     ! cold point tropopause level index
    real(kind_phys), intent(in)      :: tropP_cpp(:)       ! cold point tropopause pressure (Pa)
    real(kind_phys), intent(in)      :: tropT_cpp(:)       ! cold point tropopause temperature (K)
    real(kind_phys), intent(in)      :: tropZ_cpp(:)       ! cold point tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(kind_phys), intent(in)      :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(kind_phys), intent(in)      :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(kind_phys), intent(in)      :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local Variables
    integer       :: i
    integer       :: alg
    real(kind_phys)      :: tropFound(ncol)         ! tropopause found
    real(kind_phys)      :: tropDZ(ncol, pver)      ! relative tropopause height (m)
    real(kind_phys)      :: tropPdf(ncol, pver)     ! tropopause probability distribution

    errmsg = ' '
    errflg = 0

    ! Default algorithm output
    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue
    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ(i)
      end if
    end do

    call history_out_field('TROP_P',  tropP)
    call history_out_field('TROP_T',  tropT)
    call history_out_field('TROP_Z',  tropZ)
    call history_out_field('TROP_DZ', tropDZ)
    call history_out_field('TROP_PD', tropPdf)
    call history_out_field('TROP_FD', tropFound)

    ! Primary-only (currently TWMO) algorithm output
    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue
    do i = 1, ncol
      if (tropLev_twmo(i) /= NOTFOUND) then
        tropPdf(i, tropLev_twmo(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ(i)
      end if
    end do

    call history_out_field('TROPP_P',  tropP_twmo)
    call history_out_field('TROPP_T',  tropT_twmo)
    call history_out_field('TROPP_Z',  tropZ_twmo)
    call history_out_field('TROPP_DZ', tropDZ)
    call history_out_field('TROPP_PD', tropPdf)
    call history_out_field('TROPP_FD', tropFound)

    ! Cold point output
    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue
    do i = 1, ncol
      if (tropLev_cpp(i) /= NOTFOUND) then
        tropPdf(i, tropLev_cpp(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ_cpp(i)
      end if
    end do

    call history_out_field('TROPF_P',  tropP_cpp)
    call history_out_field('TROPF_T',  tropT_cpp)
    call history_out_field('TROPF_Z',  tropZ_cpp)
    call history_out_field('TROPF_DZ', tropDZ)
    call history_out_field('TROPF_PD', tropPdf)
    call history_out_field('TROPF_FD', tropFound)

    ! Climatology output
    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue
    do i = 1, ncol
      if (tropLev_clim(i) /= NOTFOUND) then
        tropPdf(i, tropLev_clim(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ_clim(i)
      end if
    end do

    call history_out_field('TROPC_P',  tropP_clim)
    call history_out_field('TROPC_T',  tropT_clim)
    call history_out_field('TROPC_Z',  tropZ_clim)
    call history_out_field('TROPC_DZ', tropDZ)
    call history_out_field('TROPC_PD', tropPdf)
    call history_out_field('TROPC_FD', tropFound)

    ! Hybridstobie outputs for chemistry
    call history_out_field('hstobie_trop',   hstobie_trop)
    call history_out_field('hstobie_linoz',  hstobie_linoz)
    call history_out_field('hstobie_tropop', hstobie_tropop)

  end subroutine tropopause_diagnostics_run
end module tropopause_diagnostics
