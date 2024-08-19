module tropopause_diagnostics
  ! ... output tropopause diagnostics
  ! this will be moved to cam_diagnostics when History is available. (hplin, 8/14/24)
  ! Currently stubbed out
  use ccpp_kinds,           only : kind_phys

  implicit none
  private

  ! CCPP-compliant subroutines
  public :: tropopause_diagnostics_init
  public :: tropopause_diagnostics_run

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
    call history_add_field('TROP_P',   'tropopause_air_pressure',    horiz_only,  'avg', 'Pa')
    call history_add_field('TROP_T',   'tropopause_air_temperature', horiz_only,  'avg', 'K' )
    call history_add_field('TROP_Z',   'tropopause_altitude',        horiz_only,  'avg', 'm' )
    call history_add_field('TROP_DZ',  'tropopause_altitude_relative', 'lev',       'avg', 'm')
    call history_add_field('TROP_PD',  'probability_distribution_of_model_level_number_at_tropopause', 'lev', 'avg', 'probability')
    call history_add_field('TROP_FD',  'tropopause_found', horiz_only,  'avg', 'probability')

    call history_add_field('TROPF_P',  'tropopause_air_pressure_assuming_cold_point',    horiz_only,  'avg', 'Pa')
    call history_add_field('TROPF_T',  'tropopause_air_temperature_assuming_cold_point', horiz_only,  'avg', 'K' )
    call history_add_field('TROPF_Z',  'tropopause_altitude_assuming_cold_point',        horiz_only,  'avg', 'm' )
    call history_add_field('TROPF_DZ', 'tropopause_altitude_relative_assuming_cold_point', 'lev',       'avg', 'm')
    call history_add_field('TROPF_PD', 'probability_distribution_of_model_level_number_at_tropopause_assuming_cold_point', 'lev', 'avg', 'probability')
    call history_add_field('TROPF_FD', 'tropopause_found_assuming_cold_point', horiz_only,  'avg', 'probability')

    call history_add_field('TROPC_P',  'tropopause_air_pressure_assuming_climatology',    horiz_only,  'avg', 'Pa')
    call history_add_field('TROPC_T',  'tropopause_air_temperature_assuming_climatology', horiz_only,  'avg', 'K' )
    call history_add_field('TROPC_Z',  'tropopause_altitude_assuming_climatology',        horiz_only,  'avg', 'm' )
    call history_add_field('TROPC_DZ', 'tropopause_altitude_relative_assuming_climatology', 'lev', 'avg', 'm')
    call history_add_field('TROPC_PD', 'probability_distribution_of_model_level_number_at_tropopause_assuming_climatology', 'lev', 'avg', 'probability')
    call history_add_field('TROPC_FD', 'tropopause_found_assuming_cold_point', horiz_only,  'avg', 'probability')

    ! Hybridstobie output fields
    call history_add_field('hstobie_trop', 'lower_bound_of_model_level_number_for_stratospheric_chemistry', 'lev', 'inst', 'fraction of model time')
    call history_add_field('hstobie_linoz', 'lower_bound_of_model_level_number_for_linoz_chemistry', 'lev', 'inst', 'fraction of model time')
    call history_add_field('hstobie_tropop', 'model_level_number_at_tropopause_for_chemistry', 'lev', 'inst', 'fraction of model time')

  end subroutine tropopause_diagnostics_init

  ! Output the tropopause pressure and temperature to the history files. Two sets
  ! of output will be generated, one for the default algorithm and another one
  ! using the default routine, but backed by a climatology when the default
  ! algorithm fails.
!> \section arg_table_tropopause_diagnostics_run Argument Table
!! \htmlinclude tropopause_diagnostics_run.html
  subroutine tropopause_diagnostics_run(ncol, pver, &
                                        zm, &
                                        tropLev, tropP, tropT, tropZ, & ! Default primary+backup ( twmo+climate)
                                        tropLev_clim, tropP_clim, tropT_clim, tropZ_clim, & ! Climate-only
                                        tropLev_hybstob, tropP_hybstob, tropT_hybstob, tropZ_hybstob, & !      Hybridstobie + climate backup
                                        tropLev_cpp, tropP_cpp, tropT_cpp, tropZ_cpp, & ! Cold point only
                                        hstobie_trop, hstobie_linoz, hstobie_tropop, & ! Hybridstobie only for chemistry diagnostics
                                        errmsg, errflg)
    use cam_history, only: history_out_field

    integer,         intent(in)      :: ncol               ! Number of atmospheric columns
    integer,         intent(in)      :: pver               ! Number of vertical levels

    real(kind_phys), intent(in)      :: zm(:,:)            ! Geopotential height above surface at midpoints (m), pver

    integer,         intent(out)     :: tropLev(:)         ! tropopause level index
    real(kind_phys), intent(out)     :: tropP(:)           ! tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT(:)           ! tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ(:)           ! tropopause height (m)

    integer,         intent(out)     :: tropLev_clim(:)    ! climatology-backed tropopause level index
    real(kind_phys), intent(out)     :: tropP_clim(:)      ! climatology-backed tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT_clim(:)      ! climatology-backed tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ_clim(:)      ! climatology-backed tropopause height (m)

    integer,         intent(out)     :: tropLev_hybstob(:) ! hybridstobie climatology-backed tropopause level index
    real(kind_phys), intent(out)     :: tropP_hybstob(:)   ! hybridstobie climatology-backed tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT_hybstob(:)   ! hybridstobie climatology-backed tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ_hybstob(:)   ! hybridstobie climatology-backed tropopause height (m)

    integer,         intent(out)     :: tropLev_cpp(:)     ! cold point tropopause level index
    real(kind_phys), intent(out)     :: tropP_cpp(:)       ! cold point tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT_cpp(:)       ! cold point tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ_cpp(:)       ! cold point tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(kind_phys), intent(out)   :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(kind_phys), intent(out)   :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(kind_phys), intent(out)   :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! Local Variables
    integer       :: i
    integer       :: alg
    integer       :: tropLev(ncol)           ! tropopause level index
    real(kind_phys)      :: tropP(ncol)             ! tropopause pressure (Pa)
    real(kind_phys)      :: tropT(ncol)             ! tropopause temperature (K)
    real(kind_phys)      :: tropZ(ncol)             ! tropopause height (m)
    real(kind_phys)      :: tropFound(ncol)         ! tropopause found
    real(kind_phys)      :: tropDZ(ncol, pver)      ! relative tropopause height (m)
    real(kind_phys)      :: tropPdf(ncol, pver)     ! tropopause probability distribution

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

    return
  end subroutine tropopause_diagnostics_run
end module tropopause_diagnostics