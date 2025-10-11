!-------------------------------------------------------------------------------
! This module uses the solar irradiance data 
! to provide a spectral scaling factor
! to approximate the spectral distribution of irradiance
! when the radiation scheme might use a different solar source function
!-------------------------------------------------------------------------------
module rrtmgp_sw_solar_var_setup

  use ccpp_kinds,        only : kind_phys

  implicit none
  save

  private
  public :: rrtmgp_sw_solar_var_setup_init

  real(kind_phys), public, allocatable :: irrad(:)           ! solar irradiance at model timestep in each band

  real(kind_phys), public, allocatable :: radbinmax(:)
  real(kind_phys), public, allocatable :: radbinmin(:)

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> \section arg_table_rrtmgp_sw_solar_var_setup_init Argument Table
!! \htmlinclude rrtmgp_sw_solar_var_setup_init.html
!!
  subroutine rrtmgp_sw_solar_var_setup_init(nswbands, do_spectral_scaling, has_spectrum, errmsg, errflg)
    use radiation_utils,   only : get_sw_spectral_boundaries_ccpp
    integer, intent(in) :: nswbands            ! number of shortwave bands
    logical, intent(in) :: do_spectral_scaling   ! flag to do spectral scaling
    logical, intent(in) :: has_spectrum        ! flag for whether solar input file has irradiance spectrum
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: radmax_loc
    character(len=256) :: alloc_errmsg

    if ( do_spectral_scaling ) then

       if ( .not.has_spectrum ) then
          write(errmsg, *) 'rrtmgp_sw_solar_var_setup_init: solar input fil must have irradiance spectrum'
          errflg = 1
          return
       endif

       allocate (radbinmax(nswbands),stat=errflg,errmsg=alloc_errmsg)
       if (errflg /= 0) then
          write(errmsg,*) 'rrtmgp_sw_solar_var_setup_init: Error allocating space for radbinmax - message: ', alloc_errmsg
          return
       end if

       allocate (radbinmin(nswbands),stat=errflg,errmsg=alloc_errmsg)
       if (errflg /= 0) then
          write(errmsg,*) 'rrtmgp_sw_solar_var_setup_init: Error allocating space for radbinmin - message: ', alloc_errmsg
          return
       end if

       allocate (irrad(nswbands), stat=errflg, errmsg=alloc_errmsg)
       if (errflg /= 0) then
          write(errmsg,*) 'rrtmgp_sw_solar_var_setup_init: Error allocating space for irrad - message: ', alloc_errmsg
          return
       end if

       call get_sw_spectral_boundaries_ccpp(radbinmin, radbinmax, 'nm', errmsg, errflg)
       if (errflg /= 0) then
          return
       end if

       ! Make sure that the far-IR is included, even if radiation grid does not
       ! extend that far down. 10^5 nm corresponds to a wavenumber of
       ! 100 cm^-1.
       radmax_loc = maxloc(radbinmax,1)
       radbinmax(radmax_loc) = max(100000._kind_phys,radbinmax(radmax_loc))

    endif

  end subroutine rrtmgp_sw_solar_var_setup_init

end module rrtmgp_sw_solar_var_setup
