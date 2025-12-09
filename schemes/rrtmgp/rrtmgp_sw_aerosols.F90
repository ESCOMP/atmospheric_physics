!> \file rrtmgp_sw_aerosols.F90
!!

!> This module contains the call to the RRTMGP-sw radiation routine
module rrtmgp_sw_aerosols
  implicit none
  private

  public rrtmgp_sw_aerosols_run
contains

!> \section arg_table_rrtmgp_sw_aerosols_run Argument Table
!! \htmlinclude rrtmgp_sw_aerosols_run.html
!!
   subroutine rrtmgp_sw_aerosols_run(doswrad, nday, aer_sw, errmsg, errflg)
    use ccpp_optical_props, only: ty_optical_props_2str_ccpp
    use ccpp_kinds,         only: kind_phys

    ! Inputs
    logical, intent(in) :: doswrad                                              !< Flag to perform shortwave calculation
    integer, intent(in) :: nday                                                 !< daytime points dimension

    ! Outputs
    class(ty_optical_props_2str_ccpp), intent(inout) :: aer_sw                  !< Aerosol optical properties object

    character(len=*), intent(out) :: errmsg                                     !< CCPP error message
    integer,          intent(out) :: errflg                                     !< CCPP error flag

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doswrad .or. nday == 0) return

    ! REMOVECAM: scheme is just a stub for now
    ! Set the optical properties to constants for snapshot testing
    aer_sw%optical_props%tau = 0.0_kind_phys
    aer_sw%optical_props%g = 0.0_kind_phys
    aer_sw%optical_props%ssa = 1.0_kind_phys

  end subroutine rrtmgp_sw_aerosols_run
end module rrtmgp_sw_aerosols
