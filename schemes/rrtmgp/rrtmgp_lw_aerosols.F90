!> \file rrtmgp_lw_aerosols.F90
!!

!> This module contains the call to the RRTMGP-lw radiation routine
module rrtmgp_lw_aerosols
  implicit none
  private

  public rrtmgp_lw_aerosols_run
contains

!> \section arg_table_rrtmgp_lw_aerosols_run Argument Table
!! \htmlinclude rrtmgp_lw_aerosols_run.html
!!
   subroutine rrtmgp_lw_aerosols_run(dolwrad, aer_lw, errmsg, errflg)
    use ccpp_optical_props, only: ty_optical_props_1scl_ccpp
    use ccpp_kinds,         only: kind_phys

    ! Inputs
    logical, intent(in) :: dolwrad                                              !< Flag to perform longwave calculation

    ! Outputs
    class(ty_optical_props_1scl_ccpp), intent(inout) :: aer_lw                  !< Aerosol optical properties object

    character(len=*), intent(out) :: errmsg                                     !< CCPP error message
    integer,          intent(out) :: errflg                                     !< CCPP error flag

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. dolwrad) return

    aer_lw%optical_props%tau = 0.0_kind_phys

  end subroutine rrtmgp_lw_aerosols_run
end module rrtmgp_lw_aerosols
