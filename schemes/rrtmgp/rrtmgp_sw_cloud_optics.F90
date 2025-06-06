! PEVERWHEE - dependencies = interpolate_data
!> \file rrtmgp_sw_cloud_optics.F90
!!

!> This module contains two routines: The first initializes data and functions
!! needed to compute the longwave cloud radiative properties in RRTMGP. The second routine
!! is a ccpp scheme within the "radiation loop", where the shortwave optical properties
!! (optical-depth, single-scattering albedo, asymmetry parameter) are computed for ALL
!! cloud types visible to RRTMGP.
module rrtmgp_sw_cloud_optics
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  public :: rrtmgp_sw_cloud_optics_run

contains

  ! SUBROUTINE rrtmgp_sw_cloud_optics_run()
  ! ######################################################################################
!> \section arg_table_rrtmgp_sw_cloud_optics_run Argument Table
!! \htmlinclude rrtmgp_sw_cloud_optics_run.html
!!
  subroutine rrtmgp_sw_cloud_optics_run(dosw, ncol, nlay, kdist_sw, cloud_sw, errmsg, errflg)
    use ccpp_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp_ccpp
    use ccpp_optical_props,       only: ty_optical_props_1scl_ccpp
    ! Compute combined cloud optical properties
    ! Create MCICA stochastic arrays for cloud sw optical properties
    ! Initialize optical properties object (cloud_sw) and load with MCICA columns

    ! Inputs
    integer,                           intent(in) :: ncol             ! Number of columns
    integer,                           intent(in) :: nlay             ! Number of vertical layers in radiation
    logical,                           intent(in) :: dosw             ! Flag for whether to perform longwave calculation
    class(ty_gas_optics_rrtmgp_ccpp),  intent(in) :: kdist_sw         ! Longwave gas optics object

    ! Outputs
    type(ty_optical_props_1scl_ccpp),  intent(out) :: cloud_sw        ! Longwave cloud optics object
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    character(len=*), parameter :: sub = 'rrtmgp_sw_cloud_optics_run'
    !--------------------------------------------------------------------------------

    ! Set error variables
    errmsg = ''
    errflg = 0

    ! If not doing longwave, no need to proceed
    if (.not. dosw) then
       return
    end if

    errmsg =cloud_sw%optical_props%alloc_1scl(ncol, nlay, kdist_sw%gas_props)
    if (len_trim(errmsg) > 0) then
       errflg = 1
       return
    end if

  end subroutine rrtmgp_sw_cloud_optics_run

end module rrtmgp_sw_cloud_optics
