module rrtmgp_sw_cloud_temp

   public :: rrtmgp_sw_cloud_temp_run

CONTAINS
  !> \section arg_table_rrtmgp_sw_cloud_temp_run Argument Table
  !! \htmlinclude rrtmgp_sw_cloud_temp_run.html
  subroutine rrtmgp_sw_cloud_temp_run(dosw, ncol, nlay, kdist_sw, cloud_sw, errmsg, errflg)
     use ccpp_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp_ccpp
     use ccpp_optical_props,       only: ty_optical_props_2str_ccpp
     ! Inputs
    integer,                           intent(in) :: ncol             ! Number of columns
    integer,                           intent(in) :: nlay             ! Number of vertical layers in radiation
    logical,                           intent(in) :: dosw             ! Flag for whether to perform longwave calculation
    class(ty_gas_optics_rrtmgp_ccpp),  intent(in) :: kdist_sw         ! Longwave gas optics object

    ! Outputs
    type(ty_optical_props_2str_ccpp),  intent(out) :: cloud_sw        ! Longwave cloud optics object
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errflg
    ! Set error variables
    errmsg = ''
    errflg = 0

    ! If not doing shortwave, no need to proceed
    if (.not. dosw) then
       return
    end if

    errmsg =cloud_sw%optical_props%alloc_2str(ncol, nlay, kdist_sw%gas_props)
    if (len_trim(errmsg) > 0) then
       errflg = 1
       return
    end if

  end subroutine rrtmgp_sw_cloud_temp_run

end module rrtmgp_sw_cloud_temp
