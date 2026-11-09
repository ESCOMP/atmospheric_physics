! Register-only stub that provides the O3 constituent for snapshot test suites (e.g. suite_rrtmgp)
module ozone_constituent_stub

implicit none
private

public :: ozone_constituent_stub_register
public :: ozone_constituent_stub_run

contains

!> \section arg_table_ozone_constituent_stub_register  Argument Table
!! \htmlinclude ozone_constituent_stub_register.html
  subroutine ozone_constituent_stub_register(constituents, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use ccpp_chem_utils,           only: chem_constituent_qmin

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituents(:)
    character(len=*), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

    allocate(constituents(1), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) then
      errmsg = 'ozone_constituent_stub_register: ' // trim(errmsg)
      return
    end if

    call constituents(1)%instantiate( &
         std_name          = 'O3', &
         diag_name         = 'O3', &
         long_name         = 'prescribed ozone (O3)', &
         units             = 'kg kg-1', &
         vertical_dim      = 'vertical_layer_dimension', &
         min_value         = chem_constituent_qmin('O3'), &
         advected          = .false., &
         water_species     = .false., &
         mixing_ratio_type = 'dry', &
         errcode           = errcode, &
         errmsg            = errmsg)

  end subroutine ozone_constituent_stub_register

!> \section arg_table_ozone_constituent_stub_run  Argument Table
!! \htmlinclude ozone_constituent_stub_run.html
  subroutine ozone_constituent_stub_run(errmsg, errcode)
    character(len=*), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errcode = 0
    errmsg = ''
  end subroutine ozone_constituent_stub_run

end module ozone_constituent_stub
