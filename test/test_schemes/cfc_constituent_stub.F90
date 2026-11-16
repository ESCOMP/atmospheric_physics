! Register-only stub that provides the CFC11STAR constituent for snapshot test suites (e.g. suite_rrtmgp).
! CFC11STAR is the CFC11-equivalent that CAM chemistry configurations give radiation as CFC11,
! so the snapshot's radiation inputs are reproduced with rad_climate entry 'CFC11:CFC11STAR'.
module cfc_constituent_stub

implicit none
private

public :: cfc_constituent_stub_register
public :: cfc_constituent_stub_run

contains

!> \section arg_table_cfc_constituent_stub_register  Argument Table
!! \htmlinclude cfc_constituent_stub_register.html
  subroutine cfc_constituent_stub_register(constituents, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use ccpp_chem_utils,           only: chem_constituent_qmin

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituents(:)
    character(len=*), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

    allocate(constituents(1), stat=errcode, errmsg=errmsg)
    if (errcode /= 0) then
      errmsg = 'cfc_constituent_stub_register: ' // trim(errmsg)
      return
    end if

    call constituents(1)%instantiate( &
         std_name          = 'CFC11STAR', &
         diag_name         = 'CFC11STAR', &
         long_name         = 'CFC11-equivalent of the radiatively active halocarbons (CFC11STAR)', &
         units             = 'kg kg-1', &
         vertical_dim      = 'vertical_layer_dimension', &
         min_value         = chem_constituent_qmin('CFC11STAR'), &
         advected          = .false., &
         water_species     = .false., &
         mixing_ratio_type = 'dry', &
         errcode           = errcode, &
         errmsg            = errmsg)

  end subroutine cfc_constituent_stub_register

!> \section arg_table_cfc_constituent_stub_run  Argument Table
!! \htmlinclude cfc_constituent_stub_run.html
  subroutine cfc_constituent_stub_run(errmsg, errcode)
    character(len=*), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errcode = 0
    errmsg = ''
  end subroutine cfc_constituent_stub_run

end module cfc_constituent_stub
