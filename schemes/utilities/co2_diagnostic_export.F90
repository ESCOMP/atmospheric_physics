! Export the diagnostic CO2 concentration to the surface coupler.
module co2_diagnostic_export

  implicit none
  private

  ! Public interfaces
  public :: co2_diagnostic_export_run
contains
!> \section arg_table_co2_diagnostic_export_run Argument Table
!! \htmlinclude co2_diagnostic_export_run.html
subroutine co2_diagnostic_export_run(co2_vmr, co2diag, errmsg, errcode)
   use ccpp_kinds,       only: kind_phys

   ! Input arguments
   real(kind_phys),  intent(in)  :: co2_vmr      ! Prescribed CO2 volume mixing ratio [mol mol-1]

   ! Output arguments
   real(kind_phys),  intent(out) :: co2diag(:)   ! Diagnostic CO2 concentration for coupler [ppmv]
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errcode

   errmsg = ''
   errcode = 0

   ! mol/mol -> ppmv expected by the coupler. (chem_surfvals sets CO2VMR * 1e6)
   co2diag(:) = co2_vmr * 1.0e6_kind_phys

end subroutine co2_diagnostic_export_run

end module co2_diagnostic_export
