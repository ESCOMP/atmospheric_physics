module set_deep_conv_fluxes_to_general

   use ccpp_kinds, only:  kind_phys

   implicit none

contains

!===============================================================================
!> \section arg_table_set_deep_conv_fluxes_to_general_run Argument Table
!! \htmlinclude set_deep_conv_fluxes_to_general_run.html
!!

   subroutine set_deep_conv_fluxes_to_general_run(prec_gen, prec_dp, prdprec_gen, prdprec_dp)

   real(kind_phys), intent(out) :: prec_gen(:)
   real(kind_phys), intent(in) :: prec_dp(:)
   real(kind_phys), intent(out) :: prdprec_gen(:,:)
   real(kind_phys), intent(in) :: prdprec_dp(:,:)

   prec_gen = prec_dp
   prdprec_gen = prdprec_dp

   end subroutine set_deep_conv_fluxes_to_general_run

end module set_deep_conv_fluxes_to_general
