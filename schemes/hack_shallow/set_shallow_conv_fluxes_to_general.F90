module set_shallow_conv_fluxes_to_general

   use ccpp_kinds, only:  kind_phys

   implicit none

contains

!> \section arg_table_set_shallow_conv_fluxes_to_general_run Argument Table
!! \htmlinclude set_shallow_conv_fluxes_to_general_run.html
   subroutine set_shallow_conv_fluxes_to_general_run(prec_gen, prec_sh, prdprec_gen, prdprec_sh)

   real(kind_phys), intent(out) :: prec_gen(:)
   real(kind_phys), intent(in) :: prec_sh(:)
   real(kind_phys), intent(out) :: prdprec_gen(:,:)
   real(kind_phys), intent(in) :: prdprec_sh(:,:)

   prec_gen = prec_sh
   prdprec_gen = prdprec_sh

   end subroutine set_shallow_conv_fluxes_to_general_run

end module set_shallow_conv_fluxes_to_general
