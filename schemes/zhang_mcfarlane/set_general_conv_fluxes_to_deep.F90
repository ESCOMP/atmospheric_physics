module set_general_conv_fluxes_to_deep

   use ccpp_kinds, only:  kind_phys

   implicit none

contains

!===============================================================================
!> \section arg_table_set_general_conv_fluxes_to_deep_run Argument Table
!! \htmlinclude set_general_conv_fluxes_to_deep_run.html
!!

   subroutine set_general_conv_fluxes_to_deep_run(tend_s_snwprd_gen, tend_s_snwprd_dp, tend_s_snwevmlt_gen, tend_s_snwevmlt_dp, &
          prdprec_gen, prdprec_dp, prec_gen, prec_dp, snow_gen, snow_dp, ntprprd_gen, ntprprd_dp, ntsnprd_gen, ntsnprd_dp, &
          flxprec_gen, flxprec_dp, flxsnow_gen, flxsnow_dp)

   real(kind_phys), intent(in) :: tend_s_snwprd_gen(:,:)
   real(kind_phys), intent(out) :: tend_s_snwprd_dp(:,:)
   real(kind_phys), intent(in) :: tend_s_snwevmlt_gen(:,:)
   real(kind_phys), intent(out) :: tend_s_snwevmlt_dp(:,:)
   real(kind_phys), intent(in) :: prdprec_gen(:,:)
   real(kind_phys), intent(out) :: prdprec_dp(:,:)
   real(kind_phys), intent(in) :: prec_gen(:)
   real(kind_phys), intent(out) :: prec_dp(:)
   real(kind_phys), intent(in) :: snow_gen(:)
   real(kind_phys), intent(out) :: snow_dp(:)
   real(kind_phys), intent(in) :: ntprprd_gen(:,:)
   real(kind_phys), intent(out) :: ntprprd_dp(:,:)
   real(kind_phys), intent(in) :: ntsnprd_gen(:,:)
   real(kind_phys), intent(out) :: ntsnprd_dp(:,:)
   real(kind_phys), intent(in) :: flxprec_gen(:,:)
   real(kind_phys), intent(out) :: flxprec_dp(:,:)
   real(kind_phys), intent(in) :: flxsnow_gen(:,:)
   real(kind_phys), intent(out) :: flxsnow_dp(:,:)

   tend_s_snwprd_dp = tend_s_snwprd_gen
   tend_s_snwevmlt_dp = tend_s_snwevmlt_gen
   prdprec_dp = prdprec_gen
   prec_dp = prec_gen
   snow_dp = snow_gen
   ntprprd_dp = ntprprd_gen
   ntsnprd_dp = ntsnprd_gen
   flxprec_dp = flxprec_gen
   flxsnow_dp = flxsnow_gen

   end subroutine set_general_conv_fluxes_to_deep_run

end module set_general_conv_fluxes_to_deep
