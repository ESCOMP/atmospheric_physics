! Copyright (C) 2025 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module set_general_conv_fluxes_to_shallow

   use ccpp_kinds, only:  kind_phys

   implicit none

contains
!> \section arg_table_set_general_conv_fluxes_to_shallow_run Argument Table
!! \htmlinclude set_general_conv_fluxes_to_shallow_run.html
   subroutine set_general_conv_fluxes_to_shallow_run(tend_s_snwprd_gen, tend_s_snwprd_sh, tend_s_snwevmlt_gen, tend_s_snwevmlt_sh, &
          prec_gen, prec_sh, snow_gen, snow_sh, ntprprd_gen, ntprprd_sh, ntsnprd_gen, ntsnprd_sh, &
          flxprec_gen, flxprec_sh, flxsnow_gen, flxsnow_sh)

   real(kind_phys), intent(in) :: tend_s_snwprd_gen(:,:)
   real(kind_phys), intent(out) :: tend_s_snwprd_sh(:,:)
   real(kind_phys), intent(in) :: tend_s_snwevmlt_gen(:,:)
   real(kind_phys), intent(out) :: tend_s_snwevmlt_sh(:,:)
   real(kind_phys), intent(in) :: prec_gen(:) ! precipitation rate is modified by zm_conv_evap based on computed flux
   real(kind_phys), intent(out) :: prec_sh(:) ! and thus needs to be renamed from generic to shallow.
   real(kind_phys), intent(in) :: snow_gen(:)
   real(kind_phys), intent(out) :: snow_sh(:) ! unused by diagnostics
   real(kind_phys), intent(in) :: ntprprd_gen(:,:)
   real(kind_phys), intent(out) :: ntprprd_sh(:,:)
   real(kind_phys), intent(in) :: ntsnprd_gen(:,:)
   real(kind_phys), intent(out) :: ntsnprd_sh(:,:)
   real(kind_phys), intent(in) :: flxprec_gen(:,:)
   real(kind_phys), intent(out) :: flxprec_sh(:,:)
   real(kind_phys), intent(in) :: flxsnow_gen(:,:)
   real(kind_phys), intent(out) :: flxsnow_sh(:,:)

   tend_s_snwprd_sh = tend_s_snwprd_gen
   tend_s_snwevmlt_sh = tend_s_snwevmlt_gen
   prec_sh = prec_gen
   snow_sh = snow_gen
   ntprprd_sh = ntprprd_gen
   ntsnprd_sh = ntsnprd_gen
   flxprec_sh = flxprec_gen
   flxsnow_sh = flxsnow_gen

   end subroutine set_general_conv_fluxes_to_shallow_run

end module set_general_conv_fluxes_to_shallow
