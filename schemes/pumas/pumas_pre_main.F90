!This is just a placeholder for the eventual CCPP
!interstitial schemes needed for the PUMAS cloud
!microphysics package.

!The hope that there will eventually be one set of
!portable interstitials and possibly another set
!of host-specific intersititals (if needed) above
!the portable layer.

module pumas_pre_main

   use ccpp_kinds,        only: kind_phys

   implicit none

   contains

   !Add subroutines here
   !for any pre-processing
   !steps needed before the core
   !PUMAS calls.


  !> \section arg_table_pumas_pre_main_init Argument Table
  !! \htmlinclude pumas_pre_main_init.html
   subroutine pumas_pre_main_init(spat_vary_accre_enhan_in, errmsg, errcode)


     real(kind_phys), dimension (:,:) , intent(out) :: spat_vary_accre_enhan_in
     character(len=512), intent(out) :: errmsg
     integer,            intent(out) :: errcode

     errmsg = ' '
     errcode = 0

     ! Inside CAM, the pbuf accre_enhan variable is hardwired to 1.
     spat_vary_accre_enhan_in(:,:) = 1._kind_phys

   end subroutine pumas_pre_main_init

end module pumas_pre_main
