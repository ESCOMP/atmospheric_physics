! CCPP interstitial schemes needed for the PUMAS cloud microphysics package.

!The hope is that there will eventually be one set of portable interstitials and possibly another set
!of host-specific interstitials (if needed) above the portable layer.

module pumas_pre_main

   implicit none

   private

   public :: pumas_pre_main_init

   contains

  !> \section arg_table_pumas_pre_main_init Argument Table
  !! \htmlinclude pumas_pre_main_init.html
   subroutine pumas_pre_main_init(do_clubb_sgs, remove_supersat, spat_vary_accre_enhan_in, &
                                  micro_mg_warm_rain, pumas_stochastic_tau_kernel_filename, &
                                  errmsg, errcode)


     use ccpp_kinds,        only: kind_phys
     use pumas_stochastic_collect_tau, only: pumas_stochastic_kernel_init

     logical,            intent(in)  :: do_clubb_sgs
     real(kind_phys), dimension (:,:) , intent(out) :: spat_vary_accre_enhan_in
     logical,            intent(out) :: remove_supersat
     character(len=*),   intent(in)  :: micro_mg_warm_rain                   ! warm rain method (kk2000, sb2001, tau, emulated)
     character(len=*),   intent(in)  :: pumas_stochastic_tau_kernel_filename ! TAU stochastic collection kernel file
     character(len=*),   intent(out) :: errmsg
     integer,            intent(out) :: errcode

     logical :: file_exists

     errmsg = ' '
     errcode = 0

     ! Inside CAM, the pbuf accre_enhan variable is hardwired to 1.
     spat_vary_accre_enhan_in(:,:) = 1._kind_phys

     ! set remove_supersat dependent on whether CLUBB is being run or not
     if (do_clubb_sgs) then
       remove_supersat = .false.
     else
       remove_supersat = .true.
     end if

     ! Read the collection kernel for the TAU stochastic collection code (micro_mg_warm_rain = 'tau')
     ! (this is done by stochastic_tau_init_cam in CAM's micro_pumas_cam_init)
     if (trim(micro_mg_warm_rain) == 'tau') then
       inquire(file=trim(pumas_stochastic_tau_kernel_filename), exist=file_exists)
       if (.not. file_exists) then
         errcode = 1
         errmsg = 'pumas_pre_main_init: TAU stochastic collection kernel file not found: '// &
                  trim(pumas_stochastic_tau_kernel_filename)
         return
       end if
       call pumas_stochastic_kernel_init(trim(pumas_stochastic_tau_kernel_filename))
     end if

   end subroutine pumas_pre_main_init

end module pumas_pre_main

