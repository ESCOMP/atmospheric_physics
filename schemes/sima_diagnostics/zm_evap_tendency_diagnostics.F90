module zm_evap_tendency_diagnostics
   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: zm_evap_tendency_diagnostics_init ! init routine
   public :: zm_evap_tendency_diagnostics_run  ! main routine

CONTAINS

 !> \section arg_table_zm_evap_tendency_diagnostics_init  Argument Table
 !! \htmlinclude zm_evap_tendency_diagnostics_init.html
 subroutine zm_evap_tendency_diagnostics_init(errmsg, errflg)
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables:

    errmsg = ''
    errflg = 0

    call history_add_field ('ZMEIHEAT', 'Heating by precipitation freezing/melting and evaporation in ZM convection', 'lev', 'avg', 'W kg-1')

    call history_add_field ('EVAPTZM', 'T tendency - Evaporation/snow prod from ZM convection', 'lev',  'avg', 'K s-1')
    call history_add_field ('FZSNTZM', 'T tendency - Rain to snow conversion from ZM convection', 'lev',  'avg', 'K s-1')
    call history_add_field ('EVSNTZM', 'T tendency - Snow to rain prod from ZM convection', 'lev',  'avg', 'K s-1')
    call history_add_field ('EVAPQZM', 'Q tendency - Evaporation from Zhang-McFarlane moist convection', 'lev',  'avg', &
                                       'kg kg-1 s-1')

   end subroutine zm_evap_tendency_diagnostics_init

   !> \section arg_table_zm_evap_tendency_diagnostics_run  Argument Table
   !! \htmlinclude zm_evap_tendency_diagnostics_run.html
   subroutine zm_evap_tendency_diagnostics_run(ncol, pver, cpair, tend_s_snwprd, tend_s_snwevmlt, qtnd, heat, errmsg, errflg)

      use cam_history,               only: history_out_field
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: pver
      real(kind_phys), intent(in) :: cpair
      real(kind_phys), intent(in) :: tend_s_snwprd(:,:)
      real(kind_phys), intent(in) :: tend_s_snwevmlt(:,:)
      real(kind_phys), intent(in) :: qtnd(:,:)
      real(kind_phys), intent(in) :: heat(:,:)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: lengath  ! number of columns with deep convection
      integer :: const_idx
      character(len=256) :: standard_name

      real(kind_phys) :: ftem(ncol,pver)

      errmsg = ''

      errflg = 0

   ftem(:ncol,:pver) = heat(:ncol,:pver)/cpair
   call history_out_field('EVAPTZM ',ftem)
   ftem(:ncol,:pver) = tend_s_snwprd  (:ncol,:pver)/cpair
   call history_out_field('FZSNTZM ',ftem)
   ftem(:ncol,:pver) = tend_s_snwevmlt(:ncol,:pver)/cpair
   call history_out_field('EVSNTZM ',ftem)
   call history_out_field('EVAPQZM ',qtnd)
   call history_out_field('ZMEIHEAT', heat)

   end subroutine zm_evap_tendency_diagnostics_run

   !=======================================================================

end module zm_evap_tendency_diagnostics
