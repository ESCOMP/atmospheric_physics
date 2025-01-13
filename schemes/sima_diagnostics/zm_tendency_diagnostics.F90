module zm_tendency_diagnostics
   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: zm_tendency_diagnostics_init ! init routine
   public :: zm_tendency_diagnostics_run  ! main routine

CONTAINS

 !> \section arg_table_zm_tendency_diagnostics_init  Argument Table
 !! \htmlinclude zm_tendency_diagnostics_init.html
 subroutine zm_tendency_diagnostics_init(errmsg, errflg)
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables:

    errmsg = ''
    errflg = 0

    call history_add_field ('ZMDICE', 'Cloud ice tendency - Zhang-McFarlane convection', 'lev',  'avg', 'kg kg-1 s-1')
    call history_add_field ('ZMDLIQ', 'Cloud liq tendency - Zhang-McFarlane convection', 'lev',  'avg', 'kg kg-1 s-1')

   end subroutine zm_tendency_diagnostics_init

   !> \section arg_table_zm_tendency_diagnostics_run  Argument Table
   !! \htmlinclude zm_tendency_diagnostics_run.html
   subroutine zm_tendency_diagnostics_run(ncol, pver, pverp, const_props, dqdt, errmsg, errflg)

      use cam_history,               only: history_out_field
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: pver
      integer, intent(in) :: pverp
      type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
      real(kind_phys), intent(in) :: dqdt(:,:,:)  ! Tracer tendency array  (ncol,pver,ncnst)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: lengath  ! number of columns with deep convection
      integer :: const_idx
      character(len=256) :: standard_name

      real(kind_phys) :: ftem(ncol,pver)

      errmsg = ''
      errflg = 0

      do const_idx = 1, size(const_props)
         call const_props(const_idx)%standard_name(standard_name, errflg, errmsg)

         ! Use the regular constituent names, as the location in the dqdt array match the constituent ordering in the
         ! constituent properties

         select case (trim(standard_name))

         case('cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water')
            call history_out_field('ZMDLIQ', dqdt(:,:,const_idx))

         case('cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water')
            call history_out_field('ZMDICE', dqdt(:,:,const_idx))

         end select

     end do

   end subroutine zm_tendency_diagnostics_run

   !=======================================================================

end module zm_tendency_diagnostics
