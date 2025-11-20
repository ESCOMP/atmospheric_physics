module rrtmgp_cloud_diagnostics
   implicit none
   private
   save

   public :: rrtmgp_cloud_diagnostics_init ! init routine
   public :: rrtmgp_cloud_diagnostics_run  ! main routine

CONTAINS

   !> \section arg_table_rrtmgp_cloud_diagnostics_init  Argument Table
   !! \htmlinclude rrtmgp_cloud_diagnostics_init.html
   subroutine rrtmgp_cloud_diagnostics_init(has_snow, has_graupel, graupel_in_rad, errmsg, errflg)
      use cam_history,         only: history_add_field

      logical,             intent(in) :: has_snow
      logical,             intent(in) :: has_graupel
      logical,             intent(in) :: graupel_in_rad
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:
      integer :: icall

      errmsg = ''
      errflg = 0

      ! Add diagnostic fields
      call history_add_field('TOT_CLD_VISTAU',  'Total gbx cloud extinction visible sw optical depth', 'lev', 'avg', '1', &
         flag_xyfill=.true.)
      call history_add_field('TOT_ICLD_VISTAU', 'Total in-cloud extinction visible sw optical depth',  'lev', 'avg', '1', &
         flag_xyfill=.true.)
      call history_add_field('LIQ_ICLD_VISTAU', 'Liquid in-cloud extinction visible sw optical depth', 'lev', 'avg', '1', &
         flag_xyfill=.true.)
      call history_add_field('ICE_ICLD_VISTAU', 'Ice in-cloud extinction visible sw optical depth',    'lev', 'avg', '1', &
         flag_xyfill=.true.)

      if (has_snow) then
         call history_add_field('SNOW_ICLD_VISTAU', 'Snow in-cloud extinction visible sw optical depth', 'lev', 'avg', '1', &
                flag_xyfill=.true.)
      end if
      if (has_graupel .and. graupel_in_rad) then
         call history_add_field('GRAU_ICLD_VISTAU', 'Graupel in-cloud extinction visible sw optical depth', 'lev', 'avg', '1', &
                flag_xyfill=.true.)
      end if

   end subroutine rrtmgp_cloud_diagnostics_init

   !> \section arg_table_rrtmgp_cloud_diagnostics_run  Argument Table
   !! \htmlinclude rrtmgp_cloud_diagnostics_run.html
   subroutine rrtmgp_cloud_diagnostics_run(write_output, has_snow, has_graupel, graupel_in_rad, tot_cld_vistau, tot_icld_vistau, &
                  liq_icld_vistau, ice_icld_vistau, snow_icld_vistau, grau_icld_vistau, errmsg, errflg)
      use cam_history,        only: history_out_field
      use ccpp_kinds,         only: kind_phys
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      logical,            intent(in) :: write_output
      logical,            intent(in) :: has_snow
      logical,            intent(in) :: has_graupel
      logical,            intent(in) :: graupel_in_rad
      real(kind_phys),    intent(in) :: tot_cld_vistau(:,:)
      real(kind_phys),    intent(in) :: tot_icld_vistau(:,:)
      real(kind_phys),    intent(in) :: liq_icld_vistau(:,:)
      real(kind_phys),    intent(in) :: ice_icld_vistau(:,:)
      real(kind_phys),    intent(in) :: snow_icld_vistau(:,:)
      real(kind_phys),    intent(in) :: grau_icld_vistau(:,:)
      
      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables

      errmsg = ''
      errflg = 0

      if (.not. write_output) then
         return
      end if

      if (has_snow) then
         call history_out_field('SNOW_ICLD_VISTAU', snow_icld_vistau)
      end if
      if (has_graupel .and. graupel_in_rad) then
         call history_out_field('GRAU_ICLD_VISTAU', grau_icld_vistau)
      end if

      call history_out_field('TOT_CLD_VISTAU', tot_cld_vistau)
      call history_out_field('TOT_ICLD_VISTAU', tot_icld_vistau)
      call history_out_field('LIQ_ICLD_VISTAU', liq_icld_vistau)
      call history_out_field('ICE_ICLD_VISTAU', ice_icld_vistau)
      call history_out_field('SNOW_ICLD_VISTAU', snow_icld_vistau)
      call history_out_field('GRAU_ICLD_VISTAU', grau_icld_vistau)

   end subroutine rrtmgp_cloud_diagnostics_run

end module rrtmgp_cloud_diagnostics
