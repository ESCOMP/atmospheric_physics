! Diagnostics for cloud_water_paths: the CAM cloud water path history fields.
! - The in-cloud mixing ratios, water contents and water paths come from
!   CAM cloud_diagnostics_calc (two-moment branch). The grid-box and
!   vertically integrated paths are derived here from the in-cloud paths
!   exactly as CAM does, since they are diagnostic-only quantities.
! - The snow/graupel adjusted cloud fractions come from CAM micro_pumas_cam.
!   CAM does not output the in-cloud snow/graupel water paths (icswp,
!   icgrauwp), so they are not output here either.
module cloud_water_paths_diagnostics
   implicit none
   private

   public :: cloud_water_paths_diagnostics_init
   public :: cloud_water_paths_diagnostics_run

contains

   !> \section arg_table_cloud_water_paths_diagnostics_init  Argument Table
   !! \htmlinclude cloud_water_paths_diagnostics_init.html
   subroutine cloud_water_paths_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History add field calls (CAM cloud_diagnostics field names)
      call history_add_field('ICWMR',    'in_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')
      call history_add_field('ICIMR',    'in_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')
      call history_add_field('IWC',      'gridbox_mean_ice_water_content', 'lev', 'avg', 'kg m-3')
      call history_add_field('LWC',      'gridbox_mean_liquid_water_content', 'lev', 'avg', 'kg m-3')
      call history_add_field('ICLDIWP',  'in_cloud_ice_water_path', 'lev', 'avg', 'kg m-2')
      call history_add_field('ICLDTWP',  'in_cloud_total_water_path', 'lev', 'avg', 'kg m-2')
      call history_add_field('GCLDLWP',  'gridbox_total_water_path', 'lev', 'avg', 'kg m-2')
      call history_add_field('TGCLDIWP', 'vertically_integrated_gridbox_ice_water_path', horiz_only, 'avg', 'kg m-2')
      call history_add_field('TGCLDLWP', 'vertically_integrated_gridbox_liquid_water_path', horiz_only, 'avg', 'kg m-2')
      call history_add_field('TGCLDCWP', 'vertically_integrated_gridbox_total_water_path', horiz_only, 'avg', 'kg m-2')

      ! History add field calls (CAM micro_pumas_cam field names)
      call history_add_field('CLDFSNOW', 'liquid_plus_snow_stratiform_cloud_area_fraction', 'lev', 'avg', '1')
      call history_add_field('CLDFGRAU', 'liquid_plus_graupel_stratiform_cloud_area_fraction', 'lev', 'avg', '1')

   end subroutine cloud_water_paths_diagnostics_init

   !> \section arg_table_cloud_water_paths_diagnostics_run  Argument Table
   !! \htmlinclude cloud_water_paths_diagnostics_run.html
   subroutine cloud_water_paths_diagnostics_run( &
      ncol, pver, &
      cld, iclwp, iciwp, &
      icimr, icwmr, iwc, lwc, &
      cldfsnow, cldfgrau, &
      errmsg, errflg)

      use ccpp_kinds, only: kind_phys
      use cam_history, only: history_out_field

      ! Input parameters
      integer,            intent(in)  :: ncol
      integer,            intent(in)  :: pver
      real(kind_phys),    intent(in)  :: cld(:,:)        ! Total cloud fraction [fraction]
      real(kind_phys),    intent(in)  :: iclwp(:,:)      ! In-cloud cloud liquid water path [kg m-2]
      real(kind_phys),    intent(in)  :: iciwp(:,:)      ! In-cloud cloud ice water path [kg m-2]
      real(kind_phys),    intent(in)  :: icimr(:,:)      ! In-cloud ice mixing ratio [kg kg-1]
      real(kind_phys),    intent(in)  :: icwmr(:,:)      ! In-cloud water mixing ratio [kg kg-1]
      real(kind_phys),    intent(in)  :: iwc(:,:)        ! Grid box average ice water content [kg m-3]
      real(kind_phys),    intent(in)  :: lwc(:,:)        ! Grid box average liquid water content [kg m-3]
      real(kind_phys),    intent(in)  :: cldfsnow(:,:)   ! Cloud fraction adjusted for snow [1]
      real(kind_phys),    intent(in)  :: cldfgrau(:,:)   ! Cloud fraction adjusted for graupel [1]

      ! CCPP error handling variables
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      real(kind_phys) :: gicewp(ncol,pver)  ! Grid-box cloud ice water path [kg m-2]
      real(kind_phys) :: gliqwp(ncol,pver)  ! Grid-box cloud liquid water path [kg m-2]
      real(kind_phys) :: gwp(ncol,pver)     ! Grid-box cloud (total) water path [kg m-2]
      real(kind_phys) :: cwp(ncol,pver)     ! In-cloud cloud (total) water path [kg m-2]
      real(kind_phys) :: tgicewp(ncol)      ! Vertically integrated ice water path [kg m-2]
      real(kind_phys) :: tgliqwp(ncol)      ! Vertically integrated liquid water path [kg m-2]
      real(kind_phys) :: tgwp(ncol)         ! Vertically integrated (total) cloud water path [kg m-2]

      integer :: k

      errmsg = ''
      errflg = 0

      ! In the two-moment branch the in-cloud paths are the radiation inputs
      ! themselves, so the in-cloud ice/liquid paths (CAM cicewp/cliqwp) are
      ! iciwp/iclwp and the grid-box paths are those weighted by cloud fraction.
      gicewp(:ncol,:pver) = iciwp(:ncol,:pver)*cld(:ncol,:pver)
      gliqwp(:ncol,:pver) = iclwp(:ncol,:pver)*cld(:ncol,:pver)

      tgicewp(:ncol) = 0._kind_phys
      tgliqwp(:ncol) = 0._kind_phys

      do k = 1, pver
         tgicewp(:ncol) = tgicewp(:ncol) + gicewp(:ncol,k)
         tgliqwp(:ncol) = tgliqwp(:ncol) + gliqwp(:ncol,k)
      end do

      tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
      gwp(:ncol,:pver) = gicewp(:ncol,:pver) + gliqwp(:ncol,:pver)
      cwp(:ncol,:pver) = iciwp(:ncol,:pver) + iclwp(:ncol,:pver)

      ! History out field calls
      call history_out_field('IWC',      iwc)
      call history_out_field('LWC',      lwc)
      call history_out_field('ICIMR',    icimr)
      call history_out_field('ICWMR',    icwmr)
      call history_out_field('GCLDLWP',  gwp)
      call history_out_field('TGCLDCWP', tgwp)
      call history_out_field('TGCLDLWP', tgliqwp)
      call history_out_field('TGCLDIWP', tgicewp)
      call history_out_field('ICLDTWP',  cwp)
      call history_out_field('ICLDIWP',  iciwp)
      call history_out_field('CLDFSNOW', cldfsnow)
      call history_out_field('CLDFGRAU', cldfgrau)

   end subroutine cloud_water_paths_diagnostics_run

end module cloud_water_paths_diagnostics
