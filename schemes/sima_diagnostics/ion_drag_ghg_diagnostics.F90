! Diagnostics for the ion_drag_ghg scheme (ion drag from a globally uniform tensor)
module ion_drag_ghg_diagnostics

   use ccpp_kinds, only:  kind_phys

   implicit none
   private

   public :: ion_drag_ghg_diagnostics_init ! init routine
   public :: ion_drag_ghg_diagnostics_run  ! main routine

CONTAINS

   !> \section arg_table_ion_drag_ghg_diagnostics_init  Argument Table
   !! \htmlinclude arg_table_ion_drag_ghg_diagnostics_init.html
   subroutine ion_drag_ghg_diagnostics_init(do_ion_drag, errmsg, errflg)
      use cam_history, only: history_add_field

      logical,          intent(in)  :: do_ion_drag   ! ion drag is applied [flag]
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! Fields are always present
      ! TODO: register_vector_field equivalent in CAM-SIMA?
      call history_add_field('UTEND_IONDRG', 'Zonal wind tendency by ion drag', 'lev', 'avg', 'm s-2')
      call history_add_field('VTEND_IONDRG', 'Meridional wind tendency by ion drag', 'lev', 'avg', 'm s-2')

      if (.not. do_ion_drag) return

      ! These fields exist only when ion drag is active
      call history_add_field('UIONTEND', 'u-tendency due to ion drag', 'lev', 'avg', 'm s-2')
      call history_add_field('VIONTEND', 'v-tendency due to ion drag', 'lev', 'avg', 'm s-2')

   end subroutine ion_drag_ghg_diagnostics_init

   !> \section arg_table_ion_drag_ghg_diagnostics_run  Argument Table
   !! \htmlinclude arg_table_ion_drag_ghg_diagnostics_run.html
   subroutine ion_drag_ghg_diagnostics_run(do_ion_drag, dudt, dvdt, errmsg, errflg)

      use cam_history, only: history_out_field

      logical,          intent(in)  :: do_ion_drag   ! ion drag is applied [flag]
      real(kind_phys),  intent(in)  :: dudt(:,:)     ! tendency of eastward wind due to ion drag [m s-2]
      real(kind_phys),  intent(in)  :: dvdt(:,:)     ! tendency of northward wind due to ion drag [m s-2]
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! Fields are always written; zeros when ion drag is not active
      ! Note: these fields are always added/written because they're part of the history_budget
      ! group in CAM and it could be annoying if SIMA aborted when these were part of a
      ! history configuration, but ion drag was not active. This mimics CAM behavior
      call history_out_field('UTEND_IONDRG', dudt)
      call history_out_field('VTEND_IONDRG', dvdt)

      if (.not. do_ion_drag) return

      ! These fields exist only when ion drag is active
      ! Field values are the same as above, but use WACCM naming convention.
      call history_out_field('UIONTEND', dudt)
      call history_out_field('VIONTEND', dvdt)

   end subroutine ion_drag_ghg_diagnostics_run

end module ion_drag_ghg_diagnostics
