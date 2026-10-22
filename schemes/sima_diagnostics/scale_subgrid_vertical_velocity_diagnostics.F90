! Diagnostics for subgrid vertical velocity (WSUB, WSUBI)
module scale_subgrid_vertical_velocity_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: scale_subgrid_vertical_velocity_diagnostics_init
   public :: scale_subgrid_vertical_velocity_diagnostics_run

contains

   !> \section arg_table_scale_subgrid_vertical_velocity_diagnostics_init  Argument Table
   !! \htmlinclude scale_subgrid_vertical_velocity_diagnostics_init.html
   subroutine scale_subgrid_vertical_velocity_diagnostics_init(errmsg, errflg)
      use cam_history, only: history_add_field

      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call history_add_field('WSUB',  'subgrid_scale_vertical_velocity_for_droplet_nucleation', 'lev', 'avg', 'm s-1')
      call history_add_field('WSUBI', 'subgrid_scale_vertical_velocity_for_ice_nucleation',     'lev', 'avg', 'm s-1')

   end subroutine scale_subgrid_vertical_velocity_diagnostics_init

   !> \section arg_table_scale_subgrid_vertical_velocity_diagnostics_run  Argument Table
   !! \htmlinclude scale_subgrid_vertical_velocity_diagnostics_run.html
   subroutine scale_subgrid_vertical_velocity_diagnostics_run( &
      wsub, wsubi, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      real(kind_phys),    intent(in)  :: wsub(:,:)
      real(kind_phys),    intent(in)  :: wsubi(:,:)
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call history_out_field('WSUB',  wsub)
      call history_out_field('WSUBI', wsubi)

   end subroutine scale_subgrid_vertical_velocity_diagnostics_run

end module scale_subgrid_vertical_velocity_diagnostics
