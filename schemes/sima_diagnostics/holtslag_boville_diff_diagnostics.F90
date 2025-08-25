! Diagnostics for cloud fraction
module holtslag_boville_diff_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: holtslag_boville_diff_diagnostics_init
   public :: holtslag_boville_diff_diagnostics_run

contains

   !> \section arg_table_holtslag_boville_diff_diagnostics_init  Argument Table
   !! \htmlinclude holtslag_boville_diff_diagnostics_init.html
   subroutine holtslag_boville_diff_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('HB_ri', 'Richardson Number (HB scheme)', 'lev', 'inst', '1')

      call history_add_field('USTAR', 'Surface friction velocity', horiz_only, 'inst', 'm s-1')
      call history_add_field('obklen', 'Obukhov length', horiz_only, 'avg', 'm')
      call history_add_field('PBLH', 'Planetary boundary layer height', horiz_only, 'avg', 'm')

      call history_add_field('KVH', 'Vertical diffusion diffusivity coefficient for heat/moisture', 'ilev', 'avg', 'm2 s-1')
      call history_add_field('KVM', 'Vertical diffusion diffusivity coefficient for momentum', 'ilev', 'avg', 'm2 s-1')
      call history_add_field('CGS', 'Counter-gradient coefficient on surface kinematic fluxes', 'ilev', 'avg', 's m-2')

      call history_add_field('TKE', 'Specific turbulent kinetic energy at interfaces', 'ilev', 'avg', 'm2 s-2')
      call history_add_field('TPERT', 'Perturbation temperature (eddies in PBL)', horiz_only, 'avg', 'K')
      call history_add_field('QPERT', 'Perturbation specific humidity (eddies in PBL)', horiz_only, 'avg', 'kg kg-1')

   end subroutine holtslag_boville_diff_diagnostics_init

   !> \section arg_table_holtslag_boville_diff_diagnostics_run  Argument Table
   !! \htmlinclude holtslag_boville_diff_diagnostics_run.html
   subroutine holtslag_boville_diff_diagnostics_run( &
      ri, &
      ustar, obklen, &
      pblh, &
      kvh, kvm, cgs, &
      tke, &
      tpert, qpert, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input parameters
      real(kind_phys),    intent(in)  :: ri(:,:)
      real(kind_phys),    intent(in)  :: ustar(:)
      real(kind_phys),    intent(in)  :: obklen(:)
      real(kind_phys),    intent(in)  :: pblh(:)

      real(kind_phys),    intent(in)  :: kvh(:,:)
      real(kind_phys),    intent(in)  :: kvm(:,:)
      real(kind_phys),    intent(in)  :: cgs(:,:)
      real(kind_phys),    intent(in)  :: tke(:,:)

      real(kind_phys),    intent(in)  :: tpert(:)
      real(kind_phys),    intent(in)  :: qpert(:)

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call history_out_field('HB_ri', ri)
      call history_out_field('USTAR', ustar)
      call history_out_field('obklen', obklen)
      call history_out_field('PBLH', pblh)

      call history_out_field('KVH', kvh)
      call history_out_field('KVM', kvm)
      call history_out_field('CGS', cgs)
      call history_out_field('TKE', tke)

      call history_out_field('TPERT', tpert)
      call history_out_field('QPERT', qpert)

   end subroutine holtslag_boville_diff_diagnostics_run

end module holtslag_boville_diff_diagnostics
