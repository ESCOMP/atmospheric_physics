! Diagnostics for ndrop_bam_ccpp scheme (BAM droplet activation)
!
! Requires ndrop_bam_ccpp_init to have run first (populates aername).
! This ordering is guaranteed by SDF placement: diagnostic schemes
! follow the physics schemes they diagnose.
module ndrop_bam_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: ndrop_bam_diagnostics_init
   public :: ndrop_bam_diagnostics_run

   integer :: num_aer = 0

contains

   !> \section arg_table_ndrop_bam_diagnostics_init  Argument Table
   !! \htmlinclude ndrop_bam_diagnostics_init.html
   subroutine ndrop_bam_diagnostics_init(errmsg, errflg)
      use cam_history, only: history_add_field
      use ndrop_bam,   only: aername

      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: l

      errmsg = ''
      errflg = 0

      ! CCN concentrations at fixed supersaturation levels
      call history_add_field('CCN1', 'cloud_condensation_nuclei_number_concentration_at_S_0.02pct', 'lev', 'avg', 'cm-3')
      call history_add_field('CCN2', 'cloud_condensation_nuclei_number_concentration_at_S_0.05pct', 'lev', 'avg', 'cm-3')
      call history_add_field('CCN3', 'cloud_condensation_nuclei_number_concentration_at_S_0.1pct',  'lev', 'avg', 'cm-3')
      call history_add_field('CCN4', 'cloud_condensation_nuclei_number_concentration_at_S_0.2pct',  'lev', 'avg', 'cm-3')
      call history_add_field('CCN5', 'cloud_condensation_nuclei_number_concentration_at_S_0.5pct',  'lev', 'avg', 'cm-3')
      call history_add_field('CCN6', 'cloud_condensation_nuclei_number_concentration_at_S_1.0pct',  'lev', 'avg', 'cm-3')

      ! Per-aerosol number concentration diagnostics (dynamic names from ndrop_bam)
      if (allocated(aername)) then
         num_aer = size(aername)
         do l = 1, num_aer
            call history_add_field(trim(aername(l))//'_m3', &
                 'aerosol_number_concentration_of_'//trim(aername(l)), 'lev', 'avg', 'm-3')
         end do
      end if

   end subroutine ndrop_bam_diagnostics_init

   !> \section arg_table_ndrop_bam_diagnostics_run  Argument Table
   !! \htmlinclude ndrop_bam_diagnostics_run.html
   subroutine ndrop_bam_diagnostics_run( &
      psat, naer_all, &
      ccn, naer2_diag, &
      errmsg, errflg)

      use cam_history, only: history_out_field
      use ndrop_bam,   only: aername, ccn_name

      integer,            intent(in)  :: psat
      integer,            intent(in)  :: naer_all
      real(kind_phys),    intent(in)  :: ccn(:,:,:)
      real(kind_phys),    intent(in)  :: naer2_diag(:,:,:)

      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: l

      errmsg = ''
      errflg = 0

      ! CCN concentrations (sliced from 3D array)
      do l = 1, psat
         call history_out_field(ccn_name(l), ccn(:,:,l))
      end do

      ! Per-aerosol number concentration diagnostics
      do l = 1, min(naer_all, num_aer)
         call history_out_field(trim(aername(l))//'_m3', naer2_diag(:,:,l))
      end do

   end subroutine ndrop_bam_diagnostics_run

end module ndrop_bam_diagnostics
