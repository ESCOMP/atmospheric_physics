module zm_diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! THIS IS A TEMPLATE
!   1. copy this file to a new file with the correct name
!        (SCHEME_diagnostics.F90)
!   2. do a search and replace for "SCHEME" in this file and
!        replace with your scheme name
!   3. Add desired history_add_field calls to the init phase
!   4. Add all fields that are being output as inputs to the run phase
!   5. Add desired history_out_field calls to the run phase
!   6. Run $ccpp_framework/scripts/ccpp_fortran_to_metadata.py on this .F90
!        file to generate the metadata
!   7. Complete the metadata (fill out standard names, units, dimensions)
!   8. Add this scheme to the SDF file for your suite (likely will be at end)
!   9. Delete this header section
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: zm_diagnostics_init ! init routine
   public :: zm_diagnostics_run  ! main routine

CONTAINS

   !> \section arg_table_zm_diagnostics_init  Argument Table
   !! \htmlinclude zm_diagnostics_init.html
   subroutine zm_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      call history_add_field ('PRECZ', 'total precipitation from ZM convection', horiz_only, 'avg', 'm s-1')
      call history_add_field ('ZMDT', 'T tendency - Zhang-McFarlane moist convection', 'lev', 'avg', 'K s-1')

   end subroutine zm_diagnostics_init

   !> \section arg_table_zm_diagnostics_run  Argument Table
   !! \htmlinclude zm_diagnostics_run.html
   subroutine zm_diagnostics_run(ncol, pver, cpair, heat, prec, errmsg, errflg)

      use cam_history, only: history_out_field
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: pver

      real(kind_phys), intent(in) :: cpair
      real(kind_phys), intent(in) :: heat(:,:)
      real(kind_phys), intent(in) :: prec(:)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      real(kind_phys) :: ftem(ncol,pver)

      errmsg = ''
      errflg = 0

      call history_out_field('PRECZ', prec)

      ftem(:,:) = 0._kind_phys
      ftem(:ncol,:pver) = heat(:ncol,:pver)/cpair
      call history_out_field('ZMDT', ftem)

   end subroutine zm_diagnostics_run

   !=======================================================================

end module zm_diagnostics
