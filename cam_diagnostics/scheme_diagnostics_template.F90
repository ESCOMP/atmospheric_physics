module SCHEME_diagnostics
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

   public :: SCHEME_diagnostics_init ! init routine
   public :: SCHEME_diagnostics_run  ! main routine

CONTAINS

   !> \section arg_table_SCHEME_diagnostics_init  Argument Table
   !! \htmlinclude SCHEME_diagnostics_init.html
   subroutine SCHEME_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      ! History add field calls
      ! Example:
      ! call history_add_field('TESTDIAG', 'not_a_real_diagnostic_field', horiz_only, 'avg', 'gremlin s-1')

   end subroutine SCHEME_diagnostics_init

   !> \section arg_table_SCHEME_diagnostics_run  Argument Table
   !! \htmlinclude SCHEME_diagnostics_run.html
   subroutine SCHEME_diagnostics_run(list, of, fields, errmsg, errflg)

      use cam_history, only: history_out_field
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      ! State variables
      real(kind_phys), intent(in) :: list(:,:)
      real(kind_phys), intent(in) :: of(:)
      real(kind_phys), intent(in) :: fields(:,:)
      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      ! Example:
      ! call history_out_field('TESTDIAG', of)

   end subroutine SCHEME_diagnostics_run

   !=======================================================================

end module SCHEME_diagnostics
