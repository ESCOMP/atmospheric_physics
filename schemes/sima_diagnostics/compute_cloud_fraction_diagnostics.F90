! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Diagnostics for cloud fraction
! Haipeng Lin, April 2025
module compute_cloud_fraction_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: compute_cloud_fraction_diagnostics_init
   public :: compute_cloud_fraction_diagnostics_run

contains

   !> \section arg_table_compute_cloud_fraction_diagnostics_init  Argument Table
   !! \htmlinclude compute_cloud_fraction_diagnostics_init.html
   subroutine compute_cloud_fraction_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('CLDST', 'stratiform_cloud_area_fraction', 'lev', 'avg', 'fraction')

   end subroutine compute_cloud_fraction_diagnostics_init

   !> \section arg_table_compute_cloud_fraction_diagnostics_run  Argument Table
   !! \htmlinclude compute_cloud_fraction_diagnostics_run.html
   subroutine compute_cloud_fraction_diagnostics_run( &
      cldst, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input parameters
      real(kind_phys),    intent(in)  :: cldst(:,:)


      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('CLDST', cldst)

   end subroutine compute_cloud_fraction_diagnostics_run

end module compute_cloud_fraction_diagnostics
