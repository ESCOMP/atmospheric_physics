module rrtmgp_diagnostics

   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: rrtmgp_diagnostics_init ! init routine
   public :: rrtmgp_diagnostics_run  ! main routine

CONTAINS

   !> \section arg_table_rrtmgp_diagnostics_init  Argument Table
   !! \htmlinclude rrtmgp_diagnostics_init.html
   subroutine rrtmgp_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! Heating rate needed for d(theta)/dt computation
      call history_add_field('HR',   'Heating rate needed for d(theat)/dt computation', 'lev', 'avg', 'K s-1')

   end subroutine rrtmgp_diagnostics_init

   !> \section arg_table_rrtmgp_diagnostics_run  Argument Table
   !! \htmlinclude rrtmgp_diagnostics_run.html
   subroutine rrtmgp_diagnostics_run(write_output, ncol, pver, cappa, cpair, pmid, qrs, qrl, errmsg, errflg)

      use cam_history,        only: history_out_field
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      logical,                        intent(in) :: write_output ! Flag to write output for radiation
      integer,                        intent(in) :: ncol         ! Number of horizontal points
      integer,                        intent(in) :: pver         ! Number of vertical layers
      real(kind_phys),                intent(in) :: cappa        ! Ratio of dry air gas constant to specific heat of dry air at constant pressure
      real(kind_phys),                intent(in) :: cpair        ! Specific heat of dry air [J kg-1 K-1]
      real(kind_phys),                intent(in) :: pmid(:,:)    ! Air pressure at layer midpoints [Pa]
      real(kind_phys),                intent(in) :: qrs(:,:)     ! Shortwave heating rate [J kg-1 s-1]
      real(kind_phys),                intent(in) :: qrl(:,:)     ! Longwave heating rate [J kg-1 s-1]
      
      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer :: idx, kdx
      real(kind_phys) :: ftem(ncol, pver)

      errmsg = ''
      errflg = 0

      ! Don't do anything if this subcycle is inactive or we're not configured to write radiation output
      if ((.not. active_calls(diag_index)) .or. (.not. write_output)) then
         return
      end if

      ! Compute heating rate for dtheta/dt
      do kdx = 1, pver
         do idx = 1, ncol
            ftem(idx,kdx) = (qrs(idx,kdx) + qrl(idx,kdx))/cpair * (1.e5_kind_phys/pmid(idx,kdx))**cappa
         end do
      end do

      call history_out_field('HR', ftem)

   end subroutine rrtmgp_diagnostics_run

end module rrtmgp_diagnostics
