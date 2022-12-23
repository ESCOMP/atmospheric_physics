module qneg
   ! qneg parameterization of qneg_module [qneg3 only] in CAM6
   !   -- will not be ported back to CAM6 --

   use ccpp_kinds,                only: kind_phys
   use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
!!XXgoldyXX: v Reinstate when history is implemented
#if 0
   use cam_history_support, only: max_fieldname_len
#endif
!!XXgoldyXX: ^ Reinstate when history is implemented

   implicit none
   private
   save

   ! Public interface.

   public :: qneg_init
   public :: qneg_run !qneg3 in CAM6
   public :: qneg_timestep_final
   public :: qneg_final

   ! Private module variables
   logical          :: collect_stats = .false.
   logical          :: timestep_reset = .false.
   logical          :: qneg_initialized = .false.
   integer          :: num_constituents = 0

   real(kind_phys), parameter :: tol = 1.e-12_kind_phys
   real(kind_phys), parameter :: worst_reset = 1.e35_kind_phys

   ! Diagnostic field names
   !  integer              :: num_diag_fields = 0
   !  character(len=max_fieldname_len) :: diag_names(num_diag_fields)
   !  logical             :: cnst_out_calc = .false.
   !  logical             :: cnst_outfld(num_diag_fields) = .false.

   ! Summary buffers
   integer, parameter    :: num_bins = 24 ! Max number of schemes calling qneg
   integer, parameter    :: CS = 80       ! Swiped from shr_kind_mod
   character(len=CS)     :: qneg_warn_labels(num_bins) = ''
   integer, allocatable  :: qneg_warn_num(:,:)
   real(kind_phys), allocatable :: qneg_warn_worst(:,:)

   private             :: qneg_print_summary
   !  private             :: calc_cnst_out
   private             :: find_index
   interface reset_stats
      module procedure reset_stats_scalar
      module procedure reset_stats_array
   end interface reset_stats

CONTAINS

   !> \section arg_table_qneg_init Argument Table
   !! \htmlinclude qneg_init.html
   subroutine qneg_init(print_qneg_warn, num_constituents_in, qprops,         &
        errcode, errmsg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      !use cam_history,    only: addfld, horiz_only

      character(len=*),                    intent(in)  :: print_qneg_warn
      integer,                             intent(in)  :: num_constituents_in
      ! qprops is the array of constituent properties
      type(ccpp_constituent_prop_ptr_t), intent(in)  :: qprops(:)
      integer,                             intent(out) :: errcode
      character(len=512),                  intent(out) :: errmsg

      character(len=*), parameter  :: subname = 'qneg_init'
      integer                      :: ierr

      errcode = 0
      errmsg = ''
      num_constituents = num_constituents_in

      !Check if already initialized:
      if (qneg_initialized) then
         return
      end if

      !    integer :: index

      !    do index = 1, num_constituents
      !       diag_names(index) = trim(cnst_name(index))//'_qneg'
      !       call addfld(diag_names(index), (/ 'lev' /), 'I', 'kg/kg',              &
      !            trim(cnst_longname(index))//' QNEG error (cell)')
      !       diag_names(num_constituents+index) = trim(cnst_name(index))//'_qneg_col'
      !       call addfld(diag_names(num_constituents+index), horiz_only, 'I', 'kg/kg',         &
      !            trim(cnst_longname(index))//' QNEG error (column)')
      !    end do
      !    diag_names((2*num_constituents) + 1) = 'qflux_exceeded'
      !    call addfld(diag_names((2*num_constituents) + 1), horiz_only, 'I', 'kg/m^2/s',     &
      !         'qflux excess (QNEG4)')

      !Allocate and initialize arrays whose dimensions depend on num_constituents:
      allocate(qneg_warn_num(num_constituents, num_bins), stat=ierr)
      if (ierr /= 0) then
         errcode = ierr
         errmsg = trim(subname)//': Allocate of qneg_warn_num failed'
         return
      end if
      allocate(qneg_warn_worst(num_constituents, num_bins), stat=ierr)
      if (ierr /= 0) then
         errcode = ierr
         errmsg = trim(subname)//': Allocate of qneg_warn_worst failed'
         return
      end if
      qneg_warn_num = 0
      qneg_warn_worst = worst_reset

      select case(trim(print_qneg_warn))
      case('summary')
         collect_stats = .true.
         timestep_reset = .false.
      case('timestep')
         collect_stats = .true.
         timestep_reset = .true.
      case('off')
         collect_stats = .false.
         timestep_reset =.false.
      case default
         errcode = 1
         errmsg = subname//"FATAL: '"//trim(print_qneg_warn)//                &
              "' is not a valid value for print_qneg_warn"
      end select

      !Set qneg_initialized to .true.
      qneg_initialized = .true.

   end subroutine qneg_init

!!XXgoldyXX: v Reinstate when history is implemented
#if 0
   subroutine calc_cnst_out()
      use cam_history, only: hist_fld_active, history_initialized
      integer :: index

      if (history_initialized()) then
         ! to protect against routines that call qneg too early
         do index = 1, num_diag_fields
            cnst_outfld(index) = hist_fld_active(trim(diag_names(index)))
         end do
         cnst_out_calc = .true.
      end if

   end subroutine calc_cnst_out
#endif
!!XXgoldyXX: ^ Reinstate when history is implemented

   integer function find_index(nam) result(index)
      ! Find a valid or new index for 'nam' entries
      character(len=*),  intent(in) :: nam

      integer                      :: i

      index = -1
      do i = 1, num_bins
         if (trim(nam) == trim(qneg_warn_labels(i))) then
            ! We found this entry, return its index
            index = i
            exit
         else if (len_trim(qneg_warn_labels(i)) == 0) then
            ! We have run out of known entries, use a new one and reset its stats
            qneg_warn_labels(i) = nam
            index = i
            call reset_stats(qneg_warn_num(:, index), qneg_warn_worst(:,index))
            exit
         end if
      end do
   end function find_index

   !> \section arg_table_qneg_run Argument Table
   !! \htmlinclude qneg_run.html
   subroutine qneg_run(subnam, loop_begin, loop_end, lver, num_constituents,  &
        qmin, qprops, q, errcode, errmsg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
!!XXgoldyXX: v Reinstate when history is implemented
#if 0
      use cam_history, only: outfld
#endif
!!XXgoldyXX: ^ Reinstate when history is implemented
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Check moisture and tracers for minimum value, reset any below
      ! minimum value to minimum value and return information to allow
      ! warning message to be printed. The global average is NOT preserved.
      !
      !
      !-----------------------------------------------------------------------

      !------------------------------Arguments--------------------------------
      !
      ! Input arguments
      !
      character(len=*), intent(in) :: subnam ! name of calling routine

      integer,                           intent(in)    :: loop_begin
      integer,                           intent(in)    :: loop_end
      ! lver: Number of vertical levels in column
      integer,                           intent(in)    :: lver
      ! num_constituents: Number of constituents
      integer,                           intent(in)    :: num_constituents
      ! qmin: Global minimum constituent concentration
      real(kind_phys),                   intent(in)    :: qmin(:)
      ! qprops: The array of constituent properties
      type(ccpp_constituent_prop_ptr_t), intent(in)    :: qprops(:)
      ! q: The array of constituents
      real(kind_phys),                   intent(inout) :: q(:,:,:)
      integer,                           intent(out)   :: errcode
      character(len=512),                intent(out)   :: errmsg
      !
      !---------------------------Local workspace-----------------------------
      !
      integer  :: nvals            ! number of values found < qmin
      integer  :: i, k             ! longitude, level indices
      integer  :: index            ! For storing stats
      integer  :: m                ! constituent index
      integer  :: iw,kw            ! i,k indices of worst violator

      logical  :: found            ! true => at least 1 minimum violator found

      real(kind_phys) :: badvals(loop_begin:loop_end, lver) ! Collector for outfld calls
      real(kind_phys) :: badcols(loop_begin:loop_end)  ! Column sum for outfld
      real(kind_phys) :: worst           ! biggest violator
      !
      !-----------------------------------------------------------------------
      !

      errcode = 0
      errmsg = ''

      ! The first time we call this, we need to determine whether to call outfld
      !    if (.not. cnst_out_calc) then
      !       call calc_cnst_out()
      !    end if

      if (collect_stats) then
         index = find_index(trim(subnam))
      else
         index = -1
      end if

      do m = 1, num_constituents
         nvals = 0
         found = .false.
         worst = worst_reset
         badvals(:,:) = 0.0_kind_phys
         iw = -1
         kw = -1
         !
         ! Test all field values for being less than minimum value. Set q = qmin
         ! for all such points. Trace offenders and identify worst one.
         !

         do k = 1, lver
            do i = loop_begin, loop_end
               if (q(i,k,m) < qmin(m)) then
                  found = .true.
                  nvals = nvals + 1
                  badvals(i, k) = q(i, k, m)
                  if (index > 0) then
                     qneg_warn_num(m, index) = qneg_warn_num(m, index) + 1
                  end if
                  if (q(i,k,m) < worst) then
                     worst = q(i,k,m)
                     iw = i
                     kw = k
                     if (index > 0) then
                        qneg_warn_worst(m, index) = worst
                     end if
                  end if
                  q(i,k,m) = qmin(m)
               end if
            end do
         end do
!!XXgoldyXX: v Reinstate when history is implemented
#if 0
         ! Maybe output bad values
         if ((cnst_outfld(m)) .and. (worst < worst_reset)) then
            call outfld(trim(diag_names(m)), badvals, pcols, idx)
         end if
         if ((cnst_outfld(num_constituents+m)) .and. (worst < worst_reset)) then
            do i = 1, pcols
               badcols(i) = SUM(badvals(i,:))
            end do
            call outfld(trim(diag_names(num_constituents+m)), badcols, pcols, idx)
         end if
#endif
!!XXgoldyXX: ^ Reinstate when history is implemented
      end do

   end subroutine qneg_run

   !> \section arg_table_qneg_timestep_final Argument Table
   !! \htmlinclude qneg_timestep_final.html
   subroutine qneg_timestep_final(mpi_communicator, rootprocid, isrootproc,   &
        iulog, log_output, qprops, errcode, errmsg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      integer,                           intent(in)  :: mpi_communicator
      integer,                           intent(in)  :: rootprocid
      logical,                           intent(in)  :: isrootproc
      integer,                           intent(in)  :: iulog
      logical,                           intent(in)  :: log_output
      ! qprops: The array of constituent properties
      type(ccpp_constituent_prop_ptr_t), intent(in)  :: qprops(:)
      integer,                           intent(out) :: errcode
      character(len=512),                intent(out) :: errmsg

      errcode = 0
      errmsg = ''

      if (timestep_reset .and. collect_stats) then
         call qneg_print_summary(mpi_communicator, rootprocid, isrootproc,    &
              iulog, log_output, qprops)
      end if

   end subroutine qneg_timestep_final
   !> \section arg_table_qneg_final Argument Table
   !! \htmlinclude qneg_final.html
   subroutine qneg_final(mpi_communicator, rootprocid, isrootproc,            &
        iulog, log_output, qprops, errcode, errmsg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      integer,                           intent(in)  :: mpi_communicator
      integer,                           intent(in)  :: rootprocid
      logical,                           intent(in)  :: isrootproc
      integer,                           intent(in)  :: iulog
      logical,                           intent(in)  :: log_output
      ! qprops is the array of constituent properties
      type(ccpp_constituent_prop_ptr_t), intent(in)  :: qprops(:)
      integer,                           intent(out) :: errcode
      character(len=512),                intent(out) :: errmsg

      errcode = 0
      errmsg = ''

      if (.not.timestep_reset .and. collect_stats) then
         call qneg_print_summary(mpi_communicator, rootprocid, isrootproc,    &
              iulog, log_output, qprops)
      end if
      deallocate(qneg_warn_num)
      deallocate(qneg_warn_worst)

   end subroutine qneg_final

   subroutine qneg_print_summary(mpi_communicator, rootprocid, isrootproc,    &
        iulog, log_output, qprops)
      use mpi, only: MPI_MIN, MPI_SUM, MPI_INTEGER, MPI_REAL8
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      integer,                           intent(in) :: mpi_communicator
      integer,                           intent(in) :: rootprocid
      logical,                           intent(in) :: isrootproc
      integer,                           intent(in) :: iulog
      logical,                           intent(in) :: log_output
      ! qprops is the array of constituent properties
      type(ccpp_constituent_prop_ptr_t), intent(in) :: qprops(:)

      integer            :: global_warn_num(num_constituents)
      real(kind_phys)    :: global_warn_worst(num_constituents)
      integer            :: index, m
      integer            :: ierr
      character(len=512) :: errmsg
      character(len=CS)  :: cnst_name

      do index = 1, num_bins
         ! QNEG
         call reset_stats(global_warn_num(:), global_warn_worst(:))
         call MPI_REDUCE(qneg_warn_num(:, index), global_warn_num(:),         &
              num_constituents, MPI_INTEGER, MPI_SUM, rootprocid,             &
              mpi_communicator, ierr)
         call MPI_REDUCE(qneg_warn_worst(:, index), global_warn_worst(:),     &
              num_constituents, MPI_REAL8, MPI_MIN, rootprocid,               &
              mpi_communicator, ierr)
         if (isrootproc) then
            do m = 1, num_constituents
               if ( (global_warn_num(m) > 0) .and.                            &
                    (abs(global_warn_worst(m)) > tol)) then
                  call qprops(m)%standard_name(cnst_name, ierr, errmsg)
                  write(iulog, 9100) trim(qneg_warn_labels(index)),           &
                       trim(cnst_name), global_warn_num(m),                   &
                       global_warn_worst(m)
               end if
               flush(iulog)
            end do
         end if
         call reset_stats(qneg_warn_num(:,index), qneg_warn_worst(:,index))
      end do

      return
9100  format(' QNEG from ', a, ':', a, &
           ' Min. mixing ratio violated at ', i9, ' points. Worst = ', e10.1)
   end subroutine qneg_print_summary

   subroutine reset_stats_array(num_array, worst_array)
      ! Private routine to reset statistics
      integer,  intent(inout) :: num_array(:)
      real(kind_phys), intent(inout) :: worst_array(:)

      num_array(:)    = 0
      worst_array(:)  = worst_reset
   end subroutine reset_stats_array

   subroutine reset_stats_scalar(num, worst)
      ! Private routine to reset statistics
      integer,  intent(inout) :: num
      real(kind_phys), intent(inout) :: worst

      num    = 0
      worst  = worst_reset
   end subroutine reset_stats_scalar

end module qneg
