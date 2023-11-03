module micm_wrapper
   ! Wrapper for MICM functionality

   use iso_c_binding
   use micm_core
    ! use ccpp_kinds, only:  kind_phys  !TODO(jiwon) - do we need this?
 
   implicit none
 
   private
   public :: micm_init, micm_run, micm_final
   
   type(micm_t), pointer :: micm_ptr
   
   ! type :: micm_wrapper_t
   !    private
   !    type(micm_t),     pointer :: micm_
   ! ! contains
   ! !    final :: finalize
   ! end type micm_wrapper_t

   ! type(micm_wrapper_t), allocatable :: micm_wrappers

 contains
 
   !> \section arg_table_micm_init Argument Table
   !! \htmlinclude micm_init.html
   subroutine micm_init(config_path, errmsg, errflg)
      ! Initializes MICM objects by creating solvers from configure files

      character(len=*), intent(in)     :: config_path
      character(len=512), intent(out)  :: errmsg
      integer(c_size_t), intent(out)   :: errflg
 
      errmsg = ''
      errflg = 0
      
      ! Constructs MICM object
      allocate(micm_ptr)
      micm_ptr = micm_t(config_path)
      
      write(*,*) "  * [MICM wrapper] Created MICM solver "

      ! Creates solver
      errflg = micm_ptr%create_solver()  ! TODO(jiwon)

      ! allocate(micm_wrappers)
      ! micm_wrappers%micm_ = micm_t(config_path)

   end subroutine micm_init

   !> \section arg_table_micm_run Argument Table
   !! \htmlinclude micm_run.html
   subroutine micm_run(temperature, pressure, time_step, concentrations, num_concentrations, errmsg, errflg)
      real(c_double), intent(in) :: temperature
      real(c_double), intent(in) :: pressure
      real(c_double), intent(in) :: time_step
      real(c_double), dimension(*), intent(inout) :: concentrations
      integer(c_size_t), intent(in) :: num_concentrations
      character(len=512), intent(out) :: errmsg
      integer(c_size_t), intent(out) :: errflg

      errmsg = ''
      errflg = 0

      write(*,*) "  * [MICM wrapper] Running MICM solver... "
      call micm_ptr%solve(temperature, pressure, time_step, concentrations, num_concentrations)  ! TODO(jiwon): handling error

   end subroutine micm_run
 
   subroutine micm_final()
      write(*,*) "  * [MICM wrapper] Deallocating MICM object... "
      call micm_ptr%delete()
      deallocate(micm_ptr)

   end subroutine micm_final

 end module micm_wrapper
 