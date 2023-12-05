module micm_wrapper
   ! Wrapper for MICM functionality

   use iso_c_binding
   use micm_core
   !  use ccpp_kinds, only:  kind_phys  !TODO(jiwon) - do we need this?
 
   implicit none
 
   private
   public :: micm_init, micm_run, micm_final

   type(micm_t), pointer :: micm_ptr

 contains
 
   !> \section arg_table_micm_init Argument Table
   !! \htmlinclude micm_init.html
   subroutine micm_init(config_path, iulog, errcode, errmsg)
      ! Initializes MICM objects by creating solvers from configure files

      character(len=*), intent(in)     :: config_path
      integer, intent(in)              :: iulog
      integer, intent(out)             :: errcode
      character(len=512), intent(out)  :: errmsg

      errcode = 0
      errmsg = ''

      ! Constructs MICM object
      allocate(micm_ptr)
      micm_ptr = micm_t(config_path)

      ! Creates solver
      errcode = micm_ptr%create_solver()

      if (errcode /= 0) then
         errmsg = "INIT MICM: FATAL: Failed in creating MICM solver because parsing configuration files failed. &
                  Please look over at the other file."  !TODO(jiwon) - which file?
      else 
         write(iulog,*) "INIT MICM: INFO: Created MICM solver"
      endif

   end subroutine micm_init

   !> \section arg_table_micm_run Argument Table
   !! \htmlinclude micm_run.html
   subroutine micm_run(temperature, pressure, time_step, concentrations, num_concentrations, iulog, errcode, errmsg)
      real(c_double), intent(in)                                :: temperature
      real(c_double), intent(in)                                :: pressure
      real(c_double), intent(in)                                :: time_step
      ! TODO(jiwon)
      ! real(c_double), dimension(*), intent(inout)  :: concentrations
      real(c_double), dimension(:,:,:), allocatable, intent(inout)  :: concentrations      
      integer, intent(in)                                       :: num_concentrations
      integer, intent(in)                                       :: iulog
      integer, intent(out)                                      :: errcode
      character(len=512), intent(out)                           :: errmsg

      errcode = 0
      errmsg = ''

      write(iulog,*) "RUN MICM: INFO: Running MICM solver..."
      call micm_ptr%solve(temperature, pressure, time_step, concentrations, num_concentrations)

   end subroutine micm_run

   !> \section arg_table_micm_final Argument Table
   !! \htmlinclude micm_final.html
   subroutine micm_final(iulog, errcode, errmsg)
      integer, intent(in)              :: iulog
      integer, intent(out)             :: errcode
      character(len=512), intent(out)  :: errmsg

      errcode = 0
      errmsg = ''
      
      write(iulog,*) "FINAL MICM: INFO: Deallocating MICM object..."
      deallocate(micm_ptr)

   end subroutine micm_final

 end module micm_wrapper
 