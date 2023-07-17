!< \section arg_table_micm
!! \htmlinclude micm.html
module micm
   use iso_c_binding
   use ccpp_kinds, only:  kind_phys

   implicit none

   private
   public :: micm_init, micm_run

   procedure(solver), pointer :: fsolver

   interface

      type(c_funptr) function get_solver(filepath) bind(c)
         import :: c_char, c_funptr
         character(len=1, kind=c_char), dimension(*), intent(in) :: filepath
      end function get_solver

      subroutine solver(state, state_size, time_step) bind(c)
         import :: c_ptr, c_double, c_int64_t
         real(c_double), dimension(*) :: state
         integer(c_int64_t), value :: state_size
         integer(c_int64_t), value :: time_step
      end subroutine

   end interface

contains

   !> \section arg_table_micm_init Argument Table
   !! \htmlinclude micm_init.html
   subroutine micm_init(filename_of_micm_configuration, errmsg, errflg)
      ! Arguments
      character(len=512), intent(out)  :: errmsg
      integer, intent(out)             :: errflg
      character(len=*), intent(in)     :: filename_of_micm_configuration

      ! Local variables
      type(c_funptr)                   :: csolver_func_pointer
      ! Convert Fortran character array to C character array
      character(len=len(filename_of_micm_configuration)+1, kind=c_char) :: c_filepath

      errmsg = ''
      errflg = 0

      c_filepath = transfer(filename_of_micm_configuration, c_filepath)

      csolver_func_pointer = get_solver(c_filepath)
      call c_f_procpointer(csolver_func_pointer, fsolver)
   end subroutine micm_init

   !> \section arg_table_micm_run Argument Table
   !! \htmlinclude micm_run.html
   subroutine micm_run(ccpp_num_constituents, timestep_for_physics, errmsg, errflg)
      integer, intent(in)              :: ccpp_num_constituents
      real(kind=kind_phys), intent(in) :: timestep_for_physics
      character(len=512), intent(out)  :: errmsg
      integer,            intent(out)  :: errflg

      ! Declare a temporary array of type c_double
      real(c_double), dimension(:), allocatable :: state_cdouble

      errmsg = ''
      errflg = 0

      ! Allocate and convert the state array to c_double
      allocate(state_cdouble(ccpp_num_constituents), stat=errflg)
      state_cdouble = 1

      ! Check if the allocation was successful
      if (errflg /= 0) then
         errmsg = "Failed to allocate memory to transfer the constituent concentrations to MICM."
         return
      endif

      ! call fsolver(state, state_size, time_step)
      call fsolver(state_cdouble, int(ccpp_num_constituents, c_int64_t), int(timestep_for_physics, c_int64_t))

      print *, "new state", state_cdouble
      
      ! Deallocate the temporary array
      deallocate(state_cdouble)

   end subroutine micm_run

end module micm
