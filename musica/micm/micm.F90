module micm
   ! Wrapper for MICM functionality

   use iso_c_binding

   ! Note: "micm_core" is included in an external pre-built MICM library that the host
   ! model is responsible for linking to during compilation
   use micm_core, only: micm_t
   use ccpp_kinds, only: kind_phys  ! TODO(jiwon) - temporary solution until the framework
                                    ! can handle the kind conversions automatically

   implicit none
 
   private
   public :: micm_init, micm_run, micm_final

   type(micm_t), allocatable :: micm_obj

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
      allocate(micm_obj)
      micm_obj = micm_t(config_path)

      ! Creates solver
      errcode = micm_obj%create_solver()

      if (errcode /= 0) then
         errmsg = "INIT MICM: FATAL: Failed in creating MICM solver because parsing configuration files failed. &
                  Please look over at the other file."  !TODO(jiwon) - which file?
      else 
         write(iulog,*) "INIT MICM: INFO: Created MICM solver"
      endif

   end subroutine micm_init

   !> \section arg_table_micm_run Argument Table
   !! \htmlinclude micm_run.html
   subroutine micm_run(temperature, pressure, time_step, num_concentrations, concentrations, iulog, errcode, errmsg)
      real(kind_phys), dimension(:,:), intent(in)       :: temperature
      real(kind_phys), dimension(:,:), intent(in)       :: pressure
      real(kind_phys), intent(in)                       :: time_step
      integer, intent(in)                               :: num_concentrations
      real(kind_phys), dimension(:,:,:), intent(inout)  :: concentrations      
      integer, intent(in)                               :: iulog
      integer, intent(out)                              :: errcode
      character(len=512), intent(out)                   :: errmsg

      ! local variables
      real(c_double), dimension(size(temperature, dim=1), &
                                size(temperature, dim=2))     :: c_temperature
      real(c_double), dimension(size(pressure, dim=1), &
                                size(pressure, dim=2))        :: c_pressure
      real(c_double)                                          :: c_time_step
      real(c_double), dimension(size(concentrations, dim=1), &
                                size(concentrations, dim=2), &
                                size(concentrations, dim=3))  :: c_concentrations
      integer                                                 :: horizontal_loop_extent, vertical_layer_dimension
      integer                                                 :: i_column, i_layer

      c_temperature = real(temperature, c_double)
      c_pressure = real(pressure, c_double)
      c_time_step = real(time_step, c_double)
      c_concentrations = real(concentrations, c_double)
      horizontal_loop_extent = size(concentrations, dim=1)
      vertical_layer_dimension = size(concentrations, dim=2)

      errcode = 0
      errmsg = ''

      write(iulog,*) "RUN MICM: INFO: Running MICM solver..."
      do i_column = 1, horizontal_loop_extent
         do i_layer = 1, vertical_layer_dimension
            call micm_obj%solve(c_temperature(i_column, i_layer), c_pressure(i_column, i_layer), c_time_step, &
                                num_concentrations, c_concentrations(i_column, i_layer, :))
         end do
      end do
      write(iulog,*) "RUN MICM: INFO: MICM solver has finished."

      concentrations = real(c_concentrations, kind_phys)

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
      if (allocated(micm_obj)) deallocate(micm_obj)

   end subroutine micm_final

 end module micm
 