module kessler_update

   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: kessler_update_init
   public :: kessler_update_timestep_init
   public :: kessler_update_run
   public :: kessler_update_timestep_final

   ! Private module variables
   real(kind_phys)    :: gravit

CONTAINS

   !> \section arg_table_kessler_update_init  Argument Table
   !! \htmlinclude kessler_update_init.html
   subroutine kessler_update_init(gravit_in, errmsg, errflg)
      real(kind_phys),    intent(in)  :: gravit_in
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      gravit = gravit_in

   end subroutine kessler_update_init

   !> \section arg_table_kessler_update_timestep_init  Argument Table
   !! \htmlinclude kessler_update_timestep_init.html
   subroutine kessler_update_timestep_init(temp, temp_prev, ttend_t,          &
        errmsg, errflg)

      real(kind_phys),    intent(in)  :: temp(:,:)
      real(kind_phys),    intent(out) :: temp_prev(:,:)
      real(kind_phys),    intent(out) :: ttend_t(:,:)
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      !   Initialize the previous temperature and its tendency to zero
      temp_prev(:,:)  = temp(:,:)
      ttend_t(:,:)    = 0._kind_phys

   end subroutine kessler_update_timestep_init

   !> \section arg_table_kessler_update_run  Argument Table
   !! \htmlinclude kessler_update_run.html
   subroutine kessler_update_run(nz, ncol, dt, theta, exner,                  &
        temp_prev, ttend_t, errmsg, errflg)

      integer,            intent(in)    :: nz
      integer,            intent(in)    :: ncol
      real(kind_phys),    intent(in)    :: dt             !time step
      real(kind_phys),    intent(in)    :: theta(:,:)     !potential temperature
      real(kind_phys),    intent(in)    :: exner(:,:)     !Exner function
      real(kind_phys),    intent(in)    :: temp_prev(:,:) !air temperature before kessler physics

      real(kind_phys),    intent(inout) :: ttend_t(:,:)   !total air temperature tendency due to
                                                          !kessler physics

      character(len=512), intent(out)   :: errmsg
      integer,            intent(out)   :: errflg

      !Local variables
      integer                           :: klev
      real(kind_phys)                   :: new_temp(ncol)

      errmsg = ''
      errflg = 0

      ! Back out tendencies from updated fields
      do klev = 1, nz
         new_temp(:ncol) = theta(:ncol,klev) * exner(:ncol,klev)
         ttend_t(:ncol,klev) = ttend_t(:ncol,klev) + ((new_temp(:ncol) - temp_prev(:ncol,klev)) / dt)
      end do

   end subroutine kessler_update_run

   !> \section arg_table_kessler_update_timestep_final  Argument Table
   !! \htmlinclude kessler_update_timestep_final.html
   subroutine kessler_update_timestep_final(nz, cpair, temp, zm, phis, st_energy, &
        errflg, errmsg)

      ! Dummy arguments
      integer,            intent(in)    :: nz
      real(kind_phys),    intent(in)    :: cpair(:,:) ! Specific_heat_of_dry_air_at_constant_pressure (J/kg/K)
      real(kind_phys),    intent(in)    :: temp(:,:)  ! Temperature
      real(kind_phys),    intent(in)    :: zm(:,:)
      real(kind_phys),    intent(in)    :: phis(:)
      real(kind_phys),    intent(out)   :: st_energy(:,:)

      character(len=512), intent(out)   :: errmsg
      integer,            intent(out)   :: errflg

      ! Local variable
      integer :: klev

      errmsg = ''
      errflg = 0

      do klev = 1, nz
         st_energy(:,klev) = (temp(:,klev) * cpair(:,klev)) + (gravit * zm(:,klev)) + &
              phis(:)
      end do

   end subroutine kessler_update_timestep_final


end module kessler_update
