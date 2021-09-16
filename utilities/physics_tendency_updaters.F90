module physics_tendency_updaters

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: apply_tendency_of_x_wind_run
  public :: apply_tendency_of_y_wind_run
  public :: apply_heating_rate_run
  public :: apply_tendency_of_air_temperature_run
  public :: update_dry_static_energy_init
  public :: update_dry_static_energy_run

  ! Private module variables

  real(kind_phys) :: gravit = -HUGE(1.0_kind_phys)

CONTAINS

   !> \section arg_table_apply_tendency_of_x_wind_run  Argument Table
   !! \htmlinclude apply_tendency_of_x_wind_run.html
   subroutine apply_tendency_of_x_wind_run(nz, du, u, dudt, dt,             &
        errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz        ! Num vertical  layers
      real(kind_phys),    intent(in)    :: du(:,:)   ! tendency of x wind
      real(kind_phys),    intent(inout) :: u(:,:)    ! x wind
      real(kind_phys),    intent(inout) :: dudt(:,:) ! total tendency of x wind
      real(kind_phys),    intent(in)    :: dt        ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         u(:, klev) = u(:, klev) + (du(:, klev) * dt)
         dudt(:, klev) = dudt(:, klev) + du(:, klev)
      end do

   end subroutine apply_tendency_of_x_wind_run

   !> \section arg_table_apply_tendency_of_y_wind_run  Argument Table
   !! \htmlinclude apply_tendency_of_y_wind_run.html
   subroutine apply_tendency_of_y_wind_run(nz, dv, v, dvdt, dt,             &
        errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz        ! Num vertical  layers
      real(kind_phys),    intent(in)    :: dv(:,:)   ! tendency of y wind
      real(kind_phys),    intent(inout) :: v(:,:)    ! y wind
      real(kind_phys),    intent(inout) :: dvdt(:,:) ! total tendency of y wind
      real(kind_phys),    intent(in)    :: dt        ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         v(:, klev) = v(:, klev) + (dv(:, klev) * dt)
         dvdt(:, klev) = dvdt(:, klev) + dv(:, klev)
      end do

   end subroutine apply_tendency_of_y_wind_run

   !> \section arg_table_apply_heating_rate_run  Argument Table
   !! \htmlinclude apply_heating_rate_run.html
   subroutine apply_heating_rate_run(nz, s, temp, dtdt, dt, cpair,            &
        errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz        ! Num vertical  layers
      real(kind_phys),    intent(in)    :: s(:,:)    ! heating rate
      real(kind_phys),    intent(inout) :: temp(:,:) ! air temperature
      real(kind_phys),    intent(inout) :: dtdt(:,:) ! total temperature tend.
      real(kind_phys),    intent(in)    :: dt        ! physics time step
      real(kind_phys),    intent(in)    :: cpair     ! specific heat, dry air
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         temp(:, klev) = temp(:, klev) + (s(:, klev) * dt / cpair)
         dtdt(:, klev) = dtdt(:, klev) + (s(:, klev) / cpair)
      end do

   end subroutine apply_heating_rate_run

   !> \section arg_table_apply_tendency_of_air_temperature_run  Argument Table
   !! \htmlinclude apply_tendency_of_air_temperature_run.html
   subroutine apply_tendency_of_air_temperature_run(nz, t_tend, temp, dtdt,   &
        dt, errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz          ! Num vertical  layers
      real(kind_phys),    intent(in)    :: t_tend(:,:) ! temperature tendency
      real(kind_phys),    intent(inout) :: temp(:,:)   ! air temperature
      real(kind_phys),    intent(inout) :: dtdt(:,:)   ! total temp. tendency
      real(kind_phys),    intent(in)    :: dt          ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         temp(:, klev) = temp(:, klev) + (t_tend(:, klev) * dt)
         dtdt(:, klev) = dtdt(:, klev) + t_tend(:, klev)
      end do

   end subroutine apply_tendency_of_air_temperature_run

   !> \section arg_table_update_dry_static_energy_init  Argument Table
   !! \htmlinclude update_dry_static_energy_init.html
   subroutine update_dry_static_energy_init(gravit_in)
      real(kind_phys),    intent(in)    :: gravit_in

      gravit = gravit_in

   end subroutine update_dry_static_energy_init

   !> \section arg_table_update_dry_static_energy_run  Argument Table
   !! \htmlinclude update_dry_static_energy_run.html
   subroutine update_dry_static_energy_run(nz, temp, zm, phis, st_energy,     &
        cpair, errcode, errmsg)

      ! Dummy arguments
      integer,            intent(in)  :: nz             ! Num vertical  layers
      real(kind_phys),    intent(in)  :: temp(:,:)      ! air temperature
      real(kind_phys),    intent(in)  :: zm(:,:)        ! geopotential height
      real(kind_phys),    intent(in)  :: phis(:)        ! surface geopotential
      real(kind_phys),    intent(out) :: st_energy(:,:) ! dry static energy
      real(kind_phys),    intent(in)  :: cpair          ! specific heat, dry air
      integer,            intent(out) :: errcode
      character(len=512), intent(out) :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         st_energy(:, klev) = (temp(:, klev) * cpair) +                       &
              (gravit * zm(:, klev)) + phis(:)
      end do

   end subroutine update_dry_static_energy_run

end module physics_tendency_updaters
