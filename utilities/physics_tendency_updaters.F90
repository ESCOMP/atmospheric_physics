module physics_tendency_updaters

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: apply_tendency_of_eastward_wind_run
  public :: apply_tendency_of_northward_wind_run
  public :: apply_heating_rate_run
  public :: apply_tendency_of_air_temperature_run

CONTAINS

   !> \section arg_table_apply_tendency_of_eastward_wind_run  Argument Table
   !! \htmlinclude apply_tendency_of_eastward_wind_run.html
   subroutine apply_tendency_of_eastward_wind_run(nz, dudt, u, dudt_total, dt,             &
        errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz              ! Num vertical  layers
      real(kind_phys),    intent(in)    :: dudt(:,:)       ! tendency of eastward wind
      real(kind_phys),    intent(inout) :: u(:,:)          ! eastward wind
      real(kind_phys),    intent(inout) :: dudt_total(:,:) ! total tendency of eastward wind
      real(kind_phys),    intent(in)    :: dt              ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         u(:, klev) = u(:, klev) + (dudt(:, klev) * dt)
         dudt_total(:, klev) = dudt_total(:, klev) + dudt(:, klev)
      end do

   end subroutine apply_tendency_of_eastward_wind_run

   !> \section arg_table_apply_tendency_of_northward_wind_run  Argument Table
   !! \htmlinclude apply_tendency_of_northward_wind_run.html
   subroutine apply_tendency_of_northward_wind_run(nz, dvdt, v, dvdt_total, dt,             &
        errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz              ! Num vertical layers
      real(kind_phys),    intent(in)    :: dvdt(:,:)       ! tendency of northward wind
      real(kind_phys),    intent(inout) :: v(:,:)          ! northward wind
      real(kind_phys),    intent(inout) :: dvdt_total(:,:) ! total tendency of northward wind
      real(kind_phys),    intent(in)    :: dt              ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         v(:, klev) = v(:, klev) + (dvdt(:, klev) * dt)
         dvdt_total(:, klev) = dvdt_total(:, klev) + dvdt(:, klev)
      end do

   end subroutine apply_tendency_of_northward_wind_run

   !> \section arg_table_apply_heating_rate_run  Argument Table
   !! \htmlinclude apply_heating_rate_run.html
   subroutine apply_heating_rate_run(nz, heating_rate, temp, dTdt_total, dt, cpair,  &
        errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz                ! Num vertical  layers
      real(kind_phys),    intent(in)    :: heating_rate(:,:) ! heating rate
      real(kind_phys),    intent(inout) :: temp(:,:)         ! air temperature
      real(kind_phys),    intent(inout) :: dTdt_total(:,:)   ! total temperature tend.
      real(kind_phys),    intent(in)    :: dt                ! physics time step
      real(kind_phys),    intent(in)    :: cpair(:,:)        ! specific heat, dry air
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         temp(:, klev) = temp(:, klev) + (heating_rate(:, klev) * dt / cpair(:, klev))
         dTdt_total(:, klev) = dTdt_total(:, klev) + (heating_rate(:, klev) / cpair(:,klev))
      end do

   end subroutine apply_heating_rate_run

   !> \section arg_table_apply_tendency_of_air_temperature_run  Argument Table
   !! \htmlinclude apply_tendency_of_air_temperature_run.html
   subroutine apply_tendency_of_air_temperature_run(nz, t_tend, temp, dtdT_total,   &
        dt, errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz              ! Num vertical  layers
      real(kind_phys),    intent(in)    :: t_tend(:,:)     ! temperature tendency
      real(kind_phys),    intent(inout) :: temp(:,:)       ! air temperature
      real(kind_phys),    intent(inout) :: dTdt_total(:,:) ! total temp. tendency
      real(kind_phys),    intent(in)    :: dt              ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         temp(:, klev) = temp(:, klev) + (t_tend(:, klev) * dt)
         dTdt_total(:, klev) = dTdt_total(:, klev) + t_tend(:, klev)
      end do

   end subroutine apply_tendency_of_air_temperature_run

end module physics_tendency_updaters
