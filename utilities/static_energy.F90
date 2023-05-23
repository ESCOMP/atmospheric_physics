module static_energy

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: update_dry_static_energy_init
  public :: update_dry_static_energy_run

  ! Private module variables

  real(kind_phys) :: gravit = -HUGE(1.0_kind_phys)

CONTAINS

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
      real(kind_phys),    intent(in)  :: cpair(:,:)     ! specific heat, dry air
      integer,            intent(out) :: errcode
      character(len=512), intent(out) :: errmsg

      ! Local variable
      integer :: klev

      errcode = 0
      errmsg = ''

      do klev = 1, nz
         st_energy(:, klev) = (temp(:, klev) * cpair(:,klev)) +               &
              (gravit * zm(:, klev)) + phis(:)
      end do

   end subroutine update_dry_static_energy_run

end module static_energy
