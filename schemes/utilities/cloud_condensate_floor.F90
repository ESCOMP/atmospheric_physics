! Zero cloud liquid and cloud ice, together with their number concentrations,
! where the condensate mass is below a floor. This mirrors the cloud limits
! CAM's physics_update applies after the deep convection update (the
! cldlim_names list, 1e-36 kg kg-1), so it belongs right after the constituent
! tendencies of deep convection are applied.
module cloud_condensate_floor
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: cloud_condensate_floor_init
  public :: cloud_condensate_floor_run

  ! Constituent indices, -1 when the constituent is not present
  integer :: ix_cldliq = -1
  integer :: ix_cldice = -1
  integer :: ix_numliq = -1
  integer :: ix_numice = -1

  real(kind_phys), parameter :: condensate_min = 1.e-36_kind_phys  ! mass floor [kg kg-1]

contains

!> \section arg_table_cloud_condensate_floor_init Argument Table
!! \htmlinclude cloud_condensate_floor_init.html
  subroutine cloud_condensate_floor_init(errmsg, errflg)
    use ccpp_scheme_utils, only: ccpp_constituent_index

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call ccpp_constituent_index('cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', &
                                ix_cldliq, errmsg=errmsg, errcode=errflg)
    if (errflg /= 0) return
    call ccpp_constituent_index('cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', &
                                ix_cldice, errmsg=errmsg, errcode=errflg)
    if (errflg /= 0) return
    call ccpp_constituent_index('mass_number_concentration_of_cloud_liquid_water_droplets_in_moist_air_and_condensed_water', &
                                ix_numliq, errmsg=errmsg, errcode=errflg)
    if (errflg /= 0) return
    call ccpp_constituent_index('mass_number_concentration_of_cloud_ice_water_crystals_in_moist_air_and_condensed_water', &
                                ix_numice, errmsg=errmsg, errcode=errflg)
    if (errflg /= 0) return

  end subroutine cloud_condensate_floor_init

!> \section arg_table_cloud_condensate_floor_run Argument Table
!! \htmlinclude cloud_condensate_floor_run.html
  subroutine cloud_condensate_floor_run(ncol, pver, const_q, errmsg, errflg)

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(inout) :: const_q(:,:,:)   ! constituent mixing ratios [kg kg-1 or kg-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    call apply_floor(ix_cldliq, ix_numliq)
    call apply_floor(ix_cldice, ix_numice)

  contains

    subroutine apply_floor(ix_mass, ix_num)
      integer, intent(in) :: ix_mass   ! condensate mass constituent index
      integer, intent(in) :: ix_num    ! matching number concentration index, or -1

      integer :: i, k

      if (ix_mass < 1) return
      do k = 1, pver
         do i = 1, ncol
            if (const_q(i,k,ix_mass) < condensate_min) then
               const_q(i,k,ix_mass) = 0._kind_phys
               if (ix_num > 0) const_q(i,k,ix_num) = 0._kind_phys
            end if
         end do
      end do
    end subroutine apply_floor

  end subroutine cloud_condensate_floor_run

end module cloud_condensate_floor
