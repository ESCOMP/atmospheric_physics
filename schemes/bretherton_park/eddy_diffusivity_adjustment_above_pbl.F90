! Adjusts eddy diffusivities above the planetary boundary layer
! to address potentially excessive values in the free troposphere and upper atmosphere
! in the UW PBL scheme.
module eddy_diffusivity_adjustment_above_pbl
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: eddy_diffusivity_adjustment_above_pbl_run

contains

  ! Adjust eddy diffusivities above the PBL to prevent excessive values
  ! in the free troposphere and upper atmosphere.
  !
  ! The diffusivities from diag_TKE (UW) can be much larger than from HB in the free
  ! troposphere and upper atmosphere. These seem to be larger than observations,
  ! and in WACCM the gw_drag code is already applying an eddy diffusivity in the
  ! upper atmosphere. Optionally, adjust the diffusivities in the free troposphere
  ! or the upper atmosphere.
  !
  ! NOTE: Further investigation should be done as to why the diffusivities are
  ! larger in diag_TKE.
!> \section arg_table_eddy_diffusivity_adjustment_above_pbl_run Argument Table
!! \htmlinclude eddy_diffusivity_adjustment_above_pbl_run.html
  subroutine eddy_diffusivity_adjustment_above_pbl_run( &
    ncol, pverp, &
    kv_top_pressure, &
    kv_freetrop_scale, &
    kv_top_scale, &
    zi, pint, pblh, &
    kvh, kvm, kvq, &
    errmsg, errflg)

    ! Input arguments
    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pverp
    real(kind_phys),  intent(in)    :: kv_top_pressure      ! Upper atmosphere pressure threshold [Pa]
    real(kind_phys),  intent(in)    :: kv_freetrop_scale    ! Free troposphere scale factor
    real(kind_phys),  intent(in)    :: kv_top_scale         ! Upper atmosphere scale factor
    real(kind_phys),  intent(in)    :: zi(:, :)             ! Geopotential height at interfaces [m]
    real(kind_phys),  intent(in)    :: pint(:, :)           ! Pressure at interfaces [Pa]
    real(kind_phys),  intent(in)    :: pblh(:)              ! Planetary boundary layer height [m]

    ! Input/Output arguments
    real(kind_phys),  intent(inout) :: kvh(:, :)  ! Eddy diffusivity for heat at interfaces [m2 s-1]
    real(kind_phys),  intent(inout) :: kvm(:, :)  ! Eddy diffusivity for momentum at interfaces [m2 s-1]
    real(kind_phys),  intent(inout) :: kvq(:, :)  ! Eddy diffusivity for constituents at interfaces [m2 s-1]

    ! Output arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    ! Only apply adjustment if scaling factors differ from unity
    if ((kv_freetrop_scale /= 1.0_kind_phys) .or. &
        ((kv_top_scale /= 1.0_kind_phys) .and. (kv_top_pressure > 0.0_kind_phys))) then
      do i = 1, ncol
        do k = 1, pverp
          ! Outside of the boundary layer?
          if (zi(i, k) > pblh(i)) then
            ! In the upper atmosphere?
            if (pint(i, k) <= kv_top_pressure) then
              kvh(i, k) = kvh(i, k) * kv_top_scale
              kvm(i, k) = kvm(i, k) * kv_top_scale
              kvq(i, k) = kvq(i, k) * kv_top_scale
            else
              ! In the free troposphere
              kvh(i, k) = kvh(i, k) * kv_freetrop_scale
              kvm(i, k) = kvm(i, k) * kv_freetrop_scale
              kvq(i, k) = kvq(i, k) * kv_freetrop_scale
            end if
          else ! within PBL.
            exit
          end if
        end do
      end do
    end if

  end subroutine eddy_diffusivity_adjustment_above_pbl_run

end module eddy_diffusivity_adjustment_above_pbl
