! Interstitial schemes for Beljaars drag.
!
! Put in a separate file than beljaars_drag.F90 because
! scheme files with multiple schemes cannot have a namelist file attached.
module beljaars_drag_interstitials

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: beljaars_add_wind_damping_rate_run
  public :: beljaars_add_updated_residual_stress_run
contains

  ! Add Beljaars drag to the wind damping rate for vertical diffusion.
  ! Has to run after vertical_diffusion_wind_damping_rate
!> \section arg_table_beljaars_add_wind_damping_rate_run Argument Table
!! \htmlinclude arg_table_beljaars_add_wind_damping_rate_run.html
  pure subroutine beljaars_add_wind_damping_rate_run( &
    ncol, pver, &
    dragblj, &
    tau_damp_rate, &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: dragblj(:, :)    ! Drag profile from Beljaars SGO form drag > 0. [s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: tau_damp_rate(:,:) ! Rate at which external (surface) stress damps wind speeds [s-1]

    ! Output arguments
    character(len=*),   intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    ! Beljaars et al SGO scheme incorporated here. It
    ! appears as a "3D" tau_damp_rate specification.
    tau_damp_rate(:ncol,:pver) = tau_damp_rate(:ncol,:pver) + dragblj(:ncol,:pver)

  end subroutine beljaars_add_wind_damping_rate_run

  ! Add Beljaars stress using updated provisional winds to surface stresses used
  ! for horizontal momentum diffusion - kinetic energy dissipation.
!> \section arg_table_beljaars_add_updated_residual_stress_run Argument Table
!! \htmlinclude arg_table_beljaars_add_updated_residual_stress_run.html
  subroutine beljaars_add_updated_residual_stress_run( &
    ncol, pver, &
    gravit, &
    p, &
    do_iss, &
    itaures, &
    dragblj, &
    u1, v1, & ! provisionally updated winds
    ! input/output
    tauresx, tauresy, &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys
    use coords_1d, only: Coords1D

    ! Input arguments
    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver
    real(kind_phys), intent(in)       :: gravit
    type(coords1d),  intent(in)       :: p                ! Pressure coordinates [Pa]
    logical,         intent(in)       :: do_iss           ! Flag for implicit surface stress (namelist from vdiff)
    logical,         intent(in)       :: itaures          ! Flag for updating tauresx tauresy in this subroutine.
    real(kind_phys), intent(in)       :: dragblj(:, :)    ! Drag profile from Beljaars SGO form drag > 0. [s-1]
    real(kind_phys), intent(in)       :: u1(:,:)          ! After vertical diffusion u-wind [m s-1]
    real(kind_phys), intent(in)       :: v1(:,:)          ! After vertical diffusion v-wind [m s-1]

    ! Input/output arguments
    real(kind_phys), intent(inout)    :: tauresx(:)       ! Partially updated residual surface stress using provisional winds
    real(kind_phys), intent(inout)    :: tauresy(:)
    character(len=*),   intent(out)   :: errmsg           ! error message
    integer,            intent(out)   :: errflg           ! error flag

    integer :: i, k
    real(kind_phys) :: taubljx(ncol)                      ! recomputed explicit/residual beljaars stress
    real(kind_phys) :: taubljy(ncol)                      ! recomputed explicit/residual beljaars stress

    errmsg = ''
    errflg = 0

    if(.not. do_iss) then
      ! Not added to residual stress if implicit surface stress is off
      return
    end if

    if(.not. itaures) then
      ! Not added to residual stress if residual stress is not requested to be updated.
      return
    end if

    do i = 1, ncol
      ! Add vertically-integrated Beljaars drag to residual stress
      ! these are calculated using provisionally-updated winds.
      ! These are used only locally and not written back to model state.
      taubljx(i) = 0._kind_phys
      taubljy(i) = 0._kind_phys
      do k = 1, pver
        taubljx(i) = taubljx(i) + (1._kind_phys/gravit)*dragblj(i, k)*u1(i, k)*p%del(i, k)
        taubljy(i) = taubljy(i) + (1._kind_phys/gravit)*dragblj(i, k)*v1(i, k)*p%del(i, k)
      end do
    end do

    tauresx(:ncol) = tauresx(:ncol) + taubljx(:ncol)
    tauresy(:ncol) = tauresy(:ncol) + taubljy(:ncol)

  end subroutine beljaars_add_updated_residual_stress_run

end module beljaars_drag_interstitials
