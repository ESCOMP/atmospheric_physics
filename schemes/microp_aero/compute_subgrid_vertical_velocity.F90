! Compute raw subgrid-scale vertical velocity from TKE, KVH, or CLUBB WP2.
module compute_subgrid_vertical_velocity

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: compute_subgrid_vertical_velocity_tke_run
  public :: compute_subgrid_vertical_velocity_kvh_run
  public :: compute_subgrid_vertical_velocity_clubb_run

contains

!> \section arg_table_compute_subgrid_vertical_velocity_tke_run Argument Table
!! \htmlinclude compute_subgrid_vertical_velocity_tke_run.html

  ! Note: despite the tke name, the "tke" used here is the CAM5 diag_TKE scheme TKE
  ! and is not the TKE originating from CLUBB, and thus this scheme is not used for
  ! CAM6+, it is used for CAM5 only.
  ! For CAM6+ do not use TKE, derive from WP2_nadv - see the clubb variant. hplin 4/20/26
  subroutine compute_subgrid_vertical_velocity_tke_run( &
    ncol, pver, top_lev, &
    tke,                 &
    wsub,                &
    errmsg, errflg)

    ! Input arguments
    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: top_lev ! top vertical level for cloud physics

    real(kind_phys),  intent(in)  :: tke(:, :)  ! turbulent kinetic energy at interfaces [m2 s-2]

    ! Output arguments
    real(kind_phys),  intent(out) :: wsub(:, :) ! subgrid vertical velocity [m s-1]

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    wsub(:, :top_lev-1) = 0._kind_phys

    do k = top_lev, pver
      do i = 1, ncol
        wsub(i,k) = sqrt(0.5_kind_phys * (tke(i,k) + tke(i,k+1)) * (2._kind_phys / 3._kind_phys))
        wsub(i,k) = min(wsub(i,k), 10._kind_phys)
      end do
    end do

  end subroutine compute_subgrid_vertical_velocity_tke_run

!> \section arg_table_compute_subgrid_vertical_velocity_kvh_run Argument Table
!! \htmlinclude compute_subgrid_vertical_velocity_kvh_run.html
  subroutine compute_subgrid_vertical_velocity_kvh_run( &
    ncol, pver, top_lev, &
    kvh,                 &
    wsub,                &
    errmsg, errflg)

    ! Input arguments
    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: top_lev ! top vertical level for cloud physics

    real(kind_phys),  intent(in)  :: kvh(:, :)  ! eddy diffusivity for heat at interfaces [m2 s-1]

    ! Output arguments
    real(kind_phys),  intent(out) :: wsub(:, :) ! subgrid vertical velocity [m s-1]

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    wsub(:, :top_lev-1) = 0._kind_phys

    ! Get sub-grid vertical velocity from diffusion coefficient.
    ! Following Morrison et al. 2005, JAS.
    ! Assume mixing length of 30 m.
    ! Use maximum sub-grid vertical vel of 10 m/s.
    do k = top_lev, pver
      do i = 1, ncol
        wsub(i,k) = (kvh(i,k) + kvh(i,k+1)) / 2._kind_phys / 30._kind_phys
        wsub(i,k) = min(wsub(i,k), 10._kind_phys)
      end do
    end do

  end subroutine compute_subgrid_vertical_velocity_kvh_run

!> \section arg_table_compute_subgrid_vertical_velocity_clubb_run Argument Table
!! \htmlinclude compute_subgrid_vertical_velocity_clubb_run.html
  subroutine compute_subgrid_vertical_velocity_clubb_run( &
    ncol, pver, pverp, top_lev, &
    wp2,                 &
    wsub,                &
    errmsg, errflg)

    ! Input arguments
    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: pverp
    integer,          intent(in)  :: top_lev ! top vertical level for cloud physics

    real(kind_phys),  intent(in)  :: wp2(:, :)  ! CLUBB variance of vertical velocity at interfaces [m2 s-2]

    ! Output arguments
    real(kind_phys),  intent(out) :: wsub(:, :) ! subgrid vertical velocity [m s-1]

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k
    real(kind_phys) :: tke(ncol, pver+1)

    errmsg = ''
    errflg = 0

    ! Convert wp2 to TKE: tke = (3/2) * wp2 from [pver+1, top_lev]
    ! This matches CAM microp_aero.F90 CLUBB_SGS branch exactly.
    tke(:ncol,top_lev:pverp) = (3._kind_phys/2._kind_phys)*wp2(:ncol,1:pverp-top_lev+1)
    tke(:ncol,1:top_lev-1) = 0._kind_phys

    wsub(:, :top_lev-1) = 0._kind_phys

    do k = top_lev, pver
      do i = 1, ncol
        wsub(i,k) = sqrt(0.5_kind_phys * (tke(i,k) + tke(i,k+1)) * (2._kind_phys / 3._kind_phys))
        wsub(i,k) = min(wsub(i,k), 10._kind_phys)
      end do
    end do

  end subroutine compute_subgrid_vertical_velocity_clubb_run

end module compute_subgrid_vertical_velocity
