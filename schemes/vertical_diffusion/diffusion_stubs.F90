! Stubs for schemes that have not been CCPPized but their interfaces
! to vertical diffusion have been.
! These schemes will eventually be moved to the corresponding scheme
! directory after their actual computations have been CCPPized.
!
! CCPPized: Haipeng Lin, May 2025
module diffusion_stubs
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! CCPP-compliant subroutines
  public :: zero_upper_boundary_condition_init

  public :: tms_beljaars_zero_stub_run
  public :: turbulent_mountain_stress_add_drag_coefficient_run
  public :: beljaars_add_wind_damping_rate_run

  public :: beljaars_add_updated_residual_stress_run
  public :: turbulent_mountain_stress_add_updated_surface_stress_run

  public :: vertical_diffusion_not_use_rairv_init

contains

  ! Stub for zero upper boundary conditions before UBC module is CCPPized
!> \section arg_table_zero_upper_boundary_condition_init Argument Table
!! \htmlinclude zero_upper_boundary_condition_init.html
  subroutine zero_upper_boundary_condition_init( &
    ncol, pcnst, &
    ! below output
    ubc_mmr, &
    cnst_fixed_ubc, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pcnst

    ! Output arguments
    real(kind_phys),    intent(out) :: ubc_mmr(:,:)             ! Upper boundary condition mass mixing ratios [kg kg-1]
    logical,            intent(out) :: cnst_fixed_ubc(:)         ! Flag for fixed upper boundary condition of constituents [flag]
    character(len=512), intent(out) :: errmsg                   ! Error message
    integer,            intent(out) :: errflg                   ! Error flag

    errmsg = ''
    errflg = 0

    ! Set all upper boundary condition mixing ratios to zero
    ubc_mmr(:ncol, :pcnst) = 0._kind_phys

    ! Set all fixed upper boundary condition flags to false
    cnst_fixed_ubc(:pcnst) = .false.

  end subroutine zero_upper_boundary_condition_init

  ! Stub for TMS/Beljaars to be set to zero while they are not implemented.
!> \section arg_table_tms_beljaars_zero_stub_run Argument Table
!! \htmlinclude tms_beljaars_zero_stub_run.html
  subroutine tms_beljaars_zero_stub_run( &
    ncol, pver, &
    ksrftms, &
    tautmsx, tautmsy, &
    do_beljaars, &
    dragblj, &
    taubljx, taubljy, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver

    ! Output arguments
    real(kind_phys),    intent(out) :: ksrftms(:)  ! Surface drag coefficient for turbulent mountain stress. > 0. [kg m-2 s-1]
    real(kind_phys),    intent(out) :: tautmsx(:)  ! Eastward turbulent mountain surface stress [N m-2]
    real(kind_phys),    intent(out) :: tautmsy(:)  ! Northward turbulent mountain surface stress [N m-2]
    logical,            intent(out) :: do_beljaars
    real(kind_phys),    intent(out) :: dragblj(:,:)! Drag profile from Beljaars SGO form drag > 0. [s-1]
    real(kind_phys),    intent(out) :: taubljx(:)  ! Eastward Beljaars surface stress [N m-2]
    real(kind_phys),    intent(out) :: taubljy(:)  ! Northward Beljaars surface stress [N m-2]
    character(len=512), intent(out) :: errmsg      ! Error message
    integer,            intent(out) :: errflg      ! Error flag

    errmsg = ''
    errflg = 0

    ! Set TMS drag coefficient to zero (stub implementation)
    ksrftms(:ncol) = 0._kind_phys

    ! Set all TMS and Beljaars stresses to zero (stub implementation)
    tautmsx(:ncol) = 0._kind_phys
    tautmsy(:ncol) = 0._kind_phys
    taubljx(:ncol) = 0._kind_phys
    taubljy(:ncol) = 0._kind_phys

    ! Set Beljaars 3-D drag profile to zero (stub implementation)
    dragblj(:ncol,:pver) = 0._kind_phys

    ! Set do_beljaars flag to false
    do_beljaars = .false.

  end subroutine tms_beljaars_zero_stub_run

  ! Add turbulent mountain stress to the total surface drag coefficient
!> \section arg_table_turbulent_mountain_stress_add_drag_coefficient_run Argument Table
!! \htmlinclude arg_table_turbulent_mountain_stress_add_drag_coefficient_run.html
  subroutine turbulent_mountain_stress_add_drag_coefficient_run( &
    ncol, pver, &
    ksrftms, &
    ksrf, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: ksrftms(:)       ! Surface drag coefficient for turbulent mountain stress. > 0. [kg m-2 s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: ksrf(:)          ! total surface drag coefficient [kg m-2 s-1]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    ksrf(:ncol) = ksrf(:ncol) + ksrftms(:ncol)

  end subroutine turbulent_mountain_stress_add_drag_coefficient_run

  ! Add Beljaars drag to the wind damping rate for vertical diffusion.
  ! Has to run after vertical_diffusion_wind_damping_rate
!> \section arg_table_beljaars_add_wind_damping_rate_run Argument Table
!! \htmlinclude arg_table_beljaars_add_wind_damping_rate_run.html
  subroutine beljaars_add_wind_damping_rate_run( &
    ncol, pver, &
    dragblj, &
    tau_damp_rate, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: dragblj(:, :)    ! Drag profile from Beljaars SGO form drag > 0. [s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: tau_damp_rate(:,:) ! Rate at which external (surface) stress damps wind speeds [s-1]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg  ! error message
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

    use coords_1d, only: Coords1D

    ! Input arguments
    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver
    real(kind_phys), intent(in)       :: gravit
    type(Coords1D),  intent(in)       :: p                ! Pressure coordinates [Pa]
    logical,         intent(in)       :: do_iss           ! Flag for implicit surface stress (namelist from vdiff)
    logical,         intent(in)       :: itaures          ! Flag for updating tauresx tauresy in this subroutine.
    real(kind_phys), intent(in)       :: dragblj(:, :)    ! Drag profile from Beljaars SGO form drag > 0. [s-1]
    real(kind_phys), intent(in)       :: u1(:,:)          ! After vertical diffusion u-wind [m s-1]
    real(kind_phys), intent(in)       :: v1(:,:)          ! After vertical diffusion v-wind [m s-1]

    ! Input/output arguments
    real(kind_phys), intent(inout)    :: tauresx(:)       ! Partially updated residual surface stress
    real(kind_phys), intent(inout)    :: tauresy(:)
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    integer :: i, k
    real(kind_phys) :: taubljx(ncol)             ! recomputed explicit/residual beljaars stress
    real(kind_phys) :: taubljy(ncol)             ! recomputed explicit/residual beljaars stress

    errmsg = ''
    errflg = 0

    if(.not. do_iss) then
      ! Not added to residual stress if implicit surface stress is off
      return
    endif

    if(.not. itaures) then
      ! Not added to residual stress if residual stress is not requested to be updated.
      return
    endif

    do i = 1, ncol
      ! Add vertically-integrated Beljaars drag to residual stress.
      ! note: in vertical_diffusion_tend, tautotx has taubljx added to it
      ! but this is computed separately using UPDATED u, v (and is not written back to pbuf)
      ! there is a FIXME in beljaars_drag that notes maybe updated u, v could be used there
      ! but that is too early (before vdiff). Keeping this for now. hplin 5/22/25
      taubljx(i) = 0._kind_phys
      taubljy(i) = 0._kind_phys
      do k = 1, pver
        taubljx(i) = taubljx(i) + (1._kind_phys/gravit)*dragblj(i, k)*u1(i, k)*p%del(i, k)
        taubljy(i) = taubljy(i) + (1._kind_phys/gravit)*dragblj(i, k)*v1(i, k)*p%del(i, k)
      end do
    enddo

    tauresx(:ncol) = tauresx(:ncol) + taubljx(:ncol)
    tauresy(:ncol) = tauresy(:ncol) + taubljy(:ncol)

  end subroutine beljaars_add_updated_residual_stress_run

  ! Add TMS using updated provisional winds to total surface stress in
  ! horizontal momentum diffusion - kinetic energy dissipation.
!> \section arg_table_turbulent_mountain_stress_add_updated_surface_stress_run Argument Table
!! \htmlinclude arg_table_turbulent_mountain_stress_add_updated_surface_stress_run.html
  subroutine turbulent_mountain_stress_add_updated_surface_stress_run( &
    ncol, pver, &
    do_iss, itaures, &
    ksrftms, &
    u1, v1, & ! provisionally updated winds
    ! input/output
    tauresx, tauresy, &
    tautotx, tautoty, &
    ! output
    tautmsx, tautmsy, &
    errmsg, errflg)

    ! Input Arguments
    integer,         intent(in)       :: ncol             ! Number of atmospheric columns
    integer,         intent(in)       :: pver
    logical,         intent(in)       :: do_iss           ! Flag for implicit surface stress (namelist from vdiff)
    logical,         intent(in)       :: itaures          ! Flag for updating tauresx tauresy in this subroutine.
    real(kind_phys), intent(in)       :: ksrftms(:)       ! Surface drag coefficient for turbulent mountain stress. > 0. [kg m-2 s-1]
    real(kind_phys), intent(in)       :: u1(:,:)          ! After vertical diffusion u-wind [m s-1]
    real(kind_phys), intent(in)       :: v1(:,:)          ! After vertical diffusion v-wind [m s-1]

    ! Input/output arguments
    real(kind_phys), intent(inout)    :: tauresx(:)       ! Partially updated residual surface stress
    real(kind_phys), intent(inout)    :: tauresy(:)
    real(kind_phys), intent(inout)    :: tautotx(:)       ! Total surface stress (zonal)
    real(kind_phys), intent(inout)    :: tautoty(:)       ! Total surface stress (meridional)

    ! Output arguments
    real(kind_phys), intent(out)      :: tautmsx(:)       ! Implicit zonal turbulent mountain surface stress [N m-2]
    real(kind_phys), intent(out)      :: tautmsy(:)       ! Implicit meridional turbulent mountain surface stress [N m-2]
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    ! Compute the implicit 'tms' using the updated winds.
    ! Below 'tautmsx(i),tautmsy(i)' are pure implicit mountain stresses
    ! that has been actually added into the atmosphere both for explicit
    ! and implicit approach.
    tautmsx(:ncol) = -ksrftms(:ncol)*u1(:ncol, pver)
    tautmsy(:ncol) = -ksrftms(:ncol)*v1(:ncol, pver)

    if(do_iss) then
      ! Implicit surface stress:
      if(.not. itaures) then
        ! Not added to residual stress if residual stress is not requested to be updated.
        return
      endif

      tauresx(:ncol) = tauresx(:ncol) + tautmsx(:ncol)
      tauresy(:ncol) = tauresy(:ncol) + tautmsy(:ncol)

      ! leave tautotx, tautoty unchanged (they are not updated by TMS if do_iss.)
    else ! .not. do_iss
      ! Not implicit surface stress. Add to total surface stress
      tautotx(:ncol) = tautotx(:ncol) + tautmsx(:ncol)
      tautoty(:ncol) = tautoty(:ncol) + tautmsy(:ncol)
    endif
  end subroutine turbulent_mountain_stress_add_updated_surface_stress_run

  ! Set rairv use flag for vertical_diffusion_interpolate_to_interfaces
  ! to be false since we are not in WACCM-X mode.
!> \section arg_table_vertical_diffusion_not_use_rairv_init Argument Table
!! \htmlinclude vertical_diffusion_not_use_rairv_init.html
  subroutine vertical_diffusion_not_use_rairv_init( &
    use_rairv, &
    errmsg, errflg)

    ! Output arguments
    logical,            intent(out) :: use_rairv ! Flag for constituent-dependent gas constant [flag]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Set to .false. since we are not in WACCM-X mode
    ! Only use constituent-dependent gas constant when in WACCM-X mode
    use_rairv = .false.

  end subroutine vertical_diffusion_not_use_rairv_init


end module diffusion_stubs
