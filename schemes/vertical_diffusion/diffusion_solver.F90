! Vertical diffusion module.
! Solves vertical diffusion equations using a tri-diagonal solver
! The module will also apply countergradient fluxes
! Includes interstitials for applying vertical diffusion tendencies and
! other intermediate calculations formerly in vertical_diffusion_tend.
!
! This is a reduced functionality version;
!  molecular diffusion is unsupported until WACCM-X is ported.
!
! Original authors: B. Boville and others, 1991-2004
! Modularized: J. McCaa, September 2004.
! Updates: Sungsu Park, August 2006 - January 2010
! Reduced/CCPPized: Haipeng Lin, May 2025
module diffusion_solver
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! CCPP-compliant public interfaces

  public :: vertical_diffusion_set_temperature_at_toa_default_run
  public :: vertical_diffusion_interpolate_to_interfaces_run   ! Interpolate t, rho, rair to interfaces
  public :: implicit_surface_stress_add_drag_coefficient_run   ! Add implicit turbulent surface stress (if do_iss)
  public :: vertical_diffusion_wind_damping_rate_run           ! Compute wind damping rate from drag coeff
  ! Beljaars will add wind damping rates at this point

  public :: vertical_diffusion_diffuse_horizontal_momentum_timestep_init
  public :: vertical_diffusion_diffuse_horizontal_momentum_run ! Vertical diffusion solver (horizontal momentum)

  public :: vertical_diffusion_set_dry_static_energy_at_toa_zero_run
  public :: vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run
  public :: vertical_diffusion_diffuse_dry_static_energy_run   ! Vertical diffusion solver (dry static energy)

  public :: vertical_diffusion_diffuse_tracers_init
  public :: vertical_diffusion_diffuse_tracers_run             ! Vertical diffusion solver (tracers)

  public :: vertical_diffusion_tendencies_run                  ! Using output from solver compute phys. tend

  ! Parameters for implicit surface stress treatment
  real(kind_phys), parameter :: wsmin = 1._kind_phys         ! Minimum sfc wind speed for estimating frictional
                                                             ! transfer velocity ksrf. [m s-1]
  real(kind_phys), parameter :: ksrfmin = 1.e-4_kind_phys    ! Minimum surface drag coefficient [kg s-1 m-2]
  real(kind_phys), parameter :: timeres = 7200._kind_phys    ! Relaxation time scale of residual stress (>= dt) [s]

contains

  ! Set (default) temperature used at top of atmosphere for interpolation.
  ! In the future, this can use the upper boundary condition provided for temperature,
  ! or a different extrapolation in WACCM-X mode.
!> \section arg_table_vertical_diffusion_set_temperature_at_toa_default_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_set_temperature_at_toa_default_run.html
  pure subroutine vertical_diffusion_set_temperature_at_toa_default_run( &
    ncol, pver, &
    t, &
    t_toai, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: t(:, :)     ! Temperature [K]

    ! Output arguments
    real(kind_phys),    intent(out) :: t_toai(:)   ! Temperature at interface above TOA for vertical diffusion [K]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Default
    t_toai(:ncol) = t(:ncol, 1)  ! Use temperature at top model level [K]

    ! Eventually, if ubc_fixed_temp is available, use ubc_t.

    ! Eventually, in WACCM-X, use extrapolated version:
    ! if (waccmx_mode) then
    !    ! For WACCM-X, set ubc temperature to extrapolate from next two lower interface level temperatures
    !    ! the original formulation is:
    !    ! t_toai(:ncol) = 1.5_r8*tint(:ncol,2)-.5_r8*tint(:ncol,3)
    !    ! this appears to be:
    !    !               = tint(:ncol,2) + 0.5_r8*(tint(:ncol,2) - tint(:ncol,3))
    !    ! assuming that the extrapolated gradient is 1/2 of the temperature gradient between the lower interfaces
    !    ! because the interpolation will be done later in the CCPP-ized scheme, formulate this in terms of
    !    ! the temperature (at midpoints):
    !    t_toai(:ncol) = 1.5_r8*(state%t(:ncol,2)+state%t(:ncol,1))/2._r8-.5_r8*(state%t(:ncol,3)+state%t(:ncol,2))/2._r8

  end subroutine vertical_diffusion_set_temperature_at_toa_default_run

  ! Interpolates temperature, air density (moist and dry),
  ! and sets gas constant (not constituent dependent in non-WACCM-X mode)
  ! at interfaces for vertical diffusion
!> \section arg_table_vertical_diffusion_interpolate_to_interfaces_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_interpolate_to_interfaces_run.html
  pure subroutine vertical_diffusion_interpolate_to_interfaces_run( &
    ncol, pver, pverp, &
    gravit, &
    rair, rairv, &
    flag_for_constituent_dependent_gas_constant, &
    t, &
    t_toai, &
    pint, pintdry, &
    ! below output
    ti, rairi, rhoi, rhoi_dry, &
    dpidz_sq, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol       ! Number of columns
    integer,            intent(in)  :: pver       ! Number of vertical layers
    integer,            intent(in)  :: pverp      ! Number of vertical interfaces (pver+1)
    real(kind_phys),    intent(in)  :: gravit
    real(kind_phys),    intent(in)  :: rair       ! Gas constant for dry air [J kg-1 K-1]
    real(kind_phys),    intent(in)  :: rairv(:,:) ! composition_dependent_gas_constant_of_dry_air [J kg-1 K-1]
    logical,            intent(in)  :: flag_for_constituent_dependent_gas_constant
    real(kind_phys),    intent(in)  :: t(:,:)     ! Air temperature at midpoints [K]
    real(kind_phys),    intent(in)  :: t_toai(:)  ! Air temperature to use at interface above TOA [K]
    real(kind_phys),    intent(in)  :: pint(:,:)  ! Air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)  :: pintdry(:,:) ! Dry air pressure at interfaces [Pa]

    ! Output variables
    real(kind_phys),    intent(out) :: ti(:,:)    ! Air temperature at interfaces [K]
    real(kind_phys),    intent(out) :: rairi(:,:) ! Gas constant for dry air at interfaces [J kg-1 K-1]
    real(kind_phys),    intent(out) :: rhoi(:,:)  ! Air density at interfaces [kg m-3]
    real(kind_phys),    intent(out) :: rhoi_dry(:,:) ! Dry air density at interfaces [kg m-3]
    real(kind_phys),    intent(out) :: dpidz_sq(:,:) ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4]
    character(len=512), intent(out) :: errmsg     ! Error message
    integer,            intent(out) :: errflg     ! Error flag

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    ! Interpolate temperature to interfaces
    ti(:ncol,1) = t_toai(:ncol)
    do k = 2, pver
       do i = 1, ncol
          ti(i,k) = 0.5_kind_phys * (t(i,k) + t(i,k-1))
       end do
    end do
    ! Assume temperature at surface interface is the same as the midpoint value in the lowest layer:
    ti(:ncol,pverp) = t(:ncol,pver)

    ! Set gas constant at interfaces
    if(flag_for_constituent_dependent_gas_constant) then
      rairi(:ncol,1) = rairv(:ncol,1)
      do k = 2, pver
        do i = 1, ncol
          rairi(i,k) = 0.5_kind_phys * (rairv(i,k)+rairv(i,k-1))
        end do
      end do

      ! Assume gas constant at surface interface is the same as the midpoint value in the lowest layer:
      rairi(:ncol,pverp) = rairv(:ncol,pver)
    else
      rairi(:ncol,:pverp) = rair
    endif

    ! Compute rho at interfaces.
    do k = 1, pverp
       do i = 1, ncol
          rhoi(i,k)  = pint(i,k) / (rairi(i,k)*ti(i,k))
       end do
    end do

    ! Compute rho_dry at interfaces.
    do k = 1, pverp
       do i = 1, ncol
          rhoi_dry(i,k)  = pintdry(i,k) / (rairi(i,k)*ti(i,k))
       end do
    end do

    ! Compute square of derivative of pressure with height at interfaces.
    ! Note that this equation assumes hydrostatic balance, i.e. the *derivative* dp/dz is g*rho
    dpidz_sq = gravit*rhoi(:ncol, :)
    dpidz_sq = dpidz_sq*dpidz_sq

    ! FIXME: The following two lines are kept in only to preserve answers;
    !        they really should be taken out completely.
    dpidz_sq(:, 1) = gravit*(pint(:, 1)/(rairv(:ncol, 1)*t(:ncol, 1)))
    dpidz_sq(:, 1) = dpidz_sq(:, 1)*dpidz_sq(:, 1)

  end subroutine vertical_diffusion_interpolate_to_interfaces_run

  ! Add turbulent surface stress to total surface drag coefficient, if implicit
!> \section arg_table_implicit_surface_stress_add_drag_coefficient_run Argument Table
!! \htmlinclude arg_table_implicit_surface_stress_add_drag_coefficient_run.html
  pure subroutine implicit_surface_stress_add_drag_coefficient_run( &
             ncol, pver, &
             do_iss, &
             taux, tauy, &
             u0, v0, &
             ! below output
             ksrf, &
             errmsg, errflg)

    ! Input arguments
    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver
    logical,         intent(in)       :: do_iss             ! Use implicit turbulent surface stress computation
    real(kind_phys), intent(in)       :: taux(:)            ! Surface zonal stress [N m-2]
                                                            ! Input u-momentum per unit time per unit area into the atmosphere.
    real(kind_phys), intent(in)       :: tauy(:)            ! Surface meridional stress [N m-2]
                                                            ! Input v-momentum per unit time per unit area into the atmosphere.

    real(kind_phys), intent(in)       :: u0(:,:)            ! Input u-wind [m s-1]
    real(kind_phys), intent(in)       :: v0(:,:)            ! Input v-wind [m s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: ksrf(:)            ! total surface drag coefficient [kg m-2 s-1]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    ! Local variables
    integer :: i
    real(kind_phys) :: ws(ncol)                  ! Lowest-level wind speed [m s-1]
    real(kind_phys) :: ksrfturb(ncol)            ! Surface drag coefficient of 'normal' stress. > 0.
                                                 ! Virtual mass input per unit time per unit area [kg m-2 s-1]
    real(kind_phys) :: tau(ncol)                 ! Turbulent surface stress (not including mountain stress)

    errmsg = ''
    errflg = 0

    ! Do 'Implicit Surface Stress' treatment for numerical stability
    ! in the lowest model layer.
    if (do_iss) then
      ! Add surface drag coefficient for implicit diffusion
      do i = 1, ncol
        ws(i) = max(sqrt(u0(i, pver)**2._kind_phys + v0(i, pver)**2._kind_phys), wsmin)
        tau(i) = sqrt(taux(i)**2._kind_phys + tauy(i)**2._kind_phys)
        ksrfturb(i) = max(tau(i)/ws(i), ksrfmin)
      end do
      ksrf(:ncol) = ksrf(:ncol) + ksrfturb(:ncol)  ! Do all surface stress implicitly
    endif

  end subroutine implicit_surface_stress_add_drag_coefficient_run

  ! Computes rate at which wind is exponentially damped by surface stress
  ! based on total surface stress (after all surface stress is added)
!> \section arg_table_vertical_diffusion_wind_damping_rate_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_wind_damping_rate_run.html
  pure subroutine vertical_diffusion_wind_damping_rate_run( &
    ncol, pver, &
    gravit, &
    p, &
    ksrf, &
    tau_damp_rate, &
    errmsg, errflg)

    use coords_1d, only: Coords1D

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: gravit
    type(coords1d),     intent(in)    :: p                  ! Pressure coordinates [Pa]
    real(kind_phys),    intent(in)    :: ksrf(:)            ! total surface drag coefficient [kg m-2 s-1]

    ! Output arguments
    real(kind_phys),    intent(out)   :: tau_damp_rate(:,:) ! Rate at which external (surface) stress damps wind speeds [s-1]
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    ! In most layers, no damping at all.
    tau_damp_rate(:ncol,:pver) = 0._kind_phys

    ! Compute damp rate based on total surface drag coeff
    !
    ! Physical interpretation:
    ! ksrf is stress per unit wind speed.
    ! p%del / gravit is approximately the mass in the layer per unit of
    ! surface area.
    ! Therefore, gravit*ksrf/p%del is the acceleration of wind per unit
    ! wind speed, i.e. the rate at which wind is exponentially damped by
    ! surface stress.
    tau_damp_rate(:ncol, pver) = -gravit*ksrf(:ncol)*p%rdel(:ncol, pver)

  end subroutine vertical_diffusion_wind_damping_rate_run

  ! Zero surface stresses at timestep initialization
  ! Any surface stresses, e.g., Turbulent Mountain Stress
  ! will be added to this later as they are computed.
!> \section arg_table_vertical_diffusion_diffuse_horizontal_momentum_timestep_init Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_horizontal_momentum_timestep_init.html
  pure subroutine vertical_diffusion_diffuse_horizontal_momentum_timestep_init( &
    ncol, pver, &
    ksrf, &
    tau_damp_rate, &
    errmsg, errflg)

    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver

    ! Output (initial zeroed) arguments
    real(kind_phys),    intent(out)   :: ksrf(:)            ! total surface drag coefficient [kg m-2 s-1]
    real(kind_phys),    intent(out)   :: tau_damp_rate(:,:) ! Rate at which external (surface) stress damps wind speeds [s-1]

    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    ksrf(:ncol) = 0._kind_phys

    ! In most layers, no damping at all.
    tau_damp_rate(:ncol,:pver) = 0._kind_phys

  end subroutine vertical_diffusion_diffuse_horizontal_momentum_timestep_init

  ! Driver routine to compute vertical diffusion of horizontal momentum,
  ! including the loss of momentum due to frictional (KE to heat) dissipation.
  ! Turbulent diffusivities and boundary layer nonlocal transport terms are
  ! obtained from the turbulence module.
!> \section arg_table_vertical_diffusion_diffuse_horizontal_momentum_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_horizontal_momentum_run.html
  subroutine vertical_diffusion_diffuse_horizontal_momentum_run( &
             ncol, pver, pverp, &
             dt, &
             rair, gravit, &
             do_iss, &
             am_correction, &
             itaures, &
             do_beljaars, &
             p, t, rhoi, &
             taux, tauy, &
             tau_damp_rate, &
             ksrftms, &
             dragblj, &
             kvm, &
             dpidz_sq, &
             u0, v0, dse0, &
             tauresx, tauresy, &
             dtk, &
             tautmsx, tautmsy, &
             u1, v1, dse1, &
             errmsg, errflg)

    use coords_1d, only: Coords1D
    use vertical_diffusion_solver, only: fin_vol_solve

    ! Input Arguments
    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver
    integer,         intent(in)       :: pverp
    real(kind_phys), intent(in)       :: dt               ! physics timestep [s]
    real(kind_phys), intent(in)       :: rair
    real(kind_phys), intent(in)       :: gravit
    logical,         intent(in)       :: do_iss           ! Use implicit turbulent surface stress computation
    logical,         intent(in)       :: am_correction    ! Do angular momentum conservation correction
    logical,         intent(in)       :: itaures          ! Flag for updating tauresx tauresy in this subroutine.
    logical,         intent(in)       :: do_beljaars      ! Flag indicating Beljaars drag
    type(coords1d),  intent(in)       :: p                ! Pressure coordinates [Pa]
    real(kind_phys), intent(in)       :: t(:, :)          ! Temperature [K]
    real(kind_phys), intent(in)       :: rhoi(:, :)       ! Density of air at interfaces [kg m-3]
    real(kind_phys), intent(in)       :: taux(:)          ! Surface zonal stress [N m-2]
                                                          ! Input u-momentum per unit time per unit area into the atmosphere.
    real(kind_phys), intent(in)       :: tauy(:)          ! Surface meridional stress [N m-2]
                                                          ! Input v-momentum per unit time per unit area into the atmosphere.
    real(kind_phys), intent(in)       :: tau_damp_rate(:,:) ! Rate at which external (surface) stress damps wind speeds [s-1]
    real(kind_phys), intent(in)       :: ksrftms(:)       ! Surface drag coefficient for turbulent mountain stress. > 0. [kg m-2 s-1]
    real(kind_phys), intent(in)       :: dragblj(:, :)    ! Drag profile from Beljaars SGO form drag > 0. [s-1]
    real(kind_phys), intent(in)       :: kvm(:, :)        ! Eddy viscosity (Eddy diffusivity for momentum) [m^2 s-1], interfaces
    real(kind_phys), intent(in)       :: dpidz_sq(:,:)    ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4], interfaces

    real(kind_phys), intent(in)       :: u0(:,:)          ! Input u-wind [m s-1]
    real(kind_phys), intent(in)       :: v0(:,:)          ! Input v-wind [m s-1]
    real(kind_phys), intent(in)       :: dse0(:,:)        ! Input Dry static energy [J kg-1]

    ! Input:  Reserved surface stress at previous time step
    ! Output: Reserved surface stress at current  time step
    real(kind_phys), intent(inout)    :: tauresx(:)
    real(kind_phys), intent(inout)    :: tauresy(:)

    ! Output Arguments
    real(kind_phys), intent(out)      :: dtk(:, :)        ! T tendency from KE dissipation [J kg-1]

    ! Recomputed TMS stresses using provisionally updated winds:
    real(kind_phys), intent(out)      :: tautmsx(:)       ! Implicit zonal turbulent mountain surface stress [N m-2]
    real(kind_phys), intent(out)      :: tautmsy(:)       ! Implicit meridional turbulent mountain surface stress [N m-2]

    ! Outputs from vertical diffusion will be converted to physics tendencies in separate scheme
    ! These provisional outputs are used to compute diagnostics, and thus need to be kept in non-tend form
    real(kind_phys), intent(out)      :: u1(:,:)          ! After vertical diffusion u-wind [m s-1]
    real(kind_phys), intent(out)      :: v1(:,:)          ! After vertical diffusion v-wind [m s-1]
    real(kind_phys), intent(out)      :: dse1(:,:)        ! After vertical diffusion Dry static energy [J kg-1]

    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    ! Local variables
    integer  :: i, k

    real(kind_phys) :: tmp1(ncol)                ! Temporary storage
    real(kind_phys) :: tmpi1(ncol, pverp)        ! Interface KE dissipation
    real(kind_phys) :: tmpi2(ncol, pverp)        ! dt*(g*rho)**2/dp at interfaces
    real(kind_phys) :: ramda                     ! dt/timeres [unitless]

    real(kind_phys) :: keg_in(ncol, pver)        ! KE on entry to subroutine
    real(kind_phys) :: keg_out(ncol, pver)       ! KE after U and V dissipation/diffusion

    real(kind_phys) :: tautotx(ncol)             ! Total surface stress (zonal)
    real(kind_phys) :: tautoty(ncol)             ! Total surface stress (meridional)

    real(kind_phys) :: dinp_u(ncol, pverp)       ! Vertical difference at interfaces, input u
    real(kind_phys) :: dinp_v(ncol, pverp)       ! Vertical difference at interfaces, input v
    real(kind_phys) :: dout_u                    ! Vertical difference at interfaces, output u
    real(kind_phys) :: dout_v                    ! Vertical difference at interfaces, output v

    real(kind_phys) :: usum_in(ncol)             ! Vertical integral of input u-momentum. Total zonal
                                                 ! momentum per unit area in column (sum of u*dp/g) [kg m s-1 m-2]
    real(kind_phys) :: vsum_in(ncol)             ! Vertical integral of input v-momentum. Total meridional
                                                 ! momentum per unit area in column (sum of v*dp/g) [kg m s-1 m-2]
    real(kind_phys) :: usum_out(ncol)            ! Vertical integral of u-momentum after doing implicit diffusion
    real(kind_phys) :: vsum_out(ncol)            ! Vertical integral of v-momentum after doing implicit diffusion

    real(kind_phys) :: tauimpx(ncol)             ! Actual net stress added at the current step other than mountain stress
    real(kind_phys) :: tauimpy(ncol)             ! Actual net stress added at the current step other than mountain stress
    real(kind_phys) :: taubljx(ncol)             ! recomputed explicit/residual beljaars stress
    real(kind_phys) :: taubljy(ncol)             ! recomputed explicit/residual beljaars stress

    errmsg = ''
    errflg = 0

    ! Set initial iterative outputs to input values
    u1(:ncol,:pver) = u0(:ncol,:pver)
    v1(:ncol,:pver) = v0(:ncol,:pver)
    dse1(:ncol,:pver) = dse0(:ncol,:pver)

    ! necessary temporaries used in computation
    tmp1(:ncol) = dt*gravit*p%rdel(:, pver)

    ! second term here uses dpidz_sq(:,1) as defined because
    ! the actual dpidz_sq term is "kept in only to preserve answers"
    ! but it is actually not used here ...
    tmpi2(:ncol, 1) = dt*(gravit*rhoi(:ncol,1)*gravit*rhoi(:ncol,1))/(p%mid(:, 1) - p%ifc(:, 1))
    tmpi2(:ncol, 2:pver) = dt*dpidz_sq(:, 2:pver)*p%rdst

    !---------------------------- !
    ! Diffuse Horizontal Momentum !
    !---------------------------- !
    ! Compute the vertical upward differences of the input u,v for KE dissipation
    ! at each interface.
    ! Velocity = 0 at surface, so difference at the bottom interface is -u,v(pver)
    ! These 'dinp_u, dinp_v' are computed using the non-diffused input wind.

    do i = 1, ncol
      dinp_u(i, 1) = 0._kind_phys
      dinp_v(i, 1) = 0._kind_phys
      dinp_u(i, pver + 1) = -u0(i, pver)
      dinp_v(i, pver + 1) = -v0(i, pver)
    end do
    do k = 2, pver
      do i = 1, ncol
        dinp_u(i, k) = u0(i, k) - u0(i, k - 1)
        dinp_v(i, k) = v0(i, k) - v0(i, k - 1)
      end do
    end do

    ! -------------------------------------------------------------- !
    ! Do 'Implicit Surface Stress' treatment for numerical stability !
    ! in the lowest model layer.                                     !
    ! -------------------------------------------------------------- !
    if (do_iss) then
      ! Vertical integration of input momentum.
      ! This is total horizontal momentum per unit area [ kg*m/s/m2 ] in each column.
      ! Note (u,v) are the raw input to the PBL scheme, not the
      ! provisionally-marched ones within the iteration loop of the PBL scheme.

      do i = 1, ncol
        usum_in(i) = 0._kind_phys
        vsum_in(i) = 0._kind_phys
        do k = 1, pver
          usum_in(i) = usum_in(i) + (1._kind_phys/gravit)*u0(i, k)*p%del(i, k)
          vsum_in(i) = vsum_in(i) + (1._kind_phys/gravit)*v0(i, k)*p%del(i, k)
        end do
      end do

      ! Add residual stress of previous time step explicitly into the lowest
      ! model layer with a relaxation time scale of 'timeres'.

      if (am_correction) then
        ! preserve time-mean torque
        ramda = 1._kind_phys
      else
        ramda = dt/timeres
      end if

      u1(:ncol, pver) = u1(:ncol, pver) + tmp1(:ncol)*tauresx(:ncol)*ramda
      v1(:ncol, pver) = v1(:ncol, pver) + tmp1(:ncol)*tauresy(:ncol)*ramda

    else ! .not. do_iss
      ! In this case, do 'turbulent mountain stress' implicitly,
      ! but do 'normal turbulent stress' explicitly.
      ! In this case, there is no 'residual stress' as long as 'tms' is
      ! treated in a fully implicit way, which is true.

      ! Do 'normal stress' explicitly
      u1(:ncol, pver) = u1(:ncol, pver) + tmp1(:ncol)*taux(:ncol)
      v1(:ncol, pver) = v1(:ncol, pver) + tmp1(:ncol)*tauy(:ncol)
    end if  ! End of 'do iss' (implicit surface stress)

    ! --------------------------------------------------------------------------------------- !
    ! Diffuse horizontal momentum implicitly using tri-diagnonal matrix.                      !
    ! The 'u,v' are input-output: the output 'u,v' are implicitly diffused winds.             !
    !    For implicit 'normal' stress : ksrf = ksrftms + ksrfturb,                            !
    !                                   u(pver) : explicitly include 'residual normal' stress !
    !    For explicit 'normal' stress : ksrf = ksrftms                                        !
    !                                   u(pver) : explicitly include 'normal' stress          !
    ! Note that in all the two cases above, 'tms' is fully implicitly treated.                !
    ! --------------------------------------------------------------------------------------- !

    v1(:ncol, :) = fin_vol_solve(dt, p, v1(:ncol, :), ncol, pver, &
                                 coef_q=tau_damp_rate(:ncol,:pver), &
                                 coef_q_diff=kvm(:ncol, :)*dpidz_sq(:ncol, :))

    u1(:ncol, :) = fin_vol_solve(dt, p, u1(:ncol, :), ncol, pver, &
                                 coef_q=tau_damp_rate(:ncol,:pver), &
                                 coef_q_diff=kvm(:ncol, :)*dpidz_sq(:ncol, :))

    ! ---------------------------------------------------------------------- !
    ! Calculate 'total' ( tautotx ) and 'tms' ( tautmsx ) stresses that      !
    ! have been actually added into the atmosphere at the current time step. !
    ! Also, update residual stress, if required.                             !
    ! ---------------------------------------------------------------------- !

    do i = 1, ncol
      ! Compute the implicit 'tms' using the updated winds.
      ! Below 'tautmsx(i),tautmsy(i)' are pure implicit mountain stresses
      ! that has been actually added into the atmosphere both for explicit
      ! and implicit approach.

      tautmsx(i) = -ksrftms(i)*u1(i, pver)
      tautmsy(i) = -ksrftms(i)*v1(i, pver)

      ! We want to add vertically-integrated Beljaars drag to residual stress.
      ! So this has to be calculated locally.
      ! We may want to rethink the residual drag calculation performed here on. (jtb)
      taubljx(i) = 0._kind_phys
      taubljy(i) = 0._kind_phys
      do k = 1, pver
        taubljx(i) = taubljx(i) + (1._kind_phys/gravit)*dragblj(i, k)*u1(i, k)*p%del(i, k)
        taubljy(i) = taubljy(i) + (1._kind_phys/gravit)*dragblj(i, k)*v1(i, k)*p%del(i, k)
      end do

      if (do_iss) then
        ! Compute vertical integration of final horizontal momentum
        usum_out(i) = 0._kind_phys
        vsum_out(i) = 0._kind_phys
        do k = 1, pver
          usum_out(i) = usum_out(i) + (1._kind_phys/gravit)*u1(i, k)*p%del(i, k)
          vsum_out(i) = vsum_out(i) + (1._kind_phys/gravit)*v1(i, k)*p%del(i, k)
        end do

        ! Compute net stress added into the atmosphere at the current time step.
        ! Note that the difference between 'usum_in' and 'usum_out' are induced
        ! by 'explicit residual stress + implicit total stress' for implicit case, while
        ! by 'explicit normal   stress + implicit tms   stress' for explicit case.
        ! Here, 'tautotx(i)' is net stress added into the air at the current time step.
        tauimpx(i) = (usum_out(i) - usum_in(i))/dt
        tauimpy(i) = (vsum_out(i) - vsum_in(i))/dt

        tautotx(i) = tauimpx(i)
        tautoty(i) = tauimpy(i)

        ! Compute residual stress and update if required.
        ! Note that the total stress we should have added at the current step is
        ! the sum of 'taux(i) - ksrftms(i)*u(i,pver) + tauresx(i)'.
        if (itaures) then
          tauresx(i) = taux(i) + tauresx(i) - tauimpx(i) + tautmsx(i) + taubljx(i)
          tauresy(i) = tauy(i) + tauresy(i) - tauimpy(i) + tautmsy(i) + taubljy(i)
        end if
      else
        ! .not. do_iss
        tautotx(i) = taux(i) + tautmsx(i)
        tautoty(i) = tauy(i) + tautmsy(i)
        tauresx(i) = 0._kind_phys
        tauresy(i) = 0._kind_phys
      end if  ! End of 'do_iss' if
    end do ! End of 'do i = 1, ncol' loop

    ! ------------------------------------ !
    ! Calculate kinetic energy dissipation !
    ! ------------------------------------ !

    ! Modification : In future, this should be set exactly same as
    !                the ones in the convection schemes

    if (do_beljaars) then
      ! 2. Add Kinetic Energy change across dissipation to Static Energy
      do k = 1, pver
        do i = 1, ncol
          keg_in(i, k) = 0.5_kind_phys*(u0(i, k)*u0(i, k) + v0(i, k)*v0(i, k))
          keg_out(i, k) = 0.5_kind_phys*(u1(i, k)*u1(i, k) + v1(i, k)*v1(i, k))
        end do
      end do

      do k = 1, pver
        do i = 1, ncol
          dtk(i, k) = keg_in(i, k) - keg_out(i, k)
          dse1(i, k) = dse1(i, k) + dtk(i, k) ! + dkeblj(i,k)
        end do
      end do
    else
      ! .not. do_beljaars

      ! 1. Compute dissipation term at interfaces
      !    Note that 'u,v' are already diffused wind, and 'tautotx,tautoty' are
      !    implicit stress that has been actually added. On the other hand,
      !    'dinp_u, dinp_v' were computed using non-diffused input wind.

      ! Modification : I should check whether non-consistency between 'u' and 'dinp_u'
      !                is correctly intended approach. I think so.
      k = pver + 1
      do i = 1, ncol
        tmpi1(i, 1) = 0._kind_phys
        tmpi1(i, k) = 0.5_kind_phys*dt*gravit* &
                      ((-u1(i, k - 1) + dinp_u(i, k))*tautotx(i) + (-v1(i, k - 1) + dinp_v(i, k))*tautoty(i))
      end do

      do k = 2, pver
        do i = 1, ncol
          dout_u = u1(i, k) - u1(i, k - 1)
          dout_v = v1(i, k) - v1(i, k - 1)
          tmpi1(i, k) = 0.25_kind_phys*tmpi2(i, k)*kvm(i, k)* &
                        (dout_u**2 + dout_v**2 + dout_u*dinp_u(i, k) + dout_v*dinp_v(i, k))
        end do
      end do

      ! 2. Compute dissipation term at midpoints, add to dry static energy
      do k = 1, pver
        do i = 1, ncol
          dtk(i, k) = (tmpi1(i, k + 1) + tmpi1(i, k))*p%rdel(i, k)
          dse1(i, k) = dse1(i, k) + dtk(i, k)
        end do
      end do
    end if

    ! End of kinetic energy dissipation due to horizontal momentum diffusion
  end subroutine vertical_diffusion_diffuse_horizontal_momentum_run

  ! Set dry static energy top boundary condition - zero
  ! used when (1) no molecular diffusion is active, or (2) molecular diffusion is active and not WACCM-X.
  !
  ! Note: this is an artifact to the way the vertical diffusion calculations are done;
  ! it does not represent physical reality.
!> \section arg_table_vertical_diffusion_set_dry_static_energy_at_toa_zero_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_set_dry_static_energy_at_toa_zero_run.html
  pure subroutine vertical_diffusion_set_dry_static_energy_at_toa_zero_run( &
    ncol, &
    dse_top, &
    errmsg, errflg)

    ! Input arguments
    integer,         intent(in)       :: ncol

    ! Output arguments
    real(kind_phys),    intent(out)   :: dse_top(:)       ! Dry static energy top boundary condition.
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    dse_top(:ncol) = 0._kind_phys

  end subroutine vertical_diffusion_set_dry_static_energy_at_toa_zero_run

  ! Set dry static energy top boundary condition - molecular diffusion non-WACCM-X version
!> \section arg_table_vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run.html
  pure subroutine vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run( &
    ncol, &
    gravit, &
    cpairv, &
    zi, &
    tint, &
    dse_top, &
    errmsg, errflg)

    ! Input arguments
    integer,         intent(in)       :: ncol
    real(kind_phys), intent(in)       :: gravit
    real(kind_phys), intent(in)       :: cpairv(:,:)    ! specific heat of dry air at constant pressure [J kg-1 K-1]
    real(kind_phys), intent(in)       :: zi(:,:)        ! geopotential_height_wrt_surface_at_interface [m]
    real(kind_phys), intent(in)       :: tint(:,:)      ! Air temperature at interfaces [K]

    ! Output arguments
    real(kind_phys),    intent(out)   :: dse_top(:)     ! Dry static energy top boundary condition.
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    dse_top(:ncol) = cpairv(:ncol,1) * tint(:ncol,1) + gravit * zi(:ncol,1)

  end subroutine vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run

  ! Driver routine to compute vertical diffusion of dry static energy.
  ! The new temperature is computed from the diffused dry static energy
  ! after tendencies are applied.
  ! Turbulent diffusivities and boundary layer nonlocal transport terms are
  ! obtained from the turbulence module.
!> \section arg_table_vertical_diffusion_diffuse_dry_static_energy_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_dry_static_energy_run.html
  subroutine vertical_diffusion_diffuse_dry_static_energy_run( &
    ncol, pver, &
    dt, &
    gravit, &
    p, rhoi, &
    shflx, &
    dse_top, &
    kvh, cgh, &
    dpidz_sq, &
    ! input/output (may have been updated by horizontal momentum diffusion)
    dse, &
    errmsg, errflg)

    use coords_1d, only: Coords1D
    use linear_1d_operators, only: BoundaryType, BoundaryFixedLayer, &
                                   BoundaryData
    use vertical_diffusion_solver, only: fin_vol_solve

    ! Input Arguments
    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver
    real(kind_phys), intent(in)       :: dt               ! physics timestep [s]

    real(kind_phys), intent(in)       :: gravit
    type(coords1d),  intent(in)       :: p                ! Pressure coordinates [Pa]
    real(kind_phys), intent(in)       :: rhoi(:, :)       ! Density of air at interfaces [kg m-3]
    real(kind_phys), intent(in)       :: shflx(:)         ! Surface sensible heat flux [W m-2]
    real(kind_phys), intent(in)       :: dse_top(:)       ! Dry static energy top boundary condition.
    real(kind_phys), intent(in)       :: kvh(:, :)        ! Eddy diffusivity for heat [m^2 s-1], interfaces
    real(kind_phys), intent(in)       :: cgh(:, :)        ! Counter-gradient term for heat [J kg-1 m-1], interfaces
    real(kind_phys), intent(in)       :: dpidz_sq(:,:)    ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4]

    ! Input/output arguments
    real(kind_phys), intent(inout)    :: dse(:,:)         ! After vertical diffusion Dry static energy [J kg-1]
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    ! Local variables
    integer :: k
    type(BoundaryType) :: interface_boundary     ! Boundary layer objects

    real(kind_phys) :: tmp1(ncol)                ! Temporary storage

    errmsg = ''
    errflg = 0

    ! Boundary condition for a fixed concentration directly on a boundary
    ! interface (i.e. a boundary layer of size 0).
    interface_boundary = BoundaryFixedLayer(spread(0._kind_phys, 1, ncol))

    tmp1(:ncol) = dt*gravit*p%rdel(:, pver)

    ! Modification : In future, we should diffuse the fully conservative
    !                moist static energy, not the dry static energy.

    ! Add counter-gradient to input static energy profiles
    do k = 1, pver
      dse(:ncol, k) = dse(:ncol, k) + dt*p%rdel(:, k)*gravit* &
                      (rhoi(:ncol, k + 1)*kvh(:ncol, k + 1)*cgh(:ncol, k + 1) &
                       - rhoi(:ncol, k)*kvh(:ncol, k)*cgh(:ncol, k))
    end do

    ! Add the explicit surface fluxes to the lowest layer
    dse(:ncol, pver) = dse(:ncol, pver) + tmp1(:ncol)*shflx(:ncol)

    !---------------------------------------------------
    ! Solve for temperature using thermal conductivity
    !---------------------------------------------------

    ! Boundary layer thickness of "0._kind_phys" signifies that the boundary
    ! condition is defined directly on the top interface.
    dse(:ncol, :) = fin_vol_solve(dt, p, dse(:ncol, :), ncol, pver, &
                                  coef_q_diff=kvh(:ncol, :)*dpidz_sq(:ncol, :), &
                                  upper_bndry=interface_boundary, &
                                  l_cond=BoundaryData(dse_top(:ncol)))

  end subroutine vertical_diffusion_diffuse_dry_static_energy_run

  ! Set tracers to diffuse using vertical diffusion.
!> \section arg_table_vertical_diffusion_diffuse_tracers_init Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_tracers_init.html
  pure subroutine vertical_diffusion_diffuse_tracers_init( &
    ncnst, &
    do_diffusion_const, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncnst                 ! # of constituents to diffuse. In eddy_diff, only wv. Others, all advected CCPP constituents.

    ! Output arguments
    logical,            intent(out)   :: do_diffusion_const(:) ! diffuse constituents (size ncnst) [flag]
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    ! Check if logical array is of the expected size
    if (size(do_diffusion_const) .ne. ncnst) then
      errflg = 1
      write(errmsg, *) 'compute_vdiff: do_diffusion_const size ', &
        size(do_diffusion_const), ' is not equal to ncnst ', ncnst
      return
    end if

    ! Diffuse all constituents using moist coordinates.
    do_diffusion_const(:ncnst) = .true.

    ! TODO: in the future (when dropmixnuc is implemented),
    ! handle that if constituent is treated in dropmixnuc
    ! (vertically mixed by ndrop activation process) then
    ! do not handle vertical diffusion in this module.

  end subroutine vertical_diffusion_diffuse_tracers_init

  ! Diffuse tracers (no molecular diffusion)
!> \section arg_table_vertical_diffusion_diffuse_tracers_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_tracers_run.html
  subroutine vertical_diffusion_diffuse_tracers_run( &
    ncol, pver, &
    ncnst, &
    dt, &
    rair, &
    gravit, &
    do_diffusion_const, &
    p, &
    t, &
    rhoi, &
    cflx, &
    kvh, kvq, cgs, &
    qmincg, &
    dpidz_sq, &
    ubc_mmr, &
    cnst_fixed_ubc, &
    q0, &
    ! below output
    q1, &
    errmsg, errflg)

    use coords_1d, only: Coords1D
    use linear_1d_operators, only: TriDiagDecomp
    use vdiff_lu_solver, only: fin_vol_lu_decomp

    ! Input arguments
    integer,         intent(in)       :: ncol             ! Number of atmospheric columns
    integer,         intent(in)       :: pver
    integer,         intent(in)       :: ncnst            ! # of constituents to diffuse. In eddy_diff, only wv. Others, all advected CCPP constituents.
    real(kind_phys), intent(in)       :: dt               ! physics timestep [s]
    real(kind_phys), intent(in)       :: rair
    real(kind_phys), intent(in)       :: gravit
    logical,         intent(in)       :: do_diffusion_const(:) ! diffuse constituents (size ncnst) [flag]

    type(coords1d),  intent(in)       :: p                ! Pressure coordinates [Pa]
    real(kind_phys), intent(in)       :: t(:, :)          ! Temperature [K]
    real(kind_phys), intent(in)       :: rhoi(:, :)       ! Density of air at interfaces [kg m-3]

    real(kind_phys), intent(in)       :: cflx(:, :)       ! Surface constituent flux [kg m-2 s-1] (ncol,ncnst)

    real(kind_phys), intent(in)       :: kvh(:, :)        ! Eddy diffusivity for heat [m^2 s-1], interfaces
    real(kind_phys), intent(in)       :: kvq(:, :)        ! Eddy diffusivity for constituents [m^2 s-1], interfaces
    real(kind_phys), intent(in)       :: cgs(:, :)        ! Counter-gradient star [cg/flux] [s m-2], interfaces

    real(kind_phys), intent(in)       :: qmincg(:)        ! Minimum constituent mixing ratios from cg fluxes (ncnst)
    real(kind_phys), intent(in)       :: dpidz_sq(:,:)    ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4]

    ! Upper boundary properties (for when molecular diffusion is off)
    real(kind_phys), intent(in)       :: ubc_mmr(:, :)    ! Upper boundary mixing ratios [kg kg-1]
    logical,         intent(in)       :: cnst_fixed_ubc(:)! Whether upper boundary condition is fixed

    real(kind_phys), intent(in)       :: q0(:,:,:)        ! Input Constituents [kg kg-1]

    ! Output arguments
    real(kind_phys), intent(out)      :: q1(:,:,:)        ! After vertical diffusion Constituents [kg kg-1]


    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    ! Local variables
    integer :: m, k
    type(TriDiagDecomp) :: no_molec_decomp       ! LU decomposition information.
    real(kind_phys) :: tmp1(ncol)                ! Temporary storage
    logical  :: lqtst(ncol)                      ! Adjust vertical profiles
    real(kind_phys) :: qtm(ncol, pver)           ! Temporary copy of q
    real(kind_phys) :: rrho(ncol)                ! 1./bottom level density

    errmsg = ''
    errflg = 0

    ! Set initial iterative outputs to input values
    q1(:ncol,:pver,:ncnst) = q0(:ncol,:pver,:ncnst)

    ! necessary temporaries used in computation
    rrho(:ncol) = rair*t(:ncol, pver)/p%mid(:, pver)
    tmp1(:ncol) = dt*gravit*p%rdel(:, pver)

    ! Loop through constituents
    no_molec_decomp = fin_vol_lu_decomp(dt, p, &
                                        coef_q_diff=kvq(:ncol, :)*dpidz_sq(:ncol, :))

    do m = 1, ncnst
      if (do_diffusion_const(m)) then
        ! Add the nonlocal transport terms to constituents in the PBL.
        ! Check for neg q's in each constituent and put the original vertical
        ! profile back if a neg value is found. A neg value implies that the
        ! quasi-equilibrium conditions assumed for the countergradient term are
        ! strongly violated.
        qtm(:ncol, :pver) = q1(:ncol, :pver, m)

        do k = 1, pver
          q1(:ncol, k, m) = q1(:ncol, k, m) + &
                            dt*p%rdel(:, k)*gravit*(cflx(:ncol, m)*rrho(:ncol))* &
                            (rhoi(:ncol, k + 1)*kvh(:ncol, k + 1)*cgs(:ncol, k + 1) &
                             - rhoi(:ncol, k)*kvh(:ncol, k)*cgs(:ncol, k))
        end do
        lqtst(:ncol) = all(q1(:ncol, 1:pver, m) >= qmincg(m), 2)
        do k = 1, pver
          q1(:ncol, k, m) = merge(q1(:ncol, k, m), qtm(:ncol, k), lqtst(:ncol))
        end do

        ! Add the explicit surface fluxes to the lowest layer
        q1(:ncol, pver, m) = q1(:ncol, pver, m) + tmp1(:ncol)*cflx(:ncol, m)

        ! not doing molecular diffusion
        ! explicitly set mmr in top layer for cases where molecular diffusion is not active
        if (cnst_fixed_ubc(m)) then
          ! assumes TOA is at vertical level 1 here.
          q1(:ncol, 1, m) = ubc_mmr(:ncol, m)
        end if
        call no_molec_decomp%left_div(q1(:ncol, :, m))
      end if
    end do

    call no_molec_decomp%finalize()

  end subroutine vertical_diffusion_diffuse_tracers_run

  ! Based on provisional output from vertical diffusion,
  ! calculate final physics tendencies to be applied by the tendency updater.
!> \section arg_table_vertical_diffusion_tendencies_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_tendencies_run.html
  subroutine vertical_diffusion_tendencies_run( &
    ncol, pver, pcnst, &
    const_props, &
    dt, &
    pdel, pdeldry, &
    u0, v0, s0, q0, &  ! actual values at beginning of vdiff
    u1, v1, s1, q1, &  ! provisional values after vdiff, not actual
    ! below output
    tend_s, tend_u, tend_v, tend_q, &
    scheme_name, &
    errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! Input arguments
    integer,         intent(in)     :: ncol            ! Number of atmospheric columns [count]
    integer,         intent(in)     :: pver            ! Number of vertical levels [count]
    integer,         intent(in)     :: pcnst           ! Number of ccpp constituents [count]
    type(ccpp_constituent_prop_ptr_t), &
                     intent(in)     :: const_props(:)  ! CCPP constituent properties pointer
    real(kind_phys), intent(in)     :: dt              ! Timestep [s]
    real(kind_phys), intent(in)     :: pdel(:,:)       ! Pressure thickness of layers [Pa]
    real(kind_phys), intent(in)     :: pdeldry(:,:)    ! Dry pressure thickness of layers [Pa]
    real(kind_phys), intent(in)     :: u0(:,:)         ! Initial eastward wind [m s-1]
    real(kind_phys), intent(in)     :: v0(:,:)         ! Initial northward wind [m s-1]
    real(kind_phys), intent(in)     :: s0(:,:)         ! Initial dry static energy [J kg-1]
    real(kind_phys), intent(in)     :: q0(:,:,:)       ! Initial constituent mixing ratios [kg kg-1]
    real(kind_phys), intent(in)     :: u1(:,:)          ! Provisional eastward wind after diffusion [m s-1]
    real(kind_phys), intent(in)     :: v1(:,:)          ! Provisional northward wind after diffusion [m s-1]
    real(kind_phys), intent(in)     :: s1(:,:)          ! Provisional dry static energy after diffusion [J kg-1]
    real(kind_phys), intent(in)     :: q1(:,:,:)        ! Provisional constituent mixing ratios after diffusion [kg kg-1]

    ! Output arguments
    real(kind_phys), intent(out)    :: tend_s(:,:)     ! Tendency of dry air enthalpy at constant pressure [J kg-1 s-1]
    real(kind_phys), intent(out)    :: tend_u(:,:)     ! Eastward wind tendency [m s-2]
    real(kind_phys), intent(out)    :: tend_v(:,:)     ! Northward wind tendency [m s-2]
    real(kind_phys), intent(out)    :: tend_q(:,:,:)   ! Constituent mixing ratio tendencies [kg kg-1 s-1]
    character(len=64),  intent(out) :: scheme_name    ! scheme name
    character(len=512), intent(out) :: errmsg     ! Error message
    integer,            intent(out) :: errflg     ! Error flag

    integer :: m
    logical :: const_is_dry

    ! for bit-to-bitness with CAM.
    real(kind_phys) :: rdt

    scheme_name = "vertical_diffusion_tendencies"

    rdt = 1._kind_phys/dt

    errmsg = ''
    errflg = 0

    ! calculate physics tendencies
    tend_s(:ncol,:)       = (s1(:ncol,:) - s0(:ncol,:)) * rdt
    tend_u(:ncol,:)       = (u1(:ncol,:) - u0(:ncol,:)) * rdt
    tend_v(:ncol,:)       = (v1(:ncol,:) - v0(:ncol,:)) * rdt
    tend_q(:ncol,:pver,:) = (q1(:ncol,:pver,:) - q0(:ncol,:pver,:)) * rdt

    ! convert tendencies of dry constituents to dry basis
    do m = 1, pcnst
      call const_props(m)%is_dry(const_is_dry, errflg, errmsg)
      if(const_is_dry) then
        tend_q(:ncol,:pver,m) = tend_q(:ncol,:pver,m) * pdel(:ncol,:pver) / pdeldry(:ncol,:pver)
      endif
    enddo

  end subroutine vertical_diffusion_tendencies_run

end module diffusion_solver
