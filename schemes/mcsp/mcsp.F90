! Mesoscale coherent structure parameterization (MCSP)
!
! MCSP redistributes the deep convective heating and moistening vertically with
! a prescribed vertical structure to represent the effect of mesoscale organization
! (top-heavy stratiform heating), and can also add momentum tendencies, although
! that capability has not been extensively tested.
! Only the column mean deep convective tendencies are used, so the scheme conserves
! column dry static energy, kinetic energy and water by construction.
!
! References:
!   Moncrieff, M. W., & Liu, C. (2006). Representing convective organization in
!     prediction models by a hybrid strategy. J. Atmos. Sci., 63, 3404-3420.
!     https://doi.org/10.1175/JAS3812.1
!   Moncrieff, M. W., C. Liu, and P. Bogenschutz (2017). Simulation, Modeling,
!     and Dynamically Based Parameterization of Organized Tropical Convection
!     for Global Climate Models. J. Atmos. Sci., 74, 1363-1380.
!     https://doi.org/10.1175/JAS-D-16-0166.1
!   Moncrieff, M. W. (2019). Toward a Dynamical Foundation for Organized Convection
!     Parameterization in GCMs. Geophys. Res. Lett., 46, 14103-14108.
!     https://doi.org/10.1029/2019GL085316
!   Chen, C.-C., Richter, J. H., Liu, C., Moncrieff, M. W., Tang, Q., Lin, W.,
!     et al. (2021). Effects of organized convection parameterization on the MJO
!     and precipitation in E3SMv1. Part I: Mesoscale heating. J. Adv.
!     Model. Earth Syst., 13, e2020MS002401. https://doi.org/10.1029/2020MS002401
!
! Original CAM implementation by Jack Chen (NCAR), ported from E3SM.
module mcsp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mcsp_init
  public :: mcsp_run

  ! Namelist parameters, set in mcsp_init
  real(kind_phys) :: heat_coeff        ! heating coefficient [1]
  real(kind_phys) :: moisture_coeff    ! moistening coefficient [1]
  real(kind_phys) :: uwind_coeff       ! zonal wind coefficient [1]
  real(kind_phys) :: vwind_coeff       ! meridional wind coefficient [1]
  real(kind_phys) :: storm_speed_pref  ! reference pressure of the storm-level zonal wind [Pa]
  real(kind_phys) :: conv_depth_min    ! minimum pressure depth of deep convection for activation [Pa]
  real(kind_phys) :: shear_min         ! minimum low-level zonal wind shear magnitude for activation [m s-1]

  ! Shear diagnostic value where the reference pressure level is below the surface
  real(kind_phys), parameter :: shear_undefined = -999._kind_phys

  ! E3SM additionally rejects shear magnitudes at or above this bound, which also
  ! excludes columns flagged with shear_undefined. Disabled pending evaluation in CAM;
  ! see the commented-out condition in mcsp_run.
  ! real(kind_phys), parameter :: shear_max = 200._kind_phys

contains

!> \section arg_table_mcsp_init Argument Table
!! \htmlinclude mcsp_init.html
  subroutine mcsp_init(amIRoot, iulog, &
                       mcsp_heat_coeff, mcsp_moisture_coeff, mcsp_uwind_coeff, mcsp_vwind_coeff, &
                       mcsp_storm_speed_pref, mcsp_conv_depth_min, mcsp_shear_min, &
                       errmsg, errflg)

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    real(kind_phys),    intent(in)  :: mcsp_heat_coeff        ! heating coefficient [1]
    real(kind_phys),    intent(in)  :: mcsp_moisture_coeff    ! moistening coefficient [1]
    real(kind_phys),    intent(in)  :: mcsp_uwind_coeff       ! zonal wind coefficient [1]
    real(kind_phys),    intent(in)  :: mcsp_vwind_coeff       ! meridional wind coefficient [1]
    real(kind_phys),    intent(in)  :: mcsp_storm_speed_pref  ! reference pressure of the storm-level zonal wind [Pa]
    real(kind_phys),    intent(in)  :: mcsp_conv_depth_min    ! minimum pressure depth of deep convection [Pa]
    real(kind_phys),    intent(in)  :: mcsp_shear_min         ! minimum low-level zonal wind shear magnitude [m s-1]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    heat_coeff       = mcsp_heat_coeff
    moisture_coeff   = mcsp_moisture_coeff
    uwind_coeff      = mcsp_uwind_coeff
    vwind_coeff      = mcsp_vwind_coeff
    storm_speed_pref = mcsp_storm_speed_pref
    conv_depth_min   = mcsp_conv_depth_min
    shear_min        = mcsp_shear_min

    if (amIRoot) then
       write(iulog,*) 'mcsp_init: mcsp_heat_coeff       = ', heat_coeff
       write(iulog,*) 'mcsp_init: mcsp_moisture_coeff   = ', moisture_coeff
       write(iulog,*) 'mcsp_init: mcsp_uwind_coeff      = ', uwind_coeff
       write(iulog,*) 'mcsp_init: mcsp_vwind_coeff      = ', vwind_coeff
       write(iulog,*) 'mcsp_init: mcsp_storm_speed_pref = ', storm_speed_pref
       write(iulog,*) 'mcsp_init: mcsp_conv_depth_min   = ', conv_depth_min
       write(iulog,*) 'mcsp_init: mcsp_shear_min        = ', shear_min
    end if

  end subroutine mcsp_init

!> \section arg_table_mcsp_run Argument Table
!! \htmlinclude mcsp_run.html
  subroutine mcsp_run(ncol, pver, pverp, cpair, pi, dt, jctop, &
                      pmid, pint, pdel, u, v, &
                      ttend_dp, qtend_dp, &
                      ptend_s, ptend_q, ptend_u, ptend_v, &
                      mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
                      mcsp_freq, mcsp_shear, conv_depth, mcsp_dt_max, &
                      errmsg, errflg)

    ! Arguments
    integer,            intent(in)  :: ncol              ! number of atmospheric columns
    integer,            intent(in)  :: pver              ! number of vertical layers
    integer,            intent(in)  :: pverp             ! number of vertical interfaces
    real(kind_phys),    intent(in)  :: cpair             ! specific heat of dry air at constant pressure [J kg-1 K-1]
    real(kind_phys),    intent(in)  :: pi                ! pi [1]
    real(kind_phys),    intent(in)  :: dt                ! physics time step [s]
    integer,            intent(in)  :: jctop(:)          ! vertical index at top of deep convection, pver where none [index]
    real(kind_phys),    intent(in)  :: pmid(:,:)         ! mid-point pressure [Pa]
    real(kind_phys),    intent(in)  :: pint(:,:)         ! interface pressure [Pa]
    real(kind_phys),    intent(in)  :: pdel(:,:)         ! pressure thickness [Pa]
    real(kind_phys),    intent(in)  :: u(:,:)            ! zonal wind [m s-1]
    real(kind_phys),    intent(in)  :: v(:,:)            ! meridional wind [m s-1]
    real(kind_phys),    intent(in)  :: ttend_dp(:,:)     ! core deep convective temperature tendency [K s-1]
    real(kind_phys),    intent(in)  :: qtend_dp(:,:)     ! core deep convective water vapor tendency [kg kg-1 s-1]
    real(kind_phys),    intent(out) :: ptend_s(:,:)      ! dry static energy tendency [J kg-1 s-1]
    real(kind_phys),    intent(out) :: ptend_q(:,:)      ! water vapor tendency [kg kg-1 s-1]
    real(kind_phys),    intent(out) :: ptend_u(:,:)      ! zonal wind tendency [m s-2]
    real(kind_phys),    intent(out) :: ptend_v(:,:)      ! meridional wind tendency [m s-2]
    real(kind_phys),    intent(out) :: mcsp_dt_out(:,:)  ! diagnostic temperature tendency [K s-1]
    real(kind_phys),    intent(out) :: mcsp_dq_out(:,:)  ! diagnostic water vapor tendency [kg kg-1 s-1]
    real(kind_phys),    intent(out) :: mcsp_du_out(:,:)  ! diagnostic zonal wind tendency [m s-2]
    real(kind_phys),    intent(out) :: mcsp_dv_out(:,:)  ! diagnostic meridional wind tendency [m s-2]
    real(kind_phys),    intent(out) :: mcsp_freq(:)      ! 1 where MCSP contributed a tendency, else 0 [1]
    real(kind_phys),    intent(out) :: mcsp_shear(:)     ! low-level zonal wind shear [m s-1]
    real(kind_phys),    intent(out) :: conv_depth(:)     ! pressure depth of deep convection [Pa]
    real(kind_phys),    intent(out) :: mcsp_dt_max(:)    ! heating amplitude, column mean deep heating times coefficient [K s-1]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer  :: i, k

    real(kind_phys) :: tend_k        ! kinetic energy tendency [J kg-1 s-1]
    real(kind_phys) :: pdepth_mid_k  ! pressure depth from the surface to layer k [Pa]
    real(kind_phys) :: pdepth_total  ! pressure depth from the surface to the convective top [Pa]

    real(kind_phys) :: avg_tend_s(ncol)        ! mass weighted column average deep convective DSE tendency
    real(kind_phys) :: avg_tend_q(ncol)        ! mass weighted column average deep convective qv tendency
    real(kind_phys) :: pdel_sum(ncol)          ! pressure thickness of the convective column

    real(kind_phys) :: mcsp_tend_s(ncol,pver)  ! MCSP DSE tendency before the energy fixer
    real(kind_phys) :: mcsp_tend_q(ncol,pver)  ! MCSP qv tendency before the mass fixer
    real(kind_phys) :: mcsp_tend_u(ncol,pver)  ! MCSP zonal wind tendency
    real(kind_phys) :: mcsp_tend_v(ncol,pver)  ! MCSP meridional wind tendency

    real(kind_phys) :: mcsp_avg_tend_s(ncol)   ! mass weighted column average MCSP DSE tendency
    real(kind_phys) :: mcsp_avg_tend_q(ncol)   ! mass weighted column average MCSP qv tendency
    real(kind_phys) :: mcsp_avg_tend_k(ncol)   ! mass weighted column average MCSP kinetic energy tendency

    logical :: do_mcsp_t   ! compute the temperature tendency
    logical :: do_mcsp_q   ! compute the water vapor tendency
    logical :: do_mcsp_u   ! compute the zonal wind tendency
    logical :: do_mcsp_v   ! compute the meridional wind tendency

    errmsg = ''
    errflg = 0

    do i = 1, ncol
       if (jctop(i) < 1 .or. jctop(i) > pver) then
          errflg = 1
          write(errmsg, '(a,i0,a,i0)') 'mcsp_run: deep convective top index out of range in column ', i, ': ', jctop(i)
          return
       end if
    end do

    !----------------------------------------------------------------------------
    ! initialize variables

    do_mcsp_t = heat_coeff > 0._kind_phys
    do_mcsp_q = moisture_coeff > 0._kind_phys
    do_mcsp_u = uwind_coeff > 0._kind_phys
    do_mcsp_v = vwind_coeff > 0._kind_phys

    ptend_s(:,:) = 0._kind_phys
    ptend_q(:,:) = 0._kind_phys
    ptend_u(:,:) = 0._kind_phys
    ptend_v(:,:) = 0._kind_phys

    avg_tend_s(:) = 0._kind_phys
    avg_tend_q(:) = 0._kind_phys

    pdel_sum(:) = 0._kind_phys

    mcsp_avg_tend_s(:) = 0._kind_phys
    mcsp_avg_tend_q(:) = 0._kind_phys
    mcsp_avg_tend_k(:) = 0._kind_phys

    mcsp_tend_s(:,:) = 0._kind_phys
    mcsp_tend_q(:,:) = 0._kind_phys
    mcsp_tend_u(:,:) = 0._kind_phys
    mcsp_tend_v(:,:) = 0._kind_phys

    mcsp_shear(:)  = 0._kind_phys
    conv_depth(:)  = 0._kind_phys
    mcsp_dt_max(:) = 0._kind_phys

    if (do_mcsp_t .or. do_mcsp_q .or. do_mcsp_u .or. do_mcsp_v) then

       !----------------------------------------------------------------------------
       ! calculate shear

       call mcsp_calculate_shear(ncol, pver, pmid, u, mcsp_shear)

       !----------------------------------------------------------------------------
       ! calculate mass weighted column average tendencies from deep convection

       do i = 1, ncol
          if (jctop(i) /= pver) then
             ! integrate pressure and deep convective tendencies over column
             do k = jctop(i), pver
                avg_tend_s(i) = avg_tend_s(i) + ttend_dp(i,k) * pdel(i,k) * cpair
                avg_tend_q(i) = avg_tend_q(i) + qtend_dp(i,k) * pdel(i,k)
                pdel_sum(i) = pdel_sum(i) + pdel(i,k)
             end do
             ! normalize integrated deep convective tendencies by total mass
             avg_tend_s(i) = avg_tend_s(i) / pdel_sum(i)
             avg_tend_q(i) = avg_tend_q(i) / pdel_sum(i)
             ! calculate diagnostic deep convective depth
             conv_depth(i) = pint(i,pver+1) - pmid(i,jctop(i))
          else
             avg_tend_s(i) = 0._kind_phys
             avg_tend_q(i) = 0._kind_phys
             conv_depth(i) = 0._kind_phys
          end if
       end do

       !----------------------------------------------------------------------------
       ! Note: To conserve total energy we need to account for the kinetic energy tendency
       ! which we can obtain from the velocity tendencies based on the following:
       !   KE_new = (u_new^2 + v_new^2)/2
       !          = [ (u_old+du)^2 + (v_old+dv)^2 ]/2
       !          = [ ( u_old^2 + 2*u_old*du + du^2 ) + ( v_old^2 + 2*v_old*dv + dv^2 ) ]/2
       !          = ( u_old^2 + v_old^2 )/2 + ( 2*u_old*du + du^2 + 2*v_old*dv + dv^2 )/2
       !          = KE_old + [ 2*u_old*du + du^2 + 2*v_old*dv + dv^2 ] /2

       !----------------------------------------------------------------------------
       ! calculate MCSP tendencies

       do i = 1, ncol

          ! check that deep convection produced tendencies over a depth that exceeds the threshold
          if ( conv_depth(i) >= conv_depth_min ) then
             ! check that deep convection provided a non-zero column total heating
             if ( avg_tend_s(i) > 0._kind_phys ) then
                ! check that there is sufficient wind shear to justify coherent organization
                ! E3SM also requires abs(mcsp_shear(i)) < shear_max (see module header):
                ! if ( abs(mcsp_shear(i)) >= shear_min .and. abs(mcsp_shear(i)) < shear_max ) then
                if ( abs(mcsp_shear(i)) >= shear_min ) then

                   if (do_mcsp_t) mcsp_dt_max(i) = avg_tend_s(i) * heat_coeff / cpair

                   do k = jctop(i), pver

                      ! See eq 7-8 of Moncrieff et al. (2017) - also eq (5) of Moncrieff & Liu (2006)
                      pdepth_mid_k = pint(i,pver+1) - pmid(i,k)
                      pdepth_total = pint(i,pver+1) - pmid(i,jctop(i))

                      ! specify the assumed vertical structure
                      if (do_mcsp_t) mcsp_tend_s(i,k) = -1._kind_phys*heat_coeff * sin(2.0_kind_phys*pi*(pdepth_mid_k/pdepth_total))
                      if (do_mcsp_q) mcsp_tend_q(i,k) = -1._kind_phys*moisture_coeff * sin(2.0_kind_phys*pi*(pdepth_mid_k/pdepth_total))
                      if (do_mcsp_u) mcsp_tend_u(i,k) = uwind_coeff * (cos(pi*(pdepth_mid_k/pdepth_total)))
                      if (do_mcsp_v) mcsp_tend_v(i,k) = vwind_coeff * (cos(pi*(pdepth_mid_k/pdepth_total)))

                      ! scale the vertical structure by the deep convective heating/drying tendencies
                      if (do_mcsp_t) mcsp_tend_s(i,k) = avg_tend_s(i) * mcsp_tend_s(i,k)
                      if (do_mcsp_q) mcsp_tend_q(i,k) = avg_tend_q(i) * mcsp_tend_q(i,k)

                      ! integrate the DSE/qv tendencies for energy/mass fixer
                      if (do_mcsp_t) mcsp_avg_tend_s(i) = mcsp_avg_tend_s(i) + mcsp_tend_s(i,k) * pdel(i,k) / pdel_sum(i)
                      if (do_mcsp_q) mcsp_avg_tend_q(i) = mcsp_avg_tend_q(i) + mcsp_tend_q(i,k) * pdel(i,k) / pdel_sum(i)

                      ! integrate the change in kinetic energy (KE) for energy fixer
                      if (do_mcsp_u .or. do_mcsp_v) then
                         tend_k = ( 2.0_kind_phys*mcsp_tend_u(i,k)*dt*u(i,k) + mcsp_tend_u(i,k)*mcsp_tend_u(i,k)*dt*dt &
                                   +2.0_kind_phys*mcsp_tend_v(i,k)*dt*v(i,k) + mcsp_tend_v(i,k)*mcsp_tend_v(i,k)*dt*dt &
                                  )/2.0_kind_phys/dt
                         mcsp_avg_tend_k(i) = mcsp_avg_tend_k(i) + tend_k*pdel(i,k) / pdel_sum(i)
                      end if

                   end do ! k = jctop(i), pver
                end if ! shear threshold
             end if ! avg_tend_s(i) > 0
          end if ! conv_depth(i) >= conv_depth_min
       end do

    end if

    !----------------------------------------------------------------------------
    ! calculate final output tendencies

    mcsp_dt_out(:,:) = 0._kind_phys
    mcsp_dq_out(:,:) = 0._kind_phys
    mcsp_du_out(:,:) = 0._kind_phys
    mcsp_dv_out(:,:) = 0._kind_phys

    mcsp_freq(:) = 0._kind_phys

    do i = 1, ncol
       do k = jctop(i), pver

          ! update frequency if MCSP contributes any tendency in the column
          if ( abs(mcsp_tend_s(i,k)) > 0._kind_phys .or. abs(mcsp_tend_q(i,k)) > 0._kind_phys .or. &
               abs(mcsp_tend_u(i,k)) > 0._kind_phys .or. abs(mcsp_tend_v(i,k)) > 0._kind_phys ) then
             mcsp_freq(i) = 1._kind_phys
          end if

          ! subtract mass weighted average tendencies for energy/mass conservation
          mcsp_dt_out(i,k) = mcsp_tend_s(i,k) - mcsp_avg_tend_s(i)
          mcsp_dq_out(i,k) = mcsp_tend_q(i,k) - mcsp_avg_tend_q(i)
          mcsp_du_out(i,k) = mcsp_tend_u(i,k)
          mcsp_dv_out(i,k) = mcsp_tend_v(i,k)

          ! make sure kinetic energy correction is added to DSE tendency
          ! to conserve total energy whenever momentum tendencies are calculated
          if (do_mcsp_u .or. do_mcsp_v) then
             mcsp_dt_out(i,k) = mcsp_dt_out(i,k) - mcsp_avg_tend_k(i)
          end if

          ! update output tendencies
          if (do_mcsp_t) ptend_s(i,k) = mcsp_dt_out(i,k)
          if (do_mcsp_q) ptend_q(i,k) = mcsp_dq_out(i,k)
          if (do_mcsp_u) ptend_u(i,k) = mcsp_du_out(i,k)
          if (do_mcsp_v) ptend_v(i,k) = mcsp_dv_out(i,k)

          ! adjust units for diagnostic outputs
          if (do_mcsp_t) mcsp_dt_out(i,k) = mcsp_dt_out(i,k)/cpair

       end do
    end do

  end subroutine mcsp_run

  ! Low-level zonal wind shear between the storm reference pressure level and the
  ! lowest model layer. Columns whose lowest layer is above the reference level get
  ! shear_undefined.
  subroutine mcsp_calculate_shear(ncol, pver, pmid, u, mcsp_shear)

    use interpolate_data, only: vertinterp

    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    real(kind_phys), intent(in)  :: pmid(:,:)      ! mid-point pressure [Pa]
    real(kind_phys), intent(in)  :: u(:,:)         ! zonal wind [m s-1]
    real(kind_phys), intent(out) :: mcsp_shear(:)  ! low-level zonal wind shear [m s-1]

    integer         :: i
    real(kind_phys) :: storm_u(ncol)  ! zonal wind at the storm reference pressure level [m s-1]

    ! Interpolate wind to pressure level specified by storm_speed_pref
    call vertinterp(ncol, ncol, pver, pmid, storm_speed_pref, u, storm_u)

    ! calculate low-level shear
    do i = 1, ncol
       if (pmid(i,pver) > storm_speed_pref) then
          mcsp_shear(i) = storm_u(i) - u(i,pver)
       else
          mcsp_shear(i) = shear_undefined
       end if
    end do

  end subroutine mcsp_calculate_shear

end module mcsp
