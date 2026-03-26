!-----------------------------------------------------------------------
! Module to handle data that is exchanged between the atmosphere
! model and the surface models (land, sea-ice, ocean, etc.).

! Please note that currently this is a SIMA-specific module,
! but could be made more host model-independent in the future.
!-----------------------------------------------------------------------
module set_surface_coupling_vars

  implicit none
  private

  ! Public interfaces
  public :: set_surface_coupling_vars_run

!===============================================================================
CONTAINS
!===============================================================================

!> \section arg_table_set_surface_coupling_vars_run Argument Table
!! \htmlinclude set_surface_coupling_vars_run.html
!!
subroutine set_surface_coupling_vars_run(ncol, pver, ncnst, gravit, rair, phis,     &
                                         surf_pres, air_temp, inv_ps_exner, zm,     &
                                         rho, uwnd, vwnd, air_pres, cnst_array,     &
                                         prec_dp, snow_dp, prec_sh, snow_sh,        &
                                         prec_str, snow_str, air_temp_bot,          &
                                         pot_temp_bot, zm_bot, rho_bot,             &
                                         uwnd_bot, vwnd_bot, air_pres_bot,          &
                                         topo_height, sea_lev_pres, cnst_bot,       &
                                         conv_prec, conv_snow, strat_prec,          &
                                         strat_snow, errmsg, errcode)

   ! Set surface variables needed for atmosphere-surface coupling.

   ! Use statements
   use ccpp_kinds,       only: kind_phys

   ! Input arguments
   integer, intent(in) :: ncol
   integer, intent(in) :: pver
   integer, intent(in) :: ncnst

   real(kind_phys), intent(in) :: gravit            ! Gravitational acceleration [m s-2]
   real(kind_phys), intent(in) :: rair              ! Dry air gas constant [J kg-1 K-1]
   real(kind_phys), intent(in) :: phis(:)           ! Surface geopotential [m2 s-2]
   real(kind_phys), intent(in) :: surf_pres(:)      ! Surface pressure [Pa]
   real(kind_phys), intent(in) :: air_temp(:,:)     ! Air temperature [K]
   real(kind_phys), intent(in) :: inv_ps_exner(:,:) ! Inverse Exner function w.r.t. surface pressure  [1]
   real(kind_phys), intent(in) :: zm(:,:)           ! Geopotential height wrt surface [m]
   real(kind_phys), intent(in) :: rho(:,:)          ! Dry air density [kg m-3]
   real(kind_phys), intent(in) :: uwnd(:,:)         ! Eastward wind [m s-1]
   real(kind_phys), intent(in) :: vwnd(:,:)         ! Northward wind [m s-1]
   real(kind_phys), intent(in) :: air_pres(:,:)     ! Air pressure [Pa]
   real(kind_phys), intent(in) :: cnst_array(:,:,:) ! CCPP constituents array [none]

   real(kind_phys), intent(in) :: prec_dp(:)        ! LWE deep convective precipitation (all phases) [m s-1]
   real(kind_phys), intent(in) :: snow_dp(:)        ! LWE deep convective frozen precipitation (e.g. snow) [m s-1]
   real(kind_phys), intent(in) :: prec_sh(:)        ! LWE shallow convective precipitation (all phases) [m s-1]
   real(kind_phys), intent(in) :: snow_sh(:)        ! LWE shallow convective frozen precipitation (all phases) [m s-1]
   real(kind_phys), intent(in) :: prec_str(:)       ! LWE stratiform precipitation (all phases) [m s-1]
   real(kind_phys), intent(in) :: snow_str(:)       ! LWE stratiform frozen precipitation (e.g. snow) [m s-1]

   ! Output arguments
   real(kind_phys), intent(out) :: air_temp_bot(:) ! Air temperature at bottom of atmosphere for coupling [K]
   real(kind_phys), intent(out) :: pot_temp_bot(:) ! Potential air temprature at bottom of atmosphere for coupling [K]
   real(kind_phys), intent(out) :: zm_bot(:)       ! Geopotential height wrt surface at bottom of atmosphere for coupling [m]
   real(kind_phys), intent(out) :: rho_bot(:)      ! Dry air density at bottom of atmosphere for coupling [kg m-3]
   real(kind_phys), intent(out) :: uwnd_bot(:)     ! Eastward wind at bottom of atmosphere for coupling [m s-1]
   real(kind_phys), intent(out) :: vwnd_bot(:)     ! Northward wind at bottom of atmosphere for coupling [m s-1]
   real(kind_phys), intent(out) :: air_pres_bot(:) ! Air pressure at bottom of atmosphere for coupling [Pa]
   real(kind_phys), intent(out) :: topo_height(:)  ! Geopotential height of surface topography [m]
   real(kind_phys), intent(out) :: sea_lev_pres(:) ! Sea Level Pressure [Pa]
   real(kind_phys), intent(out) :: cnst_bot(:,:)   ! Constituent mixing ratios at bottom of atmosphere for coupling [kg kg-1]

   real(kind_phys), intent(out) :: conv_prec(:)    ! LWE Convective precipitation (all phases) [m s-1]
   real(kind_phys), intent(out) :: conv_snow(:)    ! LWE Frozen convective precipitation (e.g. snow) [m s-1]
   real(kind_phys), intent(out) :: strat_prec(:)   ! LWE Stratiform (large-scale) precipitation (all phases) [m s-1]
   real(kind_phys), intent(out) :: strat_snow(:)   ! LWE Frozen stratiform (large-scale)precipitation (e.g. snow) [m s-1]

   character(len=512), intent(out) :: errmsg       ! CCPP error message
   integer,            intent(out) :: errcode      ! CCPP error code

   ! Local variables

   integer :: i              ! column index
   integer :: m              ! constituent index

   !-----------------------------------------------------------------------

   errmsg  = ''
   errcode = 0

   ! Copy over atmospheric 3-D variables to surface fields:
   do i=1,ncol
     air_temp_bot(i) = air_temp(i,pver)
     pot_temp_bot(i) = air_temp(i,pver) * inv_ps_exner(i,pver)
     zm_bot(i)       = zm(i,pver)
     rho_bot(i)      = rho(i,pver)
     uwnd_bot(i)     = uwnd(i,pver)
     vwnd_bot(i)     = vwnd(i,pver)
     air_pres_bot(i) = air_pres(i,pver)
     topo_height(i)  = phis(i)/gravit
   end do

   ! Calculate Sea Level Pressure (PSL):
   call calc_sea_level_pressure(ncol, gravit, rair, air_pres_bot, phis, &
                                surf_pres, air_temp_bot, sea_lev_pres)

   ! Set surface constituent values:
   do m = 1, ncnst
     do i = 1, ncol
       cnst_bot(i,m) = cnst_array(i,pver,m)
     end do
   end do

   !
   ! Precipation and snow rates from shallow convection, deep convection and stratiform processes.
   ! Compute total convective and stratiform precipitation and snow rates
   !
   do i=1,ncol
      conv_prec(i)  = prec_dp(i) + prec_sh(i)
      conv_snow(i)  = snow_dp(i) + snow_sh(i)

      strat_prec(i) = prec_str(i)  !Might be better to have microphysics set these directly
      strat_snow(i) = snow_str(i)
   end do

end subroutine set_surface_coupling_vars_run

!***************************************************
! private subroutine to calculate sea level pressure
!***************************************************

subroutine calc_sea_level_pressure(ncol, gravit, rair, pbot, phis, ps, tbot, psl)

!-----------------------------------------------------------------------
!
! Compute sea level pressure.
!
! Uses ECMWF formulation Algorithm: See section 3.1.b in NCAR NT-396 "Vertical
! Interpolation and Truncation of Model-Coordinate Data
!
!-----------------------------------------------------------------------

  ! Use statements
  use ccpp_kinds, only: kind_phys

  !-----------------------------Arguments---------------------------------
  integer , intent(in) :: ncol             ! longitude dimension

  real(kind_phys), intent(in) :: gravit     ! Gravitational acceleration [m s-2]
  real(kind_phys), intent(in) :: rair       ! gas constant for dry air [J kg-1 K-1]
  real(kind_phys), intent(in) :: pbot(:)    ! Bottom layer atmospheric pressure [Pa]
  real(kind_phys), intent(in) :: phis(:)    ! Surface geopotential [m2 s-2]
  real(kind_phys), intent(in) :: ps(:)      ! Surface air pressure [Pa]
  real(kind_phys), intent(in) :: Tbot(:)    ! Bottom layer Atmospheric temperature [K]

  real(kind_phys), intent(out):: psl(:)     ! Sea level pressures [Pa]

  !-----------------------------Parameters--------------------------------
  real(kind_phys), parameter :: xlapse = 6.5e-3_kind_phys   ! Temperature lapse rate [K m-1]

  !-----------------------------Local Variables---------------------------
  integer  :: i             ! Loop index
  real(kind_phys) :: alpha  ! Temperature lapse rate in terms of pressure ratio (unitless)
  real(kind_phys) :: Tstar  ! Computed surface temperature
  real(kind_phys) :: TT0    ! Computed temperature at sea-level
  real(kind_phys) :: alph   ! Power to raise P/Ps to get rate of increase of T with pressure
  real(kind_phys) :: beta   ! alpha*phis/(R*T) term used in approximation of PSL
  !-----------------------------------------------------------------------

  alpha = rair*xlapse/gravit
  do i=1,ncol
    if ( abs(phis(i)/gravit) < 1.e-4_kind_phys )then
      psl(i)=ps(i)
    else
      Tstar=Tbot(i)*(1._kind_phys+alpha*(ps(i)/pbot(i)-1._kind_phys)) ! pg 7 eq 5

      TT0=Tstar + xlapse*phis(i)/gravit                                      ! pg 8 eq 13

      if ( Tstar<=290.5_kind_phys .and. TT0>290.5_kind_phys ) then           ! pg 8 eq 14.1
        alph=rair/phis(i)*(290.5_kind_phys-Tstar)
      else if (Tstar>290.5_kind_phys  .and. TT0>290.5_kind_phys) then        ! pg 8 eq 14.2
        alph=0._kind_phys
        Tstar= 0.5_kind_phys * (290.5_kind_phys + Tstar)
      else
        alph=alpha
        if (Tstar<255._kind_phys) then
          Tstar= 0.5_kind_phys * (255._kind_phys + Tstar)                    ! pg 8 eq 14.3
        end if
      end if

      beta = phis(i)/(rair*Tstar)
      psl(i)=ps(i)*exp( beta*(1._kind_phys-alph*beta/2._kind_phys+((alph*beta)**2)/3._kind_phys))
    end if
  end do

end subroutine calc_sea_level_pressure

end module set_surface_coupling_vars
