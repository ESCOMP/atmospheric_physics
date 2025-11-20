module rrtmgp_sw_calculate_heating_rate
   public :: rrtmgp_sw_calculate_heating_rate_run

CONTAINS
   !> \section arg_table_rrtmgp_sw_calculate_heating_rate_run  Argument Table
   !! \htmlinclude rrtmgp_sw_calculate_heating_rate_run.html
   subroutine rrtmgp_sw_calculate_heating_rate_run(ktopcam, pver, gravit, rpdel, flux_net, flux_net_clrsky, &
                   hrate, hrate_clrsky, errmsg, errflg)
      use ccpp_kinds, only: kind_phys
      integer,          intent(in) :: ktopcam               ! Vertical index at top level where RRTMGP is active
      integer,          intent(in) :: pver                  ! Number of vertical layers
      real(kind_phys),  intent(in) :: gravit                ! Standard gravitiational acceleration [m s-2]
      real(kind_phys),  intent(in) :: rpdel(:,:)            ! Reciprocal of air pressure thickness [Pa-1]
      real(kind_phys),  intent(in) :: flux_net(:,:)         ! Shortwave net radiative flux [W m-2]
      real(kind_phys),  intent(in) :: flux_net_clrsky(:,:)  ! Shortwave net radiative clear-sky flux [W m-2]
      real(kind_phys), intent(out) :: hrate(:,:)            ! Tendency of dry air enthalpy due to SW radiation [J kg-1 s-1]
      real(kind_phys), intent(out) :: hrate_clrsky(:,:)     ! Tendency of dry air enthalpy due to SW clear-sky radiation [J kg-1 s-1]
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: k

      ! Set error variables
      errmsg = ''
      errflg = 0

      hrate = 0.0_kind_phys

      do k = ktopcam, pver
         ! top - bottom
         hrate(:,k) = (flux_net(:,k) - flux_net(:,k+1)) * &
                 gravit * rpdel(:,k)
         hrate_clrsky(:,k) = (flux_net_clrsky(:,k) - flux_net_clrsky(:,k+1)) * &
                 gravit * rpdel(:,k)
      end do

   end subroutine rrtmgp_sw_calculate_heating_rate_run

end module rrtmgp_sw_calculate_heating_rate
