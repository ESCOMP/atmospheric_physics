module rrtmgp_lw_calculate_heating_rate
   public :: rrtmgp_lw_calculate_heating_rate_run

CONTAINS
   !> \section arg_table_rrtmgp_lw_calculate_heating_rate_run  Argument Table
   !! \htmlinclude rrtmgp_lw_calculate_heating_rate_run.html
   subroutine rrtmgp_lw_calculate_heating_rate_run(ktopcam, pver, gravit, rpdel, flux_net, &
                   flux_net_clrsky, hrate, hrate_clrsky, errmsg, errflg)
      use ccpp_kinds, only: kind_phys
      integer,          intent(in) :: ktopcam
      integer,          intent(in) :: pver
      real(kind_phys),  intent(in) :: gravit
      real(kind_phys),  intent(in) :: rpdel(:,:)
      real(kind_phys),  intent(in) :: flux_net(:,:)
      real(kind_phys),  intent(in) :: flux_net_clrsky(:,:)
      real(kind_phys), intent(out) :: hrate(:,:)
      real(kind_phys), intent(out) :: hrate_clrsky(:,:)
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: k

      ! Set error variables
      errmsg = ''
      errflg = 0

      hrate = 0.0_kind_phys
      hrate_clrsky = 0.0_kind_phys

      do k = ktopcam, pver
         ! (flux divergence as bottom-MINUS-top) * g/dp
         hrate(:,k) = (flux_net(:,k+1) - flux_net(:,k)) * &
                              gravit * rpdel(:,k)
         hrate_clrsky(:,k) = (flux_net_clrsky(:,k+1) - flux_net_clrsky(:,k)) * &
                              gravit * rpdel(:,k)
      end do

   end subroutine rrtmgp_lw_calculate_heating_rate_run

end module rrtmgp_lw_calculate_heating_rate
