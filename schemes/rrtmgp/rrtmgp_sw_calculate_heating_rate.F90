module rrtmgp_sw_calculate_heating_rate
   public :: rrtmgp_sw_calculate_heating_rate_run

CONTAINS
   !> \section arg_table_rrtmgp_sw_calculate_heating_rate_run  Argument Table
   !! \htmlinclude rrtmgp_sw_calculate_heating_rate_run.html
   subroutine rrtmgp_sw_calculate_heating_rate_run(ktopcam, pver, gravit, rpdel, flux_net, hrate, errmsg, errflg)
      use ccpp_kinds, only: kind_phys
      integer,          intent(in) :: ktopcam
      integer,          intent(in) :: pver
      real(kind_phys),  intent(in) :: gravit
      real(kind_phys),  intent(in) :: rpdel(:,:)
      real(kind_phys),  intent(in) :: flux_net(:,:)
      real(kind_phys), intent(out) :: hrate(:,:)
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
      end do

   end subroutine rrtmgp_sw_calculate_heating_rate_run

end module rrtmgp_sw_calculate_heating_rate
