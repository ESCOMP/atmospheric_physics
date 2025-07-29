module calculate_net_heating
!-----------------------------------------------------------------------
!
! Purpose:  Provide an interface to convert shortwave and longwave
!           radiative heating terms into net heating.
!
!           This module provides a hook to allow incorporating additional
!           radiative terms (eUV heating and nonLTE longwave cooling).
!
! Original version: B.A. Boville
!-----------------------------------------------------------------------

implicit none
private

! Public interfaces
public :: calculate_net_heating_run

!===============================================================================
contains
!===============================================================================

!> \section arg_table_calculate_net_heating_run Argument Table
!! \htmlinclude calculate_net_heating_run.html
!!
subroutine calculate_net_heating_run(ncol, pver, rad_heat, qrl_prime, qrs_prime, &
                gravit, pdel, net_flx, errmsg, errflg)
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!-----------------------------------------------------------------------
   use ccpp_kinds,    only: kind_phys

    ! Arguments
   integer,                    intent(in)  :: ncol             ! horizontal dimension
   real(kind_phys),            intent(in)  :: qrl_prime(:,:)   ! longwave heating [J kg-1 s-1]
   real(kind_phys),            intent(in)  :: qrs_prime(:,:)   ! shortwave heating [J kg-1 s-1]
   real(kind_phys),            intent(in)  :: gravit           ! Standard gravitational acceleration [m s-2]
   real(kind_phys),            intent(in)  :: pdel(:,:)        ! Air pressure thickness [Pa]
   real(kind_phys),            intent(out) :: rad_heat(:,:)    ! radiative heating [J kg-1 s-1]
   real(kind_phys),            intent(out) :: net_flx(:)       ! net boundary flux [W m-2]
   character(len=*),           intent(out) :: errmsg
   integer,                    intent(out) :: errflg


   ! Local variables
   integer :: idx
   !-----------------------------------------------------------------------
   ! Set error variables
   errmsg = ''
   errflg = 0

   rad_heat(:,:) = (qrs_prime(:,:) + qrl_prime(:,:))

!   do idx = 1, ncol
!      net_flx(idx) = fsnt(idx) - fsns(idx) - flnt(idx) + flns(idx)
!   end do
   do kdx = 1, pver
      do idx = 1, ncol
         net_flx(idx) = net_flx(idx) + rad_heat(idx,kdx)*pdel(idx,kdx)/gravit
      end do
   end do

end subroutine calculate_net_heating_run

!================================================================================================
end module calculate_net_heating
