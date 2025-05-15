module rrtmgp_dry_static_energy_tendency
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
public :: rrtmgp_dry_static_energy_tendency_run

!===============================================================================
contains
!===============================================================================

!> \section arg_table_rrtmgp_dry_static_energy_tendency_run Argument Table
!! \htmlinclude rrtmgp_dry_static_energy_tendency_run.html
!!
subroutine rrtmgp_dry_static_energy_tendency_run(pdel, calc_sw_heat, calc_lw_heat, &
                qrs, qrl, qrs_prime, qrl_prime, errmsg, errflg)
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!-----------------------------------------------------------------------
   use ccpp_kinds,    only: kind_phys

   ! Arguments
   real(kind_phys), dimension(:,:), intent(in)    :: pdel          ! Layer thickness
   logical,                         intent(in)    :: calc_sw_heat  ! Flag to calculate net shortwave heating
   logical,                         intent(in)    :: calc_lw_heat  ! Flag to calculate net longwave heating
   real(kind_phys), dimension(:,:), intent(in)    :: qrs           ! shortwave heating rate adjusted by air pressure thickness (J Pa kg-1 s-1)
   real(kind_phys), dimension(:,:), intent(in)    :: qrl           ! longwave heating rate adjusted by air pressure thickness  (J Pa kg-1 s-1)
   real(kind_phys), dimension(:,:), intent(out)   :: qrs_prime     ! shortwave heating rate (J kg-1 s-1)
   real(kind_phys), dimension(:,:), intent(out)   :: qrl_prime     ! longwave heating rate  (J kg-1 s-1)
   character(len=*),                intent(out)   :: errmsg
   integer,                         intent(out)   :: errflg


   !-----------------------------------------------------------------------
   ! Set error variables
   errmsg = ''
   errflg = 0

   if (calc_sw_heat) then
      qrs_prime(:,:) = qrs(:,:) / pdel(:,:)
   end if

   if (calc_lw_heat) then
      qrl_prime(:,:) = qrl(:,:) / pdel(:,:)
   end if

end subroutine rrtmgp_dry_static_energy_tendency_run

!================================================================================================
end module rrtmgp_dry_static_energy_tendency
