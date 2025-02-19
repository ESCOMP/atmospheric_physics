module gw_front

!
! This module handles gravity waves from frontal sources, and was extracted
! from gw_drag in May 2013.
!

use ccpp_kinds, only:  kind_phys
use gw_common, only: GWBand, pi, pver

implicit none
private
save

public :: CMSourceDesc
public :: flat_cm_desc
public :: gaussian_cm_desc
public :: gw_cm_src

! Tuneable settings.
type CMSourceDesc
   ! Source level.
   integer :: ksrc
   ! Level at which to check whether the frontogenesis function is above
   ! the critical threshold.
   integer :: kfront
   ! Frontogenesis function critical threshold.
   real(kind_phys) :: frontgfc = huge(1._kind_phys)
   ! The stress spectrum to produce at the source level.
   real(kind_phys), allocatable :: src_tau(:)
end type CMSourceDesc

contains

! Create a flat profile to be launched (all wavenumbers have the same
! source strength, except that l=0 is excluded).
function flat_cm_desc(band, ksrc, kfront, frontgfc, taubgnd) result(desc)
  ! Wavelengths triggered by frontogenesis.
  type(GWBand), intent(in) :: band
  ! The following are used to set the corresponding object components.
  integer, intent(in) :: ksrc
  integer, intent(in) :: kfront
  real(kind_phys), intent(in) :: frontgfc
  ! Amount of stress to launch from each wavelength.
  real(kind_phys), intent(in) :: taubgnd

  type(CMSourceDesc) :: desc

  desc%ksrc = ksrc
  desc%kfront = kfront
  desc%frontgfc = frontgfc

  allocate(desc%src_tau(-band%ngwv:band%ngwv))
  desc%src_tau = taubgnd

  ! Prohibit wavenumber 0.
  desc%src_tau(0) = 0._kind_phys

end function flat_cm_desc

! Create a source tau profile that is a gaussian over wavenumbers (l=0 is
! excluded).
function gaussian_cm_desc(band, ksrc, kfront, frontgfc, height, width) &
     result(desc)

  use shr_spfn_mod, only: erfc => shr_spfn_erfc

  ! Wavelengths triggered by frontogenesis.
  type(GWBand), intent(in) :: band
  ! The following are used to set the corresponding object components.
  integer, intent(in) :: ksrc
  integer, intent(in) :: kfront
  real(kind_phys), intent(in) :: frontgfc
  ! Parameters of gaussian.
  real(kind_phys), intent(in) :: height
  real(kind_phys), intent(in) :: width

  type(CMSourceDesc) :: desc

  ! Bounds used to average bins of the gaussian.
  real(kind_phys) :: gaussian_bounds(2*band%ngwv+2)

  ! Wavenumber index.
  integer :: l

  desc%ksrc = ksrc
  desc%kfront = kfront
  desc%frontgfc = frontgfc

  allocate(desc%src_tau(-band%ngwv:band%ngwv))

  ! Find the boundaries of each bin.
  gaussian_bounds(:2*band%ngwv+1) = band%cref - 0.5_kind_phys*band%dc
  gaussian_bounds(2*band%ngwv+2) = band%cref(band%ngwv) + 0.5_kind_phys*band%dc

  ! Integral of the gaussian at bin interfaces (from the point to
  ! positive infinity).
  gaussian_bounds = &
       [( erfc(gaussian_bounds(l)/width)*height*width*sqrt(pi)/2._kind_phys, &
       l = 1, 2*band%ngwv+2 )]

  ! Get average in each bin using integral from right to left side.
  desc%src_tau = &
       (gaussian_bounds(:2*band%ngwv+1) - gaussian_bounds(2:)) / band%dc

  ! Prohibit wavenumber 0.
  desc%src_tau(0) = 0._kind_phys

end function gaussian_cm_desc

!==========================================================================
subroutine gw_cm_src(ncol, band, &
     desc_ksrc, &
     desc_kfront, &
     desc_frontgfc, &
     desc_src_tau, &
     u, v, frontgf, &
     src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  !-----------------------------------------------------------------------
  ! Driver for multiple gravity wave drag parameterization.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  !-----------------------------------------------------------------------

  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol

  ! Wavelengths triggered by frontogenesis.
  type(GWBand), intent(in) :: band

  ! Bandification of how to produce the gravity wave bandtrum.
  integer, intent(in)             ::  desc_ksrc
  integer, intent(in)             ::  desc_kfront
  real(kind_phys), intent(in)     ::  desc_frontgfc
  real(kind_phys), intent(in)     ::  desc_src_tau(-band%ngwv:band%ngwv)

  ! Midpoint zonal/meridional winds.
  real(kind_phys), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Frontogenesis function.
  real(kind_phys), intent(in) :: frontgf(:,:)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(kind_phys), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(kind_phys), intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(kind_phys), intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(kind_phys), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  !---------------------------Local Storage-------------------------------
  ! Column and wavenumber indices.
  integer :: k, l

  ! Whether or not to launch waves in this column.
  logical :: launch_wave(ncol)

  ! Zonal/meridional wind averaged over source region.
  real(kind_phys) :: usrc(ncol), vsrc(ncol)

  !------------------------------------------------------------------------
  ! Determine the source layer wind and unit vectors, then project winds.
  !------------------------------------------------------------------------

  ! Just use the source level interface values for the source wind speed
  ! and direction (unit vector).
  src_level = desc_ksrc
  tend_level = desc_ksrc
  usrc = 0.5_kind_phys*(u(:,desc_ksrc+1)+u(:,desc_ksrc))
  vsrc = 0.5_kind_phys*(v(:,desc_ksrc+1)+v(:,desc_ksrc))

  ! Get the unit vector components and magnitude at the surface.
  call get_unit_vector(usrc, vsrc, xv, yv, ubi(:,desc_ksrc+1))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, desc_ksrc
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:desc_ksrc) = midpoint_interp(ubm(:,1:desc_ksrc))

  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !-----------------------------------------------------------------------

  tau = 0._kind_phys

  ! GW generation depends on frontogenesis at specified level (may be below
  ! actual source level).
  launch_wave = (frontgf(:,desc_kfront) > desc_frontgfc)

  do l = -band%ngwv, band%ngwv
     where (launch_wave)
        tau(:,l,desc_ksrc+1) = desc_src_tau(l)
     end where
  end do

  ! Set phase speeds as reference speeds plus the wind speed at the source
  ! level.
  c = spread(band%cref, 1, ncol) + &
       spread(ubi(:,desc_ksrc+1), 2, 2*band%ngwv+1)

end subroutine gw_cm_src

end module gw_front
