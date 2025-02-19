module gw_convect

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

use ccpp_kinds, only:  kind_phys

implicit none
private
save

!!$public :: gw_beres_run
!!$public :: gw_beres_init
public :: BeresSourceDesc
public :: gw_beres_src

type :: BeresSourceDesc
   ! Whether wind speeds are shifted to be relative to storm cells.
   logical :: storm_shift
   ! Index for level where wind speed is used as the source speed.
   integer :: k
   ! Heating depths below this value [m] will be ignored.
   real(kind_phys) :: min_hdepth
   ! Table bounds, for convenience. (Could be inferred from shape(mfcc).)
   integer :: maxh
   integer :: maxuh
   ! Heating depths [m].
   real(kind_phys), allocatable :: hd(:)
   ! Table of source spectra.
   real(kind_phys), allocatable :: mfcc(:,:,:)
end type BeresSourceDesc

real(kind_phys), allocatable :: tau(:,:,:)  ! wave Reynolds stress

! gravity wave wind tendency for each wave
real(kind_phys), allocatable :: gwut(:,:,:)

!!$! Temperature tendencies from diffusion and kinetic energy.
!!$real(kind_phys), :: dttdf(state%ncol,pver)
!!$real(kind_phys), :: dttke(state%ncol,pver)

! Wave phase speeds for each column
real(kind_phys), allocatable :: phase_speeds(:,:)

contains

!!$!==========================================================================
!!$subroutine gw_beres_init(file_name,ttend_sh,ttend_dp,effgw)
!!$
!!$  character(len=*), intent(in) :: file_name
!!$  type(GWBand), intent(in) :: band
!!$
!!$  type(BeresSourceDesc), intent(inout) :: desc
!!$
!!$  type(file_desc_t) :: gw_file_desc
!!$
!!$  ! PIO variable ids and error code.
!!$  integer :: mfccid, hdid, stat
!!$
!!$  ! Number of wavenumbers in the input file.
!!$  integer :: ngwv_file
!!$
!!$  ! Full path to gw_drag_file.
!!$  character(len=cl) :: file_path
!!$
!!$  character(len=cl) :: msg
!!$
!!$  !----------------------------------------------------------------------
!!$  ! read in look-up table for source spectra
!!$  !-----------------------------------------------------------------------
!!$
!!$  call getfil(file_name, file_path)
!!$  call cam_pio_openfile(gw_file_desc, file_path, pio_nowrite)
!!$
!!$  ! Get HD (heating depth) dimension.
!!$
!!$  desc%maxh = get_pio_dimlen(gw_file_desc, "HD", file_path)
!!$
!!$  ! Get MW (mean wind) dimension.
!!$
!!$  desc%maxuh = get_pio_dimlen(gw_file_desc, "MW", file_path)
!!$
!!$  ! Get PS (phase speed) dimension.
!!$
!!$  ngwv_file = get_pio_dimlen(gw_file_desc, "PS", file_path)
!!$
!!$  ! Number in each direction is half of total (and minus phase speed of 0).
!!$  desc%maxuh = (desc%maxuh-1)/2
!!$  ngwv_file = (ngwv_file-1)/2
!!$
!!$  call shr_assert(ngwv_file >= band%ngwv, &
!!$       "gw_init_beres: PhaseSpeed in lookup table file does not cover the whole &
!!$       &spectrum implied by the model's ngwv. ")
!!$
!!$  ! Allocate hd and get data.
!!$
!!$  allocate(desc%hd(desc%maxh), stat=stat, errmsg=msg)
!!$
!!$  call shr_assert(stat == 0, &
!!$       "gw_init_beres: Allocation error (hd): "//msg// &
!!$       errMsg(__FILE__, __LINE__))
!!$
!!$  stat = pio_inq_varid(gw_file_desc,'HD',hdid)
!!$
!!$  call handle_pio_error(stat, &
!!$       'Error finding HD in: '//trim(file_path))
!!$
!!$  stat = pio_get_var(gw_file_desc, hdid, start=[1], count=[desc%maxh], &
!!$       ival=desc%hd)
!!$
!!$  call handle_pio_error(stat, &
!!$       'Error reading HD from: '//trim(file_path))
!!$
!!$  ! While not currently documented in the file, it uses kilometers. Convert
!!$  ! to meters.
!!$  desc%hd = desc%hd*1000._kind_phys
!!$
!!$  ! Allocate mfcc. "desc%maxh" and "desc%maxuh" are from the file, but the
!!$  ! model determines wavenumber dimension.
!!$
!!$  allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,&
!!$       -band%ngwv:band%ngwv), stat=stat, errmsg=msg)
!!$
!!$  call shr_assert(stat == 0, &
!!$       "gw_init_beres: Allocation error (mfcc): "//msg// &
!!$       errMsg(__FILE__, __LINE__))
!!$
!!$  ! Get mfcc data.
!!$
!!$  stat = pio_inq_varid(gw_file_desc,'mfcc',mfccid)
!!$
!!$  call handle_pio_error(stat, &
!!$       'Error finding mfcc in: '//trim(file_path))
!!$
!!$  stat = pio_get_var(gw_file_desc, mfccid, &
!!$       start=[1,1,ngwv_file-band%ngwv+1], count=shape(desc%mfcc), &
!!$       ival=desc%mfcc)
!!$
!!$  call handle_pio_error(stat, &
!!$       'Error reading mfcc from: '//trim(file_path))
!!$
!!$  call pio_closefile(gw_file_desc)
!!$
!!$  if (masterproc) then
!!$
!!$     write(iulog,*) "Read in source spectra from file."
!!$     write(iulog,*) "mfcc max, min = ", &
!!$          maxval(desc%mfcc),", ",minval(desc%mfcc)
!!$
!!$  endif
!!$
!!$  band = GWBand(pgwv, gw_dc, 1.0_kind_phys, wavelength_mid)
!!$
!!$! Allocate wavenumber fields.
!!$     allocate(tau(ncol,-band%ngwv:band%ngwv,pver+1),stat=istat)
!!$     call alloc_err(istat,'gw_tend','tau',ncol*(band%ngwv**2+1)*(pver+1))
!!$     allocate(gwut(ncol,pver,-band%ngwv:band%ngwv),stat=istat)
!!$     call alloc_err(istat,'gw_tend','gwut',ncol*pver*(band%ngwv**2+1))
!!$     allocate(phase_speeds(ncol,-band%ngwv:band%ngwv),stat=istat)
!!$     call alloc_err(istat,'gw_tend','phase_speeds',ncol*(band%ngwv**2+1))
!!$
!!$     ! Efficiency of gravity wave momentum transfer.
!!$     ! This is really only to remove the pole points.
!!$     where (pi/2._kind_phys - abs(lat(:ncol)) >= 4*epsilon(1._kind_phys))
!!$        effgw_dp = effgw_beres_dp
!!$        effgw_sh = effgw_beres_sh
!!$     elsewhere
!!$        effgw_dp = 0._kind_phys
!!$        effgw_sh = 0._kind_phys
!!$     end where
!!$
!!$
!!$end subroutine gw_beres_init
!!$!==========================================================================
!!$subroutine gw_beres_run(ncol, band, desc, u, v, &
!!$     netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
!!$     c, hdepth, maxq0, ptend_u, ptend_v, ptend_s, ptend_q )
!!$
!!$
!!$  !------------------------------------------------------------------
!!$  ! Convective gravity waves (Beres scheme, deep).
!!$  !------------------------------------------------------------------
!!$
!!$  tau=0._kind_phys
!!$  gwut=0._kind_phys
!!$  phase_speeds=0._kind_phys
!!$
!!$  ! Determine wave sources for Beres deep scheme
!!$  call gw_beres_src(ncol, band, desc, &
!!$       u, v, netdt(:ncol,:), zm, src_level, tend_level, tau, &
!!$       ubm, ubi, xv, yv, phase_speeds, hdepth, maxq0)
!!$
!!$  ! Solve for the drag profile with Beres source spectrum.
!!$  call gw_drag_prof(ncol, band, p, src_level, tend_level, dt, &
!!$          t, vramp,    &
!!$          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
!!$          effgw,   phase_speeds,       kvtt, q,  dse,  tau,  utgw,  vtgw, &
!!$          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
!!$          lapply_effgw_in=gw_apply_tndmax)
!!$
!!$     ! Project stress into directional components.
!!$     taucd = calc_taucd(ncol, band%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)
!!$
!!$     !  add the diffusion coefficients
!!$     do k = 1, pver+1
!!$        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
!!$     end do
!!$
!!$     ! Store constituents tendencies
!!$     do m=1, pcnst
!!$        do k = 1, pver
!!$           ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
!!$        end do
!!$     end do
!!$
!!$     ! Find momentum flux, and use it to fix the wind tendencies below
!!$     ! the gravity wave region.
!!$     call momentum_flux(tend_level, taucd, um_flux, vm_flux)
!!$     call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)
!!$
!!$     ! Add the momentum tendencies to the output tendency arrays.
!!$     do k = 1, pver
!!$        ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
!!$        ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
!!$     end do
!!$
!!$     ! Find energy change in the current state, and use fixer to apply
!!$     ! the difference in lower levels.
!!$     call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
!!$          ptend%v(:ncol,:), ptend%s(:ncol,:)+ttgw, de)
!!$     call energy_fixer(tend_level, p, de-flx_heat(:ncol), ttgw)
!!$
!!$     do k = 1, pver
!!$        ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
!!$     end do
!!$
!!$     ! Change ttgw to a temperature tendency before outputing it.
!!$     ttgw = ttgw / cpair
!!$     call gw_spec_outflds(beres_dp_pf, lchnk, ncol, band, phase_speeds, u, v, &
!!$          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
!!$          taucd)
!!$
!!$     ! Diagnostic outputs (convert hdepth to km).
!!$     call outfld('NETDT', ttend, pcols, lchnk)
!!$     call outfld('HDEPTH', hdepth/1000._kind_phys, ncol, lchnk)
!!$     call outfld('MAXQ0', maxq0, ncol, lchnk)
!!$
!!$end subroutine gw_beres_run
!!$!==========================================================================

  subroutine gw_beres_src(ncol, band, &
       desc_storm_shift, &
       desc_k, &
       desc_min_hdepth, &
       desc_maxh, &
       desc_maxuh, &
       desc_hd, &
       desc_mfcc, &
       u, v, &
     netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
     c, hdepth, maxq0)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
!
! The parameterization is assumed to operate only where water vapor
! concentrations are negligible in determining the density.
!
! Beres, J.H., M.J. Alexander, and J.R. Holton, 2004: "A method of
! specifying the gravity wave spectrum above convection based on latent
! heating properties and background wind". J. Atmos. Sci., Vol 61, No. 3,
! pp. 324-337.
!
!-----------------------------------------------------------------------
  use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  use gw_common, only: GWBand, pver, qbo_hdepth_scaling

!------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol

  ! Wavelengths triggered by convection.
  type(GWBand), intent(in) :: band

  ! Settings for convection type (e.g. deep vs shallow).
  logical, intent(in)             ::  desc_storm_shift
  integer, intent(in)             ::  desc_k
  real(kind_phys), intent(in)     ::  desc_min_hdepth
  integer, intent(in)             ::  desc_maxh
  integer, intent(in)             ::  desc_maxuh
  real(kind_phys), intent(in)     ::  desc_hd(:)
  real(kind_phys), intent(in)     ::  desc_mfcc(desc_maxh,-desc_maxuh:desc_maxuh,-band%ngwv:band%ngwv)

  ! Midpoint zonal/meridional winds.
  real(kind_phys), intent(in) :: u(:,:), v(:,:)
  ! Heating rate due to convection.
  real(kind_phys), intent(in) :: netdt(:,:)
  ! Midpoint altitudes.
  real(kind_phys), intent(in) :: zm(:,:)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(:)
  integer, intent(out) :: tend_level(:)

  ! Wave Reynolds stress.
  real(kind_phys), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(kind_phys), intent(out) :: ubm(:,:), ubi(:,:)
  ! Unit vectors of source wind (zonal and meridional components).
  real(kind_phys), intent(out) :: xv(:), yv(:)
  ! Phase speeds.
  real(kind_phys), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  ! Heating depth [m] and maximum heating in each column.
  real(kind_phys), intent(out) :: hdepth(:), maxq0(:)

!---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Zonal/meridional wind at roughly the level where the convection occurs.
  real(kind_phys) :: uconv(ncol), vconv(ncol)

  ! Maximum heating rate.
  real(kind_phys) :: q0(ncol)

  ! Bottom/top heating range index.
  integer  :: boti(ncol), topi(ncol)
  ! Index for looking up heating depth dimension in the table.
  integer  :: hd_idx(ncol)
  ! Mean wind in heating region.
  real(kind_phys) :: uh(ncol)
  ! Min/max wavenumber for critical level filtering.
  integer :: Umini(ncol), Umaxi(ncol)
  ! Source level tau for a column.
  real(kind_phys) :: tau0(-band%ngwv:band%ngwv)
  ! Speed of convective cells relative to storm.
  real(kind_phys) :: CS(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift

  ! Heating rate conversion factor.
  real(kind_phys), parameter :: CF = 20._kind_phys
  ! Averaging length.
  real(kind_phys), parameter :: AL = 1.0e5_kind_phys

  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------

  tau = 0.0_kind_phys
  hdepth = 0.0_kind_phys
  q0 = 0.0_kind_phys
  tau0 = 0.0_kind_phys

  !------------------------------------------------------------------------
  ! Determine wind and unit vectors approximately at the source level, then
  ! project winds.
  !------------------------------------------------------------------------

  ! Source wind speed and direction.
  uconv = u(:,desc_k)
  vconv = v(:,desc_k)

  ! Get the unit vector components and magnitude at the source level.
  call get_unit_vector(uconv, vconv, xv, yv, ubi(:,desc_k+1))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:pver) = midpoint_interp(ubm)

  !-----------------------------------------------------------------------
  ! Calculate heating depth.
  !
  ! Heating depth is defined as the first height range from the bottom in
  ! which heating rate is continuously positive.
  !-----------------------------------------------------------------------

  ! First find the indices for the top and bottom of the heating range.
  boti = 0
  topi = 0
  do k = pver, 1, -1
     do i = 1, ncol
        if (boti(i) == 0) then
           ! Detect if we are outside the maximum range (where z = 20 km).
           if (zm(i,k) >= 20000._kind_phys) then
              boti(i) = k
              topi(i) = k
           else
              ! First spot where heating rate is positive.
              if (netdt(i,k) > 0.0_kind_phys) boti(i) = k
           end if
        else if (topi(i) == 0) then
           ! Detect if we are outside the maximum range (z = 20 km).
           if (zm(i,k) >= 20000._kind_phys) then
              topi(i) = k
           else
              ! First spot where heating rate is no longer positive.
              if (.not. (netdt(i,k) > 0.0_kind_phys)) topi(i) = k
           end if
        end if
     end do
     ! When all done, exit.
     if (all(topi /= 0)) exit
  end do

  ! Heating depth in m.
  hdepth = [ ( (zm(i,topi(i))-zm(i,boti(i))), i = 1, ncol ) ]

  ! J. Richter: this is an effective reduction of the GW phase speeds (needed to drive the QBO)
  hdepth = hdepth*qbo_hdepth_scaling

  hd_idx = index_of_nearest(hdepth, desc_hd)

  ! hd_idx=0 signals that a heating depth is too shallow, i.e. that it is
  ! either not big enough for the lowest table entry, or it is below the
  ! minimum allowed for this convection type.
  ! Values above the max in the table still get the highest value, though.
  where (hdepth < max(desc_min_hdepth, desc_hd(1))) hd_idx = 0

  ! Maximum heating rate.
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        q0 = max(q0, netdt(:,k))
     end where
  end do

  !output max heating rate in K/day
  maxq0 = q0*24._kind_phys*3600._kind_phys

  ! Multipy by conversion factor
  q0 = q0 * CF

  if (desc_storm_shift) then

     ! Find the cell speed where the storm speed is > 10 m/s.
     ! Storm speed is taken to be the source wind speed.
     CS = sign(max(abs(ubm(:,desc_k))-10._kind_phys, 0._kind_phys), ubm(:,desc_k))

     ! Average wind in heating region, relative to storm cells.
     uh = 0._kind_phys
     do k = minval(topi), maxval(boti)
        where (k >= topi .and. k <= boti)
           uh = uh + ubm(:,k)/(boti-topi+1)
        end where
     end do

     uh = uh - CS

  else

     ! For shallow convection, wind is relative to ground, and "heating
     ! region" wind is just the source level wind.
     uh = ubm(:,desc_k)

  end if

  ! Limit uh to table range.
  uh = min(uh, real(desc_maxuh, kind_phys))
  uh = max(uh, -real(desc_maxuh, kind_phys))

  ! Speeds for critical level filtering.
  Umini =  band%ngwv
  Umaxi = -band%ngwv
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        Umini = min(Umini, nint(ubm(:,k)/band%dc))
        Umaxi = max(Umaxi, nint(ubm(:,k)/band%dc))
     end where
  end do

  Umini = max(Umini, -band%ngwv)
  Umaxi = min(Umaxi, band%ngwv)

  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !-----------------------------------------------------------------------
  ! Start loop over all columns.
  !-----------------------------------------------------------------------
  do i=1,ncol

     !---------------------------------------------------------------------
     ! Look up spectrum only if the heating depth is large enough, else set
     ! tau0 = 0.
     !---------------------------------------------------------------------

     if (hd_idx(i) > 0) then

        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------

        tau0 = desc_mfcc(hd_idx(i),nint(uh(i)),:)

        if (desc_storm_shift) then
           ! For deep convection, the wind was relative to storm cells, so
           ! shift the spectrum so that it is now relative to the ground.
           shift = -nint(CS(i)/band%dc)
           tau0 = eoshift(tau0, shift)
        end if

        ! Adjust magnitude.
        tau0 = tau0*q0(i)*q0(i)/AL

        ! Adjust for critical level filtering.
        tau0(Umini(i):Umaxi(i)) = 0.0_kind_phys

        tau(i,:,topi(i)+1) = tau0

     end if ! heating depth above min and not at the pole

  enddo

  !-----------------------------------------------------------------------
  ! End loop over all columns.
  !-----------------------------------------------------------------------

  ! Output the source level.
  src_level = topi
  tend_level = topi

  ! Set phase speeds; just use reference speeds.
  c = spread(band%cref, 1, ncol)

end subroutine gw_beres_src

! Short routine to get the indices of a set of values rounded to their
! nearest points on a grid.
function index_of_nearest(x, grid) result(idx)
  real(kind_phys), intent(in) :: x(:)
  real(kind_phys), intent(in) :: grid(:)

  integer :: idx(size(x))

  real(kind_phys) :: interfaces(size(grid)-1)
  integer :: i, n

  n = size(grid)
  interfaces = (grid(:n-1) + grid(2:))/2._kind_phys

  idx = 1
  do i = 1, n-1
     where (x > interfaces(i)) idx = i + 1
  end do

end function index_of_nearest

end module gw_convect
