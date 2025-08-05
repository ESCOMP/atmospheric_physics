module gw_convect

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

use ccpp_kinds, only:  kind_phys
use gw_common,      only: GWBand

implicit none
private
save

!jtpublic :: gw_beres_run
public :: gw_beres_init
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
   real(kind_phys), pointer :: hd(:)
   ! Table of source spectra.
   real(kind_phys), pointer :: mfcc(:,:,:)
end type BeresSourceDesc

! Beres settings and table.
type(BeresSourceDesc), public :: beres_dp_desc
type(BeresSourceDesc), public :: beres_sh_desc
type(GWBand)          :: band_mid

real(kind_phys), allocatable :: tau(:,:,:)  ! wave Reynolds stress

! gravity wave wind tendency for each wave
real(kind_phys), allocatable :: gwut(:,:,:)

!!$! Temperature tendencies from diffusion and kinetic energy.
!!$real(kind_phys), :: dttdf(state%ncol,pver)
!!$real(kind_phys), :: dttke(state%ncol,pver)

! Wave phase speeds for each column
real(kind_phys), allocatable :: phase_speeds(:,:)

contains

!==========================================================================
subroutine gw_beres_init(pver, pi, gw_drag_file_sh, gw_drag_file_dp, pref_edge, gw_dc, wavelength, pgwv, &
       use_gw_convect_dp,use_gw_convect_sh, masterproc, iulog, errmsg, errflg )

  use ccpp_io_reader, only: abstract_netcdf_reader_t, create_netcdf_reader_t


  integer, intent(in)                           :: pver
  real(kind_phys), intent(in)                   :: pi
  character(len=*), intent(in) :: gw_drag_file_sh, gw_drag_file_dp
  real(kind_phys), intent(in)                   :: pref_edge(:)
  real(kind_phys), intent(in)                   :: gw_dc
  real(kind_phys), intent(in)                   :: wavelength
  integer, intent(in)                           :: pgwv
  logical, intent(in)                           :: use_gw_convect_dp, use_gw_convect_sh

  logical, intent(in)                           :: masterproc
  integer, intent(in)                           :: iulog
  character(len=512), intent(out)               :: errmsg
  integer, intent(out)                          :: errflg

  ! Number of wavenumbers in the input file.
  integer :: ngwv_file

  ! Full path to gw_drag_file.

  class(abstract_netcdf_reader_t), allocatable :: reader
  type(BeresSourceDesc), pointer :: desc

  integer :: istat, k

  character(len=*), parameter :: sub = 'gw_beres_init'

  !----------------------------------------------------------------------
  ! read in look-up table for source spectra
  !-----------------------------------------------------------------------

  ! Initialize error variables
  errmsg =''
  errflg = 0

  if (use_gw_convect_dp .or. use_gw_convect_sh ) &
       band_mid = GWBand(pgwv, gw_dc, 1.0_kind_phys, wavelength)

  if (use_gw_convect_dp) then
     ! Set the deep scheme specification components.
     beres_dp_desc%storm_shift = .false.

     do k = 0, pver
        ! 700 hPa index
        if (pref_edge(k+1) < 70000._kind_phys) beres_dp_desc%k = k+1
     end do

     if (masterproc) then
        write (iulog,*) 'gw_beres_init: Beres deep level =',beres_dp_desc%k
     end if

     ! Don't use deep convection heating depths below this limit.
     ! This is important for QBO. Bad result if left at 2.5 km.
     beres_dp_desc%min_hdepth = 1000._kind_phys


     ! Check that deep gw file is set in namelist
     if (trim(gw_drag_file_dp) == "") then
        write(errmsg,'(a, a)') sub, "No gw_drag_file provided for Beres deep ", &
             "scheme. Set this via namelist."
             errflg = 1
        return
     end if

     call gw_init_beres_desc(gw_drag_file_dp, band_mid, beres_dp_desc, errmsg, errflg)
  end if

  if (use_gw_convect_sh) then
     ! Set the shallow scheme specification components.
     beres_sh_desc%storm_shift = .false.

     do k = 0, pver
        ! 900 hPa index
        if (pref_edge(k+1) < 90000._kind_phys) beres_sh_desc%k = k+1
     end do

     if (masterproc) then
        write (iulog,*) 'gw_beres_init: Beres shallow level =',beres_sh_desc%k
     end if

     ! Use all heating depths for shallow convection.
     beres_sh_desc%min_hdepth = 0._kind_phys

     ! Check that shallow gw file is set in namelist
     if (trim(gw_drag_file_sh) == "") then
        write(errmsg,'(a, a)') sub, "No gw_drag_file provided for Beres shallow ", &
             "scheme. Set this via namelist."
             errflg = 1
        return
     end if

     call gw_init_beres_desc(gw_drag_file_sh, band_mid, beres_sh_desc, errmsg, errflg)
  end if
  contains
    subroutine gw_init_beres_desc(file_path, band, desc, errmsg, errflg)
       type(GWBand), intent(in)                 :: band
       type(BeresSourceDesc), intent(inout)     :: desc
       character(len=*), intent(in)             :: file_path
       character(len=512), intent(out)          :: errmsg
       integer, intent(out)                     :: errflg

       integer                                  :: istat
       character(len=512)                       :: alloc_errmsg
       character(len=*), parameter :: sub = 'gw_init_beres_desc'
       real(kind_phys), pointer                 :: tmp_var1d(:)
       real(kind_phys), pointer                 :: file_mfcc(:,:,:) !is the lookup table from the file f(depth, wind, phase speed)

       ! Read Beres file.

       reader = create_netcdf_reader_t()

       ! Open file
       call reader%open_file(file_path, errmsg, errflg)
       if (errflg /= 0) then
          return !Error has occurred, so exit scheme
       end if

       ! Get HD (heating depth) dimension.

       call reader%get_var('HD',desc%hd , errmsg, errflg)
       if (errflg /= 0) then
          return !Error has occurred reading HD, so exit scheme
       end if
       desc%maxh = size(desc%hd)

       ! Get MW (mean wind) dimension.

       call reader%get_var('MW',tmp_var1d , errmsg, errflg)
       if (errflg /= 0) then
          return !Error has occurred reading MW, so exit scheme
       end if
       desc%maxuh = size(tmp_var1d)
       nullify(tmp_var1d)

       ! Get PS (phase speed) dimension.

       call reader%get_var('PS',tmp_var1d , errmsg, errflg)
       if (errflg /= 0) then
          return !Error has occurred reading PS, so exit scheme
       end if
       ngwv_file = size(tmp_var1d)
       nullify(tmp_var1d)

       ! Number in each direction is half of total (and minus phase speed of 0).
       desc%maxuh = (desc%maxuh-1)/2
       ngwv_file = (ngwv_file-1)/2

       ! Check for inconsistency between file and model variable
       if (ngwv_file < band%ngwv) then
          write(errmsg,'(a, a, i4, a, i4)') sub, "PhaseSpeed in lookup table file ", &
               ngwv_file, "does not cover the whole spectrum implied by the model ngwv. ", band%ngwv
          errflg = 1
          return
       end if

       ! While not currently documented in the file, it uses kilometers. Convert
       ! to meters.
       desc%hd = desc%hd*1000._kind_phys

       ! Allocate mfcc. "desc%maxh" and "desc%maxuh" are from the file, but the
       ! model determines wavenumber dimension.

       allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,&
            -band%ngwv:band%ngwv), stat=istat, errmsg=alloc_errmsg)

       if (istat/=0) then
          write(errmsg, '(a,a,a)') sub, ': ERROR allocating array: desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh, -band%ngwv:band%ngwv); message - ', trim(alloc_errmsg)
          errflg = 1
          return
       end if

       ! Get mfcc data.
       call reader%get_var('mfcc', file_mfcc, errmsg, errflg)
       if (errflg /= 0) then
          return !Error has occurred reading MFCC, so exit scheme
       end if

       desc%mfcc(:,-desc%maxuh:desc%maxuh,-band%ngwv:band%ngwv) = file_mfcc(:,:,ngwv_file-band%ngwv+1:)

       ! Close file
       call reader%close_file(errmsg, errflg)
       if (errflg /= 0) then
          return !Error has occurred while closing file, so exit scheme
       end if

       if (masterproc) then

          write(iulog,*) "Read in source spectra from file."
          write(iulog,*) "mfcc max, min = ", &
               maxval(desc%mfcc),", ",minval(desc%mfcc)

       endif

     end subroutine gw_init_beres_desc
   end subroutine gw_beres_init
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

  subroutine gw_beres_src(ncol, &
       desc, &
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
  use gw_common, only:  pver, qbo_hdepth_scaling

!------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  type(BeresSourceDesc) :: desc

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
  real(kind_phys), intent(out) :: tau(ncol,-band_mid%ngwv:band_mid%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(kind_phys), intent(out) :: ubm(:,:), ubi(:,:)
  ! Unit vectors of source wind (zonal and meridional components).
  real(kind_phys), intent(out) :: xv(:), yv(:)
  ! Phase speeds.
  real(kind_phys), intent(out) :: c(ncol,-band_mid%ngwv:band_mid%ngwv)

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
  real(kind_phys) :: tau0(-band_mid%ngwv:band_mid%ngwv)
  ! Speed of convective cells relative to storm.
  real(kind_phys) :: CS(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift

  ! Heating rate conversion factor.
  real(kind_phys), parameter :: CF = 20._kind_phys
  ! Averaging length.
  real(kind_phys), parameter :: AL = 1.0e5_kind_phys

  character(len=*), parameter :: sub = 'gw_beres_src'

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
  uconv = u(:,desc%k)
  vconv = v(:,desc%k)

  ! Get the unit vector components and magnitude at the source level.
  call get_unit_vector(uconv, vconv, xv, yv, ubi(:,desc%k+1))

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
           ! Detect if we are outside the top of range (where z = 20 km).
           if (zm(i,k) >= 20000._kind_phys) then
              boti(i) = k
              topi(i) = k
           else
              ! First spot where heating rate is positive.
              if (netdt(i,k) > 0.0_kind_phys) boti(i) = k
           end if
        end if
     end do
     ! When all done, exit
     if (all(boti /= 0)) exit
  end do

  do k = 1, pver
     do i = 1, ncol
        if (topi(i) == 0) then
                ! First spot where heating rate is positive.
              if ((netdt(i,k) > 0.0_kind_phys) .AND. (zm(i,k) <= 20000._kind_phys)) topi(i) = k-1
        end if
     end do
     ! When all done, exit
     if (all(topi /= 0)) exit
  end do

  ! Heating depth in m.
  hdepth = [ ( (zm(i,topi(i))-zm(i,boti(i))), i = 1, ncol ) ]

  ! J. Richter: this is an effective reduction of the GW phase speeds (needed to drive the QBO)
  hdepth = hdepth*qbo_hdepth_scaling

  hd_idx = index_of_nearest(hdepth, desc%hd)

  ! hd_idx=0 signals that a heating depth is too shallow, i.e. that it is
  ! either not big enough for the lowest table entry, or it is below the
  ! minimum allowed for this convection type.
  ! Values above the max in the table still get the highest value, though.
  where (hdepth < max(desc%min_hdepth, desc%hd(1))) hd_idx = 0

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

  if (desc%storm_shift) then

     ! Find the cell speed where the storm speed is > 10 m/s.
     ! Storm speed is taken to be the source wind speed.
     CS = sign(max(abs(ubm(:,desc%k))-10._kind_phys, 0._kind_phys), ubm(:,desc%k))

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
     uh = ubm(:,desc%k)

  end if

  ! Limit uh to table range.
  uh = min(uh, real(desc%maxuh, kind_phys))
  uh = max(uh, -real(desc%maxuh, kind_phys))

  ! Speeds for critical level filtering.
  Umini =  band_mid%ngwv
  Umaxi = -band_mid%ngwv
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        Umini = min(Umini, nint(ubm(:,k)/band_mid%dc))
        Umaxi = max(Umaxi, nint(ubm(:,k)/band_mid%dc))
     end where
  end do

  Umini = max(Umini, -band_mid%ngwv)
  Umaxi = min(Umaxi, band_mid%ngwv)

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

        tau0 = desc%mfcc(hd_idx(i),nint(uh(i)),:)

        if (desc%storm_shift) then
           ! For deep convection, the wind was relative to storm cells, so
           ! shift the spectrum so that it is now relative to the ground.
           shift = -nint(CS(i)/band_mid%dc)
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
  c = spread(band_mid%cref, 1, ncol)

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
