module gw_drag

!--------------------------------------------------------------------------
! CAM and WACCM gravity wave parameterizations were merged by Sean Patrick
! Santos in Summer 2013, and at the same time, gw_drag was split into
! various modules. This is the CAM interface and driver module. The below
! notes are for the old CAM and WACCM versions of gw_drag.
!--------------------------------------------------------------------------
! This file came from wa17 and was modified by Fabrizio: 07-02-2004
! Standard gw_drag with modification (6) of latitude profile of gw spectrum
!--------------------------------------------------------------------------
! Purpose:
!
! Module to compute the forcing due to parameterized gravity waves. Both an
! orographic and an internal source spectrum are considered.
!
! Author: Byron Boville
!
!--------------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys
  use gw_common, only: GWBand, handle_err

  implicit none
  private
  save

  ! Public CCPP-compliant interfaces.
  public :: gw_drag_init                  ! Initialization

  ! Whether or not to enforce an upper boundary condition of tau = 0.
  ! (Like many variables, this is only here to hold the value between
  ! the readnl phase and the init phase of the CAM physics; only gw_common
  ! should actually use it.)
  logical                              :: tau_0_ubc = .false.

  real(kind_phys)                      :: rearth   ! earth radius

  ! Whether or not to limit tau *before* applying any efficiency factors.
  logical                              :: gw_limit_tau_without_eff = .false.
  real(kind_phys)                      :: gw_prndl = 0.25_kind_phys

  ! Whether or not to apply tendency max
  real(kind_phys)                      :: gw_qbo_hdepth_scaling = 1._kind_phys ! heating depth scaling factor

  real(kind_phys), parameter :: unset_kind_phys = huge(1._kind_phys)

  ! Bottom level for frontal waves.
  integer :: kbot_front

  real(kind_phys)   :: gravit          ! gravitational acceleration (m s-2)
  real(kind_phys)   :: rair            ! Dry air gas constant     (J K-1 kg-1)
  real(kind_phys)   :: pi

contains

  ! Time independent initialization for multiple gravity wave parameterizations.
!> \section arg_table_gw_drag_init Argument Table
!! \htmlinclude gw_drag_init.html
  subroutine gw_drag_init( &
    iulog, &
    masterproc, &
    ktop, &
    pver, &
    gravit_in, &
    rair_in, &
    pi_in, &
    rearth_in, &
    pref_edge, &
    tau_0_ubc_nl, &
    gw_limit_tau_without_eff_nl, &
    gw_prndl_nl, &
    gw_qbo_hdepth_scaling_nl, &
    errmsg, &
    errflg)

    ! Host model dependency for interpolation.
    use interpolate_data, only: lininterp

    ! Underlying init subroutines
    use gw_common, only: gw_common_init

    use gw_common, only: wavelength_mid, wavelength_long

    integer, intent(in)             :: iulog
    logical, intent(in)             :: masterproc
    integer, intent(in)             :: ktop

    integer, intent(in)             :: pver
    real(kind_phys), intent(in)     :: gravit_in          ! gravitational acceleration (m s-2)
    real(kind_phys), intent(in)     :: rair_in            ! Dry air gas constant     (J K-1 kg-1)
    real(kind_phys), intent(in)     :: pi_in

    ! Whether or not to enforce an upper boundary condition of tau = 0.
    ! (Like many variables, this is only here to hold the value between
    ! the readnl phase and the init phase of the CAM physics; only gw_common
    ! should actually use it.)
    logical, intent(in)             :: tau_0_ubc_nl
    real(kind_phys), intent(in)     :: pref_edge(:)
    real(kind_phys), intent(in)             :: rearth_in   ! earth radius
    ! Whether or not to limit tau *before* applying any efficiency factors.
    logical, intent(in)             :: gw_limit_tau_without_eff_nl
    real(kind_phys), intent(in)             :: gw_prndl_nl
    real(kind_phys), intent(in)             :: gw_qbo_hdepth_scaling_nl

    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errflg

    ! Local variables
    character(len=*), parameter :: sub = 'gw_drag_init'

    integer :: k, l
    character(len=128) :: errstring

    ! Index for levels at specific pressures.
    integer :: topndx
    integer :: botndx
    integer :: stat

    ! Interpolated Newtonian cooling coefficients.
    real(kind_phys) :: alpha(pver + 1)

    ! Levels of pre-calculated Newtonian cooling (1/day).
    ! The following profile is digitized from:
    ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5

    integer, parameter :: nalph = 71
    real(kind_phys) :: alpha0(nalph) = [ &
                       0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, &
                       0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, &
                       0.1_kind_phys, 0.1_kind_phys, 0.10133333_kind_phys, 0.104_kind_phys, &
                       0.108_kind_phys, 0.112_kind_phys, 0.116_kind_phys, 0.12066667_kind_phys, &
                       0.126_kind_phys, 0.132_kind_phys, 0.138_kind_phys, 0.144_kind_phys, &
                       0.15133333_kind_phys, 0.16_kind_phys, 0.17_kind_phys, 0.18_kind_phys, &
                       0.19_kind_phys, 0.19933333_kind_phys, 0.208_kind_phys, 0.216_kind_phys, &
                       0.224_kind_phys, 0.232_kind_phys, 0.23466667_kind_phys, 0.232_kind_phys, &
                       0.224_kind_phys, 0.216_kind_phys, 0.208_kind_phys, 0.20133333_kind_phys, &
                       0.196_kind_phys, 0.192_kind_phys, 0.188_kind_phys, 0.184_kind_phys, &
                       0.18266667_kind_phys, 0.184_kind_phys, 0.188_kind_phys, 0.192_kind_phys, &
                       0.196_kind_phys, 0.19333333_kind_phys, 0.184_kind_phys, 0.168_kind_phys, &
                       0.152_kind_phys, 0.136_kind_phys, 0.12133333_kind_phys, 0.108_kind_phys, &
                       0.096_kind_phys, 0.084_kind_phys, 0.072_kind_phys, 0.061_kind_phys, &
                       0.051_kind_phys, 0.042_kind_phys, 0.033_kind_phys, 0.024_kind_phys, &
                       0.017666667_kind_phys, 0.014_kind_phys, 0.013_kind_phys, 0.012_kind_phys, &
                       0.011_kind_phys, 0.010333333_kind_phys, 0.01_kind_phys, 0.01_kind_phys, &
                       0.01_kind_phys, 0.01_kind_phys, 0.01_kind_phys &
                       ]

    ! Pressure levels that were used to calculate alpha0 (hPa).
    real(kind_phys) :: palph(nalph) = [ &
                       2.06115E-06_kind_phys, 2.74280E-06_kind_phys, 3.64988E-06_kind_phys, 4.85694E-06_kind_phys, &
                       6.46319E-06_kind_phys, 8.60065E-06_kind_phys, 1.14450E-05_kind_phys, 1.52300E-05_kind_phys, &
                       2.02667E-05_kind_phys, 2.69692E-05_kind_phys, 3.58882E-05_kind_phys, 4.77568E-05_kind_phys, &
                       6.35507E-05_kind_phys, 8.45676E-05_kind_phys, 0.000112535_kind_phys, 0.000149752_kind_phys, &
                       0.000199277_kind_phys, 0.000265180_kind_phys, 0.000352878_kind_phys, 0.000469579_kind_phys, &
                       0.000624875_kind_phys, 0.000831529_kind_phys, 0.00110653_kind_phys, 0.00147247_kind_phys, &
                       0.00195943_kind_phys, 0.00260744_kind_phys, 0.00346975_kind_phys, 0.00461724_kind_phys, &
                       0.00614421_kind_phys, 0.00817618_kind_phys, 0.0108801_kind_phys, 0.0144783_kind_phys, &
                       0.0192665_kind_phys, 0.0256382_kind_phys, 0.0341170_kind_phys, 0.0453999_kind_phys, &
                       0.0604142_kind_phys, 0.0803939_kind_phys, 0.106981_kind_phys, 0.142361_kind_phys, &
                       0.189442_kind_phys, 0.252093_kind_phys, 0.335463_kind_phys, 0.446404_kind_phys, &
                       0.594036_kind_phys, 0.790490_kind_phys, 1.05192_kind_phys, 1.39980_kind_phys, &
                       1.86273_kind_phys, 2.47875_kind_phys, 3.29851_kind_phys, 4.38936_kind_phys, &
                       5.84098_kind_phys, 7.77266_kind_phys, 10.3432_kind_phys, 13.7638_kind_phys, &
                       18.3156_kind_phys, 24.3728_kind_phys, 32.4332_kind_phys, 43.1593_kind_phys, &
                       57.4326_kind_phys, 76.4263_kind_phys, 101.701_kind_phys, 135.335_kind_phys, &
                       180.092_kind_phys, 239.651_kind_phys, 318.907_kind_phys, 424.373_kind_phys, &
                       564.718_kind_phys, 751.477_kind_phys, 1000._kind_phys &
                       ]

    errmsg = ''
    errflg = 0

    gravit = gravit_in
    rair = rair_in
    pi = pi_in
    tau_0_ubc = tau_0_ubc_nl
    rearth = rearth_in
    gw_limit_tau_without_eff = gw_limit_tau_without_eff_nl
    gw_prndl = gw_prndl_nl
    gw_qbo_hdepth_scaling = gw_qbo_hdepth_scaling_nl

    ! pre-calculated newtonian damping:
    !     * convert to s-1
    !     * ensure it is not smaller than 1e-6
    !     * convert palph from hpa to pa

    do k = 1, nalph
      alpha0(k) = alpha0(k)/86400._kind_phys
      alpha0(k) = max(alpha0(k), 1.e-6_kind_phys)
      palph(k) = palph(k)*1.e2_kind_phys
    end do

    call lininterp(alpha0, palph, nalph, alpha, pref_edge, pver + 1)
    if (masterproc) then
      write (iulog, *) 'gw_init: newtonian damping (1/day):'
      write (iulog, fmt='(a4,a12,a10)') ' k  ', '  pref_edge      ', &
        '  alpha   '
      do k = 1, pver + 1
        write (iulog, fmt='(i4,1e12.5,1f10.2)') k, pref_edge(k), &
          alpha(k)*86400._kind_phys
      end do
    end if

    if (masterproc) then
      write (iulog, *) 'gw_init: ktop = ', ktop
    end if

    ! Initialize subordinate modules.
    call gw_common_init(pver, &
                        tau_0_ubc, ktop, gravit, rair, alpha, gw_prndl, &
                        gw_qbo_hdepth_scaling, errstring)
  end subroutine gw_drag_init

!==========================================================================

! Add all history fields for a gravity wave spectrum source.
  subroutine gw_spec_addflds(prefix, scheme, band, history_defaults)
    use cam_history, only: addfld, add_default, register_vector_field

    !------------------------------Arguments--------------------------------

    ! One character prefix prepended to output fields.
    character(len=1), intent(in) :: prefix
    ! Gravity wave scheme name prepended to output field descriptions.
    character(len=*), intent(in) :: scheme
    ! Wave speeds.
    type(GWBand), intent(in) :: band
    ! Whether or not to call add_default for fields output by WACCM.
    logical, intent(in) :: history_defaults

    !---------------------------Local storage-------------------------------

    character(len=*), parameter :: sub = 'gw_spec_addflds'
    integer :: l
    ! 7 chars is enough for "-100.00"
    character(len=7)  :: fnum
    ! 10 chars is enough for "BTAUXSn32"
    character(len=10) :: dumc1x, dumc1y
    ! Allow 80 chars for description
    character(len=80) dumc2

    !-----------------------------------------------------------------------

    ! Overall wind tendencies.
    call addfld(trim(prefix)//'UTGWSPEC', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency - gravity wave spectrum')
    call addfld(trim(prefix)//'VTGWSPEC', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' V tendency - gravity wave spectrum')
    call register_vector_field(trim(prefix)//'UTGWSPEC', trim(prefix)//'VTGWSPEC')

    call addfld(trim(prefix)//'TTGWSPEC', (/'lev'/), 'A', 'K s-1', &
                trim(scheme)//' T tendency - gravity wave spectrum')

    ! Wind tendencies broken across five spectral bins.
    call addfld(trim(prefix)//'UTEND1', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency   c < -40')
    call addfld(trim(prefix)//'UTEND2', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency  -40 < c < -15')
    call addfld(trim(prefix)//'UTEND3', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency  -15 < c <  15')
    call addfld(trim(prefix)//'UTEND4', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency   15 < c <  40')
    call addfld(trim(prefix)//'UTEND5', (/'lev'/), 'A', 'm s-2', &
                trim(scheme)//' U tendency   40 < c ')

    ! Reynold's stress toward each cardinal direction, and net zonal stress.
    call addfld(trim(prefix)//'TAUE', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' Eastward Reynolds stress')
    call addfld(trim(prefix)//'TAUW', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' Westward Reynolds stress')
    call addfld(trim(prefix)//'TAUNET', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' E+W Reynolds stress')
    call addfld(trim(prefix)//'TAUN', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' Northward Reynolds stress')
    call addfld(trim(prefix)//'TAUS', (/'ilev'/), 'A', 'Pa', &
                trim(scheme)//' Southward Reynolds stress')

    ! Momentum flux in each direction.
    call addfld(trim(prefix)//'EMF', (/'lev'/), 'A', 'Pa', &
                trim(scheme)//' Eastward MF')
    call addfld(trim(prefix)//'WMF', (/'lev'/), 'A', 'Pa', &
                trim(scheme)//' Westward MF')
    call addfld(trim(prefix)//'NMF', (/'lev'/), 'A', 'Pa', &
                trim(scheme)//' Northward MF')
    call addfld(trim(prefix)//'SMF', (/'lev'/), 'A', 'Pa', &
                trim(scheme)//' Southward MF')

    ! Temperature tendency terms.
    call addfld(trim(prefix)//'TTGWSDF', (/'lev'/), 'A', 'K s-1', &
                trim(scheme)//' t tendency - diffusion term')
    call addfld(trim(prefix)//'TTGWSKE', (/'lev'/), 'A', 'K s-1', &
                trim(scheme)//' t tendency - kinetic energy conversion term')

    ! Gravity wave source spectra by wave number.
    do l = -band%ngwv, band%ngwv
      ! String containing reference speed.
      write (fnum, fmt='(f7.2)') band%cref(l)

      dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
      dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)
      dumc2 = trim(scheme)//" tau at c= "//trim(fnum)//" m s-1"
      call addfld(trim(dumc1x), (/'lev'/), 'A', 'Pa', dumc2)
      call addfld(trim(dumc1y), (/'lev'/), 'A', 'Pa', dumc2)

    end do

    if (history_defaults) then
      call add_default(trim(prefix)//'UTGWSPEC', 1, ' ')
      call add_default(trim(prefix)//'VTGWSPEC', 1, ' ')
      call add_default(trim(prefix)//'TTGWSPEC', 1, ' ')
      call add_default(trim(prefix)//'TAUE', 1, ' ')
      call add_default(trim(prefix)//'TAUW', 1, ' ')
      call add_default(trim(prefix)//'TAUNET', 1, ' ')
      call add_default(trim(prefix)//'TAUN', 1, ' ')
      call add_default(trim(prefix)//'TAUS', 1, ' ')
    end if

  end subroutine gw_spec_addflds

!==========================================================================

! Outputs for spectral waves.
  subroutine gw_spec_outflds(prefix, ncol, pver, band, phase_speeds, u, v, xv, yv, &
                             gwut, dttdf, dttke, tau, utgw, vtgw, ttgw, taucd)

    use gw_common, only: west, east, south, north

    ! One-character prefix prepended to output fields.
    character(len=1), intent(in) :: prefix
    ! Chunk and number of columns in the chunk.
    integer, intent(in) :: ncol
    integer, intent(in) :: pver
    ! Wave speeds.
    type(GWBand), intent(in) :: band
    ! Wave phase speeds for each column.
    real(kind_phys), intent(in) :: phase_speeds(ncol, -band%ngwv:band%ngwv)
    ! Winds at cell midpoints.
    real(kind_phys), intent(in) :: u(:, :)
    real(kind_phys), intent(in) :: v(:, :)
    ! Unit vector in the direction of wind at source level.
    real(kind_phys), intent(in) :: xv(:)
    real(kind_phys), intent(in) :: yv(:)
    ! Wind tendency for each wave.
    real(kind_phys), intent(in) :: gwut(ncol, pver, -band%ngwv:band%ngwv)
    ! Temperature tendencies from diffusion and kinetic energy.
    real(kind_phys), intent(in) :: dttdf(:, :)
    real(kind_phys), intent(in) :: dttke(:, :)
    ! Wave Reynolds stress.
    real(kind_phys), intent(in) :: tau(ncol, -band%ngwv:band%ngwv, pver)
    ! Zonal and meridional total wind tendencies.
    real(kind_phys), intent(in) :: utgw(:, :)
    real(kind_phys), intent(in) :: vtgw(:, :)
    ! Temperature tendencies.
    real(kind_phys), intent(in) :: ttgw(:, :)
    ! Reynolds stress for waves propagating in each cardinal direction.
    real(kind_phys), intent(in) :: taucd(ncol, pver + 1, 4)

    character(len=*), parameter :: sub = 'gw_spec_outflds'
    ! Indices
    integer :: i, k, l
    integer :: ix(ncol, -band%ngwv:band%ngwv), iy(ncol, -band%ngwv:band%ngwv)
    integer :: iu(ncol), iv(ncol)

    ! Zonal wind tendency, broken up into five bins.
    real(kind_phys) :: utb(ncol, pver, 5)
    ! Definition of the bin boundaries.
    real(kind_phys), parameter :: bounds(4) = (/-40._kind_phys, -15._kind_phys, &
                                                15._kind_phys, 40._kind_phys/)

    ! Momentum flux in the four cardinal directions.
    real(kind_phys) :: mf(ncol, pver, 4)

    ! Wave stress in zonal/meridional direction
    real(kind_phys) :: taux(ncol, -band%ngwv:band%ngwv, pver)
    real(kind_phys) :: tauy(ncol, -band%ngwv:band%ngwv, pver)

    ! Temporaries for output
    real(kind_phys) :: dummyx(ncol, pver)
    real(kind_phys) :: dummyy(ncol, pver)
    ! Variable names
    character(len=10) :: dumc1x, dumc1y

    ! Accumulate wind tendencies binned according to phase speed.

    utb = 0._kind_phys

    ! Find which output bin the phase speed corresponds to.
    ix = find_bin(phase_speeds)

    ! Put the wind tendency in that bin.
    do l = -band%ngwv, band%ngwv
      do k = 1, pver
        do i = 1, ncol
          utb(i, k, ix(i, l)) = utb(i, k, ix(i, l)) + gwut(i, k, l)
        end do
      end do
    end do

    ! Find just the zonal part.
    do l = 1, 5
      do k = 1, pver
        utb(:, k, l) = utb(:, k, l)*xv
      end do
    end do

!!$  call outfld(trim(prefix)//'UTEND1', utb(:,:,1), ncol, lchnk)
!!$  call outfld(trim(prefix)//'UTEND2', utb(:,:,2), ncol, lchnk)
!!$  call outfld(trim(prefix)//'UTEND3', utb(:,:,3), ncol, lchnk)
!!$  call outfld(trim(prefix)//'UTEND4', utb(:,:,4), ncol, lchnk)
!!$  call outfld(trim(prefix)//'UTEND5', utb(:,:,5), ncol, lchnk)
!!$
!!$  ! Output temperature tendencies due to diffusion and from kinetic energy.
!!$  call outfld(trim(prefix)//'TTGWSDF', dttdf / cpair, ncol, lchnk)
!!$  call outfld(trim(prefix)//'TTGWSKE', dttke / cpair, ncol, lchnk)
    ! Output tau broken down into zonal and meridional components.

    taux = 0._kind_phys
    tauy = 0._kind_phys

    ! Project phase_speeds, and convert each component to a wavenumber index.
    ! These are mappings from the wavenumber index of tau to those of taux
    ! and tauy, respectively.
    do l = -band%ngwv, band%ngwv
      ix(:, l) = c_to_l(phase_speeds(:, l)*xv)
      iy(:, l) = c_to_l(phase_speeds(:, l)*yv)
    end do

    ! Find projection of tau.
    do k = 1, pver
      do l = -band%ngwv, band%ngwv
        do i = 1, ncol
          taux(i, ix(i, l), k) = taux(i, ix(i, l), k) &
                                 + abs(tau(i, l, k)*xv(i))
          tauy(i, iy(i, l), k) = tauy(i, iy(i, l), k) &
                                 + abs(tau(i, l, k)*yv(i))
        end do
      end do
    end do

    do l = -band%ngwv, band%ngwv

      dummyx = taux(:, l, :)
      dummyy = tauy(:, l, :)

      dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
      dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)

!jt     call outfld(dumc1x,dummyx,ncol,lchnk)
!jt     call outfld(dumc1y,dummyy,ncol,lchnk)

    end do

    ! Output momentum flux in each cardinal direction.
    mf = 0._kind_phys

    do k = 1, pver

      ! Convert wind speed components to wavenumber indices.
      iu = c_to_l(u(:, k))
      iv = c_to_l(v(:, k))

      ! Sum tau components in each cardinal direction.
      ! Split west/east and north/south based on whether wave speed exceeds
      ! wind speed.
      do l = -band%ngwv, band%ngwv

        where (iu > l)
          mf(:, k, west) = mf(:, k, west) + taux(:, l, k)
        elsewhere
          mf(:, k, east) = mf(:, k, east) + taux(:, l, k)
        end where

        where (iv > l)
          mf(:, k, south) = mf(:, k, south) + tauy(:, l, k)
        elsewhere
          mf(:, k, north) = mf(:, k, north) + tauy(:, l, k)
        end where

      end do

    end do

!!$  call outfld(trim(prefix)//'WMF',mf(:,:,west),ncol,lchnk)
!!$  call outfld(trim(prefix)//'EMF',mf(:,:,east),ncol,lchnk)
!!$  call outfld(trim(prefix)//'SMF',mf(:,:,south),ncol,lchnk)
!!$  call outfld(trim(prefix)//'NMF',mf(:,:,north),ncol,lchnk)
!!$
!!$  ! Simple output fields written to history file.
!!$  ! Total wind tendencies.
!!$  call outfld (trim(prefix)//'UTGWSPEC', utgw , ncol, lchnk)
!!$  call outfld (trim(prefix)//'VTGWSPEC', vtgw , ncol, lchnk)
!!$  call outfld (trim(prefix)//'TTGWSPEC', ttgw , ncol, lchnk)
!!$
!!$  ! Tau in each direction.
!!$  call outfld (trim(prefix)//'TAUE', taucd(:,:,east), ncol, lchnk)
!!$  call outfld (trim(prefix)//'TAUW', taucd(:,:,west), ncol, lchnk)
!!$  call outfld (trim(prefix)//'TAUN', taucd(:,:,north), ncol, lchnk)
!!$  call outfld (trim(prefix)//'TAUS', taucd(:,:,south), ncol, lchnk)
!!$
!!$  call outfld (trim(prefix)//'TAUNET', taucd(:,:,east)+taucd(:,:,west), &
!!$       ncol, lchnk)

  contains

    ! Convert a speed to a wavenumber between -ngwv and ngwv.
    elemental function c_to_l(c) result(l)
      real(kind_phys), intent(in) :: c

      integer :: l

      l = min(max(int(c/band%dc), -band%ngwv), band%ngwv)

    end function c_to_l

  end subroutine gw_spec_outflds

!==========================================================================

! Generates names for tau output across the wave spectrum (e.g.
! BTAUXSn01 or TAUYSp05).
! Probably this should use a wavenumber dimension on one field rather
! than creating a ton of numbered fields.
  character(len=9) pure function tau_fld_name(l, prefix, x_not_y)
    ! Wavenumber
    integer, intent(in) :: l
    ! Single-character prefix for output
    character(len=1), intent(in) :: prefix
    ! X or Y?
    logical, intent(in) :: x_not_y

    character(len=2) :: num_str

    tau_fld_name = trim(prefix)

    tau_fld_name = trim(tau_fld_name)//"TAU"

    if (x_not_y) then
      tau_fld_name = trim(tau_fld_name)//"XS"
    else
      tau_fld_name = trim(tau_fld_name)//"YS"
    end if

    if (l < 0) then
      tau_fld_name = trim(tau_fld_name)//"n"
    else
      tau_fld_name = trim(tau_fld_name)//"p"
    end if

    write (num_str, '(I2.2)') abs(l)

    tau_fld_name = trim(tau_fld_name)//num_str

  end function tau_fld_name

end module gw_drag
